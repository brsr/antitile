# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 15:01:28 2017

@author: Bstone
"""

import argparse
from sys import stdin
import numpy as np
import pandas as pd
import warnings
from scipy import sparse
from antitile import off, tiling

def ca_step(state, rule, adj):
    x, y, _ = sparse.find(adj)
    nbx = state[x]
    px = pd.DataFrame(data=nbx, index=y)
    neighbors = px.groupby(px.index).sum()
    result = np.squeeze(neighbors.values) - 128*state
    return np.in1d(result, rule)

def ruleparse(string):
    return [int(x) for x in string.split(',')]

DESCRIPTION = """Colors tiling using semitotalistic cellular automata.
Colors indexes of the input file is used as the initial condition,
interpreting 0 as dead and anything else as alive. If colors are not
given, initial condition is assigned randomly with 50% live cells."""

def main():
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("filename", nargs='?',
                        help="Input file. Reads from stdin if not given.")
    parser.add_argument("-n", help="Number of steps", default=100, type=int)
    parser.add_argument("-b", help="Birth rule, comma-separated",
                        default=[1, 5, 6], type=ruleparse)
    parser.add_argument("-s", help="Survival rule, comma-separated",
                        default=[1, 2], type=ruleparse)
    parser.add_argument("-l", help="Lookback", type=int, default=12)
    parser.add_argument("-v", help="Use vertex adjacency instead of face",
                        action='store_true')
    parser.add_argument("-o", help="Output prefix.", default="cellular")
    args = parser.parse_args()
    file = open(args.filename) if args.filename else stdin
    with file as f:
        vertices, faces, fc, e, ec, verts, vc = off.load_off(f)
    if fc is None:
        fc = np.random.randint(2, size=len(faces))
    elif len(fc.shape) > 1:
        raise ValueError("Faces contain color values, not color indexes")
    if vc is None:
        vc = np.random.randint(2, size=len(vertices))
    elif len(vc.shape) > 1:
        raise ValueError("Vertices contain color values, not color indexes")
    if verts is None:
        verts = np.arange(len(vertices))
    init = vc if args.v else fc
    poly = tiling.Tiling(vertices, faces)
    adj = poly.vertex_adjacency if args.v else poly.face_adjacency
    rule = np.array(args.b + [i - 128 for i in args.s], dtype=np.int8)
    state = np.zeros((args.n + 1, len(init)), dtype=bool)
    state[0] = init
    lookback = args.l
    width = int(np.ceil(np.log10(args.n + 2)))
    file_template = args.o + "{:0" + str(width) + "d}" + ".off"
    for i in range(args.n):
        this_state = ca_step(state[i], rule, adj)
        state[i+1] = this_state
        if args.v:
            vc = this_state.astype(int)
        else:
            fc = this_state.astype(int)
        string = off.write_off(vertices, faces, facecolors=fc, edges=e,
                               edgecolors=ec, vertexcolors=vc)
        fn = file_template.format(i + 1)
        with open(fn, 'w') as f:
            f.write(string)
        start = max(i - lookback, 0)
        if np.any(np.all(this_state == state[start:i + 1], axis=-1)):
            string = "reached steady state at {}, breaking out early"
            warnings.warn(string.format(i+1))
            break

if __name__ == "__main__":
    main()
