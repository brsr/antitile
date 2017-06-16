#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Nslerp surfaces.
"""

import argparse
import numpy as np
from antitile import projection, flat, off

np.seterr(divide='ignore', invalid='ignore')

TRIANGLE = np.array([[0, 1],
                     [-np.sqrt(3)/2, -0.5],
                     [np.sqrt(3)/2, -0.5],
                     [0, 1]])

SQUARE = np.array([[ 1,  0],
                   [ 0,  1],
                   [-1,  0],
                   [ 0, -1],
                   [ 1,  0]])

SHAPES = {3: TRIANGLE,
          4: SQUARE}

CENTER = np.array([0, 0, 1])
#dihedral angles:
#    icosahedron: arctan(2)
#    octahedron: np.pi/2
#    tetrahedron: arccos(-1/3)
#    3-dihedron: 2*np.pi/3
#    cube: arccos(1/3)
#    4-dihedron: np.pi/2

HEIGHTS = {'ico': np.sqrt((5 + 2*np.sqrt(5))/15),
           'oct': 1/np.sqrt(3),
           'tet': 1/3,
           'cube': 1/np.sqrt(2),
           'di': 0}

def base_3d(z, shape=TRIANGLE):
    """
    Gives a base shape centered on [0,0,1].
    """
    r = np.sqrt(np.abs(1 - z**2))
    return np.append(r*shape, z*np.ones((shape.shape[0], 1)), axis=-1)

def zerotoone(string):
    if string in HEIGHTS:
        return HEIGHTS[string]
    x = float(string)
    if x < 0 or x > 1:
        msg = "z must be between 0 and 1 inclusive"
        raise argparse.ArgumentTypeError(msg)
    return x

def main():
    parser = argparse.ArgumentParser(description="Nslerp surfaces")
    parser.add_argument("factor", nargs='?', default=2, type=int,
                        help="Factor for extent of surface")
    parser.add_argument("freq", nargs='?', default=10, type=int,
                        help="Frequency of surface division")    
    parser.add_argument("-z", default=1/np.sqrt(3), type=zerotoone,
                        help="height of breakdown face, in [0, 1], default "
                             "1/sqrt(3). or a text string: " +
                             ', '.join(HEIGHTS))
    parser.add_argument("-q", action="store_true",
                        help="use quadrilateral instead of triangle")    
    parser.add_argument("-t", action="store_true",
                        help="use nslerp2 instead of nslerp (quads only)")    
        
    args = parser.parse_args()
    base = base_3d(args.z, shape=SQUARE if args.q else TRIANGLE)[:-1]
    n = (2*args.factor+1)*args.freq
    thing = flat.FlatTiling(n, n, shape=4 if args.q else 3)
    faces = thing.faces
    vx = thing.vertices/args.freq
    if args.q:
        offset = np.array([args.factor, args.factor, 0])
        coord = (vx - offset)[..., :-1]
        if args.t:
            vertices = projection.square_naive_slerp_2(coord, base)
        else:
            vertices = projection.square_naive_slerp(coord, base)
    else:
        offset = np.array([args.factor, args.factor,
                           -2*args.factor-1])
        coord = vx - offset
        
        vertices = projection.tri_naive_slerp(coord, base)
    vcolor = np.all((0 <= coord) & (coord <= 1), axis=-1)
    facecolor = np.sum(vcolor[faces], axis=-1).astype(int)
    result = off.write_off(vertices, faces, facecolors=facecolor)
    result += "#extends {}, frequency {}\n".format(args.factor, args.freq)
    result += "#base height = {}\n".format(args.z)
    if args.q and args.t:
        msg = "#nslerp2 (quad)"
    elif args.q:
        msg = "#nslerp (quad)"
    else:
        msg = "#nslerp (tri)"
    result += msg + "\n"
    print(result)
    

if __name__ == "__main__":
    main()
