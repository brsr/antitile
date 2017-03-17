#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Subdivide a tiling or polyhedra
"""
import argparse
import numpy as np
from antitile import off, tiling, breakdown, projection

def subdiv(base, freq={'t': (2,0), 'q': (2,0)}, proj='flat'):
    bkdns = {shape: breakdown.Breakdown(*freq, shape)
                for (shape, freq) in freq.items()}
    newverts = []
    #base.orient_faces()
    for face in base.faces:
        face_n = len(face)
        face_pts = base.vertices[face]
        if face_n > 4:
            raise ValueError("Tiling contains at least one face with more than"
                             " 4 sides. Try triangulating those faces first.")
        elif face_n < 2:
            #just ignore any edges or vertices that show up in the face list
            continue
        elif face_n == 3:
            bkdn = bkdns['t']
            newvert = projection.tri_bary(bkdn.coord, face_pts)
        elif face_n == 4:
            bkdn = bkdns['q']
            newvert = projection.square_to_quad(bkdn.coord[:, np.newaxis], face_pts)
        newverts.append(newvert)
    return np.concatenate(newverts, axis=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
                                     "Subdivide a tiling or polyhedron")
    parser.add_argument("filename", help="Input file")
    halp = ("Breakdown frequency to use. Default is 2. If a single integer,"
            "performs a class I subdivision with that frequency.")
    parser.add_argument("-f", default=2, help=halp)
    halp = ("Projection family to use. Default is flat, which is the only"
            "one likely to work with non-spherical polyhedra.")
    parser.add_argument("-p", help=halp)
    args = parser.parse_args()
    filename = args.filename
    vertices, faces, fc, e, ec, v, vc = off.load_off(filename)
    base = tiling.Tiling(vertices, faces)
    result = subdiv(base)
    print(off.write_off(result, []))