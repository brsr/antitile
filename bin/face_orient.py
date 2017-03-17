#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Orient faces of a tiling in counterclockwise order
"""
import argparse
from antitile import tiling, off

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
                                     "Subdivide a tiling or polyhedron")
    parser.add_argument("filename", help="Input file")
    args = parser.parse_args()
    filename = args.filename
    vertices, faces, fc, e, ec, v, vc = off.load_off(filename)
    base = tiling.Tiling(vertices, faces)
    base.orient_faces()
    print(off.write_off(base.vertices, base.faces))
