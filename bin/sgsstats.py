#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Statistics of polyhedra, with a focus on using similar grid subdivision
polyhedra to approximate the sphere 
"""
import argparse
import csv
import numpy as np
from sys import stdin
from antitile import off, tiling

def main():
    desc = "Statistics of a tiling that are relevant to use as a grid"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filename", nargs='*',
                        help="Input files. Reads from stdin if not given.")
    parser.add_argument("-s", action="store_true",
                        help="Suppress spherical measures")
    parser.add_argument("-c", help="Write to csv file")

    args = parser.parse_args()
    header = ['filename', 'n_v', 'n_f', 'energy', 'cog', 'bent_min', 'bent_max',
                  'edge_min', 'edge_max',
                  'aspect_ratio_min', 'aspect_ratio_max',
                  'faces_min', 'faces_max']
    if not args.s:
            header += ['angle_min', 'angle_max',
                       'angle_aspect_min', 'angle_aspect_max',
                       'solid_angle_min', 'solid_angle_max']

    lines = [header]
    handles = [(fn, open(fn)) for fn in args.filename]
    if not handles:
        handles = [('stdin',  stdin)]
    for fn, h in handles:
        with h as handle:
            vertices, faces, fc, e, ec, v, vc = off.load_off(handle)
        poly = tiling.Tiling(vertices, faces)
        n_v = len(vertices)
        n_f = len(faces)
        energy = tiling.energy(vertices)
        norm_cog = tiling.center_of_gravity(vertices)
        bentness = tiling.bentness(vertices, poly)
        bent_min, bent_max = bentness.min(), bentness.max()
        edges = tiling.edge_length(vertices, poly.edges)
        edge_min, edge_max = edges.min(), edges.max()
        aspect = tiling.aspect_ratio(vertices, poly)
        aspect_min, aspect_max = aspect.min(), aspect.max()
        faces = tiling.face_area(vertices, poly)
        face_min, face_max = faces.min(), faces.max()
        values = [fn, n_v, n_f, energy, norm_cog,
                  bent_min, bent_max,
                  edge_min, edge_max,
                  aspect_min, aspect_max,
                  face_min, face_max]
        print('---')
        print('File: ', fn)
        print('Vertices, faces: {}, {}'.format(n_v, n_f))
        print('Thomson energy: ', energy)
        if not np.isclose(0, norm_cog):
            print('Distance of center of gravity from center:\t', norm_cog)
        print('\t\t\t    minimum |maximum |ratio')
        print('Edge length:\t\t    {:<,F}|{:<,F}|{:<,F}'.format(
              edge_min, edge_max, edge_max/edge_min))
        print('Aspect ratio:\t\t    {:<,F}|{:<,F}|{:<,F}'.format(
              aspect_min, aspect_max, aspect_max/aspect_min))
        if np.isclose(0, bent_max, atol = 1E-6):
            print('Face area:\t\t    {:<,F}|{:<,F}|{:<,F}'.format(
                    face_min, face_max, face_max/face_min))
        else:
            print('Face bentness: \t\t    {:<,F}|{:<,F}|{:<,F}'.format(
                    bent_min, bent_max, bent_max/bent_min))
        if not args.s:
            edges = tiling.edge_length(vertices, poly.edges, spherical=True)
            edge_min, edge_max = edges.min(), edges.max()
            aspect = tiling.aspect_ratio(vertices, poly, spherical=True)
            aspect_min, aspect_max = aspect.min(), aspect.max()
            faces = tiling.face_area(vertices, poly, spherical=True)
            face_min, face_max = faces.min(), faces.max()
            values += [edge_min, edge_max,
                       aspect_min, aspect_max,
                       face_min, face_max]
            print('Central angle:\t\t    {:<,F}|{:<,F}|{:<,F}'.format(
                  edge_min, edge_max, edge_max/edge_min))
            print('Central angle aspect ratio: {:<,F}|{:<,F}|{:<,F}'.format(
                  aspect_min, aspect_max, aspect_max/aspect_min))
            print('Solid angle:\t\t    {:<,F}|{:<,F}|{:<,F}'.format(
                  face_min, face_max, face_max/face_min))
        lines.append(values)

    if args.c:
        with open(args.c, 'w') as f:
            writer = csv.writer(f)
            writer.writerows(lines)

if __name__ == "__main__":
    main()
