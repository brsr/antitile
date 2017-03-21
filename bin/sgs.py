#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Subdivide a tiling or polyhedra using a similar grid
"""
import argparse
import numpy as np
from sys import stdin
from antitile import off, tiling, breakdown, projection, xmath

def subdiv(base, freq=(2,0), proj='flat'):
    faces_dict = base.faces_by_size
    if any(x > 4 for x in faces_dict.keys()):
        raise ValueError("Tiling contains at least one face with more than"
                         " 4 sides. Try triangulating those faces first.")
    bkdn_t = breakdown.Breakdown(*freq, 't')
    bkdn_q = breakdown.Breakdown(*freq, 'q')

    newverts = []
    #base.orient_faces()
    for face in base.faces:
        face_n = len(face)
        face_pts = base.vertices[face]
        if face_n < 2:
            #just ignore any edges or vertices that show up in the face list
            continue
        elif face_n == 3:
            bkdn = bkdn_t
            newvert = projection.tri_bary(bkdn.coord, face_pts)
        elif face_n == 4:
            bkdn = bkdn_q
            newvert = projection.square_to_quad(bkdn.coord[:, np.newaxis], face_pts)
        newverts.append(newvert)
    return np.concatenate(newverts, axis=0)


if __name__ == "__main__":
    halp = """Similar grid subdivision: subdivide a tiling or polyhedron
        with a grid of similar triangles or squares."""
    parser = argparse.ArgumentParser(description=halp)
    parser.add_argument("filename", nargs='?',
                        help="Input file. Reads from stdin if not given.")
    halp = """Breakdown frequency. Can be a pair of integers
        separated by a space, like 2 1, or a single integer n, which is
        understood to be n 0. Default is 2."""
    parser.add_argument("-f", '--freq', nargs='+', type=int, help=halp,
                        default=('2','0'))
    halp = """Projection family. Default is flat, which is the only
    one likely to work with non-spherical polyhedra. """
#    May be:
#        flat: Flat subdivision of each face (Method 1 in geodesic dome jargon)
#        slerp: Spherical linear interpolation (or related method)
#        areal: Areal coordinates on the sphere (triangular faces only)
#        gc: Intersection of great circles (Method 2 in geodesic dome jargon)
#        gcv: Minor variation of gc
    parser.add_argument("-p", '--projection', default='flat', help=halp,
                        choices=['flat','slerp','areal','gc','gcv'])
    parser.add_argument("-n", "--normalize", action="store_true",
                        help="Normalize vertices onto the unit sphere")
    halp = """Projection constant for naive slerp. May be a float or a
    string from the list below. If a string is given, it will optimize k
    based on the specified measurement of the polyhedron.
    Default is 1. Ignored unless -p=slerp and the grid has triangular faces."""
#        energy: Minimizes the Thompson energy of the points.
#        fill: Maximizes the fill ratio of the polyhedron wrt the unit sphere.
#        edges: Minimizes the difference in edge length.
#        aspect: Minimizes the aspect ratio of euclidean triangles.
#        faces: Minimizes the difference in area between faces.
#        angle: Minimizes the difference in central angle between points
#            sharing edges. (On a unit sphere, same as the spherical distance.)
#        angle_aspect: Minimizes the aspect ratio of spherical triangles.
#        solid_angle: Minimizes the difference in solid angle between faces.
    parser.add_argument("-k", default=1, help=halp,
                        choices=['energy', 'fill', 'edges', 'aspect', 'faces',
                                 'angle', 'angle_aspect', 'solid_angle'])
    args = parser.parse_args()
    freq = tuple(args.freq)
    if len(freq) > 2:
        raise ValueError("Too many items in frequency option")
    elif len(freq) < 2:
        freq = (freq[0], 0)
    file = open(args.filename) if args.filename else stdin        
    if freq == (1, 0):#identity, just give input back as output
        with file as f:
            print(f.read())
    elif freq[0] <= 0 or freq[1] < 0:
        raise ValueError("Bad value " + str(freq) + " given for frequency. "
                         "First value must be greater than zero, "
                         "second must be greater or equal to zero.")
    else:
        with file as f:
            vertices, faces, fc, e, ec, v, vc = off.load_off(f)
        base = tiling.Tiling(vertices, faces)
        result = subdiv(base, freq, args.projection)
        if args.normalize:
            base.vertices = xmath.normalize(vertices)

        print(off.write_off(result, []))