#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Subdivide a tiling or polyhedra using a similar grid
"""
import argparse
import warnings
from sys import stdin
from antitile import gcopoly, off, tiling, projection

DESCRIPTION = """Goldberg-Coxeter operation: subdivide a tiling or
polyhedron with a grid of similar triangles or squares."""
EPILOG = """To use on a non-spherical polyhedron or tiling,
specify -n -p=flat"""
FILENAME = """Input file. Reads from stdin if not given. Faces should be
oriented counterclockwise (even for Class I and II)."""
FREQ_A = """First breakdown frequency. Default is 2."""
FREQ_B = """Second breakdown frequency. Default is 0."""
PROJ = """Mapping family to use. Default is (gn)omonic. The mapping for
triangles and quadrilaterals can be specified separately as a comma-separated
pair: -m ns, n2. """
ADJ = ("""Projection constant. May be a float or a string from the list
below. If a string is given, it will optimize k based on the specified
measurement of the polyhedron. Ignored unless -m=""" +
       ', '.join(projection.PARALLEL) + ". Default is 1. String values can be "
       + ', '.join(n for n in gcopoly.MEASURES))
TWEAK = """Makes a tweak to certian methods: see docs for details."""
FACTOR = """Factors the frequency and performs repeated subdivision using
those factors. Smaller subdivisions (by norm) are applied first. (If the
factors are powers of (2,0), this is Method 3 in geodesic dome parlance.)"""

def nonnegativeint(string, lowest=0):
    """Non-negative integer type for argparse"""
    x = int(string)
    if x < lowest:
        msg = "must be greater than or equal to {y}"
        raise argparse.ArgumentTypeError(msg.format(y=lowest))
    return x

def posint(string):
    """Positive integer type for argparse"""
    return nonnegativeint(string, 1)

def kparser(string):
    """Parse the k-factor argument"""
    return string if string in gcopoly.MEASURES else float(string)

def main():
    parser = argparse.ArgumentParser(description=DESCRIPTION, epilog=EPILOG)
    parser.add_argument("filename", nargs='?', help=FILENAME)
    parser.add_argument("-a", help=FREQ_A, default=2, type=posint)
    parser.add_argument("-b", help=FREQ_B, default=0,
                        type=nonnegativeint)
    parser.add_argument("-m", '--mapping', default='flat', help=PROJ,
                        choices=projection.PROJECTIONS)
    parser.add_argument("-n", "--no_normalize", action="store_true",
                        help="Don't normalize vertices onto the unit sphere")
    parser.add_argument("-k", default=1, help=ADJ, type=kparser)
    parser.add_argument("-t", "--tweak", action="store_true", help=TWEAK)
    parser.add_argument("-f", "--factor", action="store_true", help=FACTOR)

    args = parser.parse_args()
    frequency = (args.a, args.b)
    file = open(args.filename) if args.filename else stdin
    if frequency == (1, 0):#identity, just give input back as output
        with file as f:
            print(f.read())
    else:
        with file as f:
            vertices, faces, fc, _e, _ec, _v, _vc = off.load_off(f)
        base = tiling.Tiling(vertices, faces)
        classIlist = ['edges', 'aspect', 'angle', 'angle_aspect']
        is_tri_grid = all(base.face_size <= 3)
        if args.b == 0 and args.k in classIlist and is_tri_grid:
            warnings.warn(args.k + " is optimal for a large range of k for "
                          "class I subdivision on triangle grid. "
                          "using -k=energy instead.")
            k = 'energy'
        elif args.k == 'bent' and is_tri_grid:
            warnings.warn("bentness always == 0 for triangular faces. "
                          "Using -k=energy instead.")
            k = 'energy'
        else:
            k = args.k
        if args.factor:
            poly = gcopoly.build_gco_rep(base, frequency, args.projection,
                                     tweak=args.tweak, k=k,
                                     normalize=not args.no_normalize)
        else:
            poly = gcopoly.build_gco(base, frequency, args.projection,
                                 tweak=args.tweak, k=k,
                                 normalize=not args.no_normalize)
        vertcolor = poly.group.astype(int)
        facecolor = poly.face_group.astype(int)
        result = off.write_off(poly.vertices, poly.faces,
                               facecolors=facecolor,
                               vertexcolors=vertcolor)
        result += '#frequency = {}\n'.format(frequency)
        if args.filename:
            result += '#input file = {}\n'.format(args.filename)
        result += '#projection = {}\n'.format(args.projection)
        if args.projection in projection.PARALLEL:
            result += '#k = {}\n'.format(poly.k)
        if args.projection in projection.PARALLEL + ['gc']:
            result += '#tweak = {}\n'.format(args.tweak)
        result += '#normalized = {}\n'.format(not args.no_normalize)
        print(result)

if __name__ == "__main__":
    main()
