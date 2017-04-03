#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Subdivide a tiling or polyhedra using a similar grid
"""
import argparse
from sys import stdin
from antitile import sgs, off, tiling, xmath, projection

DESCRIPTION = """Similar grid subdivision: subdivide a tiling or
polyhedron with a grid of similar triangles or squares."""
EPILOG = """To use on a non-spherical polyhedron or tiling,
specify -n -p=flat"""
FILENAME = """Input file. Reads from stdin if not given. Faces should be
oriented counterclockwise (even for Class I and II)."""
FREQ_A = """First breakdown frequency. Default is 2."""
FREQ_B = """Second breakdown frequency. Default is 0."""
PROJ = """Projection family. Default is flat. areal is only valid on
    triangular faces. disk is only valid on dihedra."""
#    May be:
#        flat: Flat subdivision of each face (Method 1 in geodesic dome jargon)
#        slerp: Spherical linear interpolation (or related method)
#        areal: Areal coordinates on the sphere (triangular faces only)
#        gc: Intersection of great circles (Method 2 in geodesic dome jargon)
#        gcv: Minor variation of gc
ADJ = """Projection constant. May be a float or a string from the list
below. If a string is given, it will optimize k based on the specified
measurement of the polyhedron. Ignored unless -p=disk or -p=slerp and
the grid has triangular faces. Default is 1.
String values can be """ + ', '.join(n for n in sgs.MEASURES)
#        energy: Minimizes the Thompson energy of the points.
#        fill: Maximizes the fill ratio of the polyhedron wrt the unit sphere.
#        edges: Minimizes the difference in edge length.
#        aspect: Minimizes the aspect ratio of euclidean triangles.
#        faces: Minimizes the difference in area between faces.
#        angle: Minimizes the difference in central angle between points
#            sharing edges. (On a unit sphere, same as the spherical distance.)
#        angle_aspect: Minimizes the aspect ratio of spherical triangles.
#        solid_angle: Minimizes the difference in solid angle between faces.
TWEAK = """Makes a tweak to certian methods. For methods that use the
projection constant k, uses approximate parallels instead of exact. For
triangular gc, changes weights in the vertex calculation. May produce a
(slightly) different vertex positioning, and (very slightly) reduce
runtime."""
#COLOR = """Color vertices by base face and faces by group,
#instead of vice versa"""


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
    if string in sgs.MEASURES:
        return string
    else:
        return float(string)

def main():
    parser = argparse.ArgumentParser(description=DESCRIPTION, epilog=EPILOG)
    parser.add_argument("filename", nargs='?', help=FILENAME)
    parser.add_argument("-a", help=FREQ_A, default=2, type=posint)
    parser.add_argument("-b", help=FREQ_B, default=0,
                        type=nonnegativeint)
    parser.add_argument("-p", '--projection', default='flat', help=PROJ,
                        choices=projection.PROJECTIONS)
    parser.add_argument("-n", "--no_normalize", action="store_true",
                        help="Don't normalize vertices onto the unit sphere")
    parser.add_argument("-k", default=1, help=ADJ, type=kparser)
    parser.add_argument("-t", "--tweak", action="store_true", help=TWEAK)

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
        poly = sgs.subdiv(base, frequency, args.projection, args.tweak)
        if args.projection in ('slerp', 'disk'):
            try:
                k = float(args.k)
            except ValueError:
                k = sgs.optimize_k(poly, base, sgs.MEASURES[args.k],
                                   ~args.tweak, ~args.no_normalize)
            poly.vertices += k*sgs.parallels(poly, base, exact=True)
        if not args.no_normalize:
            poly.vertices = xmath.normalize(poly.vertices)
        facecolor = sgs.face_color_bf(poly)
        vertcolor = poly.group.astype(int)
        fx = list(poly.faces)
        fx.extend(range(len(poly.vertices)))
        colors = list(facecolor)
        colors.extend(vertcolor)
        print(off.write_off(poly.vertices, fx, colors))
        print('#frequency =', frequency)
        if args.filename:
            print('#input file =', args.filename)
        print('#projection =', args.projection)
        if args.projection == 'slerp':
            print('#k =', k)
        if args.projection in ('slerp', 'gc'):
            print('#tweak =', args.tweak)
        print('#normalized =', not args.no_normalize)

if __name__ == "__main__":
    main()
