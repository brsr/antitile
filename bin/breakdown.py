#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Demonstrative figures of the different breakdowns and methods.
"""

import argparse
#import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
#import antitile
from antitile import xmath, projection, breakdown

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
           'cube': 1/np.sqrt(3),
           'di': 0}

def base_3d(z, shape=TRIANGLE):
    """
    Gives a base shape centered on [0,0,1].
    """
    r = np.sqrt(np.abs(1 - z**2))
    return np.append(r*shape, z*np.ones((shape.shape[0], 1)), axis=-1)


def plot_sphere_outline(ax, color='grey'):
    circle = plt.Circle((0, 0), 1, color=color, fill=False, zorder=3)
    ax.add_artist(circle)


def plot_base(ax, shape, line_pts_n=50, color='blue'):
    rbp = np.roll(shape, -1, axis=0)
    t = np.linspace(0, 1, line_pts_n)
    border = xmath.slerp(shape[:, np.newaxis], rbp[:, np.newaxis],
                         t[:, np.newaxis])
    reshape_border = border.reshape((-1, 3))
    ax.plot(reshape_border[..., 0], reshape_border[..., 1], color=color)

def breakdown_into_bary(n=4, m=2, line_pts_n=50, degenerate=False):
    base_bary = np.eye(3)
    frame = breakdown.frame_triangle(base_bary, n, m)
    if not degenerate:
        frame = frame.reshape((-1, 2, 3))
        norm = np.linalg.norm(frame[:, 0]-frame[:, 1], axis=-1)
        bad = np.isclose(norm, 0)
        frame = frame[~bad]
    t = np.linspace(0, 1, line_pts_n)[:, np.newaxis]
    frame_0 = frame[..., 0, np.newaxis, :]
    frame_1 = frame[..., 1, np.newaxis, :]
    return xmath.lerp(frame_0, frame_1, t)

def breakdown_into_xy(n=4, m=2, line_pts_n=50, degenerate=False):
    frame = breakdown.frame_square(n, m)
    if not degenerate:
        frame = frame.reshape((-1, 2, 2))
        norm = np.linalg.norm(frame[:, 0]-frame[:, 1], axis=-1)
        bad = np.isclose(norm, 0)
        frame = frame[~bad]
    t = np.linspace(0, 1, line_pts_n)[:, np.newaxis]
    frame_0 = frame[..., 0, np.newaxis, :]
    frame_1 = frame[..., 1, np.newaxis, :]
    return xmath.lerp(frame_0, frame_1, t)

def start_plot(limits=(-1.05, 1.05)):
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 8)
    fig.set_tight_layout(True)
    plt.axis('equal')
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.axis('off')#why do we need both? it is a mystery
    return fig, ax

def plot_parallel(ax, in_pts, name, k=(1,), exact=True, line_pts_n=50):
    pts = in_pts + k[0]*projection.parallel(in_pts, CENTER, exact)
    sc = ax.scatter(pts[..., 0], pts[..., 1], label=name, zorder=10)
    color = list(sc.get_facecolor().flat)
    if len(k) > 1:
        ks = np.linspace(*k, num=line_pts_n)[..., np.newaxis, np.newaxis]
        in_pts = in_pts[np.newaxis]
        linpts = in_pts + ks*projection.parallel(in_pts, CENTER, exact)
        linpts = xmath.normalize(linpts)
        ax.plot(linpts[..., 0], linpts[..., 1], c=color, zorder=10)
    return pts

def plot_m2_triangles(ax, base_pts, freq, bkdn):
    triangles = projection.triangles_method2(bkdn.lindex, base_pts, freq)
    tri = xmath.normalize(triangles)
    pc = PolyCollection(tri[..., :2], facecolors='grey')
    ax.add_collection(pc)

def plot_m2_4(ax, base_pts, freq):
    frame = breakdown.frame_square(*freq)
    return projection.square_naive_slerp(frame, base_pts)


def plot_m1_bkdn(ax, base_pts, freq, shape, bkdn, proj,
                 line_pts_n=50, tweak=False):
    if shape == 3:
        c = breakdown_into_bary(*freq,
                                line_pts_n).reshape((-1, line_pts_n, 3))
    else:
        c = breakdown_into_xy(*freq, line_pts_n).reshape((-1, line_pts_n, 2))
    coords = xmath.recordify(['coord'], [c.reshape((-1, c.shape[-1]))])
    result = proj(coords, base_pts, freq, tweak)
    return result.reshape(c.shape[:-1] + (3,))

def plot_m2_bkdn(ax, base_pts, freq, shape, bkdn, line_pts_n=50, tweak=False):
    if shape == 4:
        preframe = breakdown.frame_square(*freq)
        frame = projection.square_naive_slerp(preframe, base_pts)
    else:
        plot_m2_triangles(ax, base_pts, freq, bkdn)
        frame = breakdown.frame_triangle(base_pts, n=freq[0], m=freq[1],
                                         interp=xmath.slerp)
    t = np.linspace(0, 1, num=line_pts_n)[..., np.newaxis]
    reshape = (-1, line_pts_n, 3)
    lines = xmath.slerp(frame[..., 0, np.newaxis, :],
                        frame[..., 1, np.newaxis, :], t).reshape(reshape)
    #ax.scatter(frame[..., 0], frame[..., 1], color='k')
    return lines


def inverse(pts):
    result = pts / np.linalg.norm(pts[..., :2], axis=-1, keepdims=True)**2
    result[..., 0] *= -1
    return result

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

def zerotoone(string):
    if string in HEIGHTS:
        return HEIGHTS[string]
    x = float(string)
    if x < 0 or x > 1:
        msg = "z must be between 0 and 1 inclusive"
        raise argparse.ArgumentTypeError(msg)
    return x

def parse_k(string):
    x = string.split(',')
    if len(x) > 2:
        raise argparse.ArgumentTypeError('too many values')
    #elif len(x) == 1:
        #return float(x[0])
    else:
        return [float(v) for v in x]

FREQ_A = """First breakdown frequency."""
FREQ_B = """Second breakdown frequency. Default is 0."""
PROJ = """List of projection families. Inferred from height if not given.
    disk is only valid on dihedra (z=0)."""
ADJ = ("Projection constant: either a floating point number, or two floats "
" separated with a comma representing a range of values to show. "
"Ignored unless -p=" + ', '.join(projection.PARALLEL) + ". Default is 1.")
TWEAK = """Makes a tweak to certian methods. For methods that use the
projection constant k, uses approximate parallels instead of exact. For
triangular gc, changes weights in the vertex calculation."""

def main():
    parser = argparse.ArgumentParser(description="Breakdown structures")
    parser.add_argument("a", help=FREQ_A, type=posint)
    parser.add_argument("b", help=FREQ_B, nargs='?', default=0,
                        type=nonnegativeint)
    parser.add_argument("-e", default=1.05, type=float,
                        help="extent of plot (defaults to 1.05)")
    parser.add_argument("-z", default=1/np.sqrt(3), type=zerotoone,
                        help="height of breakdown face, in [0, 1], default "
                             "1/3. or a text string: " + ', '.join(HEIGHTS))
    parser.add_argument("-p", '--projection', nargs="*",
                        help=PROJ, choices=projection.PROJECTIONS)
    parser.add_argument("-f", choices=projection.PROJECTIONS,
                        help="""Which projection to use for face figure. By
                        default uses the last one in the projection list. gc
                        (without the q option) shows the triangles created
                        by intersecting great circles instead.""")
    parser.add_argument("-q", action="store_true",
                        help="use quadrilateral instead of triangle")
    parser.add_argument("-k", default=[1], help=ADJ, type=parse_k)
    parser.add_argument("-t", "--tweak", action="store_true", help=TWEAK)
    parser.add_argument("-L", "--lindex", action="store_true",
                        help="Show linear indexes")

    args = parser.parse_args()
    freq = args.a, args.b
    shape_n = 4 if args.q else 3
    z = args.z
    if args.projection:
        names = args.projection
    elif z == 0:
        names = projection.PARALLEL
    else:
        names = list(projection.PROJECTIONS.keys())
        names.remove('disk')
    shape = SHAPES[shape_n]
    base_pts = base_3d(args.z, shape=shape)
    bkdn = breakdown.Breakdown(*freq, shape_n)
    fig, ax = start_plot([-args.e, args.e])
    basecolor = 'k' if args.b == 0 else 'blue'
    if z == 0:
        plot_sphere_outline(ax, color=basecolor)
    else:
        plot_sphere_outline(ax)
        plot_base(ax, base_pts, color=basecolor)

    faceproj = args.f if args.f else names[-1]
    if faceproj == 'gc':
        lines = plot_m2_bkdn(ax, base_pts[:-1], freq, shape_n,
                             bkdn, tweak=args.tweak, line_pts_n=100)
    else:
        proj = projection.PROJECTIONS[faceproj][shape_n]
        lines = plot_m1_bkdn(ax, base_pts[:-1], freq, shape_n,
                             bkdn, proj, tweak=args.tweak, line_pts_n=100)
        if faceproj in projection.PARALLEL:
            pl = projection.parallel
            lines += args.k[0]*pl(lines, CENTER, not args.tweak)
    lines = xmath.normalize(lines)
    lc = LineCollection(lines[..., :2], edgecolor='k')
    ax.add_collection(lc)
    for name in names:
        proj = projection.PROJECTIONS[name][shape_n]
        pts = proj(bkdn, base_pts[:-1], freq, args.tweak)
        index = bkdn.group < 200
        pts = pts[index]
        if name in projection.PARALLEL:
            pts = plot_parallel(ax, pts, name, k=args.k,
                                exact=not args.tweak)
        else:
            pts = xmath.normalize(pts)
            ax.scatter(pts[..., 0], pts[..., 1], label=name, zorder=10)

    if args.lindex:
        lindex_here = bkdn.lindex[index]
        for pt, lin in zip(pts, lindex_here):
            ax.annotate(','.join(str(x) for x in lin), xy=pt[:2],
                        textcoords='data', zorder=100)
        
    if z == 0 and len(names) == 1:
        inv = inverse(pts)
        invlines = inverse(lines)
        invlines[np.abs(invlines) > 2*args.e] = np.nan
        ax.scatter(inv[..., 0], inv[..., 1], color='grey', zorder=0,
                   label='inverse')
        lc = LineCollection(invlines[..., :2], edgecolor='grey', zorder=0)
        ax.add_collection(lc)
    ax.legend()
    plt.show()

if __name__ == "__main__":
    main()
