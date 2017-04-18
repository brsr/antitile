# -*- coding: utf-8 -*-
"""
Demonstrative figures of the different breakdowns and methods.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
import antitile
from antitile import xmath

TRIANGLE = np.array([[0, 1],
                     [-np.sqrt(3)/2, -0.5],
                     [np.sqrt(3)/2, -0.5],
                     [0, 1]])


def normalized_flat(*args):
    result = antitile.projection.tri_bary(*args)
    return xmath.normalize(result)


def naive_slerp(*args):
    #just use k = 1 here
    result = antitile.projection.tri_naive_slerp(*args)
    center = np.array([0, 0, 1])
    parallel = antitile.sgs.parallel_exact(result, center)
    return xmath.normalize(result + parallel)


def base_triangle_3d(angle):
    """
    Gives a base triangle centered on [0,0,1]
    icosahedron: arctan(2)
    octahedron: np.pi/2
    tetrahedron: arccos(-1/3)
    3-dihedron: 2*np.pi/3
    """
    z = np.sqrt(np.abs((2*np.cos(angle) + 1)/3))
    r = np.sqrt(np.abs(1 - z**2))
    return np.append(r*TRIANGLE, z*np.ones((4, 1)), axis=-1)


def breakdown_into_bary(n=4, m=2, line_pts_n=50, degenerate=False):
    base_bary = np.eye(3)
    frame = antitile.geodesic_grid.frame_method2(base_bary, n, m,
                                                interp=xmath.lerp)
    if not degenerate:
        frame = frame.reshape((-1,2,3))
        norm = np.linalg.norm(frame[:,0]-frame[:,1], axis=-1)
        bad = np.isclose(norm, 0)
        frame = frame[~bad]
    t = np.linspace(0, 1, line_pts_n)[:, np.newaxis]
    frame_0 = frame[..., 0, np.newaxis, :]
    frame_1 = frame[..., 1, np.newaxis, :]
    return xmath.lerp(frame_0, frame_1, t)



def square_explicit(frame, line_pts_n=61, degenerate=False):
    if ~degenerate:
        frame = frame.reshape((-1,2,2))
        norm = np.linalg.norm(frame[:,0]-frame[:,1], axis=-1)
        bad = np.isclose(norm, 0)
        frame = frame[~bad]
    t = np.linspace(0, 1, line_pts_n)[:, np.newaxis]
    frame_0 = frame[..., 0, np.newaxis, :]
    frame_1 = frame[..., 1, np.newaxis, :]
    return xmath.lerp(frame_0, frame_1, t)


def plot_square_circle(n=4, m=2, inverse=True, steplimit=90):
    if inverse:
        limits = [-2, 2]
    else:
        limits = [-1.05, 1.05]
    frame = antitile.breakdown.frame_square(n, m)
    euc_lines = square_explicit(frame)
    lines = antitile.projection.square_to_circle(euc_lines)
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 8)
    fig.set_tight_layout(True)
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.axis('off')#why do we need both? it is a mystery
    if inverse:
        norm = np.linalg.norm(lines, axis=-1, keepdims=True)
        inverse = lines / norm**2
        inverse[..., 0] = -inverse[..., 0]
        dx = np.linalg.norm(inverse[:, 1:] - inverse[:, :-1], axis=-1)
        if steplimit:
            index = dx > steplimit
            pad = np.zeros((index.shape[0], 1), dtype=bool)
            index = np.concatenate((index, pad), axis=-1)
            print(index.shape, inverse.shape)
            inverse = np.where(index[..., np.newaxis], np.nan, inverse)
        lc = LineCollection(inverse, color='grey', zorder=0)
        ax.add_collection(lc)
    lc = LineCollection(lines, color='k', zorder=1)
    ax.add_collection(lc)
    if m > 0:
        color = 'blue'
    else:
        color = 'black'
    circle = plt.Circle((0, 0), 1, color=color, fill=False, zorder=3)
    ax.add_artist(circle)
    #plot_m1_pts(ax, base_pts, method, n, m)
    return fig, ax

def show_breakdown(n=4, m=2, triangle=TRIANGLE):
    bkdn = antitile.Breakdown(n, m, remove_outside=True)
    pts = antitile.geodesic_grid.flat(bkdn.vertices.bary, triangle[:3])
    frame = antitile.geodesic_grid.frame_method2(triangle[:3],
                                                n, m, interp=xmath.lerp)
    reshape_frame = frame.reshape((-1, 2, 2))
    fig, ax = plt.subplots()
    plt.axis('off')

    fig.set_size_inches(8, 8)
    lc = LineCollection(reshape_frame, color='k', zorder=1)
    ax.add_collection(lc)
    ax.scatter(pts[..., 0], pts[..., 1], color='k', zorder=4)

    if m == 0:
        color = 'k'
    else:
        color = 'b'
    ax.plot(triangle[:, 0], triangle[:, 1], c=color, zorder=3)


def plot_sphere_outline(ax, color='grey'):
    circle = plt.Circle((0, 0), 1, color=color, fill=False, zorder=3)
    ax.add_artist(circle)


def plot_base_triangle(ax, triangle, line_pts_n=50, color='blue'):
    rbp = np.roll(triangle, -1, axis=0)
    t = np.linspace(0, 1, line_pts_n)
    border = xmath.slerp(triangle[:, np.newaxis], rbp[:, np.newaxis],
                         t[:, np.newaxis])
    reshape_border = border.reshape((-1, 3))
    ax.plot(reshape_border[..., 0], reshape_border[..., 1], color=color)


def plot_m1_pts(ax, base_pts, method, n=4, m=2, color='k'):
    bkdn = antitile.breakdown.Breakdown(n, m, remove_outside=True)
    pts = method(bkdn.coord, base_pts)
    ax.scatter(pts[..., 0], pts[..., 1], color=color, zorder=4)


def plot_m2_pts(ax, base_pts, n=4, m=2, normalize=True, color='k'):
    bkdn = antitile.breakdown.Breakdown(n, m, remove_outside=True)
    pts = antitile.projection.tri_intersections(bkdn.vertices, base_pts,
                                                (n, m), tweak=normalize)
    pts = antitile.xmath.normalize(pts)
    ax.scatter(pts[..., 0], pts[..., 1], color=color, zorder=4)


def m1_bkdn_lines(base_pts, method, n=4, m=2, line_pts_n=50):
    bary = breakdown_into_bary(n, m, line_pts_n).reshape((-1, line_pts_n, 3))
    lines = method(bary, base_pts)
    return lines.reshape((-1, line_pts_n, 3))


def plot_m1_bkdn(ax, base_pts, method, n=4, m=2, line_pts_n=50):
    lines = m1_bkdn_lines(base_pts, method, n, m, line_pts_n)
    lc = LineCollection(lines[..., :2], color='k', zorder=1)
    ax.add_collection(lc)


def plot_m2_bkdn_lines(ax, base_pts, n=4, m=2, line_pts_n=50):
    frame = antitile.breakdown.frame_triangle(base_pts, n, m)
    t = np.linspace(0, 1, line_pts_n)[:, np.newaxis]
    frame_0 = frame[..., 0, np.newaxis, :]
    frame_1 = frame[..., 1, np.newaxis, :]

    lines = xmath.slerp(frame_0, frame_1, t)
    reshape_lines = lines.reshape((-1, line_pts_n, 3))
    lc = LineCollection(reshape_lines[..., :2], color='k', zorder=1)
    ax.add_collection(lc)


def plot_m2_bkdn_triangles(ax, base_pts, n=4, m=2, line_pts_n=50):
    bkdn = antitile.breakdown.Breakdown(n, m, remove_outside=True)
    tm2 = antitile.projection.triangles_method2
    ptx = xmath.normalize(tm2(bkdn.vertices, base_pts, (n, m)))

    pc = PolyCollection(ptx[..., :2], color='lightgrey', zorder=0)
    ax.add_collection(pc)#these leak outside their triangle but close enough


def start_plot(limits=[-1.05, 1.05]):
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 8)
    fig.set_tight_layout(True)
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.axis('off')#why do we need both? it is a mystery

    return fig, ax


def show_order_2(n=4, m=2, inverse=True):
    if inverse:
        limits = [-1.5, 1.5]
    else:
        limits = [-1.05, 1.05]
    fig, ax = start_plot(limits)
    base_pts = base_triangle_3d(np.pi*2/3)[:3]
    method = antitile.geodesic_grid.naive_slerp
    lines = m1_bkdn_lines(base_pts, method, n=n, m=m)[..., :2]
    lc = LineCollection(lines, color='k', zorder=1)
    ax.add_collection(lc)
    if inverse:
        norm = np.linalg.norm(lines, axis=-1, keepdims=True)
        inverse = lines / norm**2
        inverse[..., 0] = -inverse[..., 0]
        lc = LineCollection(inverse, color='grey', zorder=0)
        ax.add_collection(lc)
    if m > 0:
        plot_sphere_outline(ax, color='blue')
    #plot_m1_bkdn(ax, base_pts, method, n, m)
    plot_m1_pts(ax, base_pts, method, n, m)
    return fig, ax


def show_m1_on_m1(base_pts, n=4, m=2, method = naive_slerp):
    fig, ax = start_plot()
    plot_sphere_outline(ax)
    plot_base_triangle(ax, base_pts)

    plot_m1_bkdn(ax, base_pts, method, n, m)
    plot_m1_pts(ax, base_pts, method, n, m)
    return fig, ax


def show_m2_on_m2(base_pts, n=4, m=2):
    fig, ax = start_plot()
    plot_sphere_outline(ax)
    plot_base_triangle(ax, base_pts)
    plot_m2_bkdn_lines(ax, base_pts, n, m)
    plot_m2_bkdn_triangles(ax, base_pts, n, m)
    plot_m2_pts(ax, base_pts, n, m)
    return fig, ax


def show_everything_on_m2(base_pts, n=4, m=2):
    fig, ax = show_m2_on_m2(base_pts, n, m)
    plot_m2_pts(ax, base_pts, n, m, normalize=False, color='b')
    plot_m1_pts(ax, base_pts, normalized_flat, n, m, color='green')
    plot_m1_pts(ax, base_pts, antitile.projection.tri_areal, n, m, color='red')
    plot_m1_pts(ax, base_pts, naive_slerp, n, m, color='orange')
    return fig, ax


if __name__ == "__main__":
    n = 6
    m = 0
    base_pts = base_triangle_3d(np.arccos(-1/3))[:3]
    show_everything_on_m2(base_pts, n=n, m=m)
#    plt.gca().set_position([0, 0, 1, 1])
    #plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    plt.savefig('test.svg',
                bbox_inches='tight', pad_inches = 0, facecolor='none')