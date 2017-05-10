#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualizes OFF files using matplotlib, which allows for export to png, svg,
etc.
"""
import argparse
from sys import stdin
import numpy as np
from antitile.off import load_off
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

def vis_faces(ax, vertices, faces, facecolors=None, edgecolors='k', cmap=None):
    """Plots faces
    Arguments:
        ax: Axis to plot on
        vertices: 3d coordinates of vertices (shape (n, 3))
        faces: List of lists of vertex indices to plot
        facecolors: Color of vertices (optional)
        edgecolors: Color of edges of each face (defaults to black)"""

    these_faces = [vertices[i] for i in faces]
    p3d = Poly3DCollection(these_faces, edgecolors=edgecolors)
    try:
        p3d.set_facecolor(facecolors)
    except ValueError:
        p3d.set_array(np.array(facecolors))
        p3d.set_cmap(plt.get_cmap(cmap))
    ax.add_collection(p3d)


def vis_edges(ax, vertices, edges, edgecolors=None):
    """Plots vertices
    Arguments:
        ax: Axis to plot on
        vertices: 3d coordinates of vertices (shape (n, 3))
        edges: Array of vertex indices to plot (shape (m, 2))
        edgecolors: Color of edges (optional)"""
    these_edges = vertices[edges]
    l3d = Line3DCollection(these_edges, edgecolor=edgecolors)
    ax.add_collection(l3d)


def vis_vertices(ax, vertices, verts, vertexcolors=None):
    """Plots vertices
    Arguments:
        ax: Axis to plot on
        vertices: 3d coordinates of vertices (shape (n, 3))
        verts: List of vertex indices to plot
        vertexcolors: Color of vertices (optional)"""
    these_vertices = vertices[verts]
    ax.scatter(these_vertices[..., 0], these_vertices[..., 1],
               these_vertices[..., 2], c=vertexcolors)


def vis_init(canvassize=4, viewangles=(45, 45), distance=10):
    """Initializes the visualization.
    Arguments:
        viewangles: Viewing elevation and azimuth
        distance: Viewing distance
    Returns:
        fig: Figure object of the visualization
        ax: Axis object"""
    fig = plt.figure()
    fig.set_size_inches(canvassize, canvassize)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect("equal")
    ax.view_init(*viewangles)
    ax.dist = distance
    plt.axis('off')

    return fig, ax

def vis_bounds(ax, vertices):
    """Sets the bounds on ax based on the vertices."""
    max_v = vertices.max()
    min_v = vertices.min()
    ax.set_xlim(min_v, max_v)
    ax.set_ylim(min_v, max_v)
    ax.set_zlim(min_v, max_v)
    ax.set_aspect("equal")

def impute_color(colors, default):
    """Replaces None in color listings"""
    if not colors:
        colors = default
    else:
        for i in range(len(colors)):
            if colors[i] is None:
                colors[i] = default
    return colors

def main_3d(vertices, faces, facecolors, edges, edgecolors, verts,
            vertexcolors, args):
    fig, ax = vis_init(viewangles=(args.elev, args.azim), distance=args.dist)
    vis_bounds(ax, vertices)
    if args.skeleton:
        vis_edges(ax, vertices, edges, edgecolors)
        vis_vertices(ax, vertices, verts, vertexcolors)
    else:
        vis_faces(ax, vertices, faces, facecolors,
                  edgecolors=args.color_edge, cmap=args.cmap)

def main_2d(vertices, faces, facecolors, edges, edgecolors, verts,
            vertexcolors, args, canvassize=4):
    fig = plt.figure()
    fig.set_size_inches(canvassize, canvassize)
    ax = fig.add_subplot(111)
    ax.set_aspect("equal")
    plt.axis('off')
    these_faces = [vertices[i] for i in faces]
    pd = PolyCollection(these_faces, edgecolors=edgecolors)
    try:
        pd.set_facecolor(facecolors)
    except ValueError:
        pd.set_array(np.array(facecolors))
        pd.set_cmap(plt.get_cmap(args.cmap))
    ax.add_collection(pd)
    these_edges = vertices[edges]
    ld = LineCollection(these_edges, edgecolor=edgecolors)
    ax.add_collection(ld)
    these_vertices = vertices[verts]
    ax.scatter(these_vertices[..., 0], these_vertices[..., 1], c=vertexcolors)

def main():
    parser = argparse.ArgumentParser(description="Visualize an OFF file "
                                                 "using matplotlib.")
    parser.add_argument("infile", nargs="?", help="OFF file to visualize."
                        "If not specified, reads from stdin.")
    parser.add_argument("-o", "--outfile", help="File to save image to. "
                        "If not specified, will display on screen.")
    parser.add_argument("--dimension", type=int, choices=[2, 3],
                        help="Force 2d or 3d view (normally inferred "
                             "from the file)")
    parser.add_argument("-a", "--azim", type=float, default=30,
                        help="Viewing azimuth (3d only)")
    parser.add_argument("-e", "--elev", type=float, default=30,
                        help="Viewing elevation (3d only)")
    parser.add_argument("-d", "--dist", type=float, default=9,
                        help="Viewing distance (3d only)")
    parser.add_argument("-s", "--skeleton", action="store_true",
                        help="Show edges instead of faces (3d only) "
                             "(edges must be defined in OFF file)")
    parser.add_argument("--color-face", default='y',
                        help="Color for faces (if not in OFF file)")
    parser.add_argument("--color-edge", default='k',
                        help="Color for edges (if not in OFF file)")
    parser.add_argument("--color-vertex", default='k',
                        help="Color for vertices (if not in OFF file)")
    parser.add_argument("--cmap", help="Matplotlib color map to use for" +
                        "color indexes (default: mpl default)")

    args = parser.parse_args()
    file = open(args.infile) if args.infile else stdin
    with file as f:
        x = load_off(f)
    vertices, faces, facecolors, edges, edgecolors, verts, vertexcolors = x
    facecolors = impute_color(facecolors, args.color_face)
    edgecolors = impute_color(edgecolors, args.color_edge)
    vertexcolors = impute_color(vertexcolors, args.color_vertex)

    vranges = vertices.ptp(axis=0)
    index = np.isclose(vranges, 0)
    if args.dimension == 2 and index.sum() != 1:
        index_x = np.array([False]*3)
        x = np.argwhere(index)
        if len(x) == 0:
            x = 2
        elif len(x) > 1:
            x = x[0]
        index_x[x] = True
        index = index_x
    if args.dimension == 2 or index.sum() == 1:
        main_2d(vertices[..., ~index], faces, facecolors, edges, edgecolors,
                verts, vertexcolors, args)
    else:
        main_3d(vertices, faces, facecolors, edges, edgecolors, verts,
                vertexcolors, args)

    if args.outfile:
        plt.savefig(args.outfile, bbox_inches='tight', pad_inches=0,
                    facecolor='none')
    else:
        plt.show()


if __name__ == "__main__":
    main()
