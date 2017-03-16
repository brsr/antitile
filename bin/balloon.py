#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spherical tiling of the henagonal hosohedron.
Vaguely resembles a peeled coconut.
"""
import argparse
import numpy as np
from antitile import breakdown, projection, xmath, off, tiling

_SQRT32 = np.sqrt(3)/2
_TRIANGLE = np.array([[0, 1],
                      [_SQRT32, -0.5],
                      [-_SQRT32, -0.5]])


def balloon(a, b, shape, proj=projection.lambert):
    """Perform the balloon tiling of the sphere"""
    bkdn = breakdown.Breakdown(a, b, shape)
    if shape == 'q':
        sqc = projection.square_to_circle(bkdn.coord)
    elif shape == 't':
        sqc = projection.tri_naive_slerp(bkdn.coord, _TRIANGLE)
    #find vertices outside or on the unit circle and merge them
    goodverts = bkdn.group < 0#np.linalg.norm(sqc, axis=-1) < 1-margin
    #assign one vertex on the unit circle
    base_ind = np.nonzero(~goodverts)[0].min()
    goodverts[base_ind] = True
    sqc[base_ind] = [1, 0]
    index = xmath.renumber(goodverts, base_ind)
    faces = index[bkdn.faces]
    vertices = sqc[goodverts]
    phi, theta = proj(vertices)
    sph_3d = projection.spherical_to_xyz(phi, theta)
    return sph_3d, faces

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
                                     "Balloon tiling of the sphere.")
    parser.add_argument("a", type=int,
                        help="first frequency of division")
    parser.add_argument("b", nargs="?", type=int, default=0,
                        help="second frequency of division. 0 or blank for "
                             "Class I, same as a for Class II")
    parser.add_argument("-p", action="store_true",
                        help="use lambert projection instead of equidistant")
    parser.add_argument("-q", action="store_true",
                        help="use quadrilaterals instead of triangles")
    args = parser.parse_args()
    a, b = args.a, args.b
    shape = 'q' if args.q else 't'
    proj_fun = projection.lambert if args.p else projection.equidistant
    v, f = balloon(a, b, shape, proj=proj_fun)
    f = tiling.remove_dupes(tiling.clean_triangles(tiling.strip_ev(f)))
    result = off.write_off(v, f)
    print(result)
