#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spherical tiling of the henagonal hosohedron.
"""
import argparse
import numpy as np
from antitile import breakdown, projection, xmath, off

_SQRT32 = np.sqrt(3)/2
_TRIANGLE = np.array([[0, 1],
                      [_SQRT32, -0.5],
                      [-_SQRT32, -0.5]])
_SQUARE = np.array([[1, 0],
                    [0, 1],
                    [-1, 0],
                    [0, -1]])


def balloon(base_pts, bkdn, bkdn_to_disk, disk_to_sphere=projection.lambert):
    """Perform the balloon tiling of the sphere

    Args:
        base_pts: Base points (probably either _TRIANGLE or _SQUARE)
        bkdn: Breakdown
        bkdn_to_disk: Projection from breakdown to disk
            (probably 'nslerp' or 'disk')
        disk_to_sphere: Projection from disk to sphere
            (default: ``projection.lambert``)
    """
    sqc = bkdn_to_disk(bkdn, base_pts, None, None)
    #find vertices outside or on the unit circle and merge them
    goodverts = bkdn.group < 90
    #assign one vertex on the unit circle
    base_ind = np.nonzero(~goodverts)[0].min()
    goodverts[base_ind] = True
    sqc[base_ind] = [1, 0]
    index = xmath.renumber(goodverts, base_ind)
    faces = index[bkdn.faces]
    vertices = sqc[goodverts]
    phi, theta = disk_to_sphere(vertices)
    sph_3d = projection.spherical_to_xyz(phi, theta)
    #get rid of degenerate faces
    #(unless they're all degenerate, then keep one for display)
    fc = np.array([len(set(x)) for x in faces])
    if fc.max() <= 2:
        i = np.argwhere(fc == 2)
        if len(i) > 1:
            faceindex = i[0]
        else:
            faceindex = i
    else:
        faceindex = fc >= 3
    faces = faces[faceindex]
    face_group = bkdn.face_group[faceindex]
    #remove repeated vertices but maintain orientation
    outfaces = []
    for face in faces:
        repindex = face != np.roll(face, 1)
        outfaces.append(face[repindex])
    return (sph_3d, outfaces, bkdn.group[goodverts], face_group)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
                                     "Balloon tiling of the sphere.")
    parser.add_argument("a", type=int,
                        help="first frequency of division")
    parser.add_argument("b", nargs="?", type=int, default=0,
                        help="second frequency of division. 0 or blank for "
                             "Class I, same as a for Class II")
    parser.add_argument("-p", action="store_true",
                        help="use disk projection instead of naive slerp")    
    parser.add_argument("-l", action="store_true",
                        help="use lambert projection instead of equidistant")
    parser.add_argument("-q", action="store_true",
                        help="use quadrilaterals instead of triangles")
    args = parser.parse_args()
    a, b = args.a, args.b
    shape = 4 if args.q else 3
    bkdn = breakdown.Breakdown(a, b, shape)
    bkdn_to_disk = projection.PROJECTIONS['disk' if args.p else 'nslerp'][shape]
    disk_to_sphere = projection.lambert if args.l else projection.equidistant
    base_pts = _SQUARE if args.q else _TRIANGLE                   
    v, f, vertcolor, facecolor = balloon(base_pts, bkdn, bkdn_to_disk,
                                         disk_to_sphere)
    result = off.write_off(v, f, facecolors=facecolor,
                           vertexcolors=vertcolor)
    result += '#frequency = {}\n'.format((a, b))
    result += '#projection = '
    result += 'disk' if args.p else 'nslerp'
    result += ' / '
    result += 'equidistant' if args.p else 'lambert'
    print(result)
