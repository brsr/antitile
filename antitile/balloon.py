# -*- coding: utf-8 -*-
"""
Spherical tiling of the henagonal hosohedron.
Vaguely resembles a peeled coconut.
"""

import breakdown
import projection
import xmath
import off
import numpy as np

def balloon(a, b, shape, proj = projection.lambert, margin=1E-12):
    bkdn = breakdown.Breakdown(a, b, shape)
    if shape == 'q':
        sqc = projection.square_to_circle(bkdn.coord)
    elif shape == 't':
        abc = np.array([[1,    0],
                        [-0.5,  np.sqrt(3)/2],
                        [-0.5, -np.sqrt(3)/2]])
        sqc = projection.tri_naive_slerp(bkdn.coord, abc)
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

def strip_ev(faces):
    """Strip 1- and 2-vertex "faces" from the facelist"""
    count = np.array([len(np.unique(i)) for i in faces])
    index = count >= 3
    if index.sum() == 0:
        index = count >= 2
    return faces[index]

def clean_triangles(faces):
    """Restate degenerate polygons as triangles or whatever"""
    facelist = faces.tolist()
    facelist = [list(set(x)) if len(set(x)) <= 3 else x for x in facelist]
    return facelist

if __name__ == "__main__":
    a, b = 10, 0
    shape = 'q'
    v, f = balloon(a, b, shape)
    f = clean_triangles(strip_ev(f))
    f = {tuple(x) for x in f}
    result = off.write_off(v, f)
    with open('test.off', 'w') as file:
        file.write(result)
