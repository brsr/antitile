# -*- coding: utf-8 -*-
"""
Subdivide a tiling or polyhedra using a similar grid
"""
import warnings
import numpy as np
from numpy.linalg import norm
from scipy.optimize import minimize_scalar
from scipy import sparse
from . import tiling, breakdown, projection, xmath

def stitch_3(edge, bf, bkdn, index_0, index_1, freq):
    a, b = freq
    bkdn_0 = bkdn[index_0]
    bkdn_1 = bkdn[index_1]
    #first figure out how to roll the shapes so they meet at the edge
    #TODO replace this with isin whenever that gets into numpy
    notedge = np.nonzero(~np.in1d(bf, edge).reshape(bf.shape))[1]
    roll = 2 - notedge
    offset = np.array([a+b, a, 2*a+b])#from 0,0 to a,b
    lindices = np.roll(bkdn_0.lindex, roll[0], axis=-1)
    flipped = offset - np.roll(bkdn_1.lindex, roll[-1], axis=-1)
    matches = np.nonzero(np.all(lindices[:, np.newaxis] ==
                                flipped[np.newaxis], axis=-1))
    l0 = np.nonzero(index_0)[0][matches[0]]
    l1 = np.nonzero(index_1)[0][matches[1]]
    return l0, l1

                   #0, 1,  2, 3
ROLLMAP = np.array([3, 6, 12, 9])
ROLLMASK = 2**np.arange(4)

def _rot_4(coord, n, freq):
    n = n%4
    a, b = freq
    coord = coord.astype(float)
    cco = xmath.float2d_to_complex(coord).flatten()
    rot_cco = cco * np.exp(1j*np.pi*n/2)
    if n == 0:
        offset = 0
    elif n == 1:
        offset = 1
    elif n == 2:
        offset = 1+1j
    elif n == 3:
        offset = 1j
    shift_cco = rot_cco + offset
    cxy = shift_cco * (a + b*1j) + b
    lindex = xmath.complex_to_float2d(cxy)
    return np.round(lindex).astype(int)

def stitch_4(edge, bf, bkdn, index_0, index_1, freq):
    #FIXME
    a, b = freq
    bkdn_0 = bkdn[index_0]
    bkdn_1 = bkdn[index_1]
    #first figure out how to roll the shapes so they meet at the edge
    notedge = np.in1d(bf, edge).reshape(bf.shape).dot(ROLLMASK)
    roll = -np.nonzero(ROLLMAP[np.newaxis] == notedge[..., np.newaxis])[1]
    offset = np.array([a+2*b, b])#from 0,0 to a,b
    lindices = _rot_4(bkdn_0.coord, roll[0], freq)
    flipped = offset - _rot_4(bkdn_1.coord, roll[-1], freq)
    #print(lindices[(bkdn_0.group == 0) | (bkdn_0.group == 1)])
    #print(flipped[(bkdn_1.group == 0) | (bkdn_1.group == 1)])
    #print('-')
    matches = np.argwhere(np.all(lindices[:, np.newaxis] ==
                                 flipped[np.newaxis], axis=-1)).T
    l0 = np.nonzero(index_0)[0][matches[0]]
    l1 = np.nonzero(index_1)[0][matches[1]]
    return l0, l1

def _find_dupe_verts(base, base_faces, rbkdn, freq, stitcher):
    #find redundant vertices
    base_edges = base.edges
    base_edge_corr, base_face_corr = base.faces_by_edge(base_edges)
    l0 = []
    l1 = []
    for i in range(len(base_edges)):
        edge = base_edges[i]
        index = base_edge_corr == i
        facex = base_face_corr[index]
        fn = len(facex)
        if fn > 2:
            warnings.warn("More than 2 faces meet at a single edge. "
                          "Choosing 2 faces arbitrarily...")
            facex = facex[:2]
        elif fn < 2:#external edge, ignore it
            continue
        index_0 = rbkdn.base_face == facex[0]
        index_1 = rbkdn.base_face == facex[1]
        lx0, lx1 = stitcher(edge, base_faces[facex], rbkdn,
                            index_0, index_1, freq)
        l0.extend(lx0)
        l1.extend(lx1)
    matches = np.stack([l0, l1], axis=-1)
    #TODO replace this with np.unique when I update numpy
    matches = np.array(list({tuple(sorted(t)) for t in matches}))
    #if first one lies outside the base face, swap
    index = rbkdn.group[matches[..., 0]] >= 200
    matches[index] = matches[index, ::-1]
    vno = len(rbkdn)
    conns = sparse.coo_matrix((np.ones(len(matches)),
                               (matches[:, 0], matches[:, 1])),
                              shape=(vno, vno))
    ncp, cp = sparse.csgraph.connected_components(conns)
    verts = np.arange(vno, dtype=int)
    for i in range(ncp):
        component = np.argwhere(cp == i).flatten()
        gp = rbkdn.group[component]
        order = np.argsort(gp)
        component = component[order]
        v = verts[component[0]]
        verts[component] = v
    unique_index = verts == np.arange(len(verts))
    renumbered = xmath.renumber(unique_index)
    return renumbered[verts], unique_index

def subdiv(base, freq=(2, 0), proj='flat', tweak=False):
    projections = projection.PROJECTIONS[proj]
    faces_dict = base.faces_by_size
    if any(x > 4 for x in faces_dict.keys()):
        raise ValueError("Tiling contains at least one face with more than "
                         "4 sides. Try triangulating those faces first.")
    elif 3 in faces_dict and 4 in faces_dict:
        msg = ("Subdivision on mixed triangle-quadrilateral polyhedra "
               "is not yet implemented")
        raise NotImplementedError(msg)
    elif 3 in faces_dict:
        n = 3
        stitcher = stitch_3
    elif 4 in faces_dict:
        n = 4
        stitcher = stitch_4
    else:
        raise ValueError("Polyhedron has no faces")
    #create list of basically cartesion product of bkdn with base_faces
    base_faces = np.array(faces_dict[n])
    n_bf = len(base_faces)
    bkdn = breakdown.Breakdown(*freq, n)
    n_bkdn = len(bkdn.coord)
    names = ['base_face', 'group', 'coord', 'lindex', 'vertices']
    bf = np.arange(n_bf).repeat(n_bkdn)
    group = np.tile(bkdn.group, n_bf)
    coord = np.tile(bkdn.coord, (n_bf, 1))
    lindex = np.tile(bkdn.lindex, (n_bf, 1))
    vertices = np.empty((n_bkdn*n_bf, 3))
    arrays = [bf, group, coord, lindex, vertices]
    rbkdn = xmath.recordify(names, arrays)
    faces = np.concatenate([bkdn.faces + i*n_bkdn for i in range(n_bf)])
    verts, unique_index = _find_dupe_verts(base, base_faces,
                                           rbkdn, freq, stitcher)

    rbkdn = rbkdn[unique_index]
    faces = tiling.remove_dupes(verts[faces])
    bf = bf[unique_index]
    group = group[unique_index]
    #project vertices
    proj_fun = projections[n]
    for i in range(n_bf):
        index = (rbkdn.base_face == i)
        base_face = base_faces[i]
        vbf = base.vertices[base_face]
        rbkdn.vertices[index] = proj_fun(rbkdn[index], vbf, freq, tweak)
    result = tiling.Tiling(rbkdn.vertices, faces)
    result.group = group
    result.base_face = bf
    return result

def face_color_bf(poly):
    faces = poly.faces
    bf = np.array(poly.base_face)
    result = []
    for face in faces:
        face = np.array(face)
        fbf = bf[face]
        counts = np.bincount(fbf)
        maxcount = counts.max()
        index = (counts == maxcount)
        if index.sum() > 1:
            result.append(255)
        else:
            result.append(int(np.argwhere(index)))
    return np.array(result)


def face_color_group(poly, fn=max):
    faces = poly.faces
    group = np.array(poly.group)
    result = []
    for face in faces:
        face = np.array(face)
        fbf = group[face]
        result.append(fn(fbf))
    return np.array(result)


#--Stuff having to do with the k-factor--
def parallels(poly, base, exact=True):
    """Given a subdivided polyhedron based on a base polyhedron, return
    the parallels to the base faces for each vertex in the polyhedron
    that would put the vertices onto the sphere"""
    normals = base.face_normals[poly.base_face]
    parallel = parallel_exact if exact else parallel_approx
    return parallel(poly.vertices, normals)

def parallel_sphere(xyz, pls, k=1):
    """Given vertices and parallels, return points on sphere"""
    return xyz + k*pls

def parallel_exact(pts, normal):
    """Projects points onto the sphere parallel to the normal vector.
    >>> center = np.array([0,0,1])
    >>> pts = np.array([[0.5,0.5,0],
                        [1,0,0]])
    >>> parallel_exact(pts, center)
    """
    vdotc = np.sum(pts * normal, axis=-1)
    vdotv = norm(pts, axis=-1)**2
    p = -vdotc + np.sqrt(np.fmax(1 + vdotc**2 - vdotv, 0))
    return p[..., np.newaxis] * normal


def parallel_approx(pts, normal):
    """Approximately projects points onto the sphere parallel to the
        center vector.
    >>> center = np.array([0,0,1])
    >>> pts = np.array([[0.5,0.5,0],
                        [1,0,0]])
    >>> parallel_approx(pts, center)
        """
    q = 1 - norm(pts, axis=-1)
    return q[..., np.newaxis] * normal

def optimize_k(poly, base, measure, exact=True, normalize=True):
    parallel_xyz = parallels(poly, base, exact)
    result = minimize_scalar(objective, bracket=[0, 1],
                             args=(poly, parallel_xyz, measure, normalize))
    if not result.success:
        warnings.warn('Optimization routine did not converge')
    return result.x

def objective(k, poly, parallel_xyz, measure, normalize=True):
    test_v = parallel_sphere(poly.vertices, parallel_xyz, k)
    if normalize:
        test_v = xmath.normalize(test_v)
    elif ~normalize and normalize is not None:
        #restore length of original vector
        test_v = xmath.normalize(test_v)
        test_v *= norm(poly.vertices, axis=-1, keepdims=True)
    #if it's None, don't normalize at all
    return measure(test_v, poly)

#measures
def energy(xyz, poly):
    """Energy of the vertex arrangement, as defined in
    Thomson's problem."""
    return tiling.energy(xyz)

def edge_length(xyz, poly, spherical=False):
    result = tiling.edge_length(xyz, poly.edges, spherical)
    return result.max()/result.min()

def face_area(xyz, poly, spherical=False):
    result = tiling.face_area(xyz, poly, spherical)
    return result.max()/result.min()

def aspect_ratio(xyz, poly, spherical=False):
    result = tiling.aspect_ratio(xyz, poly, spherical)
    return result.max()/result.min()

MEASURES = {'energy': energy,
            'edges': edge_length,
            'aspect': aspect_ratio,
            'faces': face_area,
            'angle': lambda xyz, poly: edge_length(xyz, poly, spherical=True),
            'angle_aspect': lambda xyz, poly: aspect_ratio(xyz, poly,
                                                           spherical=True),
            'solid_angle': lambda xyz, poly: face_area(xyz, poly,
                                                       spherical=True)}
