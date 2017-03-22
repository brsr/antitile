# -*- coding: utf-8 -*-
"""
Subdivide a tiling or polyhedra using a similar grid
"""
import warnings
import numpy as np
from numpy.linalg import norm
from scipy.optimize import minimize_scalar
from . import tiling, breakdown, projection, xmath


def subdiv(base, freq=(2, 0), proj='flat', tweak=False):
    projections = projection.PROJECTIONS[proj]
    faces_dict = base.faces_by_size
    if any(x > 4 for x in faces_dict.keys()):
        raise ValueError("Tiling contains at least one face with more than"
                         " 4 sides. Try triangulating those faces first.")
    elif 3 in faces_dict and 4 in faces_dict:
        msg = ("Subdivision on mixed triangle-quadrilateral polyhedra "
               "is not yet implemented")
        raise NotImplementedError(msg)
    elif 3 in faces_dict:
        n = 3
    elif 4 in faces_dict:
        n = 4
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
    #remove redundant vertices/faces
    #TODO

    #project vertices
    proj_fun = projections[n]
    for i in range(n_bf):
        index = (rbkdn.base_face == i)
        base_face = base_faces[i]
        vbf = base.vertices[base_face]
        rbkdn.vertices[index] = proj_fun(rbkdn[index], vbf, freq, tweak)
    #proj_fun(bkdn, abc, freq, k)
    result = tiling.Tiling(rbkdn.vertices, faces)
    result.group = group
    result.base_face = bf
    return result

#--Stuff having to do with naive slerp and the k-factor--
def parallels(poly, base, exact=True):
    #This is only ever going to be used with equilateral faces,
    #so center = normal.
    #FIXME
    face_centers = xmath.normalize(base.vertices[base.faces].sum(axis=1))
    normals = face_centers[poly.base_face]
    parallel = parallel_exact if exact else parallel_approx
    parallel_xyz = parallel(poly.vertices, normals)
    #xyz = poly.vertices
    return parallel_xyz

def parallel_sphere(xyz, pls, k=1):
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
    parallel_xyz = parallel_sphere(poly, base, exact)
    result = minimize_scalar(objective, bracket=[0, 1],
                             args=(poly, parallel_xyz, measure, normalize))
    if ~result.success:
        warnings.warn('optimization routine did not converge')
    return result.x

def objective(k, poly, parallel_xyz, measure, normalize=True):
    test_v = parallel_sphere(poly.vertices, parallel_xyz, k)
    test_v = xmath.normalize(test_v)
    if normalize:
        test_v = xmath.normalize(test_v)
    elif ~normalize and normalize is not None:
        #restore length of original vector
        test_v = xmath.normalize(test_v)
        test_v *= norm(poly.vertices, axis=-1, keepdims=True)
    #if it's None, don't normalize at all
    return measure(test_v, poly)

#In most of these, only triangle faces are considered (since naive slerp is
#the only one that optimizes k, and naive slerp only operates on triangles)
def energy(xyz, poly, exponent=1, over=None):
    """Energy of the vertex arrangement, as defined in
    Thomson's problem."""
    dists = xmath.distance(xyz[np.newaxis], xyz[:, np.newaxis])
    vertex_energy = 1/dists**exponent
    index = np.diag_indices_from(dists)
    vertex_energy[index] = 0
    return vertex_energy.sum(axis=over)/2

def fill_ratio(xyz, poly):
    """Fill ratio of the grid"""
    faces = np.array(poly.faces_by_size[3])
    face_volume = np.abs(np.linalg.det(xyz[faces]))
    return 1 - face_volume.sum()/(8*np.pi)

def edge_length(xyz, poly, spherical=False):
    """Length of each edge in the grid"""
    if spherical:
        dist = xmath.spherical_distance
    else:
        dist = xmath.distance
    edges = poly.edges
    x = xyz[edges]
    result = dist(x[:, 0], x[:, 1])
    return result.max()/result.min()

def face_area(xyz, poly, spherical=False):
    """Area of each face"""
    if spherical:
        area = xmath.spherical_triangle_area
    else:
        area = xmath.triangle_area
    faces = np.array(poly.faces_by_size[3])
    x = xyz[faces]
    result = area(x[:, 0], x[:, 1], x[:, 2])
    return result.max()/result.min()

def aspect_ratio(xyz, poly, spherical=False):
    """Ratio of longest edge to shortest edge for each triangle"""
    if spherical:
        dist = xmath.spherical_distance
    else:
        dist = xmath.distance
    faces = np.array(poly.faces_by_size[3])
    x = xyz[faces]
    y = dist(x, np.roll(x, 1, axis=1))
    result = y.max(axis=-1)/y.min(axis=-1)
    return result.max()/result.min()

MEASURES = {'energy': energy,
            'fill': fill_ratio,
            'edges': edge_length,
            'aspect': aspect_ratio,
            'faces': face_area,
            'angle': lambda xyz, poly: edge_length(xyz, poly, spherical=True),
            'angle_aspect': lambda xyz, poly: aspect_ratio(xyz, poly,
                                                           spherical=True),
            'solid_angle': lambda xyz, poly: face_area(xyz, poly,
                                                       spherical=True)}
