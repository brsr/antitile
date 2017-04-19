# -*- coding: utf-8 -*-
"""
Generic polyhedron and tiling methods
"""
import warnings
import numpy as np
from . import xmath#, off
#from scipy.linalg import circulant
from scipy import sparse


class Tiling:
    """Generic class for tilings and polyhedrons
    Attributes:
        vertices: List of vertices
        faces: List of faces
        """
    def __init__(self, vertices, faces):
        self.vertices = vertices
        self.faces = faces

#    @property
#    def off(self):
#        """Returns a string representing the OFF file for the faces
#        and vertices of this tiling. Explicit edges and vertices,
#        and colorings, are not included."""
#        return off.write_off(self.vertices, self.faces)

    #should probably cache these properties
    @property
    def edges(self):
        """Returns a list of edges in the tiling"""
        return edges_from_facelist(self.faces)

    @property
    def face_size(self):
        """Size of each face in tiling"""
        return np.array([len(x) for x in self.faces])

    @property
    def faces_by_size(self):
        """Faces in tiling, grouped by size and converted into
        numpy arrays."""
        result = dict()
        for face in self.faces:
            x = len(face)
            if x in result:
                result[x].append(face)
            else:
                result[x] = [face]
        for i in result:
            result[i] = np.array(result[i])
        return result

    @property
    def face_normals(self):
        """Normals for each face. This is ill-defined for skew faces."""
        fs = self.face_size
        fd = self.faces_by_size
        result = np.empty((len(fs), 3))
        for i in fd:
            index = fs == i
            faces = fd[i]
            a = self.vertices[faces]
            b = self.vertices[np.roll(faces, -1, axis=-1)]
            normal = np.cross(a,b).sum(axis=-2)
            result[index] = normal
        return xmath.normalize(result)


    def faces_by_edge(self, edges):
        """Given edges, lists the faces adjacent to the edges
        Arguments:
            edges: List of edges. An array of shape (n, 2). Lowest-numbered
            vertex comes first.
        Returns: two arrays:
            edges: Index of each edge
            faces: Index of each face"""
        facesize = self.face_size
        facedict = self.faces_by_size
        ex = []
        fx = []
        for i in facedict:
            face_no = np.nonzero(facesize == i)[0]
            faces = facedict[i]
            for j in range(i):
                faceedge = faces[..., [j-1, j]]
                condition = faceedge[..., 0] < faceedge[..., 1]
                a = np.where(condition[..., np.newaxis],
                             faceedge, faceedge[..., ::-1])
                edgel = np.all(edges[:, np.newaxis] == a[np.newaxis], axis=-1)
                edge_face = np.nonzero(edgel)
                ex.append(edge_face[0])
                fx.append(face_no[edge_face[1]])
        return np.concatenate(ex), np.concatenate(fx)

    @property
    def vertex_adjacency(self):
        """Vertex adjacency matrix"""
        edges = self.edges
        redges = np.concatenate([edges, edges[:, ::-1]], axis=0)
        n_r = len(redges)
        n = len(self.vertices)
        return sparse.coo_matrix((np.ones(n_r), (redges[:, 0], redges[:, 1])),
                                 shape=(n, n))

    @property
    def face_adjacency(self):
        """Face adjacency matrix"""
        edges = self.edges
        ex, fx = self.faces_by_edge(edges)
        faceadj = set()
        for i in range(len(edges)):
            index = ex == i
            these_faces = fx[index]
            for j in range(len(these_faces) - 1):
                faceadj.add((these_faces[j - 1], these_faces[j]))
        faceadj = faceadj.union( (x[1], x[0]) for x in faceadj)
        faceadj = np.array(list(faceadj))
        #rface = np.concatenate([faceadj, faceadj[:, ::-1]], axis=0)
        n_r = len(faceadj)
        n = len(self.faces)
        return sparse.coo_matrix((np.ones(n_r, dtype=np.int8),
                                  (faceadj[:, 0], faceadj[:, 1])),
                                 shape=(n, n))

def edges_from_facelist(faces):
    """Given a list of faces, returns a list of edges."""
    edges = set()
    for face in faces:
        for i in range(len(face)):
            edge = set((face[i-1], face[i]))
            if len(edge) == 2:
                edge = tuple(sorted(edge))
                edges.add(edge)
    result = np.array(sorted([edge for edge in edges]))
    return result

#def strip_ev(faces):
#    """Strip 1- and 2-vertex "faces" from the facelist"""
#    count = np.array([len(np.unique(i)) for i in faces])
#    index = count >= 3
#    if index.sum() == 0:
#        index = count >= 2
#    return faces[index]
#
#def clean_triangles(faces):
#    """Restate degenerate polygons as triangles or whatever"""
#    facelist = faces.tolist()
#    facelist = [list(set(x)) if len(set(x)) <= 3 else x for x in facelist]
#    return facelist
#
#def remove_dupes(faces):
#    """Remove duplicate entries from face list"""
#    result = set()
#    for face in faces:
#        aface = np.array(face)
#        r = aface.argmin()
#        result.add(tuple(np.roll(aface, -r)))
#    return list(result)


#Measures
def energy(xyz, exponent=1, over=None):
    """Energy of the vertex arrangement, as defined in
    Thomson's problem."""
    dists = xmath.distance(xyz[np.newaxis], xyz[:, np.newaxis])
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',
                                r'divide by zero encountered in true_divide')
        vertex_energy = 1/dists**exponent
    index = np.diag_indices_from(dists)
    vertex_energy[index] = 0
    return vertex_energy.sum(axis=over)/2

def center_of_gravity(xyz):
    return np.linalg.norm(xyz.mean(axis=0))

def bentness(xyz, poly):
    """Bentness of faces: how far away from planar they are. Faces with 3
    vertices (or less) inherently have bentness 0."""
    fs = poly.face_size
    fd = poly.faces_by_size
    result = np.zeros(len(fs))
    for i in fd:
        if i <= 3:
            continue
        index = fs == i
        faces = fd[i]
        a = xyz[faces]
        centered_a = a - np.mean(a, axis=1, keepdims=True)
        sv = np.linalg.svd(centered_a, compute_uv=False)
        result[index] = sv[..., 2]/(sv[..., :2].sum(axis=-1))
    return result

def edge_length(xyz, edges, spherical=False):
    """Length of each edge in the grid"""
    if spherical:
        dist = xmath.central_angle
    else:
        dist = xmath.distance
    x = xyz[edges]
    result = dist(x[:, 0], x[:, 1])
    return result

def face_area(xyz, poly, spherical=False):
    """Area of each face"""
    #FIXME this is inconsistent for skew faces
    if spherical:
        area = xmath.triangle_solid_angle
    else:
        area = xmath.triangle_area
    result = np.zeros(len(poly.faces))
    fs = poly.face_size
    for i, faces in poly.faces_by_size.items():
        if i <= 2:
            continue
        index = fs == i
        x = xyz[faces]
        a = np.zeros(len(x))
        for j in range(i-2):
            a += area(x[:, 0], x[:, j+1], x[:, j+2])
        result[index] = a
    return result

def aspect_ratio(xyz, poly, spherical=False):
    """Ratio of longest edge to shortest edge for each face"""
    if spherical:
        dist = xmath.central_angle
    else:
        dist = xmath.distance
    result = np.ones(len(poly.faces))
    fs = poly.face_size
    for i, faces in poly.faces_by_size.items():
        if i <= 2:
            continue
        index = fs == i
        x = xyz[faces]
        y = dist(x, np.roll(x, 1, axis=1))
        result[index] = y.max(axis=-1)/y.min(axis=-1)
    return result
