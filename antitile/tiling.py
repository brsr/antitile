# -*- coding: utf-8 -*-
"""
Generic polyhedron and tiling methods
"""
import numpy as np
#from scipy.linalg import circulant
from scipy import sparse


class Tiling:
    """Generic class for tilings and polyhedrons"""
    def __init__(self, vertices, faces):
        self.vertices = vertices
        self.faces = faces

    #should probably cache these properties
    @property
    def edges(self):
        """Returns a list of edges in the tiling"""
        return edges_from_facelist(self.faces)

    @property
    def face_size(self):
        return np.array([len(x) for x in self.faces])

    @property
    def faces_by_size(self):
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
    def vertex_adjacency(self):
        return NotImplemented

    @property
    def face_adjacency(self):
        return NotImplemented


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

def remove_dupes(faces):
    """Remove duplicate entries from face list"""
    return {tuple(x) for x in faces}
