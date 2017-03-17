# -*- coding: utf-8 -*-
"""
Generic polyhedron and tiling methods
"""
import numpy as np
from scipy.linalg import circulant

class Tiling:
    """Generic class for tilings and polyhedrons"""
    def __init__(self, vertices, faces):
        self.vertices = vertices
        self.faces = faces

    @property
    def edges(self):
        """Returns a list of edges in the tiling"""
        return edges_from_facelist(self.faces)

    @property
    def face_size(self):
        return np.array([len(x) for x in self.faces])

    @property
    def face_orientation(self):
        """Returns the orientation of the points in each face in the
        tiling with respect to (0,0,0). 1 is counterclockwise, -1 is
        clockwise, 0 means it's mixed (probably some sort of weird
        self-intersecting thing) or lies on a plane through the origin."""
        lens = self.face_size
        min_len = lens.min()
        max_len = lens.max()
        orientation = np.zeros(lens.shape, dtype=np.int8)
        for i in range(min_len, max_len+1):
            index = lens == i
            i_faces = np.array([self.faces[i] for i in np.nonzero(index)[0]])
            i_pts = self.vertices[i_faces]
            if i == 3:
                dets = np.linalg.det(i_pts)
                cclock = dets > 0
                clock = dets < 0
            else:
                det_index = circulant(np.arange(4))[:, :3]
                dets_ind = np.linalg.det(i_pts[..., det_index, :])
                cclock = np.all(dets_ind > 0, axis=-1)
                clock =  np.all(dets_ind < 0, axis=-1)
            ix = np.nonzero(index)[0]
            orientation[ix[cclock]] = 1
            orientation[ix[clock]] = -1
        return orientation

    def orient_faces(self):
        """Orient faces of tiling in counterclockwise order with
        respect (0,0,0)"""
        for face, o in zip(self.faces, self.face_orientation):
            if o == -1:
                face.reverse()

    @property
    def vertex_adjacency(self):
        return NotImplemented

    @property
    def face_adjacency(self):
        return NotImplemented

    @property
    def vertex_face_incidence(self):
        return NotImplemented

    def faces_by_edge(self, edges):
        return NotImplemented

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