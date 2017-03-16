# -*- coding: utf-8 -*-
"""
Generic polyhedron and tiling methods
"""
import numpy as np

class Tiling:
    def __init__(self, vertices, faces):
        self.vertices = vertices
        self.faces = faces

    @property
    def edges(self):
        """Returns a list of edges in the tiling"""
        return edges_from_facelist(self.faces)
        
    @property
    def face_orientation(self):
        """Returns the orientation of each face in the tiling with
        respect to (0,0,0)."""
        #only really need this for triangles and squares
        return NotImplemented

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