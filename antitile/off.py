# -*- coding: utf-8 -*-
"""
Methods to handle loading from and writing to OFF files
"""
import numpy as np

def readline_comment(file, symbol='#'):
    """Reads line from a file object, but ignores everything after the
    comment symbol (by default '#')"""
    line = file.readline()
    result = line.partition(symbol)[0]
    if len(result) == 0:
        return readline_comment(file)
    else:
        return result


def load_off(file):
    """Loads OFF files from Antiprism.
    Arguments:
        file: file handle
    Returns:
        vertices: A numpy array of shape (# of vertices , 3) representing
            the vertices of the polyhedron.
        faces: A list of list of vertex indices. All faces from the OFF file
            where Nv >= 3 are included here.
        facecolors: A list of colors. If colors are not specified, this will
            be None.
        edges: A numpy array of shape (# of edges, 2) with values of
            vertex indices. Any entry in the face list in the OFF
            file where Nv = 2 will be here.
        edgecolors: A list of colors. If colors are not specified, entry will
            be None.
        verts: A numpy array  of vertex indices. Any entry in the face
            list in the OFF file where Nv = 1 will be here.
        vertexcolors: A list of colors. If colors are not specified, will
            be None.
            """
    head = readline_comment(file)
    if head[:3] != 'OFF':
        raise ValueError('not an OFF file')
    count_line = readline_comment(file)
    nvertices, nfaces, _ = [int(x) for x in count_line.split()]
    vertices = []
    for _ in range(nvertices):
        line = readline_comment(file)
        split_line = line.split()
        if len(split_line) != 3:
            raise ValueError('bad vertex: need 3d xyz coordinates')
        vertices.append([float(x) for x in split_line])
    vertices = np.array(vertices)
    faces = []
    facecolors = []
    edges = []
    edgecolors = []
    verts = []
    vertexcolors = []
    for _ in range(nfaces):
        line = readline_comment(file)
        split_line = line.split()
        nv = int(split_line[0])
        vx = [int(x) for x in split_line[1:nv+1]]
        colorspec = split_line[nv+1:]
        lc = len(colorspec)
        if lc == 0:
            colorspec = None
        elif lc == 1:
            colorspec = int(colorspec[0])
        elif lc >= 3:
            colorspec = [float(x) for x in colorspec]
            if max(colorspec) > 1:
                colorspec = [x/255 for x in colorspec]
        if nv == 1:
            verts.append(vx[0])
            vertexcolors.append(colorspec)
        elif nv == 2:
            edges.append(vx)
            edgecolors.append(colorspec)
        else:
            faces.append(vx)
            facecolors.append(colorspec)

    edges = np.array(edges, dtype=int)
    verts = np.array(verts, dtype=int)
    facecolors = np.array(facecolors)
    edgecolors = np.array(edgecolors)
    vertexcolors = np.array(vertexcolors)
    
    return (vertices, faces, facecolors, edges, edgecolors,
            verts, vertexcolors)


def write_off(vertices, faces, facecolors=None):
    """Export a grid to an OFF file for use with Antiprism.
    Inputs:
        vertices: List of vertices
        faces: List of faces"""
    result = 'OFF\n'
    string = ' '.join([str(len(vertices)), str(len(faces)), '0'])
    result += string
    result += '\n'
    for row in vertices:
        string = ' '.join([str(i) for i in row])
        result += string
        result += '\n'
        
    if facecolors is None:
        for face in faces:
            try:
                x = len(face)
                face = list(face)
            except TypeError:
                x = 1
                face = [face]
            row = [x] + face 
            result += ' '.join(str(i) for i in row) + '\n'
    else:
        for face, color in zip(faces, facecolors):
            try:
                x = len(face)
                face = list(face)
            except TypeError:
                x = 1
                face = [face]
            try:
                color = list(color)
            except TypeError:
                color = [color]
            row = [x] + face + color
            result += ' '.join(str(i) for i in row) + '\n'
    return result
