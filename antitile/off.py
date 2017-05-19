# -*- coding: utf-8 -*-
"""
Methods to handle loading from and writing to OFF files
"""
import numpy as np

def readline_comment(file, symbol='#'):
    """Reads line from a file object, but ignores everything after the
    comment symbol (by default '#')"""
    line = file.readline()
    if not line:
        return ''
    result = line.partition(symbol)[0]
    return result if result else readline_comment(file)

def load_off(file):
    """Loads OFF files from Antiprism.

    Args:
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
            file where Nv = 2 will be here. May be None.
        edgecolors: A list of colors. If colors are not specified, entry will
            be None.
        verts: A numpy array  of vertex indices. Any entry in the face
            list in the OFF file where Nv = 1 will be here. May be None.
        vertexcolors: A list of colors. If colors are not specified, will
            be None.
            """
    head = readline_comment(file)
    if not head:
        raise ValueError('empty input')
    elif head[:3] != 'OFF':
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

    edges = np.array(edges, dtype=int) if edges else None
    verts = np.array(verts, dtype=int) if verts else None
    facecolors = facecolors if facecolors else None
    edgecolors = edgecolors if edgecolors else None
    vertexcolors = vertexcolors if vertexcolors else None

    return (vertices, faces, facecolors, edges, edgecolors,
            verts, vertexcolors)


def write_off(vertices, faces, facecolors=None, edges=None, edgecolors=None,
              verts=None, vertexcolors=None):
    """Export a grid to an OFF file for use with Antiprism.

    Args:
        vertices: List of vertex cordinates
        faces: List of faces
        facecolors: List of colors corresponding to each face (optional)
        edges: List of edges (2-vertex faces) (optional)
        edgecolors: List of colors corresponding to each edge (optional)
        verts: List of vertexes (1-vertex faces) (optional)
        vertexcolors: List of colors corresponding to each vertex (optional)

    Returns:
        A string containing the contents of the OFF file
    """
    #align the "faces" into a single list
    if (facecolors is None or len(facecolors) == 0) and faces is not None:
        facecolors = ['']*len(faces)
    if (edgecolors is None or len(edgecolors) == 0) and edges is not None:
        edgecolors = ['']*len(edges)
    if (vertexcolors is None or len(vertexcolors) == 0) and verts is not None:
        print('woo')
        vertexcolors = ['']*len(verts)
    elif vertexcolors is not None and verts is None:
        verts = list(range(len(vertices)))

    if len(faces) != len(facecolors):
        msg = "Different number of face colors than faces: {}, {}"
        raise ValueError(msg.format(len(faces), len(facecolors)))
    elif edgecolors is not None and len(edges) != len(edgecolors):
        msg = "Different number of edge colors than edges: {}, {}"
        raise ValueError(msg.format(len(edges), len(edgecolors)))
    elif vertexcolors is not None and len(verts) != len(vertexcolors):
        msg = "Different number of vertex colors than vertices: {}, {}"
        raise ValueError(msg.format(len(verts), len(vertexcolors)))

    flist = list(faces)
    if edges is not None:
        flist.extend(list(edges))
    if verts is not None:
        flist.extend([v] for v in verts)
    clist = list(facecolors)
    if edgecolors is not None:
        clist.extend(list(edgecolors))
    if vertexcolors is not None:
        clist.extend([v] for v in vertexcolors)

    lengths = [len(x) for x in flist]

    result = 'OFF\n'
    string = ' '.join([str(len(vertices)), str(len(flist)), '0'])
    result += string
    result += '\n'
    #coordinates
    for row in vertices:
        string = ' '.join([str(i) for i in row])
        result += string
        result += '\n'

    #"faces"
    for n, f, c in zip(lengths, flist, clist):
        result += str(n) + ' '
        try:
            result += ' '.join(str(i) for i in f) + ' '
        except TypeError:
            result += str(f) + ' '
        try:
            result += ' '.join(str(i) for i in c) + '\n'
        except TypeError:
            result += str(c) + '\n'
    return result
