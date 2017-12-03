# -*- coding: utf-8 -*-
"""
Generic polyhedron and tiling methods
"""
import warnings
import numpy as np
from scipy import sparse
from scipy.sparse import csgraph
from . import xmath

class Tiling:
    """Generic class for tilings and polyhedrons
    Args/Attributes:
        vertices: List of vertices
        faces: List of faces
        """
    def __init__(self, vertices, faces):
        self.vertices = vertices
        self.faces = faces

    #none of these properties are cached, that's OK for now
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
            normal = np.cross(a, b).sum(axis=-2)
            result[index] = normal
        return xmath.normalize(result)

    @property
    def true_faces(self):
        """List of true faces (those with 3 or more vertices)."""
        return [x for x in self.faces if len(x) > 2]

    def faces_by_edge(self, edges):
        """Given edges, lists the faces adjacent to the edges

        Arguments:
            edges: List of edges. An array of shape (n, 2).
                Lowest-numbered vertex comes first.

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
        """Vertex edge adjacency matrix"""
        edges = self.edges
        redges = np.concatenate([edges, edges[:, ::-1]], axis=0)
        n_r = len(redges)
        n = len(self.vertices)
        return sparse.coo_matrix((np.ones(n_r), (redges[:, 0], redges[:, 1])),
                                 shape=(n, n), dtype=np.int8)

    @property
    def face_adjacency(self):
        """Face edge adjacency matrix"""
        edges = self.edges
        ex, fx = self.faces_by_edge(edges)
        faceadj = set()
        for i in range(len(edges)):
            index = ex == i
            these_faces = fx[index]
            for j in range(len(these_faces) - 1):
                faceadj.add((these_faces[j - 1], these_faces[j]))
        faceadj = faceadj.union((x[1], x[0]) for x in faceadj)
        faceadj = np.array(list(faceadj))
        #rface = np.concatenate([faceadj, faceadj[:, ::-1]], axis=0)
        n_r = len(faceadj)
        n = len(self.faces)
        return sparse.coo_matrix((np.ones(n_r, dtype=np.int8),
                                  (faceadj[:, 0], faceadj[:, 1])),
                                 shape=(n, n))

    @property
    def face_vertex_incidence(self):
        """Face vertex incidence matrix"""
        faces = self.faces
        n_f = len(self.faces)
        n_v = len(self.vertices)
        fx = []
        vx = []
        for i in range(n_f):
            face = faces[i]
            for v in face:
                fx.append(i)
                vx.append(v)
        ones = np.ones(len(fx), dtype=np.int8)
        return sparse.coo_matrix((ones, (fx, vx)), shape=(n_f, n_v))

    @property
    def vertex_f_adjacency(self):
        """Vertex face adjacency matrix"""
        fvi = self.face_vertex_incidence
        result = fvi.T @ fvi
        result.setdiag(0)
        result.eliminate_zeros()
        return result

    @property
    def face_v_adjacency(self):
        """Face vertex adjacency matrix"""
        fvi = self.face_vertex_incidence
        result = fvi @ fvi.T
        result.setdiag(0)
        result.eliminate_zeros()
        return result

    def normalize_faces(self):
        """Normalizes faces by rotating the vertex list so that
        the lowest vertex is first."""
        faces = self.faces
        for i in range(len(faces)):
            aface = np.array(faces[i])
            r = aface.argmin()
            faces[i] = np.roll(aface, -r)

    def id_dupe_faces(self, face_group=None, oriented=False):
        """Identify duplicate faces. Should use `normalize_faces` first.

        Args:
            face_group: Optional. Of a collection with duplicate faces,
                the faces with lowest value will be preferred to those
                with higher values. By default, just uses the face that
                is earlier in the face list.
            oriented: If True, consider opposite orientations of
                the same face to be different. Default False.

        Returns:
            An array of face indexes. Non-duplicated faces will have an entry
            equal to its position in the array. If a face is duplicated, its
            antry will equal the preferred face's position
        """
        if face_group is None:
            face_group = np.zeros(len(self.faces))
        fs = self.faces_by_size
        fn = self.face_size
        n = len(fn)
        comp = np.zeros((n, n), dtype=bool)
        for i in fs:
            faces = fs[i]
            index = fn == i
            #TODO replace this with np.unique when 1.13 comes out
            reps = np.all(faces[:, np.newaxis] == faces[np.newaxis], axis=-1)
            if not oriented:
                revfaces = np.roll(faces[..., ::-1], 1, axis=-1)
                revreps = np.all(faces[:, np.newaxis] == revfaces[np.newaxis],
                                 axis=-1)
                reps |= revreps
            reps[np.diag_indices_from(reps)] = False
            comp[np.ix_(index, index)] = reps
        ncp, cp = csgraph.connected_components(comp)
        result = np.ones(n, dtype=bool)
        for i in range(ncp):
            index = cp == i
            fg = face_group[index]
            am = np.argmin(fg)
            result_i = np.zeros(fg.shape, dtype=bool)
            result_i[am] = True
            result[index] = result_i
        return result

def edges_from_facelist(faces):
    """Given a list of faces, returns a list of edges."""
    edges = set()
    for face in faces:
        for i in range(len(face)):
            edge = (face[i-1], face[i])
            if edge[0] != edge[1]:
                edge = tuple(sorted(edge))
                edges.add(edge)
    result = np.array(sorted([edge for edge in edges]))
    return result

def orient_face(face, edge):
    """Determines how to roll the list of vertices for a face so that its
    0th and 1st vertices correspond to edge.

    Args:
        face: A list of vertex numbers.
        edge: A list of vertex numbers.

    Returns:
        roll: Argument to ``np.roll`` to roll ``face`` into the right position
        flip: If the 0th and 1st indexes of ``face``, once rolled, are in
            opposite order of ``edge``, this is true."""
    indexes = np.argwhere(np.in1d(face, edge))
    if len(indexes) != 2:
        raise ValueError()
    elif ((indexes[1] + 1) % len(face)) == indexes[0]:
        roll = 1
    elif indexes[1] == indexes[0] + 1:
        roll = -indexes[0]
    else:
        raise ValueError()

    rolled_face = np.roll(face, roll)
    if np.all(rolled_face[:2] == edge):
        flip = False
    elif np.all(rolled_face[:2][::-1] == edge):
        flip = True
    else:
        raise ValueError()
    return int(roll), flip

#Measures
def energy(xyz, exponent=1, over=None):
    """Energy of the vertex arrangement, as defined in
    Thomson's problem.

    Args:
        xyz: Coordinates of vertices.
        exponent: By default, the energy is the sum of 1/r, where r is
        the distance from one point to another. If exponent is specified,
        instead 1/r^exponent is used.
        over: Which axis to sum over. By default all axes are summed over.
            To get the contribution of each individual vertex, specify
            over = 1.

    Returns:
        Energy of the vertex arrangement (or energy of individual vertexes
        if over = 1)
    """
    dists = xmath.distance(xyz[np.newaxis], xyz[:, np.newaxis])
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',
                                r'divide by zero encountered in true_divide')
        vertex_energy = 1/dists**exponent
    index = np.diag_indices_from(dists)
    vertex_energy[index] = 0
    return vertex_energy.sum(axis=over)/2

def center_of_gravity(xyz):
    """Distance between the center of gravity of the vertices xyz and
    the zero point."""
    return np.linalg.norm(xyz.mean(axis=0))

def bentness(xyz, poly):
    """Bentness of faces: how far away from planar they are. Faces with 3
    vertices (or less) inherently have bentness 0.

    Args:
        xyz: Coordinates of vertices.
        poly: Tiling object, used to get the polyhedral structure of the
            vertices. (xyz is separate so we don't have to recalculate the
            structure every time we change the vertex position.)

    Returns:
        Array containing the bentness of each face."""
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
    """Length of each edge in the grid
    Args:
        xyz: Coordinates of vertices.
        edge: Array of shape (N, 2), specifying which vertices make up an edge.
        spherical: If true, uses the spherical distance (central angle)
            along each edge instead of the Euclidean distance. Default False.

    Returns:
        Array of length N containing the length of each edge.
    """
    if spherical:
        dist = xmath.central_angle
    else:
        dist = xmath.distance
    x = xyz[edges]
    result = dist(x[:, 0], x[:, 1])
    return result

def face_area(xyz, poly, spherical=False):
    """Area of each face. If one of the faces is bent and spherical is not
    True, that face may have an inconsistent result.

    Args:
        xyz: Coordinates of vertices.
        poly: Tiling object, used to get the polyhedral structure of the
            vertices. (xyz is separate so we don't have to recalculate the
            structure every time we change the vertex position.)
        spherical: If true, uses the spherical area (solid angle)
            of each face instead of the Euclidean area. Default False.

    Returns:
        Array containing the area of each face.
    """
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
    """Ratio of longest edge to shortest edge for each face.

    Args:
        xyz: Coordinates of vertices.
        poly: Tiling object, used to get the polyhedral structure of the
            vertices. (xyz is separate so we don't have to recalculate the
            structure every time we change the vertex position.)
        spherical: If true, uses the spherical distance (central angle)
            along each edge instead of the Euclidean distance. Default False.

    Returns:
        Array containing the aspect ratio of each face.
    """
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
