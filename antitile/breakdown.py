# -*- coding: utf-8 -*-
"""
Breakdown structures
"""
import numpy as np
from . import flat, xmath

class Breakdown(flat.FlatTiling):
    """Breakdown structures.

    Args:
        a, b: Frequency of the breakdown
        shape: 3 for triangle, 4 for quadrilateral.
        remove_outside: Whether to remove features that lie outside the
            breakdown.
            True: Remove faces and vertices that are neither inside the
                breakdown nor adjacent to it,
            False: Keep everything
            By default, True.

    Attributes:
        freq: A 2-tuple (a, b) describing the breakdown structure.
        shape: Either 3 (triangular) or 4 (quadrilateral)
        vertices: Vertices of the breakdown
        coord: Barycentric coordinates if shape==3,
            or xy coordinates if shape==4
        lindex: Linear index coordinates of each vertex
        group: An integer describing where each vertex falls.
            0: Inside the breakdown
            1 - 89: Inside the breakdown, lying on a feature
            90 - 93: Base vertices
            100 - 103: Base edges
            200 - 203: Vertices that lie outside but are adjacent to
                vertices within
            255: Outside the breakdown
        faces: Array of faces in the breakdown structure
        face_group: Max of the group of each vertex within a face.
            The face_group of faces lying partially out of the breakdown
            is adjusted so faces with more points in the breakdown
            will tend to be assigned to that base face in the SGS
            routine.
            """

    def __init__(self, a, b, shape=3, remove_outside=True):
        v = a + b + 1
        super().__init__(v, v, shape)
        self.freq = (a, b)
        self.shape = shape

        vertices = self.vertices
        group = np.zeros(len(vertices), dtype=np.uint8)
        self.group = group
        #adjust the vertices so we have room for the shape
        vertices[..., 0] -= b
        vertices[..., 2] += b
        if shape == 3:
            self._t()
        elif shape == 4:
            self._q()

        cn = 222 if remove_outside else 256
        condition = group < cn
        # renumber the adjacency list too
        index = xmath.renumber(condition)
        self.vertices = self.vertices[condition]
        self.coord = self.coord[condition]
        self.lindex = self.lindex[condition]
        self.group = self.group[condition]
        faces = index[self.faces]
        badface = np.any(faces < 0, axis=-1)#remove missing faces
        self.faces = faces[~badface]
        faces = self.faces
        group = self.group
        #self.base_pts = np.nonzero(np.in1d(self.group, [90, 91, 92, 93]))[0]
#        #assign faces to groups
#        face_group = np.zeros(len(faces), dtype=np.uint8)
#        index_in = self.group < 90
        index_bv = (self.group >= 90) & (self.group < 100)
#        index_be = (self.group >= 100) & (self.group < 200)
        index_out = (self.group >= 200)
        face_out = index_out[faces]
        face_out_sum = face_out.sum(axis=-1).astype(np.uint8)
#        face_out_max = face_out.max(axis=-1)
#        face_group = face_out_max + 10*face_out_sum
        fx = group[faces]
        face_group = fx.max(axis=-1)
        face_group += face_out_sum*10
        face_bv = np.any(index_bv[faces], axis=-1)
        face_group[face_bv] = 90 + face_out_sum[face_bv]
        self.face_group = face_group


    def _t(self):
        """Initialization for triangle breakdowns"""
        a, b = self.freq
        vertices = self.vertices
        group = self.group
        anorm = a**2 + a * b + b**2
        # vertices of the triangle are (0,0), (a,b),(-b,a+b)
        x = vertices[:, 0]
        y = vertices[:, 1]
        #featural lines
        q = a - b
        line1 = x == 0
        line1_seg = y <= max(a, b)
        line2 = y == b
        line2_seg = x >= min(0, q)
        line3 = y == a - x
        line3_seg = y >= min(a, b)
        group[line1 & line1_seg] = 1
        group[line2 & line2_seg] = 2
        group[line3 & line3_seg] = 3
        if a == b:
            group[line1 & ~line1_seg] = 11
            group[line2 & ~line2_seg] = 12
            group[line3 & ~line3_seg] = 13
        side1 = a*y - b*x
        side2 = (a + b)*y + a*x
        side3 = (a + b)*x + b*y
        group[side1 == 0] = 100
        group[side2 == anorm] = 101
        group[side3 == 0] = 102
        group[(side1 < 0) | (side2 > anorm) | (side3 < 0)] = 255
        group[(x == 0) & (y == 0)] = 90
        group[(x == a) & (y == b)] = 91
        group[(x == -b) & (y == a + b)] = 92
        self._shared_group()
        group[(group == 250) & (side1 < 0)] = 200
        group[(group == 250) & (side2 > anorm)] = 201
        group[(group == 250) & (side3 < 0)] = 202

        # barycentric coordinates
        mat = np.array([[   -a, -b - a, anorm],
                        [a + b,      b, 0],
                        [   -b,      a, 0]]) / anorm

        coords = vertices.copy()
        coords[:, 2] = 1  # the ol' affine matrix trick
        self.coord = coords.dot(mat.T)

        # linear index coordinates
        u = x + y
        v = a - x
        w = a + b - y
        # u + v + w = 2a+b
        self.lindex = np.array([u, v, w]).T

    def _q(self):
        """Initialization for quad breakdowns"""
        a, b = self.freq
        vertices = self.vertices
        group = self.group
        #vertices of the square are (0,0), (a,b), (a-b, a+b), (-b,a)
        x = vertices[:, 0]
        y = vertices[:, 1]
        #featural lines
        q = a - b
        line1 = x == 0
        line1_seg = y <= max(a, b)
        line2 = y == b
        line2_seg = x >= min(0, q)
        line3 = x == a-b
        line3_seg = y >= min(a, b)
        line4 = y == a
        line4_seg = x <= max(0, q)
        group[line1 & line1_seg] = 1
        group[line2 & line2_seg] = 2
        group[line3 & line3_seg] = 3
        group[line4 & line4_seg] = 4
        if a == b:
            group[line1 & ~line1_seg] = 11
            group[line2 & ~line2_seg] = 12
            group[line3 & ~line3_seg] = 13
            group[line4 & ~line4_seg] = 14
        right = a*y - b*x
        left = a*x + b*y
        anorm = (a**2 + b**2)
        group[right == 0] = 100
        group[left == anorm] = 101
        group[right == anorm] = 102
        group[left == 0] = 103
        group[(right < 0) | (left > anorm) |
              (right > anorm) | (left < 0)] = 255
        group[(x == 0) & (y == 0)] = 90
        group[(x == a) & (y == b)] = 91
        group[(x == a - b) & (y == a + b)] = 92
        group[(x == -b) & (y == a)] = 93
        self._shared_group()
        group[(group == 250) & (right < 0)] = 200
        group[(group == 250) & (left < anorm)] = 201
        group[(group == 250) & (right > anorm)] = 202
        group[(group == 250) & (left > 0)] = 203
        # xy [0,1]^2 coordinates
        vx = a + b*1j
        cx = self.proj_complex
        xy = cx/vx
        self.coord = xmath.complex_to_float2d(xy)

        #linear index coordinates
        u = x + b
        v = y
        self.lindex = np.array([u, v]).T

    def _shared_group(self):
        """Routine for the parts of determining vertex groups
        that are shared between triangle and quad faces."""
        group = self.group
        faces = self.faces
        inside = np.nonzero(group < 90)[0]
        fv_inside = np.in1d(faces, inside).reshape(faces.shape)
        face_inside = np.any(fv_inside, axis=-1)
        shared = np.unique(faces[face_inside])
        group[shared] = np.fmin(group[shared], 250)

    def lindex_reorient(self, n, flip=False):
        """Reorient linear indexes with n and flip, as returned by
        ``tiling.orient_face``. This really only works for Class I
        and II breakdowns."""
        if self.shape == 3:
            return _reorient_3(self, n, flip)
        else:
            return _reorient_4(self, n, flip)

def _reorient_3(bkdn, n, flip=False):
    """Helper function for triangle faces"""
    if flip:
        n = n - 1
    n = n%3
    result = np.roll(bkdn.lindex, n, axis=-1)
    if flip:
        a, b = bkdn.freq
        normal = np.array([-b, a + b, -a])
        rm = xmath.reflect_through_origin(normal).T
        result = result.dot(rm)
    return np.round(result).astype(int)

_ROTMAT = np.array([[0, -1],
                    [1, 0]])

def _reorient_4(bkdn, n, flip=False):
    """Helper function for quad faces"""
    if flip:
        n = n - 2
    n = n%4
    rm = np.linalg.matrix_power(_ROTMAT, n).T
    result = bkdn.lindex.dot(rm)
    result -= result.min(axis=0)
    if flip:
        a, b = bkdn.freq
        normal = np.array([b, -a])
        offset = np.array([b, a])/2
        rm = xmath.reflect_through_origin(normal).T
        ro = result - offset
        result = ro.dot(rm)
        result += offset
    return np.round(result).astype(int)

def frame(n=4, m=2, shape=3):
    """Creates the "frame" of edge points for method 2.

    Args:
        n, m: breakdown frequency
        shape: triangle (3) or square (4) faces

    Returns a multidimensional array with dimensions:
        (x:         Which rotation of the lines
        n + m + 1:  Linear index
        2:          Which end of the line (0 is the "bottom")
        x:          barycentric or xy coordinates)
        x is 3 if shape is 3, else 2.
    """
    if shape == 3:
        return frame_triangle(n=n, m=m)
    elif shape == 4:
        return frame_square(n=n, m=m)

def frame_triangle(n=4, m=2, base_pts=np.eye(3), interp=xmath.lerp):
    """Creates the "frame" of edge points for method 2 on a triangle face.

    Args:
        n, m: breakdown frequency

    Returns a multidimensional array with dimensions:
        (3:         Which rotation of the lines
        n + m + 1:  Linear index
        2:          Which end of the line (0 is the "bottom")
        3:          barycentric coordinates)
    """
    tnm = np.linspace(0, 1, n + m + 1)[np.newaxis, :, np.newaxis]
    tn = np.linspace(0, 1, n, endpoint=False)[np.newaxis, :, np.newaxis]
    tm = np.linspace(0, 1, m + 1)[np.newaxis, :, np.newaxis]

    base_pts_0 = base_pts[:, np.newaxis]
    base_pts_1 = np.roll(base_pts, -1, axis=0)[:, np.newaxis]
    base_pts_2 = np.roll(base_pts, -2, axis=0)[:, np.newaxis]

    linterpol_nm = interp(base_pts_0, base_pts_1, tnm)
    linterpol_n = interp(base_pts_0, base_pts_2, tn)
    linterpol_m = interp(base_pts_2, base_pts_1, tm)

    to = np.concatenate((linterpol_n, linterpol_m), axis=1)
    return np.stack([linterpol_nm, to], axis=2)

def frame_square(n=4, m=2):
    """Creates the "frame" of edge points for method 2 on a quad face.

    Args:
        n, m: breakdown frequency

    Returns a multidimensional array with dimensions:
        (2:         Which rotation of the lines
        n + m + 1:  Linear index
        2:          Which end of the line (0 is the "bottom")
        2:          xy coordinates)
    """
    tn = np.linspace(0, 1, n + 1)
    tm = np.linspace(1, 0, m + 1)
    bottom = np.stack([tn, np.zeros(tn.shape)], axis=-1)
    top = np.stack([tn, np.ones(tn.shape)], axis=-1)
    right = np.stack([np.ones(tm.shape), tm], axis=-1)
    left = np.stack([np.zeros(tm.shape), tm], axis=-1)
    frm = np.concatenate((left[:-1], bottom), axis=0)
    to = np.concatenate((top[:-1], right), axis=0)
    pairs = np.stack([frm, to], axis=1)
    other = np.stack([1-pairs[..., 1], pairs[..., 0]], axis=-1)
    return np.stack([pairs, other])
