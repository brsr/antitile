# -*- coding: utf-8 -*-
"""

"""
import numpy as np
from . import flat, xmath

class Breakdown(flat.FlatTiling):
    """
    Breakdown structures
    Attributes:
        freq: A 2-tuple (n, m) describing the breakdown structure.
        shape: Either 3 (triangular) or 4 (quadrilateral)
        vertices: Vertices of the breakdown
        coord: Barycentric coordinates if shape=3,
            or xy coordinates if shape=4
        lindex: Linear index coordinates of each vertex
        group: An integer describing where each vertex falls.
            -19 through -10: Inside the breakdown, lying on a feature
            -1: Inside the breakdown, other
            0, 1, 2, (3): Base vertices
            10, 11, 12, (13): Base edges
            100, 101, 102, (103): Vertices that lie outside the triangle
                but are adjacent to vertices within the triangle
            127: Outside the breakdown triangle
        faces: Array of faces in the breakdown structure

    """
    def __init__(self, a, b, shape=3, remove_outside=True):
        """
        Constructor for the breakdown object.

        Args:
        n, m: Frequency of the breakdown
        remove_outside: Whether to remove features that lie outside the
            breakdown.
            True: Remove faces and vertices that are neither inside the
                breakdown nor adjacent to it,
            False: Keep everything
            By default, true.
        """
        v = a + b + 1
        super().__init__(v,v,shape)
        self.freq = (a, b)

        vertices = self.vertices
        group = -np.ones(len(vertices), dtype=np.int8)
        self.group = group
        #adjust the vertices so we have room for the shape
        vertices[..., 0] -= b
        vertices[..., 2] += b
        if shape==3:
            self._t()
        elif shape == 4:
            self._q()

        cn = 111 if remove_outside else 128
        condition = group < cn
        # renumber the adjacency list too
        index = xmath.renumber(condition)
        self.vertices = self.vertices[condition]
        self.coord = self.coord[condition]
        self.lindex = self.lindex[condition]
        self.group = self.group[condition]
        faces = index[self.faces]
        badface = np.any(faces < 0, axis=-1)
        self.faces = faces[~badface]
        self.base_pts = np.nonzero(np.in1d(self.group, [0,1,2,3]))[0]


    def _t(self):
        a, b = self.freq
        vertices = self.vertices
        group = self.group
        anorm = a**2 + a * b + b**2
        # vertices of the triangle are (0,0), (a,b),(-b,a+b)
        x = vertices[:, 0]
        y = vertices[:, 1]

        side1 = a*y - b*x
        side2 = (a + b)*y + a*x
        side3 = (a + b)*x + b*y
        group[side1 == 0] = 10
        group[side2 == anorm] = 11
        group[side3 == 0] = 12
        group[(side1 < 0) | (side2 > anorm) | (side3 < 0)] = 127
        group[(x == 0) & (y == 0)] = 0
        group[(x == a) & (y == b)] = 1
        group[(x == -b) & (y == a + b)] = 2
        self._shared_group()
        group[(group == 126) & (side1 < 0)] = 100
        group[(group == 126) & (side2 > anorm)] = 101
        group[(group == 126) & (side3 < 0)] = 102

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
        a, b = self.freq
        vertices = self.vertices
        group = self.group
        #vertices of the square are (0,0), (a,b), (a-b, a+b), (-b,a)
        x = vertices[:, 0]
        y = vertices[:, 1]
        right = a*y - b*x
        left = a*x + b*y
        anorm = (a**2+b**2)
        group[right == 0] = 10
        group[left == anorm] = 11
        group[right == anorm] = 12
        group[left == 0] = 13
        group[(right < 0) | (left > anorm) |
              (right > anorm) | (left < 0)] = 127
        group[(x == 0) & (y == 0)] = 0
        group[(x == a) & (y == b)] = 1
        group[(x == a - b) & (y == a + b)] = 2
        group[(x == -b) & (y == a)] = 3
        self._shared_group()
        group[(group == 126) & (right < 0)] = 100
        group[(group == 126) & (left < anorm)] = 101
        group[(group == 126) & (right > anorm)] = 102
        group[(group == 126) & (left > 0)] = 103
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
        group = self.group
        faces = self.faces
        inside = np.nonzero(group < 0)[0]
        fv_inside = np.in1d(faces, inside).reshape(faces.shape)
        face_inside = np.any(fv_inside, axis=-1)
        shared = np.unique(faces[face_inside])
        group[shared] = np.fmin(group[shared], 126)

def frame(n=4, m=2, shape=3):
    if shape == 3:
        return frame_triangle(n=n, m=m)
    elif shape == 4:
        return frame_square(n=n, m=m)

def frame_triangle(base_pts = np.eye(3), n=4, m=2, interp=xmath.lerp):
    """Creates the "frame" of edge points for method 2.
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
    """Creates the "frame" of edge points for method 2.
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
    to = np.concatenate( (top[:-1], right), axis=0)
    pairs = np.stack([frm, to], axis=1)
    other = np.stack([1-pairs[..., 1], pairs[..., 0]], axis=-1)
    return np.stack([pairs, other])

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection, LineCollection

    a, b = 2, 0
    shape = 3

    bkdn = Breakdown(a, b, shape)
    frm = frame(a, b, shape)
    if shape == 3:
        abc = np.array([[0,  0],
                        [1,  0],
                        [0.5, np.sqrt(3)/2]])
        pts_2d = bkdn.coord @ abc
    else:
        abc = np.array([[0,0],
                         [1,0],
                         [1,1],
                         [0,1]])
        pts_2d = bkdn.coord
    mx = pts_2d.max(axis=0)
    mn = pts_2d.min(axis=0)
    ptx = pts_2d[bkdn.faces]
    fig, ax = plt.subplots()
    fig.set_size_inches(6, 6)
    plt.axis('equal')
    pc = PolyCollection(ptx, edgecolors='grey')
    ax.add_collection(pc)
    ax.scatter(pts_2d[..., 0], pts_2d[..., 1], c=-bkdn.group)
    x = abc[...,0].tolist() + [abc[0,0]]
    y = abc[...,1].tolist() + [abc[0,1]]

    ax.plot(x,y, c='k')
    if shape == 't':
        frmp = frm @ abc
        pass
    else:
        frmp = frm
    lc = LineCollection(frmp.reshape((-1, frmp.shape[-2], frmp.shape[-1])), color='c')
    ax.add_collection(lc)
    lindex = bkdn.lindex
    li = xmath.line_intersection(frmp[0,lindex[:,0],0],frmp[0,lindex[:,0],1],
                                 frmp[1,lindex[:,1],0],frmp[1,lindex[:,1],1])
    ax.scatter(li[..., 0], li[..., 1], c='y')
    #x = np.stack([li[:,0],pts_2d[:,0]],axis=-1)
    #y = np.stack([li[:,1],pts_2d[:,1]],axis=-1)
    #ax.plot(x, y, c='g')

    anorm = a**2 + a*b+ b**2
    vertices = bkdn.vertices
    mat = np.array([[a+b, b],
                    [-b, a]])/anorm
    coords1 = vertices[:, :2].dot(mat.T)
    l1 = 1-coords1.sum(axis=-1, keepdims=True)
    coord1 = np.concatenate([l1, coords1], axis=1)
    mat2 = np.array([[   -a, -b - a, anorm],
                    [a + b,      b, 0],
                    [   -b,      a, 0]]) / anorm

    coords = vertices.copy()
    coords[:, 2] = 1  # the ol' affine matrix trick
    coord2 = coords.dot(mat2.T)