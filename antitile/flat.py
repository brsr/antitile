# -*- coding: utf-8 -*-
"""
Flat square {4,4} and triangular tilings {3,6} for use with breakdown.
unitile2d does these, but there's a note in its docs saying it
may change in the future, so shouldn't rely on it.
"""
import numpy as np
from . import tiling

def project_skew(coords):
    """
    Flattens cubic coordinates or 2d skew coordinates into
    2d orthogonal coordinates on the complex plane.

    >>> x = np.arange(5)
    >>> coords = np.vstack((x,np.zeros(5),-x)).T
    >>> project_skew(coords) # doctest: +NORMALIZE_WHITESPACE
    array([[ 0.,  0.],
           [ 1.,  0.],
           [ 2.,  0.],
           [ 3.,  0.],
           [ 4.,  0.]])
    """
    x = coords[..., 0] + coords[..., 1] / 2
    y = coords[..., 1] * np.sqrt(3) / 2
    return np.stack([x, y], axis=-1)

class FlatTiling(tiling.Tiling):
    """Info about the geometry of a flat square {4,4} and triangular
    tiling {3,6}
    Attributes:
        freq: A 2-tuple (n, m) describing the breakdown structure.
        shape: Either 3 (triangular) or 4 (quadrilateral)
        vertices: Vertices of the flat tiling
        faces: Faces of the flat tiling
    """

    def __init__(self, a=25, b=25, shape=3):
        """Constructor for the planar triangle/hex geometry.

        Args:
            n, m: Number of vertices along each skew-axis.
            shape: 3 for triangle, 4 for quadrilateral
        """

        if a < 1 or b < 1:
            raise ValueError("Frequency out of range")
        self.freq = (a, b)
        self.shape = shape
        n = a * b
        # use the cube coordinate method
        mesh = np.meshgrid(np.arange(0, a), np.arange(0, b), 0)
        vertices = np.array(mesh).reshape(3, n).T
        self.vertices = vertices

        if shape == 3:
            face_config = self._init_t()
        elif shape == 4:
            face_config = self._init_q()
        else:
            raise ValueError('invalid shape ' + str(shape))
        #find adjacency
        #some of these will be off of the grid
        adjacency = vertices[:, np.newaxis, np.newaxis]  + face_config
        #use clip mode so it doesn't error on us...
        adj_indexes = np.ravel_multi_index([adjacency[..., 0],
                                            adjacency[..., 1]],
                                           (a, b), mode='clip')
        #but then wipe out the bad entries with -1
        bad = ~np.all((adjacency[..., :2] >= 0) &
                      (adjacency[..., 0, np.newaxis] < a) &
                      (adjacency[..., 1, np.newaxis] < b), axis=-1)
        adj_indexes[bad] = -1
        #turn adjacency into a list of faces
        faces = adj_indexes.reshape((-1, face_config.shape[-2]))
        goodface = np.all(faces >= 0, axis=-1)#eliminate off-the-grid faces
        super().__init__(vertices, faces[goodface])

    def _init_t(self):
        vertices = self.vertices
        vertices[..., 2] = -vertices[..., 0] - vertices[..., 1]
        self.default_proj = project_skew
        return np.array([[[0, 0,  0],
                          [1, 0, -1],
                          [0, 1, -1]],
                         [[0,  0,  0],
                          [0,  1, -1],
                          [-1, 1,  0]]], dtype=int)


    def _init_q(self):
        self.default_proj = lambda x: x[..., :2]
        return np.array([[0, 0, 0],
                         [1, 0, 0],
                         [1, 1, 0],
                         [0, 1, 0]], dtype=int)

    @property
    def proj_2d(self):
        return self.default_proj(self.vertices)

    @property
    def proj_complex(self):
        twod = self.proj_2d
        return twod[..., 0] + 1j*twod[..., 1]
#
#if __name__ == "__main__":
#    import matplotlib.pyplot as plt
#    from matplotlib.collections import PolyCollection
#
#    a, b = 2, 2
#    shape = 3
#    tiling = FlatTiling(a, b, shape)
#    pts_2d = tiling.proj_2d
#    mx = pts_2d.max(axis=0)
#    mn = pts_2d.min(axis=0)
#    ptx = pts_2d[tiling.faces]
#    fig, ax = plt.subplots()
#    fig.set_size_inches(10, 6)
#    plt.axis('equal')
#    pc = PolyCollection(ptx, edgecolors='grey')
#    ax.add_collection(pc)
#    ax.scatter(pts_2d[..., 0], pts_2d[..., 1], c='k')
