# -*- coding: utf-8 -*-
"""
Transformations from the breakdowns to some surface in 3-space
(sometimes 2-space)

Attributes:
    PROJECTIONS: A map of maps. The first map is keyed on the name of a family
        of projections, and the second is keyed on shape (3 or 4). For
        instance, PROJECTIONS['nslerp'][3] gives the function to perform
        nslerp on a triangle.
    PARALLEL: A list of which projections are compatible with
        parallel projection.
"""

import numpy as np
from numpy.linalg import norm

from . import xmath, breakdown

_TEST_EQ_TRI = np.eye(3)
_TEST_EQ_Q = xmath.normalize(np.array([[1, 0, 1],
                                       [0, 1, 1],
                                       [-1, 0, 1],
                                       [0, -1, 1]]))
_TEST_SKEW_TRI = np.array([[0.8, 0.6, 0],
                           [0, 1, 0],
                           [0, 0, 1]])
_TEST_SKEW_Q = xmath.normalize(np.array([[1, 0, 1],
                                         [0, 1, 1],
                                         [-1, 0, 1],
                                         [0, -1, 0]]))
_TEST_BARY = np.array([[1, 0.8, 0.6, 0.4, 0.2, 0.0],
                       [0, 0.2, 0.2, 0.2, 0.5, 0.5],
                       [0, 0.0, 0.2, 0.4, 0.3, 0.5]]).T
_TEST_XY = np.array([[0, 1, 1, 0, 0.5, 0.5, 0.5, 0.3, 0.3, 0.3],
                     [0, 0, 1, 1, 0, 0.5, 1, 0, 0.5, 1]]).T

_TEST_FREQ = 4, 2
_TEST_TRI_LINDEX = np.array([[0, 4, 6],
                             [1, 4, 5],
                             [2, 3, 5],
                             [3, 2, 5],
                             [2, 4, 4],
                             [3, 3, 4]])
_TEST_Q_LINDEX = np.array([[2, 0],
                           [2, 1],
                           [3, 1],
                           [4, 1],
                           [1, 2],
                           [2, 2]])
_TEST_DISK_C = np.exp(np.linspace(0, 2j*np.pi, 7))*np.linspace(0, 1, 7)
_TEST_DISK_R = xmath.complex_to_float2d(_TEST_DISK_C)
_TEST_TRI_PTS = np.array([[0.8, 0.6, 0],
                          [0, 1, 0],
                          [0, 0, 1],
                          [0.4, 0.8, 0],
                          [0, 0.5, 0.5],
                          [0.4, 0.3, 0.5],
                          [4/15, 8/15, 1/3]])
_TEST_SPHERE_PTS = xmath.normalize(_TEST_TRI_PTS)

#generic methods, valid on any-d, any shape
def square_to_quad(xy, base_pts):
    """Transforms a square in [0,1]^2 to a (possibly skew) quadrilateral
    defined by base_pts.

    Args:
        xy: Array, shape [..., 2]. The 2d coordinates of the point.
        base_pts: Array, shape [4, ...]. The coordinates of the quadrilateral.
            Should be in counterclockwise order to maintain orientation.
    Returns:
        Coordinates in whatever space base_pts was defined in.

    >>> square_to_quad(_TEST_XY[:, np.newaxis], _TEST_SKEW_Q)
    array([[ 0.70710678,  0.        ,  0.70710678],
           [ 0.        ,  0.70710678,  0.70710678],
           [-0.70710678,  0.        ,  0.70710678],
           [ 0.        , -1.        ,  0.        ],
           [ 0.35355339,  0.35355339,  0.70710678],
           [ 0.        , -0.0732233 ,  0.53033009],
           [-0.35355339, -0.5       ,  0.35355339],
           [ 0.49497475,  0.21213203,  0.70710678],
           [ 0.14142136, -0.24393398,  0.45961941],
           [-0.21213203, -0.7       ,  0.21213203]])

    """
    a, b, c, d = base_pts[0], base_pts[1], base_pts[2], base_pts[3]
    x, y = xy[..., 0], xy[..., 1]
    return a + (b-a)*x + (d-a)*y + (a-b+c-d)*x*y

def tri_bary(bary, base_pts):
    """Transforms barycentric coordinates to a (Euclidean) triangle
    defined by base_pts.

    Args:
        bary: Array, shape [..., 3]. Barycentric coordinates.
        base_pts: Array, shape [3, ...]. Coordinates of the triangle.
            Should be in counterclockwise order to maintain orientation.

    Returns:
        Coordinates in whatever space base_pts was define in.

    >>> tri_bary(_TEST_BARY, _TEST_SKEW_TRI)
    array([[ 0.8 ,  0.6 ,  0.  ],
           [ 0.64,  0.68,  0.  ],
           [ 0.48,  0.56,  0.2 ],
           [ 0.32,  0.44,  0.4 ],
           [ 0.16,  0.62,  0.3 ],
           [ 0.  ,  0.5 ,  0.5 ]])

    """
    return bary.dot(base_pts)

#methods for disks
def square_to_disk(xy, rotation=1):#np.exp(1j*np.pi/4)):
    """Transforms square on [0,1]^2 to unit disk.

    >>> np.round(square_to_disk(_TEST_XY), 6)
    array([[-0.707107, -0.707107],
           [ 0.707107, -0.707107],
           [ 0.707107,  0.707107],
           [-0.707107,  0.707107],
           [ 0.      , -1.      ],
           [ 0.      ,  0.      ],
           [ 0.      ,  1.      ],
           [-0.371391, -0.928477],
           [-0.4     ,  0.      ],
           [-0.371391,  0.928477]])

    """
    pts = 2*xy - 1
    r = np.max(np.abs(pts), axis=-1)
    theta = np.arctan2(pts[..., 1], pts[..., 0])
    result = r*np.exp(1j*theta)*rotation
    return xmath.complex_to_float2d(result)

DISK_TRI_C = np.exp(2j*np.pi/3*np.arange(3))*1j
DISK_TRI_R = xmath.complex_to_float2d(DISK_TRI_C)

def tri_to_disk(bary, rotation=1, pts=DISK_TRI_C):
    """Transforms triangle in barycentric coordinates to unit disk.
    tri_naive_slerp also does this when pts are on a great circle,
    with somewhat different results.

    >>> np.round(tri_to_disk(_TEST_BARY), 6)
    array([[ 0.      ,  1.      ],
           [-0.240192,  0.970725],
           [-0.      ,  0.4     ],
           [ 0.34641 ,  0.2     ],
           [-0.261861, -0.302372],
           [-0.      , -1.      ]])
    """
    tri_pts = bary.dot(pts)
    angle = np.angle(tri_pts)
    r = 1 - 3*bary.min(axis=-1)
    result = r*np.exp(1j*angle)*rotation
    return xmath.complex_to_float2d(result)

DISK_SQ_C = np.array([1, 1j, -1, -1j])
DISK_SQ_R = xmath.complex_to_float2d(DISK_SQ_C)

def _sq_disk(bkdn, base_pts, freq, tweak):
    sc = square_to_disk(bkdn.coord)
    sc2 = sc/np.sqrt(2) + 0.5
    result = square_to_quad(sc2[:, np.newaxis], base_pts)
    return result


def _tri_disk(bkdn, base_pts, freq, tweak):
    rebary = bary_tri(tri_to_disk(bkdn.coord), DISK_TRI_R)
    return tri_bary(rebary, base_pts)

#disk -> sphere

def spherical_to_xyz(phi, theta):
    """Converts spherical coordinates to 3d xyz coordinates.

    Args:
        phi: Inclination
        theta: Azimuth

    Returns:
        An array with of shape = [..., 3].

    >>> phi = np.arccos(np.linspace(-1, 1, 7))
    >>> theta = np.arcsin(np.linspace(-1, 1, 7))
    >>> np.round(spherical_to_xyz(phi, theta), 6)
    array([[ 0.      , -0.      , -1.      ],
           [ 0.555556, -0.496904, -0.666667],
           [ 0.888889, -0.31427 , -0.333333],
           [ 1.      ,  0.      ,  0.      ],
           [ 0.888889,  0.31427 ,  0.333333],
           [ 0.555556,  0.496904,  0.666667],
           [ 0.      ,  0.      ,  1.      ]])
    """
    return np.array([np.sin(phi) * np.cos(theta), # pylint: disable=no-member
                     np.sin(phi) * np.sin(theta),
                     np.cos(phi)]).T

def lambert(disk):
    """Converts coordinates on the disk to spherical coordinates, using
    the Lambert azimuthal equal-area projection.

    Args:
        disk: Array of shape [..., 2] representing points on the disk.

    Returns:
        phi, theta: Spherical coordinates

    >>> phi, theta = lambert(_TEST_DISK_R)
    >>> np.round(phi*180/np.pi, 2)
    array([   0.  ,   19.19,   38.94,   60.  ,   83.62,  112.89,  180.  ])
    >>> np.round(theta*180/np.pi, 2)
    array([   0.,   60.,  120.,  180., -120.,  -60.,   -0.])
    """
    theta = np.arctan2(disk[..., 1], disk[..., 0])
    phi = 2*np.arcsin(np.linalg.norm(disk, axis=-1))
    return phi, theta

def equidistant(disk):
    """Converts coordinates on the disk to spherical coordinates, using
    the azimuthal equal-distance projection.

    Args:
        disk: Array of shape [..., 2] representing points on the disk.

    Returns:
        phi, theta: Spherical coordinates

    >>> phi, theta = equidistant(_TEST_DISK_R)
    >>> np.round(phi*180/np.pi, 2)
    array([   0.,   30.,   60.,   90.,  120.,  150.,  180.])
    >>> np.round(theta*180/np.pi, 2)
    array([   0.,   60.,  120.,  180., -120.,  -60.,   -0.])
    """
    theta = np.arctan2(disk[..., 1], disk[..., 0])
    phi = np.linalg.norm(disk, axis=-1)*np.pi
    return phi, theta

#methods for spheres

#triangles -> spherical triangle

def tri_naive_slerp(bary, base_pts):
    """
    Naive slerp (spherical linear interpolation) on a spherical triangle.

    Args:
        bary: Array, shape [..., 3]. Barycentric coordinates.
        base_pts: Array, shape [3, ..., 3]. Coordinates of the triangle.
            Should be in counterclockwise order to maintain orientation.

    Returns:
        An array of shape [..., 3], representing points in 3d-space.

    >>> tri_naive_slerp(_TEST_BARY, _TEST_EQ_TRI)
    array([[ 1.        ,  0.        ,  0.        ],
           [ 0.95105652,  0.30901699,  0.        ],
           [ 0.80901699,  0.30901699,  0.30901699],
           [ 0.58778525,  0.30901699,  0.58778525],
           [ 0.30901699,  0.70710678,  0.4539905 ],
           [ 0.        ,  0.70710678,  0.70710678]])

    """
    angle = xmath.central_angle_equilateral(base_pts)
    b = np.sin(angle * bary) / np.sin(angle)
    return b.dot(base_pts)

def tri_areal(bary, base_pts):
    """Given a triangle and spherical areal coordinates, returns the vectors
    cooresponding to those coordinates.

    Args:
        bary: Array, shape [..., 3]. Barycentric coordinates.
        base_pts: Array, shape [3, ..., 3]. Coordinates of the triangle.
            Should be in counterclockwise order to maintain orientation.

    Returns:
        An array of shape [..., 3], representing points in 3d-space.

    >>> tri_areal(_TEST_BARY, _TEST_SKEW_TRI)
    array([[ 0.8       ,  0.6       ,  0.        ],
           [ 0.67564273,  0.7372292 ,  0.        ],
           [ 0.59542957,  0.71123145,  0.37364883],
           [ 0.42703426,  0.59682523,  0.67929477],
           [ 0.21874817,  0.81375629,  0.53847   ],
           [ 0.        ,  0.63544017,  0.77215011]])
    """
    base_pts = xmath.normalize(base_pts)
    area = xmath.triangle_solid_angle(base_pts[0], base_pts[1], base_pts[2])
    area_i = bary * area
    base_pts_iplus1 = np.roll(base_pts, -1, axis=0)
    base_pts_iplus2 = np.roll(base_pts, 1, axis=0)
    #FIXME whytf is this commented statement not equivalent to below?
#    L = ((1 + np.cos(area_i))[:, np.newaxis]*
#         np.cross(base_pts_iplus1, base_pts_iplus2) -
#         np.sin(area_i)[:, np.newaxis]*
#         (base_pts_iplus1 + base_pts_iplus2)).transpose((0,2,1))
    L0 = ((1 + np.cos(area_i[..., 0]))[..., np.newaxis]*
          np.cross(base_pts[1], base_pts[2]) -
          np.sin(area_i[..., 0])[..., np.newaxis]*
          (base_pts[1] + base_pts[2]))
    L1 = ((1 + np.cos(area_i[..., 1]))[..., np.newaxis]*
          np.cross(base_pts[2], base_pts[0]) -
          np.sin(area_i[..., 1])[..., np.newaxis]*
          (base_pts[2] + base_pts[0]))
    L2 = ((1 + np.cos(area_i[..., 2]))[..., np.newaxis]*
          np.cross(base_pts[0], base_pts[1]) -
          np.sin(area_i[..., 2])[..., np.newaxis]*
          (base_pts[0] + base_pts[1]))
    L = np.stack([L0, L1, L2], axis=-2)
    h = np.sin(area_i)*(1 + np.sum(base_pts_iplus1*base_pts_iplus2, axis=-1))
    return np.linalg.solve(L, h)

def triangles_method2(lindex, base_pts, freq):
    """Triangles of method 2

    Args:
        lindex: Array, shape (..., 3). Linear indexes (should correspond
            to freq)
        base_pts: Array, shape (3, ..., 3). Coordinates of the triangle.
            Should be in counterclockwise order to maintain orientation.
        freq: 2-tuple. Frequency of the subdivision.

    Returns:
        Array of shape (..., 3, 3)

    >>> np.round(triangles_method2(_TEST_TRI_LINDEX[1:3], _TEST_SKEW_TRI,
    ...                            _TEST_FREQ), 6)

    array([[[ 0.205392,  0.183498,  0.051016],
            [ 0.34558 ,  0.317059,  0.10024 ],
            [ 0.098942,  0.085439,  0.03013 ]],
    <BLANKLINE>
           [[ 0.283136,  0.377856,  0.038599],
            [ 0.24816 ,  0.339397,  0.042048],
            [ 0.152133,  0.199126,  0.028172]]])
    """
    n, m = freq
    frame = breakdown.frame_triangle(base_pts, n, m, interp=xmath.slerp)
    #get the normal to the great circle corresponding to the lines
    #don't need to normalize this
    gc_normals = np.cross(frame[..., 0, :], frame[..., 1, :])
    index = np.arange(3)
    pairs = gc_normals[index, lindex[:, index]]
    #intersection of great circles = cross product of normals
    ptx = np.cross(pairs, np.roll(pairs, 2, axis=1))
    # cross product could give the point we want or its negative.
    # test to see if the points are on the correct side of the sphere
    # take the dot product of these vectors with the center of the
    # base face. if it's positive, it's right, if not, negate it
    center = np.sum(base_pts, axis=0)#don't need to normalize this
    sign_correct = np.sum(center*ptx, axis=-1, keepdims=True) >= 0
    result = np.where(sign_correct, ptx, -ptx)
    return result

def tri_intersections(lindex, base_pts, freq, tweak=False):
    """Transforms a triangle to a spherical triangle using the method of
    intersections

    Args:
        lindex: Array, shape (..., 3). Linear indexes (should correspond
            to freq)
        base_pts: Array, shape (3, ..., 3). Coordinates of the triangle.
            Should be in counterclockwise order to maintain orientation.
        freq: 2-tuple. Frequency of the subdivision.
        tweak: Whether to normalize the points of the triangle before
            taking their centroid. By default false.

    Returns:
        Array of shape (..., 3)

    >>> tri_intersections(_TEST_TRI_LINDEX, _TEST_SKEW_TRI, _TEST_FREQ)
    array([[ 0.8       ,  0.6       ,  0.        ],
           [ 0.73045891,  0.65102869,  0.20524575],
           [ 0.59608649,  0.79682621,  0.09756895],
           [ 0.4472136 ,  0.89442719,  0.        ],
           [ 0.6127549 ,  0.67069951,  0.41582597],
           [ 0.46297244,  0.82965693,  0.30963537]])
    """
    a, b = freq
    pts = triangles_method2(lindex, base_pts, freq)
    if tweak:
        pts = xmath.normalize(pts)
    result = pts.mean(axis=1)
    result[(lindex[:, 0] == 0) &
           (lindex[:, 1] == a) &
           (lindex[:, 2] == a + b)] = base_pts[0]
    result[(lindex[:, 0] == a + b) &
           (lindex[:, 1] == 0) &
           (lindex[:, 2] == a)] = base_pts[1]
    result[(lindex[:, 0] == a) &
           (lindex[:, 1] == a + b) &
           (lindex[:, 2] == 0)] = base_pts[2]
    return result


#squares -> spherical quadrilateral

def square_naive_slerp(xy, base_pts):
    """
    Naive slerp (spherical linear interpolation) on a spherical square.

    Args:
        xy: Array, shape [..., 2]. XY coordinates on the square.
        base_pts: Array, shape [4, ..., 3]. Coordinates of the square.
            Should be in counterclockwise order to maintain orientation.

    Returns:
        An array of shape [..., 3], representing points in 3d-space.

    >>> square_naive_slerp(_TEST_XY, _TEST_EQ_Q)
    array([[ 0.70710678,  0.        ,  0.70710678],
           [ 0.        ,  0.70710678,  0.70710678],
           [-0.70710678,  0.        ,  0.70710678],
           [ 0.        , -0.70710678,  0.70710678],
           [ 0.40824829,  0.40824829,  0.81649658],
           [ 0.        ,  0.        ,  0.84529946],
           [-0.40824829, -0.40824829,  0.81649658],
           [ 0.54634285,  0.25231132,  0.79865417],
           [ 0.164878  , -0.164878  ,  0.84066882],
           [-0.25231132, -0.54634285,  0.79865417]])

    """
    angle = xmath.central_angle_equilateral(base_pts)
    x, y = xy[..., 0], xy[..., 1]
    a = (1-x)*(1-y)
    b = x*(1-y)
    c = x*y
    d = (1-x)*y
    mat = np.sin(np.stack([a, b, c, d], axis=-1)*angle) / np.sin(angle)
    return mat.dot(base_pts)

def square_naive_slerp_2(xy, base_pts):
    """
    Variant naive slerp on a spherical square.

    Args:
        xy: Array, shape [..., 2]. XY coordinates on the square.
        base_pts: Array, shape [4, ..., 3]. Coordinates of the square.
            Should be in counterclockwise order to maintain orientation.

    Returns:
        An array of shape [..., 3], representing points in 3d-space.

    >>> square_naive_slerp_2(_TEST_XY, _TEST_EQ_Q)
    array([[ 0.70710678,  0.        ,  0.70710678],
           [ 0.        ,  0.70710678,  0.70710678],
           [-0.70710678,  0.        ,  0.70710678],
           [ 0.        , -0.70710678,  0.70710678],
           [ 0.40824829,  0.40824829,  0.81649658],
           [ 0.        ,  0.        ,  0.94280904],
           [-0.40824829, -0.40824829,  0.81649658],
           [ 0.54634285,  0.25231132,  0.79865417],
           [ 0.16975918, -0.16975918,  0.9222064 ],
           [-0.25231132, -0.54634285,  0.79865417]])

    """
    angle = xmath.central_angle_equilateral(base_pts)
    x, y = xy[..., 0], xy[..., 1]
    sx = np.sin(x*angle)
    sy = np.sin(y*angle)
    scx = np.sin((1-x)*angle)
    scy = np.sin((1-y)*angle)
    a = scx * scy
    b = sx * scy
    c = sx * sy
    d = scx * sy
    mat = np.stack([a, b, c, d], axis=-1) / np.sin(angle)**2
    return mat.dot(base_pts)

def _square_slerp(xy, base_pts):
    """
    Helper function for square_slerp. This does the slerp, and then
    square_slerp averages the two orientations together.

    Args:
        xy: Array, shape [..., 2]. XY coordinates on the square.
        base_pts: Array, shape [4, ..., 3]. Coordinates of the square.
            Should be in counterclockwise order to maintain orientation.

    Returns:
        An array of shape [..., 3], representing points in 3d-space.

    >>> _square_slerp(_TEST_XY[:, np.newaxis], _TEST_SKEW_Q)
    array([[ 0.70710678,  0.        ,  0.70710678],
           [ 0.        ,  0.70710678,  0.70710678],
           [-0.70710678,  0.        ,  0.70710678],
           [ 0.        , -1.        ,  0.        ],
           [ 0.40824829,  0.40824829,  0.81649658],
           [-0.06780818, -0.22086837,  0.97294358],
           [-0.5       , -0.70710678,  0.5       ],
           [ 0.54634285,  0.25231132,  0.79865417],
           [ 0.1721895 , -0.48808406,  0.85564287],
           [-0.32101976, -0.89100652,  0.32101976]])
    """
    a, b, c, d = base_pts[0], base_pts[1], base_pts[2], base_pts[3]
    x, y = xy[..., 0], xy[..., 1]
    ab = xmath.slerp(a, b, x)
    dc = xmath.slerp(d, c, x)
    result = xmath.slerp(ab, dc, y)
    return result

def square_slerp(xy, base_pts):
    """Transforms a square in [0,1]^2 to a spherical quadrilateral
    defined by base_pts, using spherical linear interpolation

    Args:
        xy: Array, shape [..., 2]. XY coordinates on the square.
        base_pts: Array, shape [4, ..., 3]. Coordinates of the square.
            Should be in counterclockwise order to maintain orientation.

    Returns:
        An array of shape [..., 3], representing points in 3d-space.

    >>> square_slerp(_TEST_XY[:, np.newaxis], _TEST_SKEW_Q)
    array([[ 0.70710678,  0.        ,  0.70710678],
           [ 0.        ,  0.70710678,  0.70710678],
           [-0.70710678,  0.        ,  0.70710678],
           [ 0.        , -1.        ,  0.        ],
           [ 0.40824829,  0.40824829,  0.81649658],
           [ 0.        , -0.2213779 ,  0.9751881 ],
           [-0.5       , -0.70710678,  0.5       ],
           [ 0.54634285,  0.25231132,  0.79865417],
           [ 0.21865594, -0.4721394 ,  0.85397539],
           [-0.32101976, -0.89100652,  0.32101976]])
    """
    #have to do this twice and average them because just doing nested
    #slerp isn't symmetric
    one = _square_slerp(xy, base_pts)
    two = _square_slerp(xy[..., ::-1], base_pts[[0, 3, 2, 1]])
    return xmath.normalize(one + two)

def square_intersections(lindex, base_pts, freq):
    """Transforms a square to a spherical quadrilateral using the method of
    intersections

    Args:
        lindex: Array, shape (..., 2). Linear indexes (should correspond
            to freq)
        base_pts: Array, shape (4, ..., 3). Coordinates of the triangle.
            Should be in counterclockwise order to maintain orientation.
        freq: 2-tuple. Frequency of the subdivision.

    Returns:
        Array of shape (..., 3)

    >>> square_intersections(_TEST_Q_LINDEX, _TEST_SKEW_Q, _TEST_FREQ)
    array([[ 0.70710678,  0.        ,  0.70710678],
           [ 0.43240784, -0.07811493,  0.54287903],
           [ 0.39952941,  0.13802823,  0.63415567],
           [ 0.28927895,  0.28927895,  0.57855789],
           [ 0.32664074, -0.46193977,  0.32664074],
           [ 0.4267767 , -0.25      ,  0.78033009]])
    """
    a, b = freq
    preframe = breakdown.frame_square(a, b)
    frame = _square_slerp(preframe[..., np.newaxis, :], base_pts)
    gc_normals = np.cross(frame[..., 0, :], frame[..., 1, :])
    index = np.arange(2)
    pairs = gc_normals[index, lindex[:, index]]
    #intersection of great circles = cross product of normals
    ptx = np.cross(pairs[:, 0], pairs[:, 1])
    # cross product could give the point we want or its negative.
    # test to see if the points are on the correct side of the sphere
    # take the dot product of these vectors with the center of the
    # base face. if it's positive, it's right, if not, negate it
    center = np.sum(base_pts, axis=0)#don't need to normalize this
    sign_correct = np.sum(center*ptx, axis=-1, keepdims=True) >= 0
    result = np.where(sign_correct, ptx, -ptx)
    result[(lindex[:, 0] == b)     & (lindex[:, 1] == 0)] = base_pts[0]
    result[(lindex[:, 0] == a + b) & (lindex[:, 1] == b)] = base_pts[1]
    result[(lindex[:, 0] == a)     & (lindex[:, 1] == a + b)] = base_pts[2]
    result[(lindex[:, 0] == 0)     & (lindex[:, 1] == a)] = base_pts[3]
    return result

#extra crap
def bary_tri(tri, vertices):
    """Transforms a triangle back into barycentric
    coordinates. Only elements [..., :2] are used.

    Args:
        tri: Triangle to calculate barycentric coordinates with respect to.
        vertices: Array of vertices

    Returns:
        Array of barycentric coordinates

    >>> bary_tri(_TEST_TRI_PTS, _TEST_SKEW_TRI)
    array([[ 1.        ,  0.        ,  0.        ],
           [ 0.        ,  1.        ,  0.        ],
           [ 0.        ,  0.        ,  1.        ],
           [ 0.5       ,  0.5       ,  0.        ],
           [ 0.        ,  0.5       ,  0.5       ],
           [ 0.5       ,  0.        ,  0.5       ],
           [ 0.33333333,  0.33333333,  0.33333333]])
    """
    afill = np.ones(vertices.shape[:-1])
    a = np.concatenate([vertices[..., :2], afill[..., np.newaxis]], axis=-1)
    bfill = np.ones(tri.shape[:-1])
    b = np.concatenate([tri[..., :2], bfill[..., np.newaxis]], axis=-1)
    #couldn't make np.linalg.solve cooperate here.
    #this should be OK numerically
    ainv = np.linalg.inv(a)
    return b.dot(ainv)

def to_sph_areal_coords(pts, triangle):
    """Given a triangle and pts within that triangle, returns the
    spherical areal coordinates of the pts with respect to the triangle.

    Args:
        tri: Array with shape (..., 3). Triangle to calculate spherical
            areal coordinates with respect to.
        vertices: Array with shape (..., 3). Array of vertices

    Returns:
        Array of spherical areal coordinates

    >>> to_sph_areal_coords(_TEST_SPHERE_PTS, _TEST_SKEW_TRI)
    array([[ 1.        ,  0.        ,  0.        ],
           [ 0.        ,  1.        ,  0.        ],
           [ 0.        ,  0.        ,  1.        ],
           [ 0.5       ,  0.5       ,  0.        ],
           [ 0.        ,  0.5595371 ,  0.4404629 ],
           [ 0.5595371 ,  0.        ,  0.4404629 ],
           [ 0.36751409,  0.36751409,  0.26497182]])
    """
    area = xmath.triangle_solid_angle(triangle[0], triangle[1], triangle[2])
    area_i = xmath.triangle_solid_angle(pts[:, np.newaxis],
                                        np.roll(triangle, 1, axis=0),
                                        np.roll(triangle, -1, axis=0))
    return area_i/area



def project_sphere(sphcoords, zfunc=np.arcsin, scale=180 / np.pi):
    """
    Projects 3d coordinates on the sphere onto a 2d rectangle.
    Corresponds to a number of rectangular map projections.

    Args:
        sphcoords: An array of shape (..., 3).
        zfunc: A function to transform the z-values on the sphere. By
               default this is np.arcsin, which makes the projection a
               "rectangular" map projection. Use zfunc = lambda x: x
               for an equal-area projection, and np.arctanh for Meractor.
        scale: A scale function, applied to both coordinates of the result.
               By default this is 180/np.pi, to
               transform radians into degrees.

    Returns:
        The 2d rectangular coordinates, in an array of shape (..., 2).
        By default, returns latitude and longitude, but if zfunc is
        specified, the second coordinate will be whatever the function
        transforms it to be.

    >>> project_sphere(_TEST_SPHERE_PTS)
    array([[ 36.86989765,   0.        ],
           [ 90.        ,   0.        ],
           [  0.        ,  90.        ],
           [ 63.43494882,   0.        ],
           [ 90.        ,  45.        ],
           [ 36.86989765,  45.        ],
           [ 63.43494882,  29.20593225]])
    """
    # specify shape of result
    newdim = list(sphcoords.shape)
    newdim[-1] = 2
    result = np.empty(newdim)
    # populate the array
    result[..., 0] = np.arctan2(sphcoords[..., 1],
                                sphcoords[..., 0]) * scale
    result[..., 1] = zfunc(sphcoords[..., 2]) * scale
    return result

#parallel projections
def parallel_exact(pts, normal):
    """Projects points exactly onto the sphere parallel to the normal vector.
    Args:
        pts: Points to project
        normal: Normal vector

    >>> center = np.array([0, 0, 1])
    >>> pts = np.array([[0.5, 0.5, 0],
    ...                 [1, 0, 0]])
    >>> parallel_exact(pts, center)
    array([[ 0.        ,  0.        ,  0.70710678],
           [ 0.        ,  0.        ,  0.        ]])
    """
    vdotc = np.sum(pts * normal, axis=-1)
    vdotv = norm(pts, axis=-1)**2
    p = -vdotc + np.sqrt(np.fmax(1 + vdotc**2 - vdotv, 0))
    return p[..., np.newaxis] * normal


def parallel_approx(pts, normal):
    """Approximately projects points onto the sphere parallel to the
        normal vector.
    Args:
        pts: Points to project
        normal: Normal vector

    >>> center = np.array([0, 0, 1])
    >>> pts = np.array([[0.5, 0.5, 0],
    ...                 [1, 0, 0]])
    >>> parallel_approx(pts, center)
    array([[ 0.        ,  0.        ,  0.29289322],
           [ 0.        ,  0.        ,  0.        ]])
        """
    q = 1 - norm(pts, axis=-1)
    return q[..., np.newaxis] * normal

def parallel(pts, normal, exact=True):
    """Projects points onto the sphere parallel to the normal vector.
    Args:
        pts: Points to project
        normal: Normal vector
        exact: Whether to project exactly or approximately.
            Defaults to exact (True).
    """
    if exact:
        result = parallel_exact(pts, normal)
    else:
        result = parallel_approx(pts, normal)
    return result

_FLAT = {3: lambda bkdn, base_pts, freq, tweak: tri_bary(bkdn.coord, base_pts),
         4: lambda bkdn, base_pts, freq, tweak:
            square_to_quad(bkdn.coord[:, np.newaxis], base_pts)}

_SLERP = {3: lambda bkdn, base_pts, freq, tweak: tri_naive_slerp(bkdn.coord, base_pts),
          4: lambda bkdn, base_pts, freq, tweak:
             square_naive_slerp(bkdn.coord, base_pts)}

_SLERP2 = {3: lambda bkdn, base_pts, freq, tweak: tri_naive_slerp(bkdn.coord, base_pts),
           4: lambda bkdn, base_pts, freq, tweak:
              square_naive_slerp_2(bkdn.coord, base_pts)}

_OTHER = {3: lambda bkdn, base_pts, freq, tweak: tri_areal(bkdn.coord, base_pts),
          4: lambda bkdn, base_pts, freq, tweak:
             square_slerp(bkdn.coord[:, np.newaxis], base_pts)}

_GC = {3: lambda bkdn, base_pts, freq, tweak:
          tri_intersections(bkdn.lindex, base_pts, freq, tweak),
       4: lambda bkdn, base_pts, freq, tweak:
          square_intersections(bkdn.lindex, base_pts, freq)}

_DISK = {3: _tri_disk,
         4: _sq_disk}

PROJECTIONS = {'flat':  _FLAT,
               'nslerp': _SLERP,
               'nslerp2': _SLERP2,
               'other': _OTHER,
               'gc':    _GC,
               'disk':  _DISK}

PARALLEL = ['nslerp', 'nslerp2', 'disk']
