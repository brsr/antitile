# -*- coding: utf-8 -*-
"""
Extra math and array functions.
"""

import numpy as np
from numpy.linalg import norm
import pandas as pd

def line_intersection(a1, a2, b1, b2):
    a1x, a1y = a1[..., 0], a1[..., 1]
    a2x, a2y = a2[..., 0], a2[..., 1]
    b1x, b1y = b1[..., 0], b1[..., 1]
    b2x, b2y = b2[..., 0], b2[..., 1]
    
    
    numx = (a1x*a2y - a1y*a2x)*(b1x - b2x) - (b1x*b2y - b1y*b2x)*(a1x - a2x)
    numy = (a1x*a2y - a1y*a2x)*(b1y - b2y) - (b1x*b2y - b1y*b2x)*(a1y - a2y)
    denom = (a1x-a2x)*(b1y-b2y) - (a1y-a2y)*(b1x-b2x)
    num = np.stack([numx, numy], axis=-1)
    return num/denom[..., np.newaxis]

#def square_to_quad(xy, abcd):
#    """Transforms a square in [0,1]^2 to a (possibly skew) quadrilateral
#    defined by abcd.
#
#    Args:
#        xy: The 2d coordinates of the point. Last element of shape should be 2.
#            If not, treated like xy[..., :2]: the values in xy[..., 2:] are
#            ignored.
#        abcd: The coordinates of the quadrilateral. Should be in
#            counterclockwise order to maintain orientation. First element of
#            shape should be 4.
#    Returns:
#        Coordinates in whatever space abcd was defined in."""
#    a, b, c, d = abcd[0], abcd[1], abcd[2], abcd[3]
#    x, y = xy[..., 0], xy[..., 1]
#    return a + (b-a)*x + (d-a)*y + (a-b+c-d)*x*y
#
#def square_to_sphquad(xy, abcd):
#    """Transforms a square in [0,1]^2 to a spherical quadrilateral
#    defined by abcd.
#
#    Args:
#        xy: The 2d coordinates of the point. Last element of shape should be 2.
#            If not, treated like xy[..., :2]: the values in xy[..., 2:] are
#            ignored.
#        abcd: The coordinates of the quadrilateral. Should be in
#            counterclockwise order to maintain orientation. Shape should be
#            (4, ..., 3)
#    Returns:
#        Coordinates on the sphere."""
#    a, b, c, d = abcd[0], abcd[1], abcd[2], abcd[3]
#    x, y = xy[..., 0], xy[..., 1]
#    ab = slerp(a, b, x)
#    dc = slerp(d, c, x)
#    return slerp(ab, dc, y)

def spherical_to_xyz(phi, theta):
    """Converts spherical coordinates to 3d xyz coordinates"""
    return np.array([np.sin(phi) * np.cos(theta), # pylint: disable=no-member
                     np.sin(phi) * np.sin(theta),
                     np.cos(phi)]).T

def lambert(disk):
    """Converts coordinates on the disk to spherical coordinates, using
    the Lambert azimuthal equal-area projection."""
    theta = np.arctan2(disk[..., 1], disk[..., 0])
    phi = 2*np.arcsin(norm(disk, axis=-1))
    return phi, theta

def equidistant(disk):
    """Converts coordinates on the disk to spherical coordinates, using
    the azimuthal equal-distance projection."""
    theta = np.arctan2(disk[..., 1], disk[..., 0])
    phi = norm(disk, axis=-1)*np.pi
    return phi, theta


def record_initialize(shape, dtype, default_bool=False,
                      default_int=-1,
                      default_float=np.nan):
    """
    Creates and initializes a record array with specified defaults
    instead of whatever garbage happened to already be there.

    >>> dtype = np.dtype([('a', bool), ('b', int), ('c', float)])
    >>> record_initialize(1, dtype)  #doctest: +ELLIPSIS
    rec.array([(False, -1, nan)], ...])
    """
    sctypes = np.sctypes
    result = np.recarray(shape, dtype=dtype)
    fields = dtype.fields
    for field, t in fields.items():
        base = t[0].base
        if base in sctypes['int'] + sctypes['uint']:
            result[field] = default_int
        elif base in sctypes['float'] + sctypes['complex']:
            result[field] = default_float
        elif base == np.bool:
            result[field] = default_bool
    return result


def transpose(long_array, filler=-1):
    """
    Transposes an array of shape N x 2 into a wide 2d array, using the first
    column as an index. Blanks in the array are filled in with the argument
    filler (by default -1)
    >>> x = np.tile(np.arange(2),5)[:9]
    >>> y = np.arange(9)
    >>> transpose(np.stack((x,y),axis=-1))
    array([[ 0,  2,  4,  6,  8],
           [ 1,  3,  5,  7, -1]])

    """
    df = pd.DataFrame(long_array, columns=["a", "b"])
    df["rank"] = (df.groupby('a').b.rank() - 1)
    step2 = df[["a", 'b', "rank"]].set_index(["a", "rank"]).unstack()
    return step2.fillna(filler).values.astype(long_array.dtype)


def recordify(names, arrays):
    """Take a bunch of arrays of the same first dimension and combine them
    into a record array

    >>> x = np.linspace(0,1,4)
    >>> y = np.eye(4)
    >>> z = np.arange(4**3).reshape(4,4,4)
    >>> r = recordify(['x','y','z'],[x,y,z])
    >>> r.dtype
    dtype((numpy.record, [('x', '<f8'), ('y', '<f8', (4,)), ('z', '<i4', (4, 4))]))"""
    type_list = [(name, array.dtype, array.shape[1:])
                 for (name, array) in zip(names, arrays)]
    wtype = np.dtype(type_list)
    result = np.recarray(len(arrays[0]), dtype=wtype)
    for name, array in zip(names, arrays):
        result[name] = array
    return result


def renumber(in_array, fill=-1):
    """Renumbers an index of a subarray given by the boolean array c.
    (This will probably make more sense if you look at the doctest.)

    Arguments:
        inarray: A 1d boolean array
        fill: An integer to fill in elements that don't appear in the subarray
            By default, -1.
    Returns:
        A 1d numeric array. Takes the value of fill if element does not
            exist in the subarray.
    >>> x = np.arange(5)
    >>> condition = (x > 2)
    >>> y = x[condition]
    >>> i = renumber(condition)
    >>> np.where(i >= 0, y[i], -1)
    array([-1, -1, -1,  3,  4])
    """
    count = np.sum(in_array)
    index = np.repeat(fill, in_array.shape)
    index[in_array] = np.arange(count)
    return index


def triple_product(a, b, c):
    """The scalar triple product of 3 vectors
    a dot b cross c = determinant of [a b c]
    a,b, and c must have last dimension = 3
    >>> a = np.zeros(3)
    >>> a[0] = 1
    >>> b = np.arange(4,16).reshape(4,3)
    >>> c = np.arange(12)[::-1].reshape(4,3)
    >>> triple_product(a,b,c)
    array([-15., -15., -15., -15.])
    """
    return np.sum(a * np.cross(b, c), axis=-1)


def normalize(vectors, axis=-1):
    """Normalizes vectors in n-space. The zero vector stays the zero vector.

    >>> x = np.stack((np.ones(5),np.arange(5)),axis=-1)
    >>> normalize(x) #doctest: +NORMALIZE_WHITESPACE
    array([[ 1.        ,  0.        ],
           [ 0.70710678,  0.70710678],
           [ 0.4472136 ,  0.89442719],
           [ 0.31622777,  0.9486833 ],
           [ 0.24253563,  0.9701425 ]])
    """
    n = norm(vectors, axis=axis, keepdims=True)
    return np.where(n <= 0, 0, vectors / n)


def slerp(pt1, pt2, intervals, axis=-1):
    """Spherical linear interpolation.
    >>> x = np.array([1,0,0])
    >>> y = np.array([0,0,1])
    >>> t = np.linspace(0, 1, 4)[:, np.newaxis]
    >>> slerp(x, y, t) #doctest: +NORMALIZE_WHITESPACE
    array([[ 1.       ,  0.       ,  0.       ],
           [ 0.8660254,  0.       ,  0.5      ],
           [ 0.5      ,  0.       ,  0.8660254],
           [ 0.       ,  0.       ,  1.       ]])
    """
    t = intervals
    a = np.arccos(np.sum(pt1 * pt2, axis=axis, keepdims=True))    
    return (np.sin((1 - t)*a)*pt1 + np.sin((t)*a)*pt2)/np.sin(a)


def lerp(pt1, pt2, intervals):
    """Linear interpolation.
    >>> x = np.array([1,0,0])
    >>> y = np.array([0,0,1])
    >>> t = np.linspace(0, 1, 4)[:, np.newaxis]
    >>> lerp(x, y, t) #doctest: +NORMALIZE_WHITESPACE
    array([[ 1.        ,  0.        ,  0.        ],
           [ 0.66666667,  0.        ,  0.33333333],
           [ 0.33333333,  0.        ,  0.66666667],
           [ 0.        ,  0.        ,  1.        ]])
    """
    t = intervals
    return (1 - t)*pt1 + t*pt2


def nlerp(*args):
    """Normalized linear interpolation.
    >>> x = np.array([1,0,0])
    >>> y = np.array([0,0,1])
    >>> t = np.linspace(0, 1, 4)[:, np.newaxis]
    >>> nlerp(x, y, t) #doctest: +NORMALIZE_WHITESPACE
    array([[ 1.        ,  0.        ,  0.        ],
           [ 0.89442719,  0.        ,  0.4472136 ],
           [ 0.4472136 ,  0.        ,  0.89442719],
           [ 0.        ,  0.        ,  1.        ]])
           """
    return normalize(lerp(*args))


def distance(x, y, p=2, axis=-1, scale=1):
    """Distance between points in Euclidean space.
    Args:
        x, y: Coordinates of points.
        p: Which norm to use. 2 = euclidean, 1 = taxicab.
        axis: Which axis the vectors lie along. By default, -1.
        scale: Scale factor for the distance.
    Returns: Array of distances.

    >>> t = np.linspace(0,1,5)[:,np.newaxis]
    >>> x = np.array([[0,0,0]])*t+np.array([[0,10,-10]])*(1-t)
    >>> y = np.array([[0,0,0]])*t+np.array([[10,0,-10]])*(1-t)
    >>> np.round(distance(x, y), 2) #doctest: +NORMALIZE_WHITESPACE
    array([ 14.14,  10.61,   7.07,   3.54,   0.  ])
    """
    return norm((x - y) / scale, ord=p, axis=axis)


def triangle_area(a, b, c):
    """Area of Euclidean triangle given by a, b, and c.

    Args:
        a, b, c: Coordinates of points.

    Returns: Array of areas.
    >>> t = np.linspace(0,np.pi,5)
    >>> a = np.stack([np.cos(t), np.sin(t), np.zeros(5)],axis=-1)
    >>> b = np.array([0,1,1])/np.sqrt(2)
    >>> c = np.array([0,-1,1])/np.sqrt(2)
    >>> np.round(triangle_area(a, b, c), 4)
    array([ 0.866 ,  0.7071,  0.5   ,  0.7071,  0.866 ])
    """
    ab = a - b
    ac = a - c
    return norm(np.cross(ab, ac), axis=-1) / 2


def bearing(origin, destination, pole=np.array([0, 0, 1])):
    """ Returns the bearing (angle) between points. By default,
        the bearing is calculated with respect to the +y direction.

    Args:
        origin: Origin points
        destination: Destination points
        direction: A vector giving the direction the bearing is
            calculated with respect to. By default, [0, 1].
    Returns: Array of bearings.
    >>> a = np.array([0,0])
    >>> b = np.array([1,0])
    >>> c = np.array([0,1])

    >>> bearing(a,b,c)/np.pi*180 #doctest: +ELLIPSIS
    90...
    """
    direction = origin - pole
    pv = destination - origin
    d = np.sum(pv * direction, axis=-1)
    x = norm(np.cross(pv, direction), axis=-1)
    return np.arctan2(x, d)


def spherical_distance(x, y):
    """Spherical distance, i.e. central angle, between vectors.
    Args:
        x, y: Coordinates of points on the sphere.
        axis: Which axis the vectors lie along. By default, -1.
    Returns: Array of spherical distances.

    >>> t = np.linspace(0,np.pi,5)
    >>> c = np.cos(t)
    >>> s = np.sin(t)
    >>> z = np.zeros(t.shape)
    >>> x = np.stack((c,s,z),axis=-1)
    >>> y = np.stack((c,z,s),axis=-1)
    >>> spherical_distance(x,y)/np.pi*180 # doctest: +NORMALIZE_WHITESPACE
    array([  0.,  60.,  90.,  60.,   0.])
    """
    # technically this is not the most numerically sound way to do this.
    # if issues arise, can change it.
    return np.arccos(np.clip(np.sum(x * y, axis=-1), -1, 1))


def spherical_triangle_area(a, b, c):
    """Spherical area, i.e. solid angle, of a triangle. Note there are two
    areas defined by three points on a sphere: inside the triangle and
    outside it. This will always return the smaller of the two. (The other
    is 4*pi minus what this function returns.)

    Args:
        a, b, c: Coordinates of points on the sphere.

    Returns: Array of spherical areas.

    >>> t = np.linspace(0,np.pi,5)
    >>> a = np.stack([np.cos(t), np.sin(t), np.zeros(5)],axis=-1)
    >>> b = np.array([0,1,1])/np.sqrt(2)
    >>> c = np.array([0,-1,1])/np.sqrt(2)
    >>> np.round(spherical_triangle_area(a, b, c), 4)
    array([ 1.5708,  1.231 ,  0.    ,  1.231 ,  1.5708])
    """
    top = np.abs(triple_product(a, b, c))
    bottom = (1 + np.sum(a * b, axis=-1)
                + np.sum(b * c, axis=-1)
                + np.sum(c * a, axis=-1))
    return 2 * (np.arctan(top / bottom) % np.pi)


def spherical_bearing(origin, destination, pole=np.array([0, 0, 1])):
    """ Returns the bearing (angle) between points. By default,
        the bearing is calculated with respect to the north pole.
        Can also be considered as the angle adjacent to origin in the
        triangle formed by origin, destination, and pole.

    Args:
        origin: Origin points
        destination: Destination points
        pole: Point bearing is calculated with respect to.
            By default, the north pole.

    Returns: Array of bearings.

    >>> x = np.array([1,0,0])
    >>> spherical_bearing(x,np.roll(x,1))/np.pi
    0.5
    """
    c_1 = np.cross(origin, destination)
    c_2 = np.cross(origin, pole)
    cos_theta = np.sum(c_1 * c_2, axis=-1)
    sin_theta = triple_product(origin, destination, pole)
    return np.arctan2(sin_theta, cos_theta)

def to_sph_areal_coords(pts, triangle):
    """Given a triangle and pts within that triangle, returns the
    spherical areal coordinates of the pts with respect to the triangle.
    >>> triangle = np.eye(3)
    >>> triangle[2,1] = 1
    >>> triangle = normalize(triangle)
    >>> x = np.linspace(0,1,5)
    >>> z = np.linspace(0.7,0,5)
    >>> y = np.sqrt(1-x**2-z**2)
    >>> pts = normalize(np.stack([x,y,z],axis=-1))
    >>> to_sph_areal_coords(pts, triangle) # doctest: +NORMALIZE_WHITESPACE
    array([[ 0.        ,  0.01273324,  0.98726676],
           [ 0.12972226,  0.2358742 ,  0.63440354],
           [ 0.27122545,  0.34291999,  0.38585456],
           [ 0.45754141,  0.35616745,  0.18629114],
           [ 1.        ,  0.        ,  0.        ]])
    """
    area = spherical_triangle_area(triangle[0], triangle[1], triangle[2])
    area_i = spherical_triangle_area(pts[:, np.newaxis],
                                     np.roll(triangle, 1, axis=0),
                                     np.roll(triangle, -1, axis=0))
    return area_i/area

def from_sph_areal_coords(beta, triangle):
    """Given a triangle and spherical areal coordinates, returns the vectors
    cooresponding to those coordinates.

    Args:
        beta: spherical areal coordinates
        triangle: vertices of the spherical triangle in a 3x3 array

    Returns: Points on the sphere

    >>> beta = np.array([[ 1, 0.8, 0.6, 0.4, 0.2, 0.0 ],
    ...                  [ 0, 0.2, 0.2, 0.2, 0.5, 0.5 ],
    ...                  [ 0, 0.0, 0.2, 0.4, 0.3, 0.5 ]]).T
    >>> triangle = np.eye(3)
    >>> triangle[2,1] = 1
    >>> triangle = normalize(triangle)
    >>> from_sph_areal_coords(beta, triangle) # doctest: +NORMALIZE_WHITESPACE
    array([[ 1.        ,  0.        ,  0.        ],
           [ 0.97123031,  0.23814214,  0.        ],
           [ 0.87887868,  0.44073244,  0.18255735],
           [ 0.68229421,  0.63250197,  0.3666277 ],
           [ 0.37935999,  0.88556406,  0.26807143],
           [ 0.        ,  0.92387953,  0.38268343]])
    """
    area = spherical_triangle_area(triangle[0], triangle[1], triangle[2])
    area_i = beta * area
    triangle_iplus1 = np.roll(triangle, -1, axis=0)
    triangle_iplus2 = np.roll(triangle, 1, axis=0)
    #FIXME whytf is this commented statement not equivalent to below?
#    L = ((1 + np.cos(area_i))[:, np.newaxis]*
#         np.cross(triangle_iplus1, triangle_iplus2) -
#         np.sin(area_i)[:, np.newaxis]*
#         (triangle_iplus1 + triangle_iplus2)).transpose((0,2,1))
    L0 = ((1 + np.cos(area_i[..., 0]))[..., np.newaxis]*
          np.cross(triangle[1], triangle[2]) -
          np.sin(area_i[..., 0])[..., np.newaxis]*
          (triangle[1] + triangle[2]))
    L1 = ((1 + np.cos(area_i[..., 1]))[..., np.newaxis]*
          np.cross(triangle[2], triangle[0]) -
          np.sin(area_i[..., 1])[..., np.newaxis]*
          (triangle[2] + triangle[0]))
    L2 = ((1 + np.cos(area_i[..., 2]))[..., np.newaxis]*
          np.cross(triangle[0], triangle[1]) -
          np.sin(area_i[..., 2])[..., np.newaxis]*
          (triangle[0] + triangle[1]))
    L = np.stack([L0, L1, L2], axis=-2)
    h = np.sin(area_i)*(1 + np.sum(triangle_iplus1*triangle_iplus2, axis=-1))
    return np.linalg.solve(L, h)

def project_sphere(sphcoords, zfunc=np.arcsin, scale=180 / np.pi):
    """
    Projects 3d coordinates on the sphere onto a 2d rectangle.

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

    >>> project_sphere(np.eye(3)) # doctest: +NORMALIZE_WHITESPACE
    array([[  0.,   0.],
           [ 90.,   0.],
           [  0.,  90.]])
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
