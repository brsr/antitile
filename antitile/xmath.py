# -*- coding: utf-8 -*-
"""
Extra math and array functions.
"""

import numpy as np
from numpy.linalg import norm
import pandas as pd

def reflect_through_origin(normal):
    """Reflection matrix for reflecting through a plane through the origin
    specified by its normal

    >>> x = np.array([1, 0, 0])
    >>> reflect_through_origin(x)  #doctest: +NORMALIZE_WHITESPACE
    array([[-1.,  0.,  0.],
           [ 0.,  1.,  0.],
           [ 0.,  0.,  1.]])    
    """
    return (np.eye(len(normal)) -
            2 * np.outer(normal, normal) / np.inner(normal, normal))

def complex_to_float2d(arr):
    """Converts a complex array to a multidimensional float array.
    >>> x = np.exp(2j*np.pi*np.linspace(0,1,5)).round()
    >>> complex_to_float2d(x.round())
    array([[ 1.,  0.],
           [ 0.,  1.],
           [-1.,  0.],
           [-0., -1.],
           [ 1., -0.]])
    """
    return arr.view(float).reshape(list(arr.shape) + [-1])

def float2d_to_complex(arr):
    """Converts a multidimensional float array to a complex array.
    >>> y = np.arange(8, dtype=float).reshape((-1, 2))
    >>> float2d_to_complex(y)
    array([[ 0.+1.j],
           [ 2.+3.j],
           [ 4.+5.j],
           [ 6.+7.j]])
    """
    return arr.view(complex)

def line_intersection(a1, a2, b1, b2):
    """Finds the point in the plane that lies at the intersection of the
    line from a1 to a2 and the line from b1 to b2.
    >>> a1 = np.array([0, 0])
    >>> a2 = np.array([1, 1])
    >>> b1 = np.array([1, 0])
    >>> b2 = np.array([0, 1])
    >>> line_intersection(a1, a2, b1, b2)
    array([ 0.5,  0.5])
    """
    a1x, a1y = a1[..., 0], a1[..., 1]
    a2x, a2y = a2[..., 0], a2[..., 1]
    b1x, b1y = b1[..., 0], b1[..., 1]
    b2x, b2y = b2[..., 0], b2[..., 1]

    numx = (a1x*a2y - a1y*a2x)*(b1x - b2x) - (b1x*b2y - b1y*b2x)*(a1x - a2x)
    numy = (a1x*a2y - a1y*a2x)*(b1y - b2y) - (b1x*b2y - b1y*b2x)*(a1y - a2y)
    denom = (a1x-a2x)*(b1y-b2y) - (a1y-a2y)*(b1x-b2x)
    return np.stack([numx, numy], axis=-1)/denom[..., np.newaxis]


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
    """Renumbers an index of a subarray given by the boolean array inarray.
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


def slerp(pt1, pt2, intervals):
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
    angle = central_angle(pt1, pt2)[..., np.newaxis]
    return (np.sin((1 - t)*angle)*pt1 + np.sin((t)*angle)*pt2)/np.sin(angle)


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
    >>> a = np.array([0,0,0])
    >>> b = np.array([1,0,0])
    >>> c = np.array([0,0,1])

    >>> bearing(a,b,c)/np.pi*180 #doctest: +ELLIPSIS
    90...
    """
    direction = origin - pole
    pv = destination - origin
    d = np.sum(pv * direction, axis=-1)
    x = norm(np.cross(pv, direction), axis=-1)
    return np.arctan2(x, d)


def central_angle(x, y, signed=False):
    """Central angle between vectors with respect to 0. If vectors have norm
    1, this is the spherical distance between them.
    Args:
        x, y: Coordinates of points on the sphere.
        axis: Which axis the vectors lie along. By default, -1.
    Returns: Array of central angles.

    >>> t = np.linspace(0,np.pi,5)
    >>> c = np.cos(t)
    >>> s = np.sin(t)
    >>> z = np.zeros(t.shape)
    >>> x = np.stack((c,s,z),axis=-1)
    >>> y = np.stack((c,z,s),axis=-1)
    >>> np.round(central_angle(x,y)/np.pi*180) # doctest: +NORMALIZE_WHITESPACE
    array([  0.,  60.,  90.,  60.,   0.])
    """
    cos = np.sum(x*y, axis=-1)
    sin = norm(np.cross(x, y), axis=-1)
    result = np.arctan2(sin, cos)
    return result if signed else abs(result)


def central_angle_equilateral(pts):
    """For use with the naive slerp methods. Takes the central angle between
    each of the points in pts. If they are close, returns the central angle.
    If not, raises an error.
    >>> x = np.eye(3)
    >>> central_angle_equilateral(x)/np.pi*180 #doctest: +ELLIPSIS
    90...
    >>> y = x[[0,0,2]]
    >>> central_angle_equilateral(y)/np.pi*180 #doctest: +ELLIPSIS
    Traceback (most recent call last):
        ...
    ValueError: naive_slerp used with non-equilateral face. Difference is 1.57... radians.

    """
    omegas = central_angle(pts, np.roll(pts, 1, axis=0))
    max_diff = np.abs(omegas - np.roll(omegas, 1)).max()
    if not np.isclose(max_diff, 0):
        raise ValueError("naive_slerp used with non-equilateral face. " +
                         "Difference is " + str(max_diff) + " radians.")
    return omegas[0]


def triangle_solid_angle(a, b, c):
    """Solid angle of a triangle with respect to 0. If vectors have norm 1,
    this is the spherical area. Note there are two solid angles defined by
    three points: this will always return the smaller of the two. (The other
    is 4*pi minus what this function returns.)

    Formula is from Van Oosterom, A; Strackee, J (1983).
    "The Solid Angle of a Plane Triangle". IEEE Trans. Biom. Eng.
    BME-30 (2): 125–126. doi:10.1109/TBME.1983.325207.

    Args:
        a, b, c: Coordinates of points on the sphere.

    Returns: Array of solid angles.

    >>> t = np.linspace(0,np.pi,5)
    >>> a = np.stack([np.cos(t), np.sin(t), np.zeros(5)],axis=-1)
    >>> b = np.array([0,1,1])/np.sqrt(2)
    >>> c = np.array([0,-1,1])/np.sqrt(2)
    >>> np.round(triangle_solid_angle(a, b, c), 4)
    array([ 1.5708,  1.231 ,  0.    ,  1.231 ,  1.5708])
    """

    top = np.abs(triple_product(a, b, c))
    na = norm(a, axis=-1)
    nb = norm(b, axis=-1)
    nc = norm(c, axis=-1)
    bottom = (na*nb*nc + np.sum(a * b, axis=-1)*nc
              + np.sum(b * c, axis=-1)*na
              + np.sum(c * a, axis=-1)*nb)
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
    >>> spherical_bearing(x,np.roll(x,1))/np.pi*180 #doctest: +ELLIPSIS
    90...
    """
    c_1 = np.cross(origin, destination)
    c_2 = np.cross(origin, pole)
    cos_theta = np.sum(c_1 * c_2, axis=-1)
    sin_theta = triple_product(origin, destination, pole)
    return np.arctan2(sin_theta, cos_theta)
