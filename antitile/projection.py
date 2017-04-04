# -*- coding: utf-8 -*-
"""
Transformations from the breakdowns to 3-space (sometimes 2-space)
"""

import numpy as np
from . import xmath, breakdown

#generic methods, valid on any-d, any shape
def square_to_quad(xy, abcd):
    """Transforms a square in [0,1]^2 to a (possibly skew) quadrilateral
    defined by abcd.

    Args:
        xy: The 2d coordinates of the point. Last element of shape should be 2.
            If not, treated like xy[..., :2]: the values in xy[..., 2:] are
            ignored.
        abcd: The coordinates of the quadrilateral. Should be in
            counterclockwise order to maintain orientation. First element of
            shape should be 4.
    Returns:
        Coordinates in whatever space abcd was defined in."""
    a, b, c, d = abcd[0], abcd[1], abcd[2], abcd[3]
    x, y = xy[..., 0], xy[..., 1]
    return a + (b-a)*x + (d-a)*y + (a-b+c-d)*x*y

def tri_bary(bary, abc):
    """Transforms barycentric coordinates to a triangle defined by abc.

    Args:
        bary: Barycentric coordinates. Last element of shape should be 3.
        abc: Coordinates of the triangle. Should be in counterclockwise order
            to maintain orientation. First element of shape should be 3.

    Returns:
        Coordinates in whatever space abc was define in.
    """
    return bary.dot(abc)
#methods for disks

def bary_tri(tri, vertices):
    """Transforms a triangle back into barycentric
    coordinates. Only elements [..., :2] are used."""
    afill = np.ones(vertices.shape[:-1])
    a = np.concatenate([vertices[..., :2], afill[..., np.newaxis]], axis=-1)
    bfill = np.ones(tri.shape[:-1])
    b = np.concatenate([tri[..., :2], bfill[..., np.newaxis]], axis=-1)
    #couldn't make np.linalg.solve cooperate here.
    #this should be OK numerically
    ainv = np.linalg.inv(a)
    return b.dot(ainv)

def square_to_circle(xy, rotation=1):#np.exp(1j*np.pi/4)):
    """Transforms square on [0,1]^2 to unit circle"""
    pts = 2*xy - 1
    r = np.max(np.abs(pts), axis=-1)
    theta = np.arctan2(pts[..., 1], pts[..., 0])
    result = r*np.exp(1j*theta)*rotation
    return xmath.complex_to_float2d(result)

TRI_C = np.exp(2j*np.pi/3*np.arange(3))*1j
TRI_R = xmath.complex_to_float2d(TRI_C)

def tri_to_circle(beta, rotation=1, pts=TRI_C):
    """Transforms triangle in barycentric coordinates to unit circle.
    tri_naive_slerp also does this when the pts are on a great circle,
    with somewhat different results."""
    tri_pts = beta.dot(pts)
    angle = np.angle(tri_pts)
    r = 1 - 3*beta.min(axis=-1)
    result = r*np.exp(1j*angle)
    return xmath.complex_to_float2d(result)

def _sq_cir(bkdn, abc, freq, tweak):
    #FIXME this is hella broken
    return square_to_quad(square_to_circle(bkdn.coord)[:, np.newaxis], abc)


def _tri_cir(bkdn, abc, freq, tweak):
    rebary = bary_tri(tri_to_circle(bkdn.coord), TRI_R)
    return tri_bary(rebary, abc)

#disk -> sphere

def spherical_to_xyz(phi, theta):
    """Converts spherical coordinates to 3d xyz coordinates"""
    return np.array([np.sin(phi) * np.cos(theta), # pylint: disable=no-member
                     np.sin(phi) * np.sin(theta),
                     np.cos(phi)]).T

def lambert(disk):
    """Converts coordinates on the disk to spherical coordinates, using
    the Lambert azimuthal equal-area projection."""
    theta = np.arctan2(disk[..., 1], disk[..., 0])
    phi = 2*np.arcsin(np.linalg.norm(disk, axis=-1))
    return phi, theta

def equidistant(disk):
    """Converts coordinates on the disk to spherical coordinates, using
    the azimuthal equal-distance projection."""
    theta = np.arctan2(disk[..., 1], disk[..., 0])
    phi = np.linalg.norm(disk, axis=-1)*np.pi
    return phi, theta

#methods for spheres

#triangles -> spherical triangle

def tri_naive_slerp(bary, base_pts):
    omegas = xmath.spherical_distance(base_pts, np.roll(base_pts, 1, axis=0))
    max_diff = np.abs(omegas - np.roll(omegas, 1)).max()
    if not np.isclose(max_diff, 0):
        raise ValueError("naive_slerp used with non-equilateral face. " +
                         "Difference is " + str(max_diff) + " radians.")
    omega = omegas[0]
    b = np.sin(omega * bary) / np.sin(omega)
    return b.dot(base_pts)

def tri_areal(beta, triangle):
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
    triangle = xmath.normalize(triangle)
    area = xmath.spherical_triangle_area(triangle[0], triangle[1], triangle[2])
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

def triangles_method2(lindex, base_pts, freq):
    n, m = freq
    frame = breakdown.frame_triangle(base_pts, n, m)
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
    #result may be zero for base vertices, so fill those in
#    group = vertices.group
#    for i in range(3):
#        result[group == i] = base_pts[i]

    return result

def tri_intersections(lindex, base_pts, freq, tweak=False):
    """Transforms a triangle to a spherical triangle using the method of
    intersections"""
    a, b = freq
    pts = triangles_method2(lindex, base_pts, freq)
    if ~tweak:
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
def square_slerp(xy, abcd):
    """Transforms a square in [0,1]^2 to a spherical quadrilateral
    defined by abcd, using spherical linear interpolation

    Args:
        xy: The 2d coordinates of the point. Last element of shape should be 2.
            If not, treated like xy[..., :2]: the values in xy[..., 2:] are
            ignored.
        abcd: The coordinates of the quadrilateral. Should be in
            counterclockwise order to maintain orientation. Shape should be
            (4, ..., 3)
    Returns:
        Coordinates on the sphere."""
    a, b, c, d = abcd[0], abcd[1], abcd[2], abcd[3]
    x, y = xy[..., 0], xy[..., 1]
    ab = xmath.slerp(a, b, x)
    dc = xmath.slerp(d, c, x)
    result = xmath.slerp(ab, dc, y)
    return result

def square_intersections(lindex, abcd, freq):
    """Transforms a square to a spherical quadrilateral using the method of
    intersections"""
    n, m = freq
    preframe = breakdown.frame_square(n, m)
    frame = square_slerp(preframe[..., np.newaxis, :],
                         abcd[:, np.newaxis, np.newaxis, np.newaxis])
    gc_normals = np.cross(frame[..., 0, :], frame[..., 1, :])
    index = np.arange(2)
    pairs = gc_normals[index, lindex[:, index]]
    #intersection of great circles = cross product of normals
    ptx = np.cross(pairs[:, 0], pairs[:, 1])
    # cross product could give the point we want or its negative.
    # test to see if the points are on the correct side of the sphere
    # take the dot product of these vectors with the center of the
    # base face. if it's positive, it's right, if not, negate it
    center = np.sum(abcd, axis=0)#don't need to normalize this
    sign_correct = np.sum(center*ptx, axis=-1, keepdims=True) >= 0
    result = np.where(sign_correct, ptx, -ptx)
    result[(lindex[:, 0] == b)     & (lindex[:, 1] == 0)] = abcd[0]
    result[(lindex[:, 0] == a + b) & (lindex[:, 1] == b)] = abcd[1]
    result[(lindex[:, 0] == a)     & (lindex[:, 1] == a + b)] = abcd[2]
    result[(lindex[:, 0] == 0)     & (lindex[:, 1] == a)] = abcd[3]

    return result

#extra crap
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
    area = xmath.spherical_triangle_area(triangle[0], triangle[1], triangle[2])
    area_i = xmath.spherical_triangle_area(pts[:, np.newaxis],
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


FLAT = {3: lambda bkdn, abc, freq, tweak: tri_bary(bkdn.coord, abc),
        4: lambda bkdn, abc, freq, tweak:
            square_to_quad(bkdn.coord[:, np.newaxis], abc)}

SLERP = {3: lambda bkdn, abc, freq, tweak: tri_naive_slerp(bkdn.coord, abc),
         4: lambda bkdn, abc, freq, tweak:
             square_slerp(bkdn.coord[:, np.newaxis], abc)}

AREAL = {3: lambda bkdn, abc, freq, tweak: tri_areal(bkdn.coord, abc)}

GC = {3: lambda bkdn, abc, freq, tweak:
            tri_intersections(bkdn.lindex, abc, freq, tweak),
      4: lambda bkdn, abc, freq, tweak:
          square_intersections(bkdn.lindex, abc, freq)}

DISK = {3: _tri_cir,
        4: _sq_cir}

PROJECTIONS = {'flat':  FLAT,
               'slerp': SLERP,
               'areal': AREAL,
               'gc':    GC,
               'disk':  DISK}

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection, LineCollection

    a, b = 2, 1
    bkdn = breakdown.Breakdown(a, b, 'q')
    abcd = xmath.normalize(np.array([[1, 1, 1],
                                     [1, -1, 1],
                                     [-1, -1, 1],
                                     [-1, 1, 1]]))
    borders = xmath.slerp(abcd[:, np.newaxis],
                          np.roll(abcd, -1, axis=0)[:, np.newaxis],
                          np.linspace(0, 1)[:, np.newaxis]).reshape((-1, 3))

    sq1 = xmath.normalize(square_to_quad(bkdn.coord[:, np.newaxis], abcd))
    sqslerp = xmath.normalize(square_slerp(bkdn.coord[:, np.newaxis], abcd))
    sq2 = xmath.normalize(square_intersections(bkdn.lindex, abcd, (a, b)))
    badsq2 = np.linalg.norm(sq2, axis=-1) < 1E-6
    sq2[badsq2] = sq1[badsq2]
    #x = np.stack([sqslerp[:,0],sq2[:,0]],axis=-1)
    #y = np.stack([sqslerp[:,1],sq2[:,1]],axis=-1)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharex=True, sharey=True)
    fig.set_size_inches(8, 2)
    plt.axis('equal')
    #ax1.plot(x, y, c='g')
    ax1.plot(borders[..., 0], borders[..., 1], c='b')
    for ax in (ax1, ax2, ax3, ax4):
        circle = plt.Circle((0, 0), 1, color='k', fill=False, zorder=3)
        ax.add_artist(circle)
    ax1.scatter(sq1[..., 0], sq1[..., 1])
    ax1.scatter(sqslerp[..., 0], sqslerp[..., 1])
    ax1.scatter(sq2[..., 0], sq2[..., 1])
    for ax, pts in [(ax2, sq1), (ax3, sqslerp), (ax4, sq2)]:
        ptx = pts[bkdn.faces][..., :2]
        pc = PolyCollection(ptx, edgecolors='grey')
        ax.add_collection(pc)

    bkdn = breakdown.Breakdown(a, b, 't')
    abc = xmath.normalize(np.array([[0, 1, 1],
                                    [ np.sqrt(3)/2, -0.5, 1],
                                    [-np.sqrt(3)/2, -0.5, 1]]))
    borders = xmath.slerp(abc[:, np.newaxis],
                          np.roll(abc, -1, axis=0)[:, np.newaxis],
                          np.linspace(0, 1)[:, np.newaxis]).reshape((-1, 3))

    tr1 = xmath.normalize(tri_bary(bkdn.coord, abc))
    trslerp = xmath.normalize(tri_naive_slerp(bkdn.coord, abc))
    tr2 = xmath.normalize(tri_intersections(bkdn.lindex, abc, (a, b)))
    badtr2 = np.linalg.norm(tr2, axis=-1) < 1E-6
    tr2[badtr2] = tr1[badtr2]
    #x = np.stack([sqslerp[:,0],sq2[:,0]],axis=-1)
    #y = np.stack([sqslerp[:,1],sq2[:,1]],axis=-1)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharex=True, sharey=True)
    fig.set_size_inches(8, 2)
    plt.axis('equal')
    #ax1.plot(x, y, c='g')
    ax1.plot(borders[..., 0], borders[..., 1], c='b')
    for ax in (ax1, ax2, ax3, ax4):
        circle = plt.Circle((0, 0), 1, color='k', fill=False, zorder=3)
        ax.add_artist(circle)
    ax1.scatter(tr1[..., 0], tr1[..., 1])
    ax1.scatter(trslerp[..., 0], trslerp[..., 1])
    ax1.scatter(tr2[..., 0], tr2[..., 1])
    for ax, pts in [(ax2, tr1), (ax3, trslerp), (ax4, tr2)]:
        ptx = pts[bkdn.faces][..., :2]
        pc = PolyCollection(ptx, edgecolors='grey')
        ax.add_collection(pc)
