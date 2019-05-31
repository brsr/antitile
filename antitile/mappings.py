# -*- coding: utf-8 -*-
"""
Mappings between various spaces.

All functions are named X2Y_blahblah(), indicating that the function
maps from the space X to Y. The function may be valid for inputs outside of X
(which may produce outputs outside of Y). X and Y are as in this list:

c: Complex plane. Array of complex type with shape (...), viewable as (n, 2).
q: Standard square, complex. As c,
    with real and imaginary parts in the interval [-1,1].
t: Standard triangle, complex. As c, inside the triangle defined by the points
    np.exp(2j*np.pi/3*np.arange(3)).
d: Unit disk, complex. As c, with abs(d) <= 1.
h: Upper half-plane, complex. As c, with h.imag >= 0.
r: 3-space. Array of float type with shape (..., 3).
b: Barycentric coordinates. As r, but with the restriction that
    np.sum(b, axis=-1) == 1 for all entries.
s: Unit sphere. As r, but with the added restriction that
    s[..., 0]**2 + s[..., 1]**2 + s[..., 2]**2 = 1.

Some functions also have this argument:
base: Float array, with shape (3,3) if triangular, (4,3) if quadrilateral.
If not, the base polygon is assumed to be the standard square or triangle.

Other arguments are defined as in the functions where they are used.

References:
http://marc-b-reynolds.github.io/math/2017/01/08/SquareDisc.html
Snyder, John P. (1987). Map Projections − A Working Manual,
    Professional Paper 1395. US Geological Survey.
"""

import numpy as np
from numpy.linalg import norm
from scipy.special import hyp2f1, gamma, ellipj, ellipk, ellipkinc
from scipy.optimize import minimize
from . import xmath

#from . import xmath, breakdown

#mappings that can be higher-dimensional if needed

def b2r(bary, base):
    """Transforms barycentric coordinates to a (Euclidean) triangle
    (or simplex) defined by base.
    """
    return bary.dot(base)

def r2b(pts, base):
    """Transforms a triangle (or simplex) back into barycentric coordinates.
    Add rows of 1s if v and triangle have different final dimensions.
    """
    #couldn't make np.linalg.lstsq cooperate here.
    #this should be OK numerically
    inv = np.linalg.pinv(base)
    return pts.dot(inv)

def q2r(z, base):
    """Transforms the standard square to a (possibly skew) quadrilateral
    defined by base.
    """
    a, b, c, d = base[0], base[1], base[2], base[3]
    x = z.real[:, np.newaxis]
    y = z.imag[:, np.newaxis]
    return ((1-x)*(1-y)*a + (1+x)*(1-y)*b + (1+x)*(1+y)*c + (1-x)*(1+y)*d)/4

def r2s_normalize(pts, axis=-1, eps=0):
    """Normalizes vectors, mapping n-space to the n-1-sphere,
    except for the zero vector, which remains the zero vector.
    Args:
        pts: Array of vectors
        axis: Which axis to take the norm over (by default the last axis, -1)
        eps: If the norm of the vector is below or equal to eps, it is
            considered to be zero. By default, zero.
    """
    n = norm(pts, axis=axis, keepdims=True)
    return np.where(n <= eps, 0, pts / n)

def r2s_parallel_exact(pts, normal):
    """Projects points exactly onto the unit sphere parallel to the normal vector.
    Args:
        pts: Points to project
        normal: Normal vector
    """
    vdotc = np.sum(pts * normal, axis=-1)
    vdotv = norm(pts, axis=-1)**2
    p = -vdotc + xmath.sqrt(1 + vdotc**2 - vdotv)
    return p[..., np.newaxis] * normal


def r2s_parallel_approx(pts, normal):
    """Approximately projects points onto the unit sphere parallel to the
        normal vector.
    Args:
        pts: Points to project
        normal: Normal vector
        """
    q = 1 - norm(pts, axis=-1)
    return q[..., np.newaxis] * normal

def r2s_parallel(pts, normal, exact=True):
    """Projects points onto the unit sphere parallel to the normal vector.
    Args:
        pts: Points to project
        normal: Normal vector
        exact: Whether to project exactly or approximately.
            Defaults to exact (True).
    """
    if exact:
        return r2s_parallel_exact(pts, normal)
    else:
        return r2s_parallel_approx(pts, normal)

#Mappings from the sphere to the plane (map projections). Except for gnomonic,
#these all are scaled to map the equator to the unit circle.

def c2s_stereographic(z):
    """
    Stereographic projection from the plane to the sphere.
    """
    x = z.real
    y = z.imag
    u = 2*x
    v = 2*y
    w = 1 - x**2 - y**2
    return np.stack([u,v,w], axis=-1)/(1+x**2+y**2)[..., np.newaxis]

def s2c_stereographic(sph):
    """
    Stereographic projection from the sphere to the plane.
    """
    u = sph[..., 0]
    v = sph[..., 1]
    w = sph[..., 2]
    return (u + 1j*v)/(1+w)

def c2s_gnomonic(z):
    """Converts coordinates on the plane to spherical coordinates, using
    the gnomonic projection.
    """
    u = z.real
    v = z.imag
    w = np.ones(z.shape)
    res = np.stack([u,v,w], axis=-1)
    return r2s_normalize(res)

def s2c_gnomonic(sph):
    """Converts coordinates on the unit sphere to the plane, using
    the gnomonic projection.
    """
    u = sph[..., 0]
    v = sph[..., 1]
    w = sph[..., 2]
    return (u + 1j*v)/w

def d2s_orthographic(z):
    """Converts coordinates on the disk to spherical coordinates, using
    the orthographic projection.
    """
    u = z.real
    v = z.imag
    w = xmath.sqrt(1-abs(z)**2)
    return np.stack([u,v,w], axis=-1)

def s2d_orthographic(sph):
    """Converts coordinates on the unit sphere to the disk, using
    the orthographic projection.
    """
    return sph[..., 0] + 1j*sph[..., 1]

def d2s_ea(z, scale=1):
    """Converts coordinates on the disk to spherical coordinates, using
    the Lambert azimuthal equal-area projection. By default (scale=1), the unit
    circle is mapped to the equator.
    """
    az = scale*abs(z)**2
    factor = xmath.sqrt(1-az/2)
    u = factor*z.real
    v = factor*z.imag
    w = 1 - az
    return np.stack([u,v,w], axis=-1)

def s2d_ea(sph, scale=1):
    """Converts coordinates on the unit sphere to the disk, using
    the Lambert azimuthal equal-area projection. By default (scale=1), the
    equator is mapped to the unit circle.
    """
    u = sph[..., 0]
    v = sph[..., 1]
    w = sph[..., 2]
    return scale*(u + 1j*v)/xmath.sqrt(1+w)

def d2s_equidistant(z, scale=1/2):
    """Converts coordinates on the disk to spherical coordinates, using
    the azimuthal equal-distance projection.
    """
    z = z*scale
    az = abs(z)
    factor = np.where(az>0, np.sin(az*np.pi)/az,np.pi)
    u = factor*z.real
    v = factor*z.imag
    w = np.cos(az*np.pi)
    return np.stack([u,v,w], axis=-1)

def s2d_equidistant(sph, scale=2):
    """Converts coordinates on the unit sphere to the disk, using
    the azimuthal equal-distance projection.
    """
    u = sph[..., 0]
    v = sph[..., 1]
    w = sph[..., 2]
    return scale*(u + 1j*v)/(np.pi*xmath.sqrt(1-w**2))*np.arccos(w)

#Disk - polygon mappings
#Conformal
K = ellipk(1/2)

def d2q_conformal(z):
    """Conformal map from the unit disk to the standard square"""
    diskp = z*(1+1j)/np.sqrt(2)
    up = diskp.real
    vp = diskp.imag
    A = up**2 + vp**2
    B = up**2 - vp**2
    T = np.sqrt((1+A**2)**2 - 4*B**2)
    U = 1 + 2*B - A**2
    alpha = np.arccos(np.clip((2*A-T)/U, -1, 1))
    beta = np.arccos(np.clip(U/(2*A+T), -1, 1))
    xp = np.sign(up)*(2*K - ellipkinc(alpha, 1/2))
    yp = np.sign(vp)*ellipkinc(beta, 1/2)
    return (xp + 1j*yp)*(1-1j)/2/K

def _complexcn(z):
    """Jacobi elliptical function cn with parameter 1/2,
    valid for complex numbers."""
    x = z.real
    y = z.imag
    snx, cnx, dnx, phx = ellipj(x, 1/2)
    sny, cny, dny, phy = ellipj(y, 1/2)
    numer = cnx * cny - 1j * snx * dnx * sny * dny
    denom = (1 - dnx**2 * sny**2)
    return numer / denom

def q2d_conformal(z):
    """Conformal map from the standard square to the unit disk.
    """
    csp = z * (1+1j)/2
    w = _complexcn(K*(1-csp))
    return w * (1 - 1j)/np.sqrt(2)

def d2t_conformal(z):
    """Conformal map from the unit disk to the standard triangle. Note that
    there is a bug in scipy's hyp2f1 (#8054 on github) that may cause incorrect
    results in a small area of the disk."""
    hp1 = (z+1)/(z-1)
    result = -3*gamma(5/6)/gamma(1/3)/np.sqrt(np.pi)*hp1*hyp2f1(1/2, 2/3, 3/2, -3*hp1**2)
    #sometimes values on a branch cut wind up with the wrong sign or undefined, fix those
    result.real = abs(result.real) - 1/2
    result[z == 1] = 1
    #result[abs(result) > 1] = np.nan # to lessen the effect of the SciPy bug
    return result

def t2d_conformal(z):
    """Conformal map from the standard triangle to the unit disk. Note that this
    function numerically inverts d2t_conformal() and may be slow.
    """
    conf_tridisk = z.copy()
    for i in range(len(z)):
        zi = z[i]
        xy = np.array([zi.real, zi.imag])
        def objective(t):
            return abs(d2t_conformal(t.view(dtype=np.complex128)) - zi)
        res = minimize(objective, x0=xy, method='Nelder-Mead')
        conf_tridisk[i] = res.x.view(dtype=np.complex128)
    return conf_tridisk

def b2d_conformal(bary, base=np.exp(2j/3*np.arange(3))):
    """Wrapper for t2d_conformal to put it in terms of barycentric coordinates.
    """
    return t2d_conformal(b2r(bary,base))

#Radial stretch
def q2d_radial(z):
    """Radial stretch, square to disk"""
    sfactor = np.maximum(abs(z.real),abs(z.imag))/abs(z)
    sfactor[np.isnan(sfactor)] = 1
    return sfactor*z

def d2q_radial(z):
    """Radial stretch, disk to square"""
    dfactor = abs(z)/np.maximum(abs(z.real),abs(z.imag))
    dfactor[np.isnan(dfactor)] = 1
    return dfactor*z

def b2d_radial(bary):
    """Radial stretch, triangle (barycentric) to disk"""
    base = np.exp(2j/3*np.pi*np.arange(3))
    tri = b2r(bary, base)
    c = np.angle(tri)
    r = 1 - 3*np.min(bary, axis=-1)
    return r * np.exp(1j*c)

def d2t_radial(z):
    """Radial stretch, disk to triangle"""
    c = np.angle(z)
    r = abs(z)/(2*np.cos(c % (np.pi*2/3) - np.pi/3))
    return r * np.exp(1j*c)

#Equal-area
def q2d_equalarea(z):
    """Equal-area, square to disk"""
    x = z.real
    y = z.imag
    ea_sqdisk = np.where(abs(x) >= abs(y),
                     x*np.exp(1j*np.pi/4*y/x),
                     1j*y*np.exp(-1j*np.pi/4*x/y))
    ea_sqdisk[np.isnan(ea_sqdisk)] = 0
    return ea_sqdisk

def d2q_equalarea(z):
    """Equal-area, disk to square"""
    u = z.real
    v = z.imag
    ea_square = np.where(abs(u) >= abs(v),
                np.sign(u)*abs(z)*(1 + 4j/np.pi*np.arctan(v/u)),
                np.sign(v)*abs(z)*(4/np.pi*np.arctan(u/v) + 1j))
    ea_square[np.isnan(ea_square)] = 0
    return ea_square

def b2d_equalarea(bary):
    """Equal area, triangle (barycentric) to disk"""
    a = bary[..., 0]
    b = bary[..., 1]
    c = bary[..., 2]
    eatri0 = (3*a-1)*np.exp(1j*np.pi/3 * (b - c)/(3*a-1))
    eatri1 = (3*b-1)*np.exp(1j*np.pi/3 * ((c - a)/(3*b-1) + 2) )
    eatri2 = (3*c-1)*np.exp(1j*np.pi/3 * ((a - b)/(3*c-1) - 2) )

    mb = np.min(bary, axis=0)

    ea_tridisk = eatri0.copy()
    ea_tridisk[bary[1] == mb] = eatri1[bary[1] == mb]
    ea_tridisk[bary[2] == mb] = eatri2[bary[2] == mb]
    ea_tridisk[np.isnan(ea_tridisk)] = 0
    return ea_tridisk

def d2b_equalarea(z, eps=0):
    """Equal area, disk to triangle (barycentric)"""
    ea_beta = np.zeros((3, z.size))
    for i in range(3):
        this_disk = z.flatten() * np.exp(i*np.pi*2j/3)
        mask = abs(np.angle(this_disk)) >= 2/3*np.pi
        uz = this_disk.real
        vz = this_disk.imag
        absz = abs(this_disk)
        beta0 = (1 - absz)/3
        b = 6/np.pi*absz*np.arctan(vz/(absz-uz))
        beta1 = b/2 - beta0/2 + 1/2
        beta2 = -b/2 - beta0/2 + 1/2
        ea_b = np.roll(np.vstack([beta0, beta1, beta2]),-i,axis=0)
        ea_beta[:,mask] = ea_b[:,mask]
    ea_beta[:,abs(z.flatten()) <=eps] = 1/3
    return ea_beta.reshape([3] + list(z.shape))

#Approximate equal-area
def q2d_aea(z):
    """Approximate equal-area, square to disk"""
    x = z.real
    y = z.imag
    aea_sqdisk = np.where(abs(x) >= abs(y),
                  x*np.sqrt(1-y**2/2/x**2) + 1j*y/np.sqrt(2),
                  x/np.sqrt(2) + 1j*y*np.sqrt(1-x**2/2/y**2) )
    aea_sqdisk[np.isnan(aea_sqdisk)] = 0
    return aea_sqdisk

def d2q_aea(z):
    """Approximate equal-area, disk to square"""
    u = z.real
    v = z.imag
    aea_square = np.where(abs(u) >= abs(v),
                np.sign(u)*abs(z) + 1j*np.sqrt(2)*v,
                np.sqrt(2)*u + 1j*np.sign(v)*abs(z))
    aea_square[np.isnan(aea_square)] = 0
    return aea_square

def b2d_aea(bary):
    """Approximate equal-area, triangle (barycentric) to disk"""
    a = bary[..., 0]
    b = bary[..., 1]
    c = bary[..., 2]
    aeatri0 = (3*a-1)*np.sqrt(1- ((b - c)/(3*a-1))**2*3/4) + 1j*np.sqrt(3)/2*(b - c)
    aeatri1 = (3*b-1)*np.sqrt(1- ((c - a)/(3*b-1))**2*3/4) + 1j*np.sqrt(3)/2*(c - a)
    aeatri2 = (3*c-1)*np.sqrt(1- ((a - b)/(3*c-1))**2*3/4) + 1j*np.sqrt(3)/2*(a - b)

    mb = np.min(bary, axis=0)

    aea_tridisk = aeatri0.copy()
    aea_tridisk[bary[1] == mb] = aeatri1[bary[1] == mb]*np.exp(2j*np.pi/3)
    aea_tridisk[bary[2] == mb] = aeatri2[bary[2] == mb]*np.exp(-2j*np.pi/3)
    aea_tridisk[np.isnan(aea_tridisk)] = 0
    return aea_tridisk

def d2b_aea(z):
    """Approximate equal-area, disk to triangle (barycentric)"""
    aea_beta = np.zeros((3, z.size))
    for i in range(3):
        this_disk = (z * np.exp(i*np.pi*2j/3)).flatten()
        mask = abs(np.angle(this_disk)) >= 2/3*np.pi
        #uz = this_disk.real
        vz = this_disk.imag
        absz = abs(this_disk)
        beta0 = (1 - absz)/3
        beta1 = (2 + 2*np.sqrt(3)*vz + absz)/6
        beta2 = (2 - 2*np.sqrt(3)*vz + absz)/6
        aea_b = np.roll(np.vstack([beta0, beta1, beta2]),-i,axis=0)
        aea_beta[:,mask] = aea_b[:,mask]
    return aea_beta.reshape([3] + list(z.shape))

#*ircle
def q2d_squircle(z):
    """Squircle, square to disk"""
    x = z.real
    y = z.imag
    sfactor = np.sqrt(x**2+y**2-x**2*y**2)/abs(z)
    sfactor[np.isnan(sfactor)] = 1
    return sfactor*z

def d2q_squircle(z, eps=0):
    """Squircle, disk to square"""
    u = z.real
    v = z.imag
    n = u**2 + v**2
    s = np.sign(u*v)/np.sqrt(2) * np.sqrt(n - np.sqrt(n*(n-4*u**2*v**2)))
    sr_square = s * (1/v + 1j/u)
    index = (abs(z.real) <= 0) | (abs(z.imag) <= 0)
    sr_square[index] = z[index]
    return sr_square

def t2d_squircle(z):
    """Squircle (trircle?), triangle (barycentric) to disk"""
    x = z.real
    y = z.imag
    sfactor = (6*x*y**2 -2*x**3 + 3*y**2 + 3*x**2) / (x**2 + y**2)
    sr_tridisk = np.sqrt(sfactor)*z
    return sr_tridisk

def d2t_squircle(z, eps=1E-6):
    """Squircle (trircle?), disk to triangle"""
    u = z.real
    v = z.imag
    us = u**2
    vs = v**2
    k = (us - 3*vs)/(us + vs)
    rots = np.where(u*k < 0, np.pi*2/3, 0)
    x = 1/k*(np.sin(1/3*np.arcsin(np.clip(2*us*k**2 - 1, -1, 1)) - rots) + 1/2)
    sr_tri = x *(1 + 1j*v/u)
    index = (abs(u) < eps) | (abs(k) < eps)
    sr_tri[index] = z[index]/np.sqrt(3)
    return sr_tri

def b2d_squircle(bary, base=np.exp(2j/3*np.arange(3))):
    """Wrapper for t2d_squircle(z) to put it in terms of barycentric coordinates.
    """
    return t2d_squircle(b2r(bary,base))

#elliptic
def q2d_el(z):
    """Elliptic, square to disk"""
    x = z.real
    y = z.imag
    return x*np.sqrt(1-y**2/2) + 1j*y*np.sqrt(1-x**2/2)

def d2q_el(z):
    """Elliptic, disk to square"""
    u = z.real
    v = z.imag
    t = u**2 - v**2
    return ((np.sqrt(2+t+2*np.sqrt(2)*u) - np.sqrt(2+t-2*np.sqrt(2)*u))
            +1j*(np.sqrt(2-t+2*np.sqrt(2)*v) -np.sqrt(2-t-2*np.sqrt(2)*v)))/2

#naive slerp
def q2d_naive_slerp(z):
    """Naive Slerp, Square to disk"""
    x = z.real
    y = z.imag
    return np.sqrt(2)*(np.sin(np.pi/4*x)*np.cos(np.pi/4*y) +
                        1j*np.cos(np.pi/4*x)*np.sin(np.pi/4*y))

def d2q_naive_slerp(z):
    """Naive Slerp, disk to square"""
    us = z.real**2
    vs = z.imag**2
    m = vs - us
    k = sqrt_safe((m+2)**2 - 8 * vs)
    return (np.sign(z.real)*np.arccos(sqrt_safe(2+m+k)/2)
            + 1j*np.sign(z.imag)*np.arccos(sqrt_safe(2-m+k)/2) )*4/np.pi

def b2d_naive_slerp(bary):
    """Naive slerp, triangle (barycentric) to disk"""
    a = bary[..., 0]
    b = bary[..., 1]
    c = bary[..., 2]
    #or tri_pts @ np.sin(2*np.pi/3*bary) * 2 / np.sqrt(3)
    return ( (2*np.sin(np.pi*2/3*a)-np.sin(np.pi*2/3*b)-
                np.sin(np.pi*2/3*c))/np.sqrt(3)
               + 1j*(np.sin(np.pi*2/3*b)-np.sin(np.pi*2/3*c)) )

def q2d_naive_slerp_2(z):
    """Naive slerp 2, square to disk"""
    x = z.real
    y = z.imag
    return np.sqrt(2)*(
        +np.sqrt(2+np.sqrt(2))*np.sin(np.pi/8*x)*np.cos(np.pi/8*y)*
            np.cos(np.pi/8*x*y)
        -np.sqrt(2-np.sqrt(2))*np.cos(np.pi/8*x)*np.sin(np.pi/8*y)*
            np.sin(np.pi/8*x*y)
        +1j*np.sqrt(2+np.sqrt(2))*np.cos(np.pi/8*x)*np.sin(np.pi/8*y)*
            np.cos(np.pi/8*x*y)
        -1j*np.sqrt(2-np.sqrt(2))*np.sin(np.pi/8*x)*np.cos(np.pi/8*y)*
            np.sin(np.pi/8*x*y))

#Polygons to spherical polygons
def _b2s_fix_corners(bary, base, result):
    result[np.all(bary == np.array([1,0,0]), axis=-1)] = base[0]
    result[np.all(bary == np.array([0,1,0]), axis=-1)] = base[1]
    result[np.all(bary == np.array([0,0,1]), axis=-1)] = base[2]
    return result

def _q2s_fix_corners(z, base, result):
    result[z == -1-1j] = base[0]
    result[z == +1-1j] = base[1]
    result[z == +1+1j] = base[2]
    result[z == -1+1j] = base[3]
    return result

#gnomonic
def b2s_gnomonic(bary, base):
    """
    Gnomonic map from barycentric coordinates to a spherical triangle.
    """
    return r2s_normalize(b2r(bary, base))

def s2b_gnomonic(sph, base):
    """
    Gnomonic map from a spherical triangle to barycentric coordinates.
    """
    normal = np.sum(np.cross(base, np.roll(base, 1, axis=0)), axis=0)
    r = np.mean(base.dot(normal), axis=0)
    v = sph * r /(normal.dot(sph))
    return r2b(v, base)

def q2s_gnomonic(z, base):
    """
    Gnomonic map from standard square to a spherical triangle.
    """
    return r2s_normalize(q2r(z[:,np.newaxis], base))

#areal
def b2s_areal(bary, base):
    """Given a triangle and spherical areal coordinates, returns the vectors
    cooresponding to those coordinates.
    """
    base = xmath.normalize(base)
    area = xmath.triangle_solid_angle(base[0], base[1], base[2])
    area_i = bary * area
    base_iplus1 = np.roll(base, -1, axis=0)
    base_iplus2 = np.roll(base, 1, axis=0)
    L0 = ((1 + np.cos(area_i[..., 0]))[..., np.newaxis]*
          np.cross(base[1], base[2]) -
          np.sin(area_i[..., 0])[..., np.newaxis]*
          (base[1] + base[2]))
    L1 = ((1 + np.cos(area_i[..., 1]))[..., np.newaxis]*
          np.cross(base[2], base[0]) -
          np.sin(area_i[..., 1])[..., np.newaxis]*
          (base[2] + base[0]))
    L2 = ((1 + np.cos(area_i[..., 2]))[..., np.newaxis]*
          np.cross(base[0], base[1]) -
          np.sin(area_i[..., 2])[..., np.newaxis]*
          (base[0] + base[1]))
    L = np.stack([L0, L1, L2], axis=-2)
    detL = np.linalg.det(L)
    #to avoid errors when the matrix is singular
    #only seems to happen on the corners when area = 2pi
    L[detL == 0] = np.nan
    h = np.sin(area_i)*(1 + np.sum(base_iplus1*base_iplus2, axis=-1))
    result = np.linalg.solve(L, h)
    return _b2s_fix_corners(bary, base, result)

def s2b_areal(pts, triangle):
    """Given a triangle and pts within that triangle, returns the
    spherical areal coordinates of the pts with respect to the triangle.
    """
    area = xmath.triangle_solid_angle(triangle[0], triangle[1], triangle[2])
    area_i = xmath.triangle_solid_angle(pts[:, np.newaxis],
                                        np.roll(triangle, 1, axis=0),
                                        np.roll(triangle, -1, axis=0))
    return area_i/area

#conformal
def c2c_mobius_01inf(z, z0=0, z1=1, zinf=1j ):
    """Mobius transformation defined by mapping 3 points to 0, 1, infinity"""
    if ~np.isfinite(zinf):
        return (z-z0)/(z1-z0)
    elif ~np.isfinite(z1):
        return (z-z0)/(z-zinf)
    elif ~np.isfinite(z0):
        return (z1-zinf)/(z-zinf)
    else:
        return (z-z0)*(z1-zinf)/((z-zinf)*(z1-z0))


def c2c_mobius_finite(z,zi,wi):
    """Mobius transformation defined by mapping the points in zi to the points
    in wi."""
    ones = np.ones(zi.shape)
    a = np.linalg.det(np.stack([zi*wi,wi,ones]))
    b = np.linalg.det(np.stack([zi*wi,zi,wi]))
    c = np.linalg.det(np.stack([zi,wi,ones]))
    d = np.linalg.det(np.stack([zi*wi,zi,ones]))
    return (a*z+b)/(c*z+d)

def schwarz_fp(alpha, beta, gam):
    """Parameters of the Schwarz triangle map.
    Args:
        alpha, beta, gamma: Equal to pi times an angle of the triangle.
    Returns:
        s1: Value of the Schwarz triangle map at z=1.
        sinf: Value of the Schwarz triangle map at z=infinity.
        scale: Scale factor for spherical triangles. Will be zero or undefined
        if alpha + beta + gamma <= 1.
    """
    a = (1 - alpha - beta - gam)/2
    b = (1 - alpha + beta - gam)/2
    c = 1 - alpha
    palpha = np.pi*alpha
    pbeta = np.pi*beta
    pgam = np.pi*gam
    gfact = gamma(2-c)/(gamma(1-a)*gamma(c))
    s1 = gamma(c-a)*gamma(c-b)/gamma(1-b)*gfact
    sinf = np.exp(1j*palpha)*gamma(b)*gamma(c-a)*gfact/gamma(b-c+1)
    scale = np.sqrt(abs((np.cos(palpha+pbeta)+np.cos(pgam))/
                 (np.cos(palpha-pbeta)+np.cos(pgam))))
    return s1, sinf, scale

def h2c_schwarz(alpha, beta, gam, z):
    """Schwarz triangle map, from the upper half plane to a triangle with
    circular arcs for edges. These points on the half-plane correspond to these
    angles in the triangle:
    0: alpha*pi
    1: beta*pi
    inf: gam*pi"""
    a = (1 - alpha - beta - gam)/2
    b = (1 - alpha + beta - gam)/2
    c = 1 - alpha
    ap = (1 + alpha - beta - gam)/2#a - c + 1
    bp = (1 + alpha + beta - gam)/2#b - c + 1
    cp = 1 + alpha#2-c
    if alpha + beta + gam == 1:
        sch_res = z**(1-c)*hyp2f1(ap,bp,cp,z)
    else:
        sch_res = z**(1-c)*hyp2f1(ap,bp,cp,z)/hyp2f1(a,b,c,z)
    #fix signs near branch cuts
    sch_res.imag = np.where(z.imag >= 0, abs(sch_res.imag), -abs(sch_res.imag))
    return sch_res

def b2s_conformal(bary, base):
    """Conformal transformation from the standard triangle to a spherical
    triangle."""
    bends = xmath.spherical_bearing(base, np.roll(base,2,axis=0),
                                     np.roll(base,1,axis=0))/np.pi
    alpha, gam, beta = bends
    bendsum = np.sum(bends, axis=0, keepdims=True)
    euc_bends = bends/bendsum
    alphae, game, betae = euc_bends
    s1, sinf, _ = schwarz_fp(alphae, betae, game)
    z = b2r(bary, np.array([0,s1,sinf]))
    initial = c2c_mobius_01inf(z, z1=s1, zinf=sinf)
    conf_h = np.empty_like(z)
    for i in range(len(z)):
        zi = z[i]
        def objective(t):
            result = h2c_schwarz(alphae, betae, game, t.view(dtype=np.complex128))
            return abs(result - zi)
        res = minimize(objective, x0=[initial[i].real, initial[i].imag], method='Nelder-Mead')
        conf_h[i] = res.x.view(dtype=np.complex128)
    res_riemann_sph = h2c_schwarz(alpha, beta, gam, conf_h)
    s1, sinf, _ = schwarz_fp(alpha, beta, gam)
    base_riemann_sph = s2c_stereographic(base)
    zi = np.array([0,s1,sinf])
    moved_riemann_sph = c2c_mobius_finite(res_riemann_sph, zi=zi, wi=base_riemann_sph)
    return _b2s_fix_corners(bary, base, c2s_stereographic(moved_riemann_sph))

def q2s_conformal(z, base, eps=1E-6):
    """Conformal transformation from the standard square to a spherical
    quadrilateral. This is only conformal over the entire square if base
    defines a rhombus or square."""
    xy_base = np.array([[-1,1,1,-1],[-1,-1,1,1],np.ones(4)]).T
    results = []
    angles = []
    for i in range(len(base)):
        sph_tri_base = np.roll(base, i, axis=0)[0:3]
        angles.append(xmath.triangle_solid_angle(sph_tri_base[0],
                                           sph_tri_base[1], sph_tri_base[2]))
        tri = np.roll(xy_base, i, axis=0)[0:3]
        v = np.stack([z.real, z.imag, np.ones(z.shape)], axis=-1)
        bary_rt_tri = r2b(v, tri)
        index1 = np.all(np.round(bary_rt_tri,1) > -eps, axis=-1)
        pts_conf = b2s_conformal(bary_rt_tri, sph_tri_base)
        pts_conf[~index1] = np.nan
        results.append(pts_conf)
    i1 = abs(angles[0] - angles[2]) < eps
    i2 = abs(angles[1] - angles[3]) < eps
    if i1 and not i2:
        results = [results[0], results[2]]
    elif not i1 and i2:
        results = [results[1], results[3]]
    elif not i1 and not i2:
        print('#not an appropriate quadrilateral for q2s_conformal, good luck')
    return np.nanmedian(np.stack(results, axis=0), axis=0)

#great circle
def _b2s_greatcircle_points(bary, base):
    """Triangles of the great circle mapping.    """
    a0 = xmath.slerp(base[2], base[0], bary[:,0,np.newaxis])
    b0 = xmath.slerp(base[1], base[0], bary[:,0,np.newaxis])
    h0 = np.cross(a0, b0)
    a1 = xmath.slerp(base[0], base[1], bary[:,1,np.newaxis])
    b1 = xmath.slerp(base[2], base[1], bary[:,1,np.newaxis])
    h1 = np.cross(a1, b1)
    a2 = xmath.slerp(base[1], base[2], bary[:,2,np.newaxis])
    b2 = xmath.slerp(base[0], base[2], bary[:,2,np.newaxis])
    h2 = np.cross(a2, b2)
    pts1 = np.cross(h0, h1)
    pts2 = np.cross(h1, h2)
    pts3 = np.cross(h2, h0)
    return np.stack([pts1, pts2, pts3], axis=1)

def b2r_greatcircle(bary, base, tweak=False):
    """Transforms a triangle to a spherical triangle using the great circle
    method (Method 2).
    """
    pts = _b2s_greatcircle_points(bary, base)
    if tweak:
        pts = xmath.normalize(pts)
    result = pts.sum(axis=1)
    return result#_b2s_fix_corners(bary, base, result)

def q2r_greatcircle(z, base):
    """Transforms a square to a spherical quadrilateral using the method of
    intersections
    """
    a = base[0]
    b = base[1]
    c = base[2]
    d = base[3]
    x = (z.real + 1)/2
    y = (z.imag + 1)/2
    f = xmath.slerp(a,b,x[:, np.newaxis])
    g = xmath.slerp(d,c,x[:, np.newaxis])
    h = xmath.slerp(b,c,y[:, np.newaxis])
    k = xmath.slerp(a,d,y[:, np.newaxis])
    return np.cross(np.cross(f, g), np.cross(h, k))/4

#double slerp
def _doubleslerp_points(bary, base):
    """Triangles of double slerp
    """
    x = (np.roll(bary, -1, axis=-1)/(1-bary))
    a0 = xmath.slerp(base[2], base[0], bary[:,0,np.newaxis])
    b0 = xmath.slerp(base[1], base[0], bary[:,0,np.newaxis])
    a1 = xmath.slerp(base[0], base[1], bary[:,1,np.newaxis])
    b1 = xmath.slerp(base[2], base[1], bary[:,1,np.newaxis])
    a2 = xmath.slerp(base[1], base[2], bary[:,2,np.newaxis])
    b2 = xmath.slerp(base[0], base[2], bary[:,2,np.newaxis])
    pts0 = xmath.slerp(a0, b0, x[:,0,np.newaxis])
    pts1 = xmath.slerp(a1, b1, x[:,1,np.newaxis])
    pts2 = xmath.slerp(a2, b2, x[:,2,np.newaxis])
    return np.stack([pts0, pts1, pts2], axis=1)

def b2r_doubleslerp(bary, base):
    """Transforms a triangle to a spherical triangle using the method of
    intersections
    """
    pts = _doubleslerp_points(bary, base)
    result = pts.mean(axis=1)
    return _b2s_fix_corners(bary, base, result)

def _square_doubleslerp(z, base):
    """
    Helper function for square_slerp. This does the slerp, and then
    square_slerp averages the two orientations together.
    Args:
        xy: Array, shape [..., 2]. XY coordinates on the square.
        base: Array, shape [4, ..., 3]. Coordinates of the square.
            Should be in counterclockwise order to maintain orientation.
    Returns:
        An array of shape [..., 3], representing points in 3d-space.
    """
    a, b, c, d = base[0], base[1], base[2], base[3]
    x = (z.real+1)/2
    y = (z.imag+1)/2
    ab = xmath.slerp(a, b, x[:, np.newaxis])
    dc = xmath.slerp(d, c, x[:, np.newaxis])
    result = xmath.slerp(ab, dc, y[:, np.newaxis])
    return result

def q2r_doubleslerp(z, base):
    """Transforms a square in [0,1]^2 to a spherical quadrilateral
    defined by base, using spherical linear interpolation
    """
    #have to do this twice and average them because just doing double
    #slerp isn't symmetric
    one = _square_doubleslerp(z, base)
    two = _square_doubleslerp(1j*z, np.roll(base, 1, axis=0))
    return (one + two)/2

#naive slerp

def _tri_naive_slerp_angles(bary, base, pow=1, eps=0):
    """Interpolates the angle factor so that it's equal to the
    angle between pts 1 and 2 when beta_3=0, etc.
    """
    angles = xmath.central_angle(base,
                                 np.roll(base, -1, axis=0))
    if np.max(angles) - np.min(angles) <= eps:
        return np.mean(angles)
    #if np.allclose(angles):
    #    return angles[..., 0]
    a = bary[...,0]
    b = bary[...,1]
    c = bary[...,2]
    ab = (a*b)**pow
    bc = (b*c)**pow
    ca = (c*a)**pow
    denom = ab + bc + ca
    numer = ab*angles[...,0] + bc*angles[...,1] + ca*angles[...,2]
    return numer/denom

def b2r_naive_slerp(bary, base, pow=1):
    """
    Naive slerp on a spherical triangle.
    """
    angles = _tri_naive_slerp_angles(bary, base, pow)[..., np.newaxis]
    b = np.sin(angles * bary) / np.sin(angles)
    result = b.dot(base)
    return _b2s_fix_corners(bary, base, result)


def _q_naive_slerp_angles(z, base, pow=1, eps=0):
    """Interpolates the angle factors separately so that it's equal to the
    angle between pts 1 and 2 when y=-1, etc.
    """
    x = z.real
    y = z.imag
    angles = xmath.central_angle(base, np.roll(base, -1, axis=0))
    ax = angles[0]
    bx = angles[2]
    ay = angles[3]
    by = angles[1]
    result1 = (ax*(1-y)**pow + bx*(1+y)**pow)/((1-y)**pow + (1+y)**pow)
    result2 = (ay*(1-x)**pow + by*(1+x)**pow)/((1-x)**pow + (1+x)**pow)
    return result1, result2

def q2r_naive_slerp(z, base, pow=1):
    """
    Naive slerp on a spherical quadrilateral.
    """
    anglex, angley = _q_naive_slerp_angles(z, base, pow)
    x = z.real
    y = z.imag
    sx = np.sin((1+x)*anglex/2)
    sy = np.sin((1+y)*angley/2)
    scx = np.sin((1-x)*anglex/2)
    scy = np.sin((1-y)*angley/2)
    a = scx * scy
    b = sx * scy
    c = sx * sy
    d = scx * sy
    mat = (np.stack([a, b, c, d], axis=-1) /
        (np.sin(anglex)* np.sin(angley))[:, np.newaxis] )
    result = mat.dot(base)
    return _q2s_fix_corners(z, base, result)

def _q_naive_slerp_2_angles(z, base, pow=1, eps=0):
    """Interpolates the angle factor together that it's equal to the
    angle between pts 1 and 2 when y=-1, etc.
    """
    angles = xmath.central_angle(base, np.roll(base, -1, axis=0))
    if np.max(angles) - np.min(angles) <= eps:
        return np.mean(angles)
    x = z.real
    y = z.imag
    a = ((1-x)*(1-y)*(1+x))**pow
    b = ((1-y)*(1+x)*(1+y))**pow
    c = ((1-x)*(1+x)*(1+y))**pow
    d = ((1-x)*(1-y)*(1+y))**pow
    numer = a*angles[0] + b*angles[1] + c*angles[2] + d*angles[3]
    denom = a + b + c + d
    return numer/denom

def q2r_naive_slerp_2(z, base, pow=1):
    """
    Variant naive slerp on a spherical quadrilateral.
    """
    angle =  _q_naive_slerp_2_angles(z, base, pow)[...,np.newaxis]
    x = z.real
    y = z.imag
    a = (1-x)*(1-y)
    b = (1+x)*(1-y)
    c = (1+x)*(1+y)
    d = (1-x)*(1+y)
    mat = np.sin(np.stack([a, b, c, d], axis=-1)*angle/4) / np.sin(angle)
    result = mat.dot(base)
    return _q2s_fix_corners(z, base, result)

#
def q2r_elliptical(z, base):
    """An extension of the elliptical map.
    """
    #FIXME needs rotations
    rot_base = base
    a = rot_base[0,0]
    b = rot_base[0,1]
    c = rot_base[0,2]
    x = z.real
    y = z.imag
    u = a * x * xmath.sqrt((1 - b**2*y**2)/(1-b**2))
    v = b * y * xmath.sqrt((1 - a**2*x**2)/(1-a**2))
    w = c * xmath.sqrt((1 - a**2*x**2)*(1 - b**2*y**2)/
                            ((1-a**2)*(1-b**2)))
    return np.stack([u,v,w], axis=-1)

#one weird cube-sphere mapping (topologists hate it)
def r2r_cube_sphere(pts):
    """
    Cube-to-sphere mapping, as per
    http://mathproofs.blogspot.com/2005/07/mapping-cube-to-sphere.html
    """
    x = pts[..., 0]
    y = pts[..., 1]
    z = pts[..., 2]
    xx = x*xmath.sqrt(1-y**2/2-z**2/2+y**2*z**2/3)
    yy = y*xmath.sqrt(1-z**2/2-x**2/2+z**2*x**2/3)
    zz = z*xmath.sqrt(1-x**2/2-y**2/2+x**2*y**2/3)
    return np.stack([xx,yy,zz], axis=-1)

_cube_face = np.array([[-1, -1, 1],
                       [1,  -1, 1],
                       [1,  1, 1],
                       [-1, 1, 1]])

def q2s_cube_sphere(xy, base=_cube_face):
    pts = q2r(xy, base)
    return r2r_cube_sphere(pts)

def q2s_qsc(xy, base=None):
    """
    http://www.sai.msu.su/~megera/wiki/SphereCube
    """
    x = xy.real
    y = xy.imag
    #x >= y?
    s_over_rx = np.sin(np.pi/12*(y/x))/(np.cos(np.pi/12*(y/x))-1/np.sqrt(2))
    qx = 1 - x*x*(1-1/xmath.sqrt(2+(s_over_rx)**2))
    rx = np.sign(x)*xmath.sqrt((1 - qx**2)/(1+(s_over_rx)**2))#
    sx = np.sign(y)*xmath.sqrt(1 - qx**2 - rx**2)#
    # x < y?
    s_over_ry = np.sin(np.pi/12*(x/y))/(np.cos(np.pi/12*(x/y))-1/np.sqrt(2))
    qy = 1 - y*y*(1-1/xmath.sqrt(2+(s_over_ry)**2))
    sy = np.sign(y)*xmath.sqrt((1 - qy**2)/(1+(s_over_ry)**2))#
    ry = np.sign(x)*xmath.sqrt(1 - qy**2 - sy**2)#
    result = np.stack([np.where(abs(x) >= abs(y),rx,ry),
                       np.where(abs(x) >= abs(y),sx,sy),
                       np.where(abs(x) >= abs(y),qx,qy)], axis=-1)
    result[xy == 0] = np.array([0,0,1])
    return result

def s2q_qsc(pts, base=None):
    """
    inverse of above
    """
    q = pts[..., 0]
    r = pts[..., 1]
    s = pts[..., 2]
    x = np.sign(r)*xmath.sqrt((1-q)/(1-1/xmath.sqrt(2+(s/r)**2)))
    y_over_x = (12/np.pi)*(np.arctan(s/r)- np.arcsin(s/xmath.sqrt(2*r*r+2*s*s)))
    return x + y_over_x * x * 1j

MAPPINGS_B2R = {
    'gn': lambda bary, base, tweak: b2r(bary, base),
    'ar': lambda bary, base, tweak: b2s_areal(bary, base),
    'co': lambda bary, base, tweak: b2s_conformal(bary, base),
    'gc': lambda bary, base, tweak: b2r_greatcircle(bary, base, tweak),
    'ds': lambda bary, base, tweak: b2r_doubleslerp(bary, base),
    'ns': lambda bary, base, tweak: b2r_naive_slerp(bary, base, pow=tweak)}

MAPPINGS_Q2R = {
    'gn': lambda z, base, tweak: q2r(z, base),
    'co': lambda z, base, tweak: q2s_conformal(z, base),
    'gc': lambda z, base, tweak: q2r_greatcircle(z, base),
    'ds': lambda z, base, tweak: q2r_doubleslerp(z, base),
    'ns': lambda z, base, tweak: q2r_naive_slerp(z, base, pow=tweak),
    'n2': lambda z, base, tweak: q2r_naive_slerp_2(z, base, pow=tweak)}

MAPPINGS_B2D = {
    'sq': lambda bary, base, tweak: b2d_squircle(bary),
    'ae': lambda bary, base, tweak: b2d_aea(bary),
    'ea': lambda bary, base, tweak: b2d_equalarea(bary),
    'rs': lambda bary, base, tweak: b2d_radial(bary)}

MAPPINGS_Q2D = {
    'sq': lambda z, base, tweak: q2d_squircle(z),
    'el': lambda z, base, tweak: q2d_el(z),
    'ae': lambda z, base, tweak: q2d_aea(z),
    'ea': lambda z, base, tweak: q2d_equalarea(z),
    'rs': lambda z, base, tweak: q2d_radial(z)}

MAPPINGS_B = {**MAPPINGS_B2D, **MAPPINGS_B2R}
MAPPINGS_Q = {**MAPPINGS_Q2D, **MAPPINGS_Q2R}

MAPPINGS_ALL = {**MAPPINGS_B, **MAPPINGS_Q}
MAPPINGS_D2S = {
    #'gn': c2s_gnomonic,
    'st': c2s_stereographic,
    'or': d2s_orthographic,
    'ea': d2s_ea,
    'ed': d2s_equidistant}

PARALLEL = ['ns', 'ns2', 'sq', 'ae', 'ea', 'rs', 'el']
