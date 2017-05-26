# -*- coding: utf-8 -*-
"""
Subdivide a tiling or polyhedra using a similar grid

Attributes:
    MEASURES: A map from the measure's name to a function that can be used
        in the optimization routine for k.
"""
import warnings
import numpy as np
from numpy.linalg import norm
from scipy.optimize import minimize_scalar
from scipy import sparse
from . import tiling, breakdown, projection, xmath, factor

def _stitch(edge, faces, bkdns, freq):
    """
    Stitch together two faces sharing an edge.

    Args:
        edge: Edge shared by the faces in faces.
        faces: Two lists of vertices in each face.
        bkdns: A map from {3, 4} to the breakdown structures for
            triangle and quadrilaterals, respectively.
        freq: Frequency of subdivision (a 2-tuple)

    Returns:
        An array of shape (..., 2) containing vertices that
        should be merged together.
    """
    a, b = freq
    face_0, face_1 = faces
    shape_0 = len(face_0)
    shape_1 = len(face_1)
    bkdn_0 = bkdns[shape_0]
    bkdn_1 = bkdns[shape_1]
    improper = (a == b or b == 0)
    roll_0, flip_0 = tiling.orient_face(face_0, edge)
    roll_1, flip_1 = tiling.orient_face(face_1, edge)
    flip = (flip_0 == flip_1) and improper
    lindices = bkdn_0.lindex_reorient(roll_0, flip=flip)
    if shape_1 == 3:
        offset = np.array([a+b, a, 2*a+b])
    else:
        offset = np.array([a+2*b, b])
    flipped = offset - bkdn_1.lindex_reorient(roll_1)
    if shape_0 == shape_1:
        comparison = np.all(lindices[:, np.newaxis] == flipped[np.newaxis],
                            axis=-1)
    elif b == 0:#class I, mixed grid
        cmp_1 = lindices[:, 0, np.newaxis] == flipped[np.newaxis, ..., 0]
        cmp_2 = lindices[:, -1, np.newaxis] == (0 if shape_0 == 4 else a)
        cmp_3 = flipped[np.newaxis, ..., -1] == (a if shape_0 == 4 else 0)
        comparison = (cmp_1 & cmp_2 & cmp_3)
    else:
        msg = """Class II and III subdivision on mixed triangle-
                 quadrilateral polyhedra is not implemented"""
        raise NotImplementedError(msg)
    matches = np.argwhere(comparison)
    return matches

def _find_dupe_verts(base, bf, group, freq, bkdns):
    """
    Find duplicate vertices in a subdivided polyhedron.

    Args:
        base: Base polyhedron
        bf: Corresponding base face of each vertex
        group: Group number of each vertex
        freq: Frequency of subdivision (a 2-tuple)
        bkdns: A map from {3, 4} to the breakdown structures for
            triangle and quadrilaterals, respectively.
    """
    #find redundant vertices
    base_edges = base.edges
    base_edge_corr, base_face_corr = base.faces_by_edge(base_edges)
    l0 = []
    l1 = []
    for i in range(len(base_edges)):
        edge = base_edges[i]
        index = base_edge_corr == i
        facex = base_face_corr[index]
        fn = len(facex)
        if fn > 2:
            warnings.warn("More than 2 faces meet at a single edge. "
                          "Choosing 2 faces arbitrarily...")
            facex = facex[:2]
        elif fn < 2:#external edge, skip it
            continue
        index_0 = bf == facex[0]
        index_1 = bf == facex[1]
        faces = [base.faces[facex[0]], base.faces[facex[1]]]
        match = _stitch(edge, faces, bkdns, freq)
        lx0 = np.argwhere(index_0)[match[..., 0]].flatten()
        lx1 = np.argwhere(index_1)[match[..., 1]].flatten()
        l0.extend(lx0)
        l1.extend(lx1)
    matches = np.stack([l0, l1], axis=-1)
    #TODO replace this with np.unique when 1.13 comes out
    matches = np.array(sorted(list({tuple(sorted(t)) for t in matches})))
    vno = len(group)
    conns = sparse.coo_matrix((np.ones(len(matches)),
                               (matches[:, 0], matches[:, 1])),
                              shape=(vno, vno))
    ncp, cp = sparse.csgraph.connected_components(conns)
    verts = np.arange(vno, dtype=int)
    for i in range(ncp):
        component = np.argwhere(cp == i).flatten()
        gp = group[component]
        order = np.argsort(gp)
        component = component[order]
        v = verts[component[0]]
        verts[component] = v
    unique_index = verts == np.arange(len(verts))
    renumbered = xmath.renumber(unique_index)
    return renumbered[verts], unique_index

class SGS(tiling.Tiling):
    """Similar grid subdivision of a tiling.
    Args:
        base: Tiling/polyhedron to base the subdivision on. If the tiling has
            faces with 5 or more sides, ValueError will be raised.
        freq: Subdivision frequency.
        proj: Projection family to use
        tweak: Tweaks to certain calculations

    Attributes:
        vertices: List of vertices
        faces: List of faces
        base: Pointer to the base polyhedron that was subdivided
        freq: Frequency of the subdivision
        base_face: Which base face each vertex falls into. If vertex is
            on the border of more than one base face, will be the lower
            numbered base face.
        group: Vertex grouping (see ``Breakdown`` for definition)
        face_bf: Which base face each face lies on. If face lies over
            more than one face, will be whichever base face contains more
            points of the face; if evenly split, it will be the lower
            numbered base face.
        face_group: Face grouping (see ``Breakdown`` for definition)
        """

    def __init__(self, base, freq=(2, 0), proj='flat', tweak=False):
        self.base = base
        self.freq = freq
        a, b = freq
        projections = projection.PROJECTIONS[proj]
        base_faces = base.faces
        base_verts = base.vertices
        fs = base.face_size
        if np.any(fs > 4):
            raise ValueError("Tiling contains at least one face with more than "
                             "4 sides. Try triangulating those faces first.")
        elif np.all(fs < 2):
            raise ValueError("Polyhedron has no faces")
        bkdns = dict()
        if 3 in fs:
            bkdns[3] = breakdown.Breakdown(*freq, 3)
        if 4 in fs:
            bkdns[4] = breakdown.Breakdown(*freq, 4)
        n_bf = len(fs)
        v_f_list = []
        v_f_bf = []
        v_f_group = []
        faces = []
        f_faces_bf = []
        f_face_group = []
        for i in range(n_bf):
            n = fs[i]
            proj_fun = projections[n]
            bkdn = bkdns[n]
            n_pts = len(bkdn.group)
            face = base_verts[base_faces[i]]
            n_v_i = sum(len(x) for x in v_f_list)
            with np.errstate(divide='ignore', invalid='ignore'):
                v_f_list.append(proj_fun(bkdn, face, freq, tweak))
            v_f_bf.append(i*np.ones(n_pts, dtype=int))
            v_f_group.append(bkdn.group)
            faces.extend(x + n_v_i for x in bkdn.faces)
            f_faces_bf.append(i*np.ones(len(bkdn.faces)))
            f_face_group.append(bkdn.face_group)
        vertices = np.concatenate(v_f_list, axis=0)
        bf = np.concatenate(v_f_bf, axis=0)
        group = np.concatenate(v_f_group, axis=0)
        face_bf = np.concatenate(f_faces_bf, axis=0)
        face_group = np.concatenate(f_face_group, axis=0)
        vert_index, unique_v = _find_dupe_verts(base, bf, group, freq, bkdns)
        vertices = vertices[unique_v]
        bf = bf[unique_v]
        group = group[unique_v]
        faces = [vert_index[face] for face in faces]
        super().__init__(vertices, faces)
        self.normalize_faces()
        self.group = group
        self.base_face = bf
        if b == 0:
            #faces don't overlap if b == 0, so just take them as they are
            self.face_bf = face_bf
            self.face_group = face_group
        else:
            unique_f = self.id_dupe_faces(face_group)
            for i in np.argwhere(~unique_f):
                self.faces[int(i)] = None
            self.faces = [face for face in self.faces if face is not None]
            self.face_bf = face_bf[unique_f]
            self.face_group = face_group[unique_f]

def build_sgs(base, frequency, proj, k=1, tweak=False, normalize=True):
    """Wrapper around the SGS constructor that performs parallel projection
    and normalization.

    Args:
        base: Base polyhedron/tiling.
        frequency: Frequency of subdivision
        proj: Projection family
        k: Parallel projection constant. Can be a float, or an element of
        ``MEASURES``, in which case it will choose a k that optimizes that
        measure.
        tweak: Tweak to certain calculations. False by default.
        normalize: Whether to normalize the polyhedron's vertices.
            True by default.
    """
    poly = SGS(base, frequency, proj, tweak)
    if proj in projection.PARALLEL:
        if k in MEASURES:
            measure = MEASURES[k]
            k = optimize_k(poly, base, measure, not tweak, normalize)
        else:
            k = float(k)
        poly.vertices += k*parallels(poly, base, exact=not tweak)
    if normalize:
        poly.vertices = xmath.normalize(poly.vertices)
    poly.k = k
    return poly

def build_sgs_rep(base, frequency, proj, k=1, tweak=False,
                  normalize=True):
    """Wrapper around the SGS constructor that factors the frequency
    and uses the factors to repeatedly subdivide the grid. Uses
    lower-norm factors first. Also does everything `build_sgs` does.

    Args:
        base: Base polyhedron/tiling.
        frequency: Frequency of subdivision
        proj: Projection family
        k: Parallel projection constant. Can be a float, or an element of
        ``MEASURES``, in which case it will choose a k that optimizes that
        measure.
        tweak: Tweak to certain calculations. False by default.
        normalize: Whether to normalize the polyhedron's vertices.
            True by default.
    """
    fs = base.face_size
    fsu = np.unique(fs)
    fsu = fsu[fsu > 2]
    if len(fsu) != 1:
        raise ValueError('cannot perform repeated subdivision on mixed grid')
    elif fsu[0] == 3:
        element = factor.Nietsnesie
    elif fsu[0] == 4:
        element = factor.Gaussian
    else:
        msg = ("cannot perform repeated subdivision on grid of faces with "
               + str(fsu[0]) + " sides")
        raise ValueError(msg)
    q = element(*frequency)
    fq = q.factor()

    result = base
    for f in fq:
        if f.anorm() <= 1:
            continue
        freq = f.tuple
        result = build_sgs(result, freq, proj, k=k, tweak=tweak,
                           normalize=normalize)
    return result


#--Stuff having to do with the k-factor--
def parallels(poly, base, exact=True):
    """Given a subdivided polyhedron based on a base polyhedron, return
    the parallels to the base faces for each vertex in the polyhedron
    that would put the vertices onto the sphere

    Args:
        poly: Subdivided polyhedron
        base: Base polyhedron
        exact: Whether to use exact or approximate parallels. By default, True.

    Returns:
        Parallel vectors for each vertex in poly.
    """
    normals = base.face_normals[poly.base_face]
    return projection.parallel(poly.vertices, normals, exact)

def parallel_sphere(xyz, pls, k=1):
    """Given vertices and parallels, return points on (or near) the sphere.

    Args:
        xyz: Vertices
        pls: Parallel vectors
        k: Projection constant. Defaults to 1.

    Returns:
        Vertices, parallel projected."""
    return xyz + k*pls

def optimize_k(poly, base, measure, exact=True, normalize=True):
    """Routine to optimize the factor k

    Args:
        poly: Subdivided polyhedron
        base: Base polyhedron
        measure: Which measure to optimize on.
        exact: Whether to use exact or approximate parallels. By default, True.
        normalize: Whether to normalize the vectors. Default: True:
            True: Normalize the vectors
            False: Preserve vector norms
            None: Don't normalize at all"""
    parallel_xyz = parallels(poly, base, exact)
    result = minimize_scalar(_objective, bracket=[0, 2],
                             args=(poly, parallel_xyz, measure, normalize))
    if not result.success:
        warnings.warn('Optimization routine did not converge')
    return result.x

def _objective(k, poly, parallel_xyz, measure, normalize=True):
    """Objective function for the optimization routine

    Args:
        k: Projection constant value
        poly: Subdivided polyhedron
        parallel_xyz: Parallel vectors
        measure: Which measure to optimize on.
        normalize: Whether to normalize the vectors. Default: True:
            True: Normalize the vectors
            False: Preserve vector norms
            None: Don't normalize at all

    Return:
        Measure value for optimization rountine
    """
    test_v = parallel_sphere(poly.vertices, parallel_xyz, k)
    if normalize:
        test_v = xmath.normalize(test_v)
    elif not normalize and normalize is not None:
        #restore length of original vector
        test_v = xmath.normalize(test_v)
        test_v *= norm(poly.vertices, axis=-1, keepdims=True)
    #if it's None, don't normalize at all
    return measure(test_v, poly)

#measures
def energy(xyz, poly=None):
    """Energy of the vertex arrangement, as defined in
    Thomson's problem.

    Args:
        xyz: Coordinates of vertices.
        poly: Ignored, only present to give this function
            the same signature as the other measures. Defaults to None.

    Returns:
        Energy of the vertex arrangement
    """
    return tiling.energy(xyz)

def bentness(xyz, poly):
    """Measure: face bentness.

    Args:
        xyz: Coordinates of vertices.
        poly: Tiling object, used to get the polyhedral structure of the
            vertices. (xyz is separate so we don't have to recalculate the
            structure every time we change the vertex position.)

    Returns:
        Maximum face bentness. (Not ratio of maximum to
        minimum face bentness.)
    """
    return tiling.bentness(xyz, poly).max()

def edge_length(xyz, poly, spherical=False):
    """Measure: edge length.

    Args:
        xyz: Coordinates of vertices.
        poly: Tiling object, used to get the polyhedral structure of the
            vertices. (xyz is separate so we don't have to recalculate the
            structure every time we change the vertex position.)
        spherical: If true, uses the spherical distance (central angle)
            along each edge instead of the Euclidean distance. Default False.

    Returns:
        Ratio of maximum to minimum edge length.
    """
    result = tiling.edge_length(xyz, poly.edges, spherical)
    return result.max()/result.min()

def face_area(xyz, poly, spherical=False):
    """Measure: face area.

    Args:
        xyz: Coordinates of vertices.
        poly: Tiling object, used to get the polyhedral structure of the
            vertices. (xyz is separate so we don't have to recalculate the
            structure every time we change the vertex position.)
        spherical: If true, uses the spherical area (solid angle)
            of each face instead of the Euclidean area. Default False.

    Returns:
        Ratio of maximum to minimum face area.
    """
    result = tiling.face_area(xyz, poly, spherical)
    return result.max()/result.min()

def aspect_ratio(xyz, poly, spherical=False):
    """Measure: aspect ratio.

    Args:
        xyz: Coordinates of vertices.
        poly: Tiling object, used to get the polyhedral structure of the
            vertices. (xyz is separate so we don't have to recalculate the
            structure every time we change the vertex position.)
        spherical: If true, uses the spherical distance (central angle)
            along each edge instead of the Euclidean distance. Default False.

    Returns:
        Maximum aspect ratio of any one face.
        (Not ratio of max to min aspect ratio.)
    """
    result = tiling.aspect_ratio(xyz, poly, spherical)
    return result.max()

MEASURES = {'energy': energy,
            'bent': bentness,
            'edges': edge_length,
            'aspect': aspect_ratio,
            'faces': face_area,
            'angle': lambda xyz, poly: edge_length(xyz, poly, spherical=True),
            'angle_aspect': lambda xyz, poly: aspect_ratio(xyz, poly,
                                                           spherical=True),
            'solid_angle': lambda xyz, poly: face_area(xyz, poly,
                                                       spherical=True)}
