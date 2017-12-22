Appendix: Goldberg-Coxeter Operation on improper spherical polyhedra
====================================================================
There are some improper spherical polyhedra that cannot be tesselated with
the traditional geodesic dome methodology but can with some methods developed
earlier. The Naive Slerp methods are able to project a triangle or square
onto a full hemisphere, including the boundary.

Transformations to the disk
---------------------------
When the base vertices all lie on the same great circle (i.e. in a plane
through 0), Naive Slerp before projection or normalization
transforms the triangle or square to the disk.

There are also some disk-specific transformations, accessible with
``gcopoly -p=disk``. These amount to converting the vertex to polar form
and scaling the `r` variable. Unlike Naive Slerp, there are points where the
transformation is not smooth.

Triangular:

.. math::
    r = 1 - 3 \min(\beta_i)

Quadrilateral:

.. math::
    r = \max(|2x - 1|, |2y - 1|)

These disk operations can be projected with a k factor the same
as naive slerp, although here we require k to be `>0`: otherwise,
all points would be projected onto the great circle.

Dihedra
-------
With the transformations defined above, performing the GC operation on a
dihedron is fairly straightforward. The only caveat has to do with the
presence of vertices with valence 2. On the quadrilateral, this produces
a degenerate edge: a valence-2 vertex and its edges can be replaced with
a single edge. On the triangle, this produces dangling faces, which can be
removed. In both cases, one can also add or switch edges so that the
vertex is no longer degenerate.

The process of doing the GC operation on a dihedron resembles stitching
together two flat sheets and inflating them, like the children's novelty
the whoopi cushion. The dual of that could be called a Whoopi Goldberg
polyhedron. (Thank you, I'll be here all night.)

Henahedron
----------
The henahedron is the "polyhedron" with one face, one vertex, and no edges:
the entire sphere is one big face. Topologically, we can take the disk defined
earlier and collapse the boundary to a single point. This amounts to
projecting the disk in the manner of a map projection. Given :math:`\theta =
\arctan2(v_x, v_y)` and :math:`\mathbf{\hat{v}} = (\sin(\phi) \cos(\theta),
\sin(\phi) \sin(\theta), \cos(\phi))`, two useful projections used in
cartograph can be described as such:

Lambert azimuthal equal-area

.. math::
    \phi = 2 arcsin(\|\mathbf v\|)

Azimuthal equal-distance

.. math::
    \phi = \pi \|\mathbf v\|

Because of the radical change to the boundary, many small-parameter GC
operations produce the same polyhedron. Quad faces may be reduced to triangle
faces, and triangle faces may be reduced to degenerate faces. Many produced
polyhedra are not convex. Class III polyhedra tend to have dangling faces.
