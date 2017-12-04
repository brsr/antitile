Appendix: Spherical geometry with vectors
=========================================
In the geometric dome world, geometry is usually synthetic rather than
analytic. That is, geometric constructions are usually expressed in terms of
step-by-step geometric constructions rather than in equations and numeric
values. See, for instance, [Kenner]_. Implementing a synthetic geometric 
construction on a computer can be a pain; analytic specifications are more
amenable to programming idioms.

There are a number of ways to numerically specify points on a sphere. By far
the most common is by latitude and longitude, which appear on every modern map
of the Earth and probably some other planets. Latitude and longitude are
familiar and convenient, but performing extended geometry calculations using 
latitude and longitude is a complicated task.

Doing spherical geometry using a 3-component unit vector is more convenient
in a number of ways: the equations are often simpler, and there are no
singularities at the poles. This appendix isn't a full course in spherical
geometry with vectors, but is here to explain the notation and methods
that are used in this text.

A 3d unit sphere can be defined as the set of all unit vectors in 3-space;
i.e., vectors :math:`\mathbf v = [v_x, v_y, v_z]` such that the vector norm
:math:`\|\mathbf v \|=1`. Unit vectors are often denoted using a hat:
:math:`\mathbf \hat{v}`. We'll often find ourselves normalizing vectors, so
if the numerator and denominator are the same, we'll supress the denominator
with an ellipsis like so:

.. math::
   \mathbf \hat{v} = \frac{\mathbf{some+really+long+statement}}{\|\dots\|}

Basics
------
The shortest distance between two points in Euclidean space is a straight
line. On the sphere, the shortest distance is an arc of the great circle
between those points. The great circle is the intersection of the sphere and a
plane passing through the origin. A plane through the origin can be specified
as :math:`\mathbf \hat{n} \cdot \mathbf v = 0`, where
:math:`\mathbf \hat{n}` is a unit vector normal to the plane; this vector
:math:`\mathbf \hat{n}` can be used to specify a great circle. Given two
points :math:`\mathbf{\hat{v}_1, \hat{v}_2}` on the sphere, the
:math:`\mathbf \hat{n}` of the great circle between those two points is
(up to normalization) their cross product:
:math:`\mathbf \hat{n} = \frac{\mathbf \hat{v}_1 \times \mathbf \hat{v}_2}{\|\dots\|}`

Small circles are the intersection of the sphere with a plane not through
the origin. All planes may be specified in Hessian normal form as
:math:`\mathbf \hat{n} \cdot \mathbf v = r`, where `r` is the
minimum distance between the plane and the origin.

=========== =================================================
     r      Intersection with unit sphere
=========== =================================================
0           Great circle
> 0 and < 1 Small circle
1           Point
> 1         None
=========== =================================================

Measurements
------------
.. list-table::
   :header-rows: 1

   * - Measurement
     - Euclidean
     - Spherical
   * - Length
     - :math:`\eta = \|\mathbf v_1-\mathbf v_2\|`
     - Central angle: :math:`\theta = \arctan\left(
       \frac{\|\mathbf \hat{v}_1\times \mathbf \hat{v}_2\|}
       {\mathbf \hat{v}_1 \cdot \mathbf \hat{v}_2}\right)`
   * - Angle
     - :math:`\cos \phi_1 = \frac{(\mathbf v_1 - \mathbf v_2) \cdot
       (\mathbf v_1 - \mathbf v_3)}
       {\|\mathbf v_1 - \mathbf v_2\|\|\mathbf v_1 - \mathbf v_3\|}`
     - Dihedral angle:
       :math:`\mathbf c_{12} = \mathbf \hat v_1 \times \mathbf \hat v_2`,
       :math:`\mathbf c_{13} = \mathbf \hat v_1 \times \mathbf \hat v_3`,
       :math:`\phi_1 = \mathrm{arctan2}\left(\mathbf \hat{v}_1
       \cdot \mathbf \hat{v}_2 \times \mathbf \hat{v}_3,
       \mathbf c_{12} \cdot \mathbf c_{13}\right)`
   * - Triangle area
     - :math:`A = \frac{\|(\mathbf v_1-\mathbf v_3)\times
       (\mathbf v_2-\mathbf v_3)\|}{2}`
     - Solid angle: :math:`\tan(\Omega/2) = \frac{|\mathbf \hat v_1 \cdot
       \mathbf \hat v_2 \times \mathbf \hat v_3|}
       {1+\mathbf \hat v_1\cdot \mathbf \hat v_2+\mathbf \hat v_2
       \cdot \mathbf \hat v_3+\mathbf \hat v_3\cdot \mathbf \hat v_1}`
       [Oosterom]_

In general, the spherical measures approach the euclidean measures when the
measures are small. (A sphere looks like a plane if you get close enough.)

Constructions
-------------
.. list-table::
   :header-rows: 1

   * - Construction
     - Euclidean
     - Spherical
   * - Mean (midpoint when n=2, centroid when n=3)
     - :math:`\mathbf v_\mu = \frac{\sum\mathbf v_i}{n}`
     - :math:`\mathbf \hat v_\mu = \frac{\sum\mathbf \hat v_i}{\|\dots\|}`
   * - Interpolation
     - :math:`\mathrm{Lerp}(\mathbf{v_1}, \mathbf{v_2}; t) =
       (1-t) \mathbf{v_1} + t \mathbf{v_2}`
     - :math:`\mathrm{Slerp}(\mathbf{\hat{v}_1}, \mathbf{\hat{v}_2}; t) =
       \frac{\sin {((1-t)w)}}{\sin (w)} \mathbf{\hat{v}_1} +
       \frac{\sin (tw)}{\sin (w)} \mathbf{\hat{v}_2}`
       where :math:`\cos(w) = \mathbf{\hat{v}_1} \cdot \mathbf{\hat{v}_2}`

Face normals
------------
This program uses this definition for the normal to a Euclidean polygon:

.. math::
   \mathbf \hat n = \frac{\sum v_i \times v_{i+1}}{\|\dots\|}`

`i` should be treated as if it's :math:`\mod n`: it loops around. This 
definition allows for a somewhat sensible extension to skew polygons:
the normal points in the general direction that's expected.

The normal will be outward-facing if the points are ordered counterclockwise,
and inward-facing if the points are ordered clockwise.
