Goldberg-Coxeter Operation on spherical polyhedra
=================================================
Spherical (genus 0) polyhedra are the standard use case for the GC operation.
There are many methods in the geodesic dome literature for placing the
vertices on the surface of a sphere. Many of these are synthetic (step-by-step
geometric contstructions) rather than analytic (expressed in equations). For
a computer application, analytic forms can be more convenient. This section
might not cover every method ever written, but it expresses the major
methods in analytic form and introduces a couple new ones.

In this section, we'll use :math:`\mathbf v^*` to denote a pre-normalized
vector: :math:`\mathbf \hat{v} = \frac{\mathbf v^*}{\|\mathbf v^*\|}`

Gnomonic
--------
``gcopoly -p=flat``

The gnomonic projection was known to the ancient Greeks, and is the simplest
of the transformations listed here. It has the nice property that all lines in
Euclidean space are transformed into great circles on the sphere: that is,
geodesics stay geodesics. This is in fact the motivation for the name
"geodesic dome": Buckminster Fuller used this projection to project triangles
on the sphere. This is referred to as Method 1 in geodesic dome circles.
The main downside is that the transformation causes shapes
near the corners to appear bunched up;
this is particularly bad for larger faces e.g. on the tetrahedron.

In general, the gnomonic projection is defined as:

* To sphere: :math:`\mathbf \hat{v} = \frac{\mathbf p}{\|\mathbf p\|}`
* From sphere: :math:`\mathbf p = \frac{r\mathbf \hat{v}}
  {\mathbf \hat{n} \cdot \mathbf\hat{v}}`

where :math:`\mathbf p` is a point on a plane given in Hessian normal
form by :math:`\mathbf \hat{n} \cdot \mathbf p = r`. Projection from Euclidean
space to the sphere is literally just normalizing the vector.

For the Goldberg-Coxeter operation, this amounts to just normalizing
the vectors produced by the coordinate form:

.. math::
   \mathbf v^* =
   \beta_1 \mathbf v_1 + \beta_2 \mathbf v_2 + \beta_3 \mathbf v_3

.. math::
   \mathbf v^* = \mathbf v_1 + (\mathbf v_2-\mathbf v_1) x +
   (\mathbf v_4-\mathbf v_1) y +
   (\mathbf v_1-\mathbf v_2+\mathbf v_3-\mathbf v_4)xy

where :math:`\beta_i` are (planar) barycentric coordinates and :math:`x,y` are
x-y quadrilateral coordinates.

The case of barycentric coordinates can be thought of as generalized
barycentric coordinates, where :math:`\mathbf v = \sum\beta_i\mathbf v_i`
still holds but :math:`\sum \beta_i` is not necessarily `=1`. If the
generalized coordinates are :math:`\beta^\prime_i`, then
:math:`\beta^\prime_i = \frac{\beta_1}
{\beta_1 \mathbf v_1 + \beta_2 \mathbf v_2 + \beta_3 \mathbf v_3}`. On the
interior of the triangle, :math:`\sum \beta^\prime_i > 1`.

Spherical areal
---------------
``gcopoly -p=other``

This method applies to triangles only. The above section discussed generalized
barycentric coordinates, which removes the restriction that
:math:`\sum \beta_i = 1`. Instead we look to the relation between barycentric
coordinates and area; :math:`\beta_i` is the proportion of area in the
triangle that is opposite the vertex :math:`v_i`. Let :math:`\Omega` be the
spherical area (solid angle) of the spherical triangle and
:math:`\Omega_i = \beta_i\Omega` be the area of the smaller triangle
opposite vertex :math:`v_i`.

While this method introduces less bunching-up than the gnomonic method,
the equation to find :math:`\mathbf \hat{v}` given :math:`\beta_i` is much
more complicated than for barycentric coordinates. Let
:math:`h_i = \sin\Omega_i\left(1+\mathbf v_{i-1}\cdot\mathbf v_{i+1}\right)`,
and
:math:`\mathbf g_{i} = \left(1+\cos \Omega_{i}\right) \mathbf v_{i-1} \times
\mathbf v_{i+1} - \sin\Omega_{i}\left(\mathbf v_{i-1} + \mathbf v_{i+1}\right)`
where the subscripts loop around: 0 should be interpreted as 3 and 4 should be
interpreted as 1. Then

.. math::
   \mathbf G = \begin{bmatrix} \mathbf g_1 & \mathbf g_2 & \mathbf g_3 \end{bmatrix}

.. math::
   \mathbf h = \begin{bmatrix} h_1  & h_2 & h_3  \end{bmatrix}^T

such that :math:`\mathbf G \mathbf \hat{v} = \mathbf h` To clarify,
:math:`\mathbf G` is the 3x3 matrix where the `i`th column is
:math:`\mathbf g_i`, and :math:`\mathbf h` is the column vector where the
`i`th element is :math:`h_i`. The vector :math:`\mathbf \hat{v}` can be solved for
using standard matrix methods. (This mess can be derived from the formula
for solid angle in the appendix.)

Naive Slerp
-----------
``gcopoly -p=nslerp``, ``gcopoly -p=nslerp2``

This method shares with the gnomonic method an analytic form for the
transformation from Euclidean space to the sphere, but doesn't have the
problem with bunching up. It has two downsides: there is no analytic form
for the reverse transformation, and it can only be used on equilateral faces.

The Naive Slerp method on a triangular face resembles a naive extension of
spherical linear interpolation (Slerp) to barycentric coordinates, thus the
name. The Naive Slerp methods reduce to slerp on the edges of the face.

Let :math:`\cos(w) = \mathbf v_i \cdot \mathbf v_{i+1}` for all :math:`i`. (As
usual, the subscripts loop around.)

Triangle:
:math:`\mathbf v^* =
\sum_{i=1}^3\frac{\sin(w\beta_i)}{\sin(w)}  \mathbf v_i`

Quadrilateral 1:
:math:`\mathbf v^* =
\sum_{i=1}^4\frac{\sin(w\gamma_i)}{\sin(w)}  \mathbf v_i`
where
:math:`\gamma_1 = (1-x)(1-y)`,
:math:`\gamma_2 = x(1-y)`,
:math:`\gamma_3 = xy`,
:math:`\gamma_4 = (1-x)y`

Quadrilateral 2:
:math:`\mathbf v^* = \sum_{i=1}^4\frac{s_i}{\sin^2(w)}  \mathbf v_i`
where
:math:`s_1 = \sin (w(1-x))\sin (w(1-y))`,
:math:`s_2 = \sin (wx)\sin (w(1-y))`,
:math:`s_3 = \sin (wx)\sin (wy)`,
:math:`s_4 = \sin (w(1-x))\sin (wy)`

Because the projected edges already lie on the sphere, there is a lot of
freedom in how to adjust :math:`\mathbf v^*` to lie on the sphere.
The easiest is just to centrally project the vertices, that is, to normalize
:math:`\mathbf v^*` like we have been. Another option is to perform a parallel
projection along the face normal. (See appendix for a formula for the "normal"
to a skew face.) We need the parallel distance `p` from the vertex to the
sphere surface in the direction of the face normal \mathbf \hat{n}, such that
:math:`\mathbf \hat v = \mathbf v^* + p\mathbf \hat{n}`. That is given by:

.. math::
   p = -\mathbf v^* \cdot \mathbf \hat{n} +
   \sqrt{1+\mathbf v^* \cdot \mathbf \hat{n}-\mathbf v^* \cdot \mathbf v^*}

`p` can also be approximated as :math:`\widetilde{p} = 1 - \|\mathbf v^*\|
\leq p`, which is somewhat fewer operations and doesn't require
calculation of the face normal.

Sometimes the best polyhedron comes from a compromise of the central and
parallel projections. Choose a constant `k`, typically between 0 and 1, then:
.. math::
   \mathbf \hat v = \frac{\mathbf v^* + kp\mathbf c}{\|\dots\|}

`p` may be replaced by :math:`\widetilde{p}`. 

If our goal is to optimize a measurement of the polyhedron, we can do a 
1-variable optimization on `k`, which is usually more tractable than 
the multivariate optimization of the location of every vertex.

Great circle intersection
-------------------------
``gcopoly -p=gc``

This method draws great circles between points on the edges of the polygon,
and uses the points of intersections of those great circles to determine
the vertices. This is the only method that directly uses the linear indexes
defined in the previous chapter. In the geodesic dome world, this is called
Method 2, although this description is considerably more complicated to
accomodate Class III grids and quad faces.

Specify the linear index as :math:`\ell = (e,f,g)` for triangular faces or
:math:`\ell = (e,f)` for quad faces. Using slerp, calculate the points
:math:`\mathbf{\hat{b}}_{i,j,k}` where each line from the breakdown structure
crosses the face edge. `i` is the coordinate of the linear index, `j` is the
linear index (:math:`\ell_i = j`), and `k` is which point of intersection with
the polygon. The line corresponds to a great circle normal given by
:math:`\mathbf{\hat{n}}_{i,j}
= \frac{\mathbf{\hat{b}}_{i,j,0} \times \mathbf{\hat{b}}_{i,j,1}}{\|\dots\|}`.

We are going to calculate the intersection of these great circle normals.
The intersection of two planes is a line: the intersection of two great
circles is two antipodal points. We need to choose the point on the correct
side of the sphere. Let :math:`\mathbf{c}` be the centroid of the face: then
:math:`\mathbf{v}` is on the right side of the sphere if :math:`\mathbf{v}
\dot \mathbf{c} >0`. If not, just multiply :math:`\mathbf{v}` times -1 to put
it on the right side.

For quad faces, there are only two intersecting great
circles, so the new vertices are :math:`\mathbf{v^*}_{\ell} =
\mathbf{\hat{n}}_{1,\ell_{1}} \times \mathbf{\hat{n}}_{2,\ell_{2}}` (possibly
times -1).

For triangular faces, there are three intersecting great circles, and unlike
on the plane, on the sphere they do not necessarily intersect in the same
place. Each pair of great circles forms a vertex of a triangle as
:math:`\mathbf{\hat{v}}_{\ell, m} = \frac{\mathbf{\hat{n}}_{m,\ell_{m}} \times
\mathbf{\hat{n}}_{m+1,\ell_{m+1}}}{\|\dots\|}`,
Make sure all the :math:`\mathbf{\hat{v}}_{\ell, m}` lie on the correct side
of the sphere, and then take the centroid of that triangle to get the vertex:
:math:`\mathbf{v^*}_{\ell} = \sum_m \mathbf{v}_{\ell, m}`. (It is not
strictly necessary to take the centroid: the sum of the unnormalized
:math:`\mathbf{v}_{\ell, m}` will also be a point somewhere within the
triangle.)

Summary of methods
------------------
.. list-table::
   :header-rows: 1

   * - Method
     - Gnomonic
     - Spherical areal
     - Naive slerp
     - Great circle intersection
   * - Geodesic dome name
     - Method 1
     - New
     - New
     - Method 2
   * - Input
     - Coordinates (barycentric or xy)
     - Barycentric coordinates
     - Coordinates (barycentric or xy)
     - Linear index (triangular or quadrilateral)
   * - Adjustment to sphere
     - Central projection
     - Not needed
     - Any projection
     - Central projection
   * - Face
     - :math:`\Omega < 2\pi`
     - :math:`\Omega < 2\pi`
     - :math:`\Omega <= 2\pi`, must be equilateral
     - :math:`\Omega < 2\pi`


Multi-step methods
------------------
As mentioned earlier, the operators :math:`\Delta(a,b)` and :math:`\Box(a,b)`
may be able to be decomposed into a series of smaller operators. Many of the
smaller operators are constrained by symmetry: in particular,
:math:`\Delta(2,0)` adds vertices at the midpoints of each edge, independent
of the method used. Method 3 in geodesic dome terminology is simply repeated
application of :math:`\Delta(2,0)`.

In a more general sense, an operator can be factored into a series of "prime"
operators, and applied in order from small to large. The faces of the
polyhedron will become progressively smaller and therefore progressively
flatter, and as the faces get flatter, the differences between methods becomes
smaller. As an example, :math:`\Box(16,4) = \Box^4(1,1)\Box(4,1)`, so apply
the highly-symmetric operator :math:`\Box(1,1)` (which creates one vertex at
the centroid of a face) four times and then :math:`\Box(4,1)` once with
a simple method like Gnomonic.

``geodesic`` in Antitile performs class II and III subdivision by finding the
smallest class I operator that can be decomposed into the desired operator
and some other factor. Effectively, given :math:`\Delta(a,b)`, it finds the
smallest `n` such that :math:`\Delta(a,b)\Delta(c,d) = \Delta(n,0)` for some
`c` and `d`. That is given by :math:`n = \frac{t(a,b)}{\hcf{a,b}}`,
where hcf is the highest common factor. It calculates :math:`\Delta(n,0)` using
method 2 and then uses the vertices from that that are shared with
:math:`\Delta(a,b)`.

Skew faces
----------
Skew faces are impossible on a polyhedra with triangular faces. On a polyhedron
with quadrilateral faces, however, all of the above methods produce skew
faces. There are basically two solutions to the issue. The first is to treat
the polyhedron purely as a spherical polyhedron: all the faces are curved tiles
on the surface of a sphere, and we can ignore whether they're skewed in
Euclidean space. The second is to canonicalize the polyhedron. As per
[Hart1997]_, all convex polyhedra can be put into a unique
`canonical form` such that:

* All the edges are tangent to the unit sphere,
* The origin is the average of the points at which the edges touch the sphere,
  and
* The faces are flat (not skew)

The ``canonical`` program in Antiprism performs canonicalization via a simple
iterative process. The vertices of the faces probably do not lie on the
unit sphere. If a polyhedron created by Goldberg-Coxeter
operations is to be canonicalized, the choice of method does not matter except
as a starting point.

Choosing a method
-----------------
Which method is better depends on your criteria. If your only criteria is
speed, gnomonic is the simplest and fastest to run, and produces reasonable
results on typical cases like the icosahedron and cube. 

A common criteria is to have a polyhedron that minimizes some measure of 
deviation. That could be a number of things:

* Thomson energy	
* Edge length (spherical or euclidean)
* Aspect ratio for each face (spherical or euclidean)
* Face area (spherical or euclidean, only well-defined for non-skew faces)
* Face bentness (only applies to quad faces)

Often the objective will be to minimize either the maximum value or the ratio
of maximum to minimum values. This package includes, in the ``data`` directory,
a csv file containing measurements of polyhedra produced by GC operations
on the icosahedron, octahedron, tetrahedron, and cube. No method always comes 
out the winner, but some patterns emerge:

* Gnomonic doesn't ourperform any other method (except on speed)
* Naive slerp often outperforms all other methods. This is true when the 
  k-factor is fixed at 1 and moreso when k-factor is optimized
* On quad faces, the first Naive Slerp generally outperforms the second.
* ``geodesic`` in Antiprism performs well with regards to edge length ratio.
* Spherical areal performs well for certain measures on large faces, 
  e.g. edge length ratio and face area for the tetrahedron
* No method (aside from canonicalization) consistently creates all 
  quadrilateral faces with 0 bentness.
  
Due to symmetry, all methods produce the same results for parameters (1,0), 
(1,1), and (2,0). With :math:`\Delta(3,0)`, gc and naive slerp (with any k 
value) produce the same results.
  
[Altschuler]_ suggests (although doesn't prove) that the closer a geodesic
sphere is to being class I, the lower its Thomson energy will be. This seems
to be the case for most other measurements: Class I operations produces lower 
Thomson energy, edge length ratios, differences in aspect ratio, 
face area ratios, and so forth.

As a side note, often the geometries produced by the GC operation on a 
triangle-faced seed are convex polyhedra. In that case, the stitching 
operation outlined in the last section can be replaced with the operation of 
taking the convex hull of vertices. Almost all GC operations on an icosahedron 
or octahedron under all methods and reasonable values for `k` produce a 
convex polyhedron. 
