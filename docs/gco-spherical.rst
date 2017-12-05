Goldberg-Coxeter Operation on spherical polyhedra
=================================================
In this section, we'll use :math:`\mathbf v^*` to denote a pre-normalized 
vector: :math:`\hat{\mathbf v} = \frac{\mathbf v^*}{\|\mathbf v^*\|}`

Gnomonic
--------
The gnomonic projection was known to the ancient Greeks, and is the simplest 
of the transformations listed here. It has the nice property that all lines in 
Euclidean space are transformed into great circles on the sphere: that is, 
geodesics stay geodesics. This is in fact the motivation for the name 
"geodesic dome": Buckminster Fuller used this projection to 
project triangles on the sphere.

In general, the gnomonic projection is defined as:

* To sphere: :math:`\hat{\mathbf v} = \frac{\mathbf p}{\|\mathbf p\|}`
* From sphere: :math:`\mathbf p = \frac{r\hat{\mathbf v}}
  {\hat{\mathbf n} \cdot \hat{\mathbf v}}`
  
where :math:`\mathbf p` is a point on a plane given in Hessian normal
form by :math:`\hat{\mathbf n} \cdot \mathbf p = r`. Projection from Euclidean 
space to the sphere is literally just normalizing the vector. 

For the Goldberg-Coxeter operation, this amounts to just normalizing 
the vectors produced by the coordinate form:

.. math::
   \mathbf v^* = \beta_1 \mathbf v_1 + \beta_2 \mathbf v_2 + \beta_3 \mathbf v_3 
   
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
:math:`\beta^\prime_i = \frac{\beta_i}
{\beta_1 \mathbf v_1 + \beta_2 \mathbf v_2 + \beta_3 \mathbf v_3}`. On the 
interior of the triangle, :math:`\sum \beta^\prime_i > 1`.

Spherical areal
---------------
This method applies to triangles only. The above section discussed generalized
barycentric coordinates, which removes the restriction that 
:math:`\sum \beta_i = 1`. Instead we look to the relation between barycentric
coordinates and area; :math:`\beta_i` is the proportion of area in the 
triangle that is opposite the vertex :math:`v_i`. Let :math:`\Omega` be the 
spherical area (solid angle) of the spherical triangle and 
:math:`\Omega_i = \beta_i\Omega` be the area of the smaller triangle 
opposite vertex :math:`v_i`.

The equation to find :math:`\hat{\mathbf v}` is much more complicated than
that for barycentric coordinates. Let 
:math:`h_i = \sin\Omega_i
\left(1+\mathbf{\hat{v}}_{i-1} \cdot \mathbf{\hat{v}}_{i+1}\right)`, 
and
:math:`\mathbf g_{i} = \left(1+\cos \Omega_{i}\right) 
\mathbf{\hat{v}}_{i-1} \times \mathbf{\hat{v}}_{i+1} - 
\sin\Omega_{i}\left(\mathbf{\hat{v}}_{i-1} + \mathbf{\hat{v}}_{i+1}\right)`
where the subscripts loop around: 0 should be interpreted as 3 and 4 should be 
interpreted as 1. Then 

.. math::
   \mathbf G = \begin{bmatrix} \mathbf g_1 & \mathbf g_2 & \mathbf g_3 \end{bmatrix}

.. math::
   \mathbf h = \begin{bmatrix} h_1  & h_2 & h_3  \end{bmatrix}^T
   
such that :math:`\mathbf G \hat{\mathbf v} = \mathbf h`. To clarify, 
:math:`\mathbf G` is the 3x3 matrix where the ith column is 
:math:`\mathbf g_i`, and :math:`\mathbf h` is the column vector where the 
ith element is :math:`h_i`. The vector :math:`\hat{\mathbf v}` 
can be solved for using standard matrix methods.

Naive Slerp
-----------
Only applicable to equilateral faces

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

where :math:`\cos(w) = \mathbf v_i \cdot \mathbf v_{i+1}` for all :math:`i`

Projection

Great circle intersection
-------------------------
Method 2

Summary of methods
------------------

==================================== ======================== ==============
Method                               Input                    Face size
==================================== ======================== ==============
Gnomonic (method 1)                  Coordinates (either)     :math:`< 2\pi`
Spherical areal                      Barycentric coordinates  :math:`< 2\pi`
Naive Slerp                          Coordinates (either)     :math:`<=2\pi`
Great circle intersection (method 2) Linear index (either)    :math:`< 2\pi`
==================================== ======================== ==============

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
smallest n such that :math:`\Delta(a,b)\Delta(c,d) = \Delta(n,0)` for some c 
and d.

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
the closer it is to Class I, the more even it is: [Altschuler]_.
