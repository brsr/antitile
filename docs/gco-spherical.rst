Goldberg-Coxeter Operation on spherical polyhedra
=================================================

Canonical form
--------------
[Hart1997]_

Gnomonic
--------
The gnomonic projection was known to the ancient Greeks, and is the simplest 
of the transformations listed here. It has the nice property that all lines in 
Euclidean space are transformed into great circles on the sphere: that is, 
geodesics stay geodesics. This is in fact the motivation for the name 
"geodesic dome": Buckminster Fuller was originally using this projection to 
project triangles on the sphere.

In general, the gnomonic projection is defined as:

* To sphere: :math:`\mathbf \hat{v} = \frac{\mathbf p}{\|\mathbf p\|}`
* From sphere: :math:`\mathbf p = \frac{r\mathbf \hat{v}}
  {\mathbf \hat{n} \cdot \mathbf\hat{v}}`
  
where :math:`\mathbf p` is a point on a plane given in Hessian normal
form by :math:`\mathbf \hat{n} \cdot \mathbf p = r`. Projection from Euclidean 
space to the sphere is literally just normalizing the vector. 

where :math:`\beta_i` are (planar) barycentric coordinates

.. math::
   \mathbf v^* = \beta_1 \mathbf v_1 + \beta_2 \mathbf v_2 + \beta_3 \mathbf v_3, 
   

Spherical areal
---------------
Triangles only

:math:`\Omega_i = \beta_i\Omega_{total}`, 
:math:`h_i = \sin\Omega_i\left(1+\mathbf v_{i-1}\cdot\mathbf v_{i+1}\right)`, 
:math:`\mathbf g_{i} = \left(1+\cos \Omega_{i}\right) \mathbf v_{i-1} \times 
\mathbf v_{i+1} - \sin\Omega_{i}\left(\mathbf v_{i-1} + \mathbf v_{i+1}\right)`
where the subscripts loop around: 0 should be interpreted as 3 and 4 should be 
interpreted as 1. Then 

.. math::
   \mathbf G = \begin{bmatrix} \mathbf g_1 & \mathbf g_2 & \mathbf g_3 \end{bmatrix}
   \mathbf h = \begin{bmatrix} h_1  & h_2 & h_3  \end{bmatrix}^T
   
such that :math:`\mathbf G \mathbf v = \mathbf h` To clarify, 
:math:`\mathbf G` is the 3x3 matrix where the `i`th column is 
:math:`\mathbf g_i`, and :math:`\mathbf h` is the column vector where the 
`i`th element is :math:`h_i`. The vector :math:`\mathbf v` can be solved for 
using standard matrix methods.

Naive Slerp
-----------
Triangle:
:math:`\mathbf v^* = \sum_{i=1}^3\frac{\sin(w\beta_i)}{\sin(w)}  \mathbf v_i`

Quadrilateral 1:
:math:`\mathbf v^* = \sum_{i=1}^4\frac{\sin(w\gamma_i)}{\sin(w)}  \mathbf v_i`
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

Optimization
------------
the closer it is to Class I, the more even it is: [Altschuler]_.
