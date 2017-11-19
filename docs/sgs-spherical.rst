Transformations between Euclidean and spherical triangles and quadrilaterals
============================================================================

Review of spherical geometry
----------------------------

Gnomonic
--------

Spherical areal
---------------

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
