Goldberg-Coxeter Operation on arbitrary tilings and polyhedra
=============================================================

Coordinate form
---------------
For our purposes, triangles are specified by barycentric coordinates 
:math:`\beta_i` where :math:`\beta_1 + \beta_2 + \beta_3 = 1`.

:math:\mathbf v = \beta_1\mathbf v_1+\beta_2\mathbf v_2+\beta_3\mathbf v_3`

Quadrilaterals are instead specified by what we'll call "xy coordinates"
where :math:`x` and :math:`y` are between 0 and 1. 

:math:`\mathbf v_1 + (\mathbf v_2-\mathbf v_1) x + (\mathbf v_4-\mathbf v_1) y
+ (\mathbf v_1-\mathbf v_2+\mathbf v_3-\mathbf v_4)xy`

If the quadrilateral is a skew quadrilateral, x and y smoothly parameterize a 
surface over that skew quadrilateral.

:math:`\mathbf \hat v = \frac{\mathbf v^*}{\|\mathbf v^*\|}`

Linear index form
-----------------
Another way to think about the breakdown structure is to count lines that 
cross into the breakdown polygon.

For :math:`\Delta(a,b)`, draw an equilateral triangle in the plane. Along one 
edge, mark `a+b+1` evenly spaced vertices (including the base vertices). Going 
clockwise, mark `a+1` points along the next edge, then `b+1` along the last. 
Number the vertices along the first side, then along the remaining two sides 
sequentially (counting the point shared between the other two sides only 
once). Draw lines between points with the same number. Then take those lines 
and make copies rotated 120° and 240° about the center. The vertices of the 
breakdown structure correspond to the points where 3 lines intersect (not 
necessarily the marked points).

For :math`\Box(a,b)`, draw a square in the plane.