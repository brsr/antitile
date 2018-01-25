Appendix: Towards Goldberg-Coxeter operation on mixed polyhedra
===============================================================

The methods discussed up to now only handle polyhedra with a single type of
face, either 3- or 4-sided. This section discusses the extension of the
Goldberg-Coxeter operation to polyhedra with mixed triangle and quad faces.
It's possible to divide a polyhedron by applying 
a certain :math:`\Delta` operator to the triangle faces,
a certain :math:`\Box` operator to the quadrilateral faces,
and dealing with what happens over the edges of the original polyhedron.

Two operators are ''compatible'' if:

#. The same number of vertices lie along the base edge.
#. The same number of edges cross the base edge, excluding ones that meet a
   vertex at the edge.
#. The same number of edges lie on the base edge.
#. Edges that cross the base edge do not create digons (faces with 2 edges).

For the operators :math:`\Delta(a,b)` and :math:`\Box(c,d)`, the first 3
requirements can be distilled to these equations:

* :math:`\gcd(a, b) = \gcd(c, d) = g`
* :math:`c + d + g = 2(a + b)`. This equation has no solution if
  `c` and `d` are both odd numbers.

The 4th requirement is more complicated to put in equational form.

Consequences of these rules include:

#. :math:`\Box(1,0)` is compatible with :math:`\Delta(1,0)`
#. Iff :math:`\Box(a,b)` is compatible with :math:`\Delta(c,d)`, then
   :math:`\Box(na,nb)` is compatible with :math:`\Delta(nc,nd)`
   for all positive integers n.
#. :math:`\Box(n,0)` is compatible with :math:`\Delta(n,0)` (Class I)
#. Iff :math:`\Box(a,b)` is compatible with :math:`\Delta(c,d)`, then
   :math:`\Box(b,a)` is compatible with :math:`\Delta(d,c)`.
#. :math:`\Delta(n,n)` is compatible with :math:`\Box(n,2n)`
   and :math:`\Box(2n,n)` (Class II triangles)
#. :math:`\Box(n,n)` is not compatible with any :math:`\Delta`.
   (Class II squares)
#. If :math:`\Box(a_1,b_1)` is compatible with :math:`\Delta(c_1,d_1)`, and
   :math:`\Box(a_2,b_2)` is compatible with :math:`\Delta(c_2,d_2)`, it does
   not follow that the composed subdivision :math:`\Box(a_1,b_1)\Box(a_2,b_2)`
   is compatible with :math:`\Delta(c_1,d_1)\Delta(c_2,d_2)`.
   :math:`\Delta(1,1)\Delta(1,1) = \Delta(3,0)`;
   :math:`\Box(1,2)\Box(1,2) = \Box(4,3)`;
   :math:`\Box(1,2)\Box(2,1) = \Box(5,0)`.
