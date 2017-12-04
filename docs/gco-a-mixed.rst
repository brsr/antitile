Appendix: Towards Goldberg-Coxeter operation on mixed polyhedra
===============================================================

The methods discussed up to now only handle polyhedra with a single type of
face, either 3- or 4-sided. This section discusses the extension
of the Goldberg-Coxeter operation to polyhedra with arbitrary faces.

Triangulation
-------------
The easy way to apply GC operations to non-triangular faces is to triangulate
the faces. This has the drawback that it can destroy symmetry if done naively:
for instance, a square can easily be split into two triangles, but it no longer
has D4 symmetry. Often triangulation is good enough for many applications.

Conway polyhedron operations
----------------------------
    "If I could remember the names of all these particles, I'd be a botanist."
    -Enrico Fermi, allegedly

John Conway developed a system of operators for dividing and transforming
polyhedra into other polyhedra: George Hart expanded and promoted it, and
antiprism implements even more beyond that. While useful, there are quite a
few operators (Wikipedia lists almost 40, not including composite operators),
and their definitions are often arbitrary without much apparent structure
among them. [Brinkmann]_ has commented on this situation.

The Goldberg-Coxeter operators correspond to some Conway operators
or their duals.

=================================== ===== ===== ===== =====
SGS operator                        ?     d?d   d?    ?d
=================================== ===== ===== ===== =====
:math:`\Box(1,1)`                   j     a     a     j
:math:`\Box(2,0) = \Box(1,1)`       o=jj  e=aa  aj    ja
:math:`\Delta(1,1)`                 n     z     t     k
:math:`\Delta(2,0)`                 u     c     du    dc
:math:`\Delta(2,1)`                 v     w     dv    dw
:math:`\Delta(3,0) = \Delta^2(1,1)` kt=nn tk=zz tn    nk
=================================== ===== ===== ===== =====

The :math:`\Delta` operators as defined only operate on triangles, and the
:math:`\Box` operators on quadrilaterals, while the Conway operators work on
all faces. The :math:`\Delta` and :math:`\Box` operators can be extended to
arbitrary faces by exploiting their rotational symmetry.

Simultaneous Delta and Box operators on a mixed triangle-quad grid
------------------------------------------------------------------
Rather than extending a :math:`\Delta` or :math:`\Box` operator to work on
other faces, it's possible to divide a polyhedron by applying a certain
:math:`\Delta` operator to the triangle faces, a certain :math:`\Box`
operator to the quadrilateral faces, and dealing with what happens over the
edges of the original polyhedron. 

Two Conway operators are ''compatible'' if:

#. The same number of vertices lie along the base edge.
#. The same number of edges cross the base edge, excluding ones that meet a
   vertex at the edge.
#. The same number of edges lie on the base edge.
#. Edges that cross the base edge do not create digons (faces with 2 edges).

For the operators :math:`\Delta(a,b)` and :math:`\Box(c,d)`, the first 3 
requirements can be distilled to these equations:

* :math:`\gcd(a, b) = \gcd(c, d) = g`
* :math:`c + d + g = 2(a + b)`. This equation has no solution if `c=d` or if
  `c` and `d` are both odd numbers.

The 4th requirement is more complicated.

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
