Towards subdivision of mixed grids
==================================

Conway polyhedron operations
----------------------------
John Conway developed a system of operators for dividing and transforming polyhedra into other polyhedra: George Hart expanded and promoted it. While useful, there are quite a few operators (Wikipedia lists almost 40, not including composite operators), and their definitions are often arbitrary without much apparent structure among them. One is reminded of the (possibly apocryphal) quote from Enrico Fermi: "If I could remember the names of all these particles, I'd be a botanist." The SGS system doesn't succeed in subsuming all the Conway operators, but it does provide structure to some of them.

* SGS operator - Conway operator / dual Conway operator
* :math:`\Box(1,1)` - j (join) / a (ambo)
* :math:`\Box(2,0)` - o (ortho) / e (expand)
* :math:`\Delta(1,1)` - n (needle) / z (zip) [neither in Hart's list]
* :math:`\Delta(2,0)` - u (subdivide) / c (chamfer) [neither in Hart's list]
* :math:`\Delta(2,1)` - v (volute) / w (whirl) [neither in Hart's list]

Simultaneous Delta and Box operators on a mixed triangle-quad grid
------------------------------------------------------------------
A triangular subdivision :math:`\Delta(a,b)`  and a quad subdivision
:math:`\Box(c,d)` are ''compatible'' if:

* The same number of vertices lie along the base edge. This is equivalent to
  :math:`\gcd(a, b) - 1 = \gcd(c, d) - 1`
* The same number of edges cross the base edge, excluding ones that meet a
  vertex at the edge. This is equivalent to :math:`2a + 2b - 3gcd(a,b) =
  c + d - 2gcd(c, d)`. Since the above requires the `gcd` s to be equal,
  :math:`c + d + g = 2(a + b)` where :math:`g = gcd(a, b) = gcd(c, d)`. As a
  consequence, for there to be a compatible :math:`\Delta` for
  :math:`\Box(c,d)`, :math:`c + d + g`
  must be even. This is true except when c and d are both odd.
* The same number of edges lie on the base edge. (only affects Class I. Under
  the usual definition of gcd, the above equation covers this.)
* Something about edges that cross the base edge not creating repeated edges

Therefore:

#. :math:`\Box(1,0)` is compatible with :math:`\Delta(1,0)`
#. Iff :math:`\Box(a,b)` is compatible with :math:`\Delta(c,d)`, then
   :math:`\Box(na,nb)` is compatible with :math:`\Delta(nc,nd)`
   for all positive integers n.
#. By 1 and 2, :math:`\Box(n,0)` is compatible with :math:`\Delta(n,0)` 
   (Class I)
#. Iff :math:`\Box(a,b)` is compatible with :math:`\Delta(c,d)`, then
   :math:`\Box(b,a)` is compatible with :math:`\Delta(d,c)` (symmetry).
#. :math:`\Box(1,2)` is compatible with :math:`\Delta(1,1)` (by inspection)
#. By 5 and 2, :math:`\Delta(n,n)` is compatible with :math:`\Box(n,2n)`
   and :math:`\Box(2n,n)` (Class II triangles)
#. :math:`\Box(1,1)` has 2 vertices lying on the base edge, no edges crossing
   the base edge, and no edges lying on the base edge. There is no compatible
   :math:`\Delta`.
#. By 7 and 2, :math:`\Box(n,n)` is not compatible with any :math:`\Delta`.
   (Class II squares)
#. If :math:`\Box(a_1,b_1)` is compatible with :math:`\Delta(c_1,d_1)`, and
   :math:`\Box(a_2,b_2)` is compatible with :math:`\Delta(c_2,d_2)`, it does 
   not follow that the composed subdivision :math:`\Box(a_1,b_1)\Box(a_2,b_2)` 
   is compatible with :math:`\Delta(c_1,d_1)\Delta(c_2,d_2)`.
   :math:`\Delta(1,1)\Delta(1,1) = \Delta(3,0)`;
   :math:`\Box(1,2)\Box(1,2) = \Box(4,3)`;
   :math:`\Box(1,2)\Box(2,1) = \Box(5,0)`.
