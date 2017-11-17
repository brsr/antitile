Towards subdivision of mixed grids
==================================

Conway polyhedron operations
----------------------------
John Conway developed a system of operators for dividing and transforming polyhedra into other polyhedra: George Hart expanded and promoted it. While useful, there are quite a few operators (Wikipedia lists almost 40, not including composite operators), and their definitions are often arbitrary without much apparent structure among them. One is reminded of the (possibly apocryphal) quote from Enrico Fermi: "If I could remember the names of all these particles, I'd be a botanist." The SGS system doesn't succeed in subsuming all the Conway operators, but it does provide structure to some of them.

* SGS operator - Conway operator / dual Conway operator
* □(1,1) - j (join) / a (ambo)
* □(2,0) - o (ortho) / e (expand)
* △(1,1) - n (needle) / z (zip) [neither in Hart's list]
* △(2,0) - u (subdivide) / c (chamfer) [neither in Hart's list]
* △(2,1) - v (volute) / w (whirl) [neither in Hart's list]

Simultaneous Delta and Box operators on a mixed triangle-quad grid
------------------------------------------------------------------
A triangular subdivision △(a,b)  and a quad subdivision □(c,d) are ''compatible'' if:
* The same number of vertices lie along the base edge. This is equivalent to $gcd(a, b) - 1 = gcd(c, d) - 1$
* The same number of edges cross the base edge, excluding ones that meet a vertex at the edge. This is equivalent to $2a + 2b - 3gcd(a,b) = c + d - 2gcd(c, d)$. Since the above requires the `gcd` s to be equal, $c + d + g = 2(a + b)$ where $g = gcd(a, b) = gcd(c, d)$. As a consequence, for there to be a compatible △ for □(c,d), $c + d + g$ must be even. This is true except when c and d are both odd.
* The same number of edges lie on the base edge. (only affects Class I. Under the usual definition of gcd, the above equation covers this.)
* Something about edges that cross the base edge not creating repeated edges

Therefore:
# □(1,0) is compatible with △(1,0)
# Iff □(a,b) is compatible with △(c,d), then □(na,nb) is compatible with △(nc,nd) for all positive integers n.
# By 1 and 2, □(n,0) is compatible with △(n,0) (Class I)
# Iff □(a,b) is compatible with △(c,d), then □(b,a) is compatible with △(d,c) (symmetry).
# □(1,2) is compatible with △(1,1) (by inspection)
# By 5 and 2, △(n,n) is compatible with □(n,2n) and □(2n,n) (Class II triangles)
# □(1,1) has 2 vertices lying on the base edge, no edges crossing the base edge, and no edges lying on the base edge. There is no compatible △.
# By 7 and 2, □(n,n) is not compatible with any △. (Class II squares)
# If □(a_1,b_1) is compatible with △(c_1,d_1), and □(a_2,b_2) is compatible with △(c_2,d_2), it does not follow that the composed subdivision □(a_1,b_1)□(a_2,b_2) is compatible with △(c_1,d_1)△(c_2,d_2). △(1,1)△(1,1) = △(3,0); □(1,2)□(1,2) = □(4,3); □(1,2)□(2,1) = □(5,0)