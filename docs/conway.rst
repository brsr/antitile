Notes on operations on polyhedra
================================

Operations on polyhedra to produce other polyhedra date back as far as Kepler.
Conway defined a set of operations that could be performed on the Platonic
solids to obtain the Archimedean and Catalan solids, and others added operators
after him. Initially there was not much theory supporting operations on
polyhedra, but [Brinkmann]_'s paper provides a framework. This text
is an attempt to use Brinkmann's work to find ways to quantify, analyze,
and expand these operators.

Preliminaries
-------------
This assumes some familiarity with basic graph theory, solid geometry, and
Conway operators. See [HartConway]_ for a basic overview of Conway operators, or
better still, spend some time playing with [Polyhedronisme]_ (a web app) or
``conway`` in [Antiprism]_. Some paper to doodle on is helpful too. In general,
this text uses the same terms as ``conway``. Also beware that the term "Conway
operator" is not well-defined; it can refer to any operation on a polyhedron,
Conway's original set, operations that retain the symmetry of the seed
polyhedron like Conway's operators, etc. depending on your source.

Until we get to the extensions section, we'll limit ourselves to operations on
convex polyhedra. As a sphere looks like a plane if you zoom in enough, it can
be useful to visualize operators on a planar grid,
which can be thought of as a zoomed-in section of a large convex polyhedra.

Faces with `k` sides may be called `k`-degree faces, by analogy with `k`-degree
vertices.

Chamber structure
-----------------
.. _facechambers:
.. figure:: triangle_chambers.svg
   :align: right
   :figwidth: 25%

   Chambers of a triangular face.

.. _edgechambers:
.. figure:: edge_chambers.svg
   :align: right
   :figwidth: 25%

   Chambers adjacent to an edge: what this text calls "chamber structure"

[Brinkmann]_ et al. observed that Conway's operators, and operators like it,
can be described in terms of chambers. Each face may be divided into chambers
by identifying the face center and drawing lines from there to each vertex and
edge midpoint, as in :numref:`facechambers`. Similarly, each vertex of degree
`n` is surrounded by `n` white and `n` grey chambers. Each edge has a white and
grey chamber on each side of the edge, as shown in :numref:`edgechambers`. The
operator may then be specified by a structure of vertices and edges within
those chambers, possibly with edges crossing from one chamber to another.

[Brinkmann]_ et al. note that for all operators that can be expressed in terms
of these chambers, the number of edges in the result polyhedron are an integer
multiple of those in the seed polyhedron. They call this the inflation rate,
and we'll denote it :math:`g`. It turns out that an edge-focused view of
these operators is fruitful: we can view it as replacing each edge and its
surroundings with a structure like that in :numref:`edgechambers`, possibly
rotated or stretched, but maintaining orientation with respect to the
polyhedron. Therefore (and since Brinkmann et al. don't actually introduce an
overarching term for these operators) we'll call them edge replacement
operators, or EROs. If an operator's grey chamber is a reflection of the white
chamber, we call it achiral: otherwise the operator is chiral. (Brinkmann et al.
call these local symmetry-preserving operations (LSP) and local operations
that preserve orientation-preserving symmetries (LOPSP), respectively.)

The chamber structure of the composed operator `xy` can be drawn by applying `x`
to the edges of the chamber structure of `y`. In particular, for a given
ERO `x`, the chamber structure of `xd` is simply the chamber structure
of `x` rotated one quarter turn.

There is some freedom in where vertices are placed within the chambers.
This is more apparent with chiral EROs. Often an ERO is drawn
so that most of the vertices lie on the seed edge, but this is not necessary.
For instance, see :numref:`propeller` for two topologically equivalent ways to
draw George Hart's propeller operator (see [HartPropeller]_).

.. _propeller:
.. list-table:: Chambers for the operator `p` (Propeller)

   * - .. image:: edge_chambers_propeller.svg
     - .. image:: edge_chambers_propeller-square_grid.svg
   * - George Hart's original drawing
     - Drawing emphasizing relationship with a square grid

Particular sets of edge-replacement operators
---------------------------------------------
Conway's original set of operations is denoted with the letters `abdegjkmost`.
Some of these are reducible: `e=aa`, `o=jj`, `m=kj`, and `b=ta`.

The Goldberg-Coxeter operations :math:`\Box_{a,b}` and :math:`\Delta_{a,b}`
described in :ref:`Goldberg-Coxeter Operations on Polyhedra and Tilings` can be
fairly simply extended to a ERO. In terms of the complex plane used in
:ref:`Master Polygons`, the chamber structure of :math:`\Box_{a,b}` is the
section contained in the quadrilateral :math:`0, x(1-i)/2, x, x(1+i)/2` of a
square grid on the Gaussian integers, where :math:`x=a+bi`. For
:math:`\Delta_{a,b}`, the chamber structure is the quadrilateral section
:math:`0, x(2-u)/3, x, x(1+u)/3` of a triangular grid on the Eisenstein
integers, where :math:`x=a+bu` and :math:`u=\exp(i \pi /3)`.
GC operators have an invariant `T`, the "trianglation number",
which is identical to the inflation factor `g`.

* :math:`\Box_{a,b}`: :math:`g = T = a^2 + b^2`
* :math:`\Delta_{a,b}`: :math:`g = T = a^2 + ab + b^2`

All of the nice qualities of GC operators carry over to this extension; for
instance, :math:`\Box_{a,b}` operators commute with each other, as do
:math:`\Delta_{a,b}` operators, and the operators can be decomposed in relation
to the Gaussian or Eisenstein integers respectively. Except for `g` and `s`,
all of Conway's original operators are GC operations,
related by duality, or compositions of GC operators or their duals.

The simplest operators (aside from the identity) are :math:`\Box_{1,1} = j` and
:math:`\Delta_{1,1} = n = kd`. One useful relation is that if
:math:`a=b \mod 3`, :math:`\Delta_{a,b} = n \Delta_{(2a+b)/3, (b-a)/3}`, and if
:math:`a=b \mod 2`, :math:`\Box_{a,b} = j \Box_{(a+b)/2,(b-a)/2}`.
(These formula may result in negative values, which should be interpreted as
per :ref:`Master Polygons`.)

Alternating operators
---------------------
.. _facealtchambers:
.. figure:: square_alternating_chambers.svg
   :align: right
   :figwidth: 25%

   Alternating chambers of a quadrilateral face.

.. _edgealtchambers:
.. figure:: edge_chambers_alternating.svg
   :align: right
   :figwidth: 25%

   Alternating chambers adjacent to an edge.

.. _semi:
.. figure:: edge_chambers_alternating_semi.svg
   :align: right
   :figwidth: 25%

   Alternating chambers of the Coxeter semi operator (without digon reduction)

In [Coxeter8]_ (specifically section 8.6), Coxeter defines an alternation
operation `h` on regular polyhedra with only even-sided faces. (He actually
defines it on general polytopes, but let's not complicate things by considering
higher dimensions.) Each face is replaced
with a face with half as many sides, and alternate vertices are either retained
as part of the faces or converted into vertices with number of sides equal to
the degree of the seed vertex. (He also defines a snub operation in section 8.4,
different from the `s` snub Conway defined, that is equivalent to `ht`.) The
alternation operation converts quadrilateral faces into digons. Usually the
digons are converted into edges, but for now, let digons be digons.

This motivates the definition of "alternating operators" and an "alternating
chamber" structure, as depicted in :numref:`facealtchambers` and
:numref:`edgealtchambers`. Like earlier, we can think of this as replacing each
edge with :numref:`edgealtchambers`, stretched or rotated but maintaining
orientation with respect to the polyhedron, so we can call these operators AEROs
(alternating EROs) for short. This structure is only applicable to polyhedra
with even-sided faces. The dual operators of those are applicable to polyhedra
with even-degree vertices, and should be visualized as having chambers on the
left and right rather than top and bottom. Like EROs, the chamber
structure of `xd` is that of `x` rotated a quarter turn; but now, the direction
of rotation matters, and depends on how the alternating vertices (or faces) of
the underlying polyhedron are specified. For the sake of simplicity, we'll only
look at AEROs on even-sided faces (vertex-AEROs, or VAEROs) instead of on
even-degree vertices (face-AEROs, or FAEROs).

VAEROs depend on the ability to partition vertices into two disjoint sets, none
of which are adjacent to a vertex in the same set; i.e. it applies to bipartite
graphs. We'll denote those sets as :math:`+` and :math:`-`. By basic graph
theory, planar bipartite graphs have faces of even degree. However, this does
not mean that the two sets of vertexes have the same size, let alone that the
sets of vertices of a given degree will have a convenient partition. The cube
and many other small even-faced polyhedra do partition into two equal sets of
vertices, so beware that examining simple, highly-symmetric polyhedra can be
misleading. (A section on AEROs briefly appeared on the Wikipedia page for
Conway operators. It made some errors that seemed to result from assuming
that the partitions were of equal size.)

Strictly, since AEROs map polyhedra with even-sided faces to arbitrary
polyhedra, they are not operators in the strict mathematical sense. (In
particular, since AEROs do not necessarily produce even-sided faces or
even-degree vertices, they cannot be composed together arbitrarily.) However,
calling them "transformations" instead felt awkward, since the term "operator"
is so commonly used. You can call them AERTs, VAERTs, and FAERTs instead if
you like.

Digons and degree-2 vertices are an unavoidable fact of certain VAEROs,
particularly on quadrilateral faces. Two important special cases are where
the seed polyhedron has only quadrilateral faces, and when it has only faces of
degree 6 or more (although the latter case only appears in infinite tilings).
In the former case, the degree-2 features can be uniformly smoothed out.
In the latter, degree-2 features are not created.

Other Operators
---------------
There are some important operations on polyhedra that don't fix into the
edge-replacement schema.

* `r`, the reflection operator. This produces the mirror image of the
  polyhedron. If an operator `x` is chiral, `rxr` is its chiral pair.
* `$`, the smoothing operator (newly defined here). This operator smooths
  degree-2 vertices and digons, as produced by some AEROs. This operator is
  recursive, and will smooth features until there are no degree-2 features
  left to smooth. For instance, two vertices may be
  connected by one edge and another edge split by a degree-2 vertex; one
  smoothing iteration would smooth that degree-2 vertex into a single edge,
  creating a digon, and the next would reduce the digon into a single edge.
* `@`, the alternation operator (newly defined here).
  This operator just exchanges the :math:`+` and :math:`-` partitions.
  Applied to an operator, it reflects its chamber structure horizontally.

Representations of operators
----------------------------
In abstract algebraic terms, EROs form a monoid: a group without an inverse, or
a semigroup with an identity element. Let :math:`[v,e,f]` be the count of
vertices, edges, and faces of the seed,
and :math:`v_i` and :math:`f_i` be the count of vertices/faces of degree
:math:`i` such that :math:`\sum v_i = v` and :math:`\sum f_i = f`.
There is a series of monoids and homomorphisms between the monoids, as so:

* ERO `x` (acts on polyhedra)
* Infinite-dimensional linear operator :math:`L_x` (acts on :math:`v_i, e, f_i`)
* 3x3 matrix :math:`M_x` (acts on :math:`[v,e,f]`)
* Inflation factor `g` (acts on :math:`e`) and operator outline

AEROs do not form a monoid (since in general they cannot be composed together)
but do admit a similar representation. For VAEROs, the count of vertices of
degree :math:`i` in the :math:`+` partition are denoted :math:`v^+_i` and those
in the :math:`-` partition as :math:`v^-_i`. :math:`\sum v^+_i = v^+`, and
similarly for :math:`-`. :math:`v^+_i + v^-_i = v_i`, and :math:`v^+ + v^- = v`.
Partitions of :math:`f` for FAEROs are denoted similarly.

Each bullet will be handled in turn.

The action of an ERO on the vertices of degree :math:`i`, edges, and faces with
:math:`i` sides can be described with an infinite linear operator :math:`L_x`.
This operator can be determined by counting elements off the chamber structure.
Step by step:

* Seed vertices are either retained or converted into faces centered on that
  vertex. (Other options are precluded by symmetry). Let :math:`a = 1` if the
  seed vertices are retained, and 0 otherwise. Also, the degree of the vertex
  or face is either the same as the seed vertex, or a multiple of it;
  let :math:`k` be that multiple.
* Seed face centers are either retained (possibly of in a smaller face) or
  converted into vertices. (Again, other options are precluded by symmetry).
  Let :math:`c = 0` if the seed faces are retained, and 1 otherwise. Let
  :math:`\ell` serve a similar role as :math:`k` above: the degree of the vertex
  or face corresponding to the seed face center is :math:`k` times the degree of
  the seed vertex.
* Except for the faces or vertices corresponding to the seed vertices and face
  centers, the added elements are in proportion to to the number of edges in the
  seed. :math:`g` is the count of added edges (the edge multiplier or inflation
  rate), :math:`b_i` is the number of vertices of degree :math:`i` added, and
  :math:`b'_i` is the number of faces of degree :math:`i` added.

Count elements lying on or crossing the outer edge of the chamber structure as
half. It may help to draw an adjacent chamber, particularly when determining
the number of sides on a face. The result of the counting process can be
described in the following operator form;
variables in capital letters are the result of the operator.

.. math::
   E &= ge

   V_i &= a v_{i/k} + e b_i + c f_{i/\ell}

   F_i &= a' v_{i/k} + e b'_i + c' f_{i/\ell}

where :math:`a`, :math:`a'`, `c`, and :math:`c'` are either 0 or 1, `g` is a
positive integer, all :math:`b_i` and :math:`b'_i` are nonnegative integers, and
:math:`k` and :math:`\ell` are positive integers. The subscripted values like
:math:`v_{i/k}` should be interpreted as 0 if :math:`i/k` is not an integer.

The only alteration needed to accommodate VAEROs is that the action on seed
vertices may be different depending on which partition they are in. (Counting
elements may be more complicated: it's possible to have an edge pass through
one chamber without meeting any vertices.)

.. math::
   E &= ge

   V_i &= a^+ v^+_{i/k^+} + a^- v^-_{i/k^-} + e b_i + c f_{i/\ell}

   F_i &= a'^+ v^+_{i/k^+} + a'^- v^-_{i/k^-} + e b'_i + c' f_{i/\ell}

:math:`a^+`, :math:`a^-`, :math:`a'^+`,  and :math:`a'^-` are either 0 or 1.
:math:`k^+`, :math:`k^-` are positive integers and :math:`\ell` may take values
in :math:`\mathbb{N}/2 = \{1/2, 1, 3/2, 2, ...\}`.

We'll refer to :math:`g, a, a', b_i, b'_i, c, c', k, \ell` as the invariants of
an ERO, and :math:`g, a^+, a'^+, b_i, b'_i, c, c', k^+, k^-, \ell` as the
invariants of a VAERO. If :math:`a^+ = a^-` both may be written as :math:`a`,
and similarly for :math:`a'` and :math:`k`. :math:`a'` and :math:`c'` may be
omitted since they can be calculated from :math:`a` and :math:`c`.
FAEROs would be described correspondingly.

Explicitly the composition of two EROs `xy` can be described as so.
Let :math:`g, a, a', b_i, b'_i, c, c' k, \ell` be the invariants for :math:`L_y`;
:math:`G, A, A', B_i, B'_i, C, C', K, L` for :math:`L_x`; and
:math:`\gamma, \alpha, \alpha', \beta_i, \beta'_i, \sigma, \sigma',
\kappa, \lambda` for :math:`L_{xy}`:

.. math::
   \gamma &= Gg

   \alpha &= Aa + Ca'

   \beta_i &= A b_{i/K} + g B_i + C b'_{i/L}

   \beta'_i &= A' b_{i/K} + g B'_i + C' b'_{i/L}

   \sigma &= Ac + Cc'

.. math::
   \kappa &= \left\{
    \begin{array}{ll}
      Kk & if a=1\\
      Lk & if a=0
    \end{array}
   \right.

   \lambda &= \left\{
    \begin{array}{ll}
      K \ell & if c=1\\
      L \ell & if c=0
    \end{array}
   \right.

Under the constraint that an ERO preserves the Euler characteristic,
it can be shown that :math:`a + a' = 1`, :math:`c + c' = 1`, and
:math:`g= b + b' + 1` where :math:`\sum b_i = b` and :math:`\sum b'_i = b'`.
For VAEROs, :math:`a^+ + a'^+ = 1` and :math:`a^- + a'^- = 1`.
Also, since :math:`b_i` and :math:`b'_i` are nonnegative integers, only a
finite number of their values can be non-zero. This makes the operator form
more manageable than the term "infinite linear operator" may suggest; in
reality, nearly all applications will only use a finite number of different
vertex and face degrees.

Applying the handshake lemma gives relations between the values for EROs:

.. math::
   2g &= 2ak + 2c\ell + \sum i b_i

   2g &= 2a'k + 2c'\ell + \sum i b'_i

or for VAEROs:

.. math::
   2g &= a^+ k^+ + a^- k^- + 2c\ell + \sum i b_i

   2g &= a'^+ k^+ + a'^- k^- + 2c'\ell + \sum i b'_i

For EROs, these relations can be manipulated into the form

.. math::
   2k + 2\ell - 4 = \sum (4-i) (b_i + b'_i),

which is interesting because it eliminates `g`, `a` and `c`,
and because it suggests that features with degree 5 or more exist
in balance with features of degree 3 (triangles and degree-3 vertices),
and that in some sense degree 4 features come "for free". The relationship
for VAEROs is the same except replace :math:`2k` with :math:`k^+ + k^-`.
(For FAEROs, replace :math:`2\ell` with :math:`\ell^+ + \ell^-`.)

With these relations, and the assumption that there are no degree 2 features
and therefore :math:`i \ge 3`, a series of inequalities can be derived for EROs:

.. math::
   g + 1 \le 2a + 3b + 2c \le 2g

   2k + 2\ell \le g + 3

   0 \le 2k + 2\ell - 4 \le b_3 + b'_3

and for VAEROs:

.. math::
   1 \le a^+ + a^- + 2b + c \le 2g

   k^+ + k^- + 2\ell \le 2g + 2

The dual ERO :math:`L_d` has the form :math:`E = e, V_i = f_i, F_i = v_i`.
With a little manipulation, it is easy to see that if :math:`L_x` has invariants
`a`, :math:`b_i`, `c`, etc, then applications of the dual operator have related
forms. :math:`L_x L_d`'s invariants exchange `a` with `c`, :math:`a'` with
:math:`c'`, and `k` with :math:`\ell`. :math:`L_d L_x`'s invariants exchange `a`
with :math:`a'`, `c` with :math:`c'`, and each :math:`b_i` with each
:math:`b'_i`. Finally, :math:`L_d L_x L_d`'s invariants exchange `a` with
:math:`c'`, and :math:`a'` with `c`, `k` with :math:`\ell`,
and each :math:`b_i` with each :math:`b'_i`.

For EROs, the matrix form :math:`M_x` can be obtained from :math:`L_x` by
summing :math:`\sum v_i = v` and :math:`\sum f_i = f`, or from counting elements
directly from the chamber structure without distinguishing between vertices and
faces of different degrees. (The conversion from :math:`L_x` to :math:`M_x` is
itself a linear operator.) The matrix takes the form:

.. math::
   \mathbf{M}_x = \begin{bmatrix}
   a & b & c \\
   0 & g & 0 \\
   a' & b' & c' \end{bmatrix}

The matrix for the identity operator `S` is just the 3x3 identity matrix.
The matrix for the dual operator is the reverse of that:

.. math::
   \mathbf{M}_d = \begin{bmatrix}
   0 & 0 & 1 \\
   0 & 1 & 0 \\
   1 & 0 & 0 \end{bmatrix}

The dual matrix operates on other matrices by mirroring the values either
horizontally or vertically.

.. math::
   \mathbf{M}_x \mathbf{M}_d = \begin{bmatrix}
   c & b & a \\
   0 & g & 0 \\
   c' & b' & a' \end{bmatrix}, \mathbf{M}_d \mathbf{M}_x  = \begin{bmatrix}
   a' & b' & c' \\
   0 & g & 0 \\
   a & b & c \end{bmatrix},
   \mathbf{M}_d \mathbf{M}_x \mathbf{M}_d = \begin{bmatrix}
   c' & b' & a' \\
   0 & g & 0 \\
   c & b & a \end{bmatrix}

VAEROs with :math:`a^+ = a^-` can also be written as a 3x3 matrix. In general,
VAEROs can be written as a 4x3 matrix mapping :math:`[v^+,v^-,e,f]` to
:math:`[v,e,f]`. FAEROs can be written as a 4x3 matrix as well, but that one
mapping :math:`[v,e,f^+,f^-]` to :math:`[v,e,f]`. Since the :math:`e` row
is zero except for the value :math:`g` in the :math:`e` column, there shouldn't
be much ambiguity.

.. math::
   \mathbf{M}_x = \begin{bmatrix}
   a^+ & a^- & b & c \\
   0 & 0 & g & 0 \\
   a'^+ & a'^- & b' & c' \end{bmatrix}

It can be seen from the composition equations that for an ERO `xy`, the
expansion factor g is the product of the g invariants for operators `x` and `y`.
It can also be seen that :math:`a, a', c, c'` form their own linear system,
a submatrix of :math:`M_x`: let
:math:`\Lambda_x = \begin{bmatrix} a & c \\ a' & c' \end{bmatrix}`,
then :math:`\Lambda_{xy} = \Lambda_x \Lambda_y`. :math:`\Lambda_x` represents
the effect of the operator on the seed faces and vertices: this can also be
represented as a drawing of those seed faces and vertices, called the "outline"
of the operator. By cofactor
expansion, :math:`\det (M_x) = g \det (\Lambda_x)`. :math:`\Lambda_x` has a
determinant of -1, 0, or 1. (In fact, :math:`\Lambda_x` has two eigenvalues, one
of which is always 1, and one of which may be -1, 0, or 1. :math:`M_x` has three
eigenvalues: two it shares with :math:`\Lambda_x`, and one is `g`.) The dual
operator has :math:`\det (M_x) = \det (\Lambda_x) = -1`, and it is easy to see
that of the four possible :math:`\Lambda_x`, the first two and last two in the
table below are related by the dual operator. With that motivation, we define the
"Type" of the operator as the absolute value of the determinant of
:math:`\Lambda_x`.

Like earlier, VAEROs with :math:`a^+ = a^-` are also associated with a 2x2
matrix :math:`\Lambda_x`. All VAEROs are associated with a 3x2 matrix
:math:`\Lambda_x = \left[\begin{array}{cc|c}a^+ & a^- & c \\ a'^+ & a'^- & c'\end{array}\right]`.
FAEROs are associated with a 3x2 matrix
:math:`\Lambda_x = \left[\begin{array}{c|cc}a & c^+ & c^- \\ a' & c'^+ & c'^-\end{array}\right]`.
To reduce ambiguity, a vertical bar is included to separate the :math:`a` values
from the :math:`c` values. VAEROs and FAEROs with :math:`a^+ \ne a^-`
can be shoehorned into the 2x2 matrix form if the matrix is allowed to have
undefined values for its entries, treated like NaN in floating-point numbers,
which is denoted :math:`?`. 3x2 matrixes don't have determinants, so the
type of a VAERO with :math:`a^+ \ne a^-` is not defined.

.. list-table:: Outlines and their matrix representation
   :header-rows: 1
   :widths: 1 3 3 3

   * - Outline
     - Kind & Type
     - 2x2 Matrix
     - 3x2 Matrix
   * - .. image:: outline_1_0.svg
     - Any - 1
     - :math:`\begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix}`
     - :math:`\left[\begin{array}{cc|c}1 & 1 & 0 \\ 0 & 0 & 1\end{array}\right]` or
       :math:`\left[\begin{array}{c|cc}1 & 0 & 0 \\ 0 & 1 & 1\end{array}\right]`
   * - .. image:: outline_0_1.svg
     - Any - 1
     - :math:`\begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix}`
     - :math:`\left[\begin{array}{cc|c}0 & 0 & 1 \\ 1 & 1 & 0\end{array}\right]` or
       :math:`\left[\begin{array}{c|cc}0 & 1 & 1 \\ 1 & 0 & 0\end{array}\right]`
   * - .. image:: outline_1_1.svg
     - Any - 0
     - :math:`\begin{bmatrix} 1 & 1 \\ 0 & 0 \end{bmatrix}`
     - :math:`\begin{bmatrix} 1 & 1 & 1 \\ 0 & 0 & 0 \end{bmatrix}`
   * - .. image:: outline_0_0.svg
     - Any - 0
     - :math:`\begin{bmatrix} 0 & 0 \\ 1 & 1 \end{bmatrix}`
     - :math:`\begin{bmatrix} 0 & 0 & 0 \\ 1 & 1 & 1 \end{bmatrix}`
   * - .. image:: outline_+_0.svg
     - VAERO
     - :math:`\begin{bmatrix} ? & 0 \\ ? & 1 \end{bmatrix}`
     - :math:`\left[\begin{array}{cc|c}1 & 0 & 0 \\ 0 & 1 & 1\end{array}\right]`
   * - .. image:: outline_-_1.svg
     - VAERO
     - :math:`\begin{bmatrix} ? & 1 \\ ? & 0 \end{bmatrix}`
     - :math:`\left[\begin{array}{cc|c}0 & 1 & 1 \\ 1 & 0 & 0\end{array}\right]`
   * - .. image:: outline_+_1.svg
     - VAERO
     - :math:`\begin{bmatrix} ? & 1 \\ ? & 0 \end{bmatrix}`
     - :math:`\left[\begin{array}{cc|c}1 & 0 & 1 \\ 0 & 1 & 0\end{array}\right]`
   * - .. image:: outline_-_0.svg
     - VAERO
     - :math:`\begin{bmatrix} ? & 0 \\ ? & 1 \end{bmatrix}`
     - :math:`\left[\begin{array}{cc|c}0 & 1 & 0 \\ 1 & 0 & 1\end{array}\right]`
   * - .. image:: outline_0_+.svg
     - FAERO
     - :math:`\begin{bmatrix} 0 & ? \\ 1 & ? \end{bmatrix}`
     - :math:`\left[\begin{array}{c|cc}0 & 1 & 0 \\ 1 & 0 & 1\end{array}\right]`
   * - .. image:: outline_1_-.svg
     - FAERO
     - :math:`\begin{bmatrix} 1 & ? \\ 0 & ? \end{bmatrix}`
     - :math:`\left[\begin{array}{c|cc}1 & 0 & 1 \\ 0 & 1 & 0\end{array}\right]`
   * - .. image:: outline_1_+.svg
     - FAERO
     - :math:`\begin{bmatrix} 1 & ? \\ 0 & ? \end{bmatrix}`
     - :math:`\left[\begin{array}{c|cc}1 & 1 & 0 \\ 0 & 0 & 1\end{array}\right]`
   * - .. image:: outline_0_-.svg
     - FAERO
     - :math:`\begin{bmatrix} 0 & ? \\ 1 & ? \end{bmatrix}`
     - :math:`\left[\begin{array}{c|cc}0 & 0 & 1 \\ 1 & 1 & 0\end{array}\right]`

The composition of EROs affects their outlines like so:

.. list-table:: ERO outline composition table
   :header-rows: 1
   :stub-columns: 1

   * -
     - .. image:: outline_1_0.svg
     - .. image:: outline_0_1.svg
     - .. image:: outline_1_1.svg
     - .. image:: outline_0_0.svg
   * - .. image:: outline_1_0.svg
     - .. image:: outline_1_0.svg
     - .. image:: outline_0_1.svg
     - .. image:: outline_1_1.svg
     - .. image:: outline_0_0.svg
   * - .. image:: outline_0_1.svg
     - .. image:: outline_0_1.svg
     - .. image:: outline_1_0.svg
     - .. image:: outline_0_0.svg
     - .. image:: outline_1_1.svg
   * - .. image:: outline_1_1.svg
     - .. image:: outline_1_1.svg
     - .. image:: outline_1_1.svg
     - .. image:: outline_1_1.svg
     - .. image:: outline_1_1.svg
   * - .. image:: outline_0_0.svg
     - .. image:: outline_0_0.svg
     - .. image:: outline_0_0.svg
     - .. image:: outline_0_0.svg
     - .. image:: outline_0_0.svg

In general, AEROs cannot be composed together, but the result of an AERO is just
another polyhedron, so any AERO can be composed with an ERO on the left.

.. list-table:: VAERO outline composition table
   :header-rows: 1
   :stub-columns: 1

   * -
     - .. image:: outline_+_0.svg
     - .. image:: outline_-_1.svg
     - .. image:: outline_-_0.svg
     - .. image:: outline_+_1.svg
   * - .. image:: outline_1_0.svg
     - .. image:: outline_+_0.svg
     - .. image:: outline_-_1.svg
     - .. image:: outline_-_0.svg
     - .. image:: outline_+_1.svg
   * - .. image:: outline_0_1.svg
     - .. image:: outline_-_1.svg
     - .. image:: outline_+_0.svg
     - .. image:: outline_+_1.svg
     - .. image:: outline_-_0.svg
   * - .. image:: outline_1_1.svg
     - .. image:: outline_1_1.svg
     - .. image:: outline_1_1.svg
     - .. image:: outline_1_1.svg
     - .. image:: outline_1_1.svg
   * - .. image:: outline_0_0.svg
     - .. image:: outline_0_0.svg
     - .. image:: outline_0_0.svg
     - .. image:: outline_0_0.svg
     - .. image:: outline_0_0.svg

.. list-table:: FAERO composition table
   :header-rows: 1
   :stub-columns: 1

   * -
     - .. image:: outline_0_+.svg
     - .. image:: outline_1_-.svg
     - .. image:: outline_0_-.svg
     - .. image:: outline_1_+.svg
   * - .. image:: outline_1_0.svg
     - .. image:: outline_0_+.svg
     - .. image:: outline_1_-.svg
     - .. image:: outline_0_-.svg
     - .. image:: outline_1_+.svg
   * - .. image:: outline_0_1.svg
     - .. image:: outline_1_-.svg
     - .. image:: outline_0_+.svg
     - .. image:: outline_1_+.svg
     - .. image:: outline_0_-.svg
   * - .. image:: outline_1_1.svg
     - .. image:: outline_1_1.svg
     - .. image:: outline_1_1.svg
     - .. image:: outline_1_1.svg
     - .. image:: outline_1_1.svg
   * - .. image:: outline_0_0.svg
     - .. image:: outline_0_0.svg
     - .. image:: outline_0_0.svg
     - .. image:: outline_0_0.svg
     - .. image:: outline_0_0.svg

For EROs, the parity of the invariants :math:`g` and :math:`b` also describe the center of the chamber structure.
In particular, an ERO with both :math:`g` and :math:`b` odd is not possible.
(This does not apply to AEROs, which have different symmetry structure.)

.. list-table:: Chamber center
   :header-rows: 1
   :stub-columns: 2

   * - :math:`g`
     - :math:`b`
     - Description
   * - Even
     - Even
     - A face with even degree lies at the center
   * - Even
     - Odd
     - A vertex with even degree lies at the center
   * - Odd
     - Even
     - An edge crosses the center
   * - Odd
     - Odd
     - Excluded by symmetry

Decomposition
-------------
An operator that cannot be expressed in terms of operators aside from `d` and
`r` is "irreducible". For instance, `k` (Kis) and `j` (Join) are irreducible
in terms of EROs, but `m` (Meta) is not (it is equal to `kj`). A polyhedron that
cannot be expressed in terms of another polyhedron and one or more EROs other
than `S` and `d` is an irreducible polyhedron. An interesting fact: the only
platonic solid that is irreducible is the tetrahedron; the others can be
expressed as some operation on the tetrahedron (`O = aT`, `C = jT`, `I = sT`,
`D = gT`). Consequently, all of the Archimedean and Catalan solids can be
expressed as some series of operators and T.


.. _waffle:
.. figure:: edge_chambers_waffle.svg
   :align: right
   :figwidth: 25%

   The waffle operator (W)

The relations defined above can be used to help reduce an operator, with some
caveats. We haven't proven that the relations given in the previous section are
sufficient to discern invariants that do or do not correspond to an actual ERO.
Furthermore, none of these homomorphisms are injections: there are certain
:math:`L_x` or :math:`M_x` that correspond to more than one EROs.
Examples for :math:`M_x` are easy to come by: where `n = kd`, :math:`M_k = M_n`.
For an example where the operators are not related by duality,
:math:`M_l = M_p`. For :math:`L_x`, :math:`L_{prp} = L_{pp}` but `prp` is not
the same as `pp` (one's chiral, one's not). For the operator depicted in
:numref:`waffle`, :math:`W \ne Wd`, but :math:`L_W = L_{Wd}`.
(This is a newly named operator, introduced in this text.) A general
counterexample would be operators with sufficiently large `g` based on
:math:`\Box_{a,b}`, with a single square face (not touching the seed vertices
or face centers) divided into two triangles:
the counts of vertices of each degree, faces of each degree, and edges would be
the same no matter which faces was chosen, but the operators would be different.

The above representations do not give us a 100% reliable way to decompose an
arbitrary operator into a sequence of operators, it does suggest a (clunky,
trial-and-error filled) heuristic to reduce an operator into two operators by
starting at the bottom of the homomorphism chain and going up.

* Determine the :math:`g` of the two operators from the factors of the
  :math:`g` of the operator to be factored.
* Determine the outline (:math:`a, a', c, c'`) of the two operators.
* Determine :math:`b, b'` for the two operators.
* Determine :math:`k, \ell, b_i, b'_i`. for the two operators.
* Figure out if the representations you've produced actually corresponds to
  an ERO.

Some facts relating to decomposition that can be derived from what we have
so far:

* If a polyhedron has a prime number of edges, it is irreducible.
* Operators where `g` is a prime number are irreducible.
* If `x=xd` or `rxr=xd`, `x` has type 0.
* If `x=dxd` or `rxr=dxd`, `x` has type 1, :math:`g` is odd, and :math:`b=b'` is even.
* If an ERO has type 1, its decomposition cannot contain any EROs of
  type 0. Correspondingly, if an ERO has type 0,
  its decomposition must contain at least one type 0 ERO.
* There are no type 1 EROs with :math:`g=2`, so therefore type 1 EROs
  with :math:`g=2p`, where p is prime, are irreducible in terms of EROs.
  (However, see the section below,
  :ref:`All EROs can be expressed with smoothing, an AERO, and the join operator`.)
* :math:`\Box_{a,b}` that correspond to the Gaussian primes, and :math:`\Delta_{a,b}`
  that correspond to the Eisenstein primes, are irreducible in terms of EROs. (Proof below.)
  As a consequence of this, there are an infinite number of irreducible EROs.

Proof of the last statement: A Gaussian integer :math:`a + bi` is prime if its square norm
:math:`a^2 + b^2` is prime or the square of a prime. In the first case, that prime has
the form :math:`p=4n+1`; in the latter, :math:`p=4n+3`. Remember that the squared norm of
the integer is just the inflation factor `g` for the corresponding operator. If `g` is prime,
the operator is irreducible. If `g` is the square of a prime, the operator :math:`\Box_{a,b}` is
type 1, specifically, :math:`\det(\Lambda_{\Box_{a,b}}) = 1`. Suppose the operator can be
decomposed into :math:`\Box_{a,b} = xy`, where `x` and `y` both have inflation factor :math:`g' = \sqrt(g)`.
Without loss of generality, assume :math:`\det(\Lambda_x) = \det(\Lambda_y) = 1`. Their matrix forms are:

.. math::
   \mathbf{M}_x \mathbf{M}_y = \begin{bmatrix}
   1 & b & 0 \\
   0 & g' & 0 \\
   0 & b' & 1 \end{bmatrix} \begin{bmatrix}
   1 & B & 0 \\
   0 & g' & 0 \\
   0 & B' & 1 \end{bmatrix}
   = \begin{bmatrix}
   1 & B+bg' & 0 \\
   0 & g & 0 \\
   0 & B'+b'g' & 1 \end{bmatrix}
   = \mathbf{M}_{\Box_{a,b}} = \begin{bmatrix}
   1 & (T-1)/2 & 0 \\
   0 & T & 0 \\
   0 & (T-1)/2 & 1 \end{bmatrix}

therefore, :math:`B+bg' = B'+b'g'`. It can be demonstrated using the ERO invariant inequalities from earlier that the
only solution to this that could correspond to an actual ERO is :math:`b=b'` and :math:`B=B'`.
:math:`g' = p = 4n + 3`, so :math:`b, b', B, B'` must all be odd. As mentioned earlier, there are no EROs with both `b` and `g` odd, so we have a contradiction, and :math:`\Box_{a,b}` is irreducible.

The proof for :math:`\Delta_{a,b}` is analogous. An Eisenstein integer :math:`a + bu`, :math:`u=\exp(\pi i/3)`, is prime
if its square norm :math:`a^2 + ab + b^2` is prime or the square of a prime. The prior (except for :math:`(1 + u)`, which
we corresponds to the ERO `n` which we already know is irreducible) have the form :math:`p=3n+1`; the latter, :math:`p=3n+2`. When the prime is of the latter form, the ERO is type 1 with :math:`\det(\Lambda_{\Delta_{a,b}}) = 1` and its matrix form is:

.. math::
   \mathbf{M}_{\Delta_{a,b}} = \begin{bmatrix}
          1 & (T-1)/3 & 0 \\
          0 & T & 0 \\
          0 & 2(T-1)/3 & 1 \end{bmatrix}.

Define `x` and `y` as before: then :math:`2(B+bg') = B'+b'g'`. Using the inequalities to exclude other choices, :math:`B' = 2B` and :math:`b' = 2b`. `g = 3n + 2`, but `g = b+ b' + 1 = 3b+1`: there is no simultaneous integer solution to both equations, so we have a contradiction, and :math:`\Delta_{a,b}` is irreducible.

Chirality
---------
.. _bowtie:
.. figure:: edge_chambers_bowtie.svg
   :align: right
   :figwidth: 25%

   The bowtie operator (B)

It may be possible to introduce another invariant into these operators and
distinguish operators not discerned by :math:`L_x` or :math:`M_x`. The most
desirable may be a measure for chirality; in theory that would distinguish,
e.g. `pp` vs `prp`. However, this does not appear as simple as assigning
achiral operators to 0 and :math:`\pm 1` to chiral operators. The composition
of a chiral operator and an achiral operator is always chiral, but:

* Two chiral operators can produce an achiral operator: `prp`
* Two chiral operators can produce another chiral operator:
  `pp`, `pg`, `prg`, `gg`, `grg`

Further confusing things are chiral EROs where r and d interact. Some
chiral EROs have `xd = x`, while some others have `xd = rxr`. (Some have
`x = dxd`, but none with `rxr = dxd` have been observed or proven/disproven to exist.)
The `gyro` operator is one example of the latter, and the bowtie operator
in :numref:`bowtie` is another, maybe easier-to-visualize example.
(Bowtie is a newly named operator, introduced in this text.)

Operators that produce alternating polyhedra
--------------------------------------------

The alternation operator `@` just exchanges :math:`+` and :math:`-`, so its
matrix form is a simple permutation matrix.

.. list-table:: Alternation operator `@` on bipartite structures

   * - :math:`[v^+,v^-,e,f]` to :math:`[v^+,v^-,e,f]`
     - :math:`[v,e,f^+,f^-]` to :math:`[v,e,f^+,f^-]`
   * - .. math:: \mathbf{M}_@ = \begin{bmatrix}
          0 & 1 & 0 & 0 \\
          1 & 0 & 0 & 0 \\
          0 & 0 & 1 & 0 \\
          0 & 0 & 0 & 1 \end{bmatrix}
     - .. math:: \mathbf{M}_@ = \begin{bmatrix}
          1 & 0 & 0 & 0 \\
          0 & 1 & 0 & 0 \\
          0 & 0 & 0 & 1 \\
          0 & 0 & 1 & 0 \end{bmatrix}

When considered with the bipartite structure, the dual operator `d` can be
considered to transform polyhedra with bipartite vertices into polyhedra with
bipartite faces and vice versa. On operators, it converts VAEROs to FAEROs (and
vice versa). Its matrix is also a simple permutation matrix.

.. list-table:: Dual operator `d` on bipartite structures

   * - :math:`[v^+,v^-,e,f]` to  :math:`[v,e,f^+,f^-]`
     - :math:`[v,e,f^+,f^-]` to :math:`[v^+,v^-,e,f]`
   * - .. math:: \mathbf{M}_d = \begin{bmatrix}
          0 & 0 & 0 & 1 \\
          0 & 0 & 1 & 0 \\
          1 & 0 & 0 & 0 \\
          0 & 1 & 0 & 0 \end{bmatrix}
     - .. math:: \mathbf{M}_d = \begin{bmatrix}
          0 & 0 & 1 & 0 \\
          0 & 0 & 0 & 1 \\
          0 & 1 & 0 & 0 \\
          1 & 0 & 0 & 0 \end{bmatrix}

The join operator `j` produces quadrilateral faces only. In fact, all type 0
:math:`\Box_{a,b}` operators produce quadrilateral faces, but those can be
reduced into :math:`j\Box_{c,d}` for some :math:`c, d`, so it's enough to look
at `j` for those operators. One way to assign a bipartite structure to the
vertices of `j` is to mark the seed vertices as :math:`+` and the vertices corresponding
to the seed faces as :math:`-`. Expressed as a matrix from :math:`[v,e,f]` to
:math:`[v^+,v^-,e,f]`:

.. math::
   \mathbf{M}_j = \begin{bmatrix}
   1 & 0 & 0 \\
   0 & 0 & 1 \\
   0 & 2 & 0 \\
   0 & 1 & 0 \end{bmatrix}

The opposite bipartite structure would simply be the same matrix, flipped from
left to right. This corresponds to applying the dual operator on the right:
`jd = @j`, so the relation gets a little more complicated when considering
alternating operators. The ambo operator produces bipartite faces,
and since `a=dj`, it can be expressed in terms of `j`, `d`, and `@`.

There are some tilings where an bipartite structure can be defined on both
the vertices and the faces. The square grid is one, as well as some regular
hyperbolic tilings (in general, any regular tiling with Schläfli symbol {n,m}
where n and m are both even). However, we haven't defined any operators
that require both vertices and faces to have an bipartite structure, so it's
enough to consider one at a time.

All EROs can be expressed with smoothing, an AERO, and the join operator
------------------------------------------------------------------------
The operator `$xj`, where `x` is a VAERO, is an ERO. If `x` is type 0 or 1
VAERO, then `$xj` is a type 0 operator. If `x` has undefined type, then `$xj`
is a type 1 operator. Although `$` does not in general have a :math:`M_x` form,
in the expression `$xj` it either does nothing, removes an edge and a vertex,
or removes an edge and a face. These operations can be represented by taking
the matrix form of `xj` and subtracting the zero matrix or these two following
matrices, respectively:

.. math::
   \begin{bmatrix}
   0 & 1 & 0 \\
   0 & 1 & 0 \\
   0 & 0 & 0 \end{bmatrix} ,
   \begin{bmatrix}
    0 & 0 & 0 \\
    0 & 1 & 0 \\
    0 & 1 & 0 \end{bmatrix} .

In fact, all EROs `y` can be expressed as `y = $xj`, where `x` is some VAERO or
ERO. This is easier to see by going backwards from the operator. As mentioned
earlier, if `g` is odd, there is an edge that lies on or crosses the center
point of the seed edge in the chamber structure. Otherwise `g` is even and
either a vertex lies there or a face contains the center point. If `g` is odd,
either split the edge with a degree-2 vertex at the center point, or replace the
edge with a digon. Then the alternating chamber structure of `x` is just the
white and grey chambers of `y`, stacked along their long edge. More
specifically, given an ERO `y`, if `g` is even, then `y = xj` for an ERO or
VAERO `x`: if `g` is odd, then `y = $xj` for (at least) two VAEROs `x`
corresponding to splitting the edge with a vertex or replacing an edge with a
digon. (Even though it can be reduced further in a larger set of operators, the
ERO form is usually preferable because including all those `$` and `j`
operators would get tedious.) A VAERO `x` may be named "pre-(Name)" where
(Name) is the name of `y`.

Note that since `xjd = x@j`, the ERO of the dual corresponds to the
opposite-partition VAERO. EROs may also be decomposed into FAEROs with the form
`y = $xa`, but since `a = dj` and `xd` has the chamber structure of `x` rotated,
it's simpler to just look at VAEROs.

Extension - Operations on different polyhedra
----------------------------------------------
With some care, operators can be applied to any polyhedron or tiling; toruses,
polyhedra with multiple holes, planar tilings, hyperbolic tilings, and even
non-orientable polyhedra, although the latter is restricted to the achiral
operators. The main restriction is that the graph must be embeddable on a
certain surface. Planar tilings may be easier to analyze by taking a finite
section and treating it as a torus.
It's worth noting that applying :math:`\Delta` to the regular triangular grid
on the plane, or :math:`\Box` to the regular square grid on the plane,
just creates a topologically equivalent grid on the plane.

.. _lozenge:
.. figure:: edge_chambers_lozenge.svg
   :align: right
   :figwidth: 25%

   The lozenge operator

Convex polyhedra may be put into "canonical form" such that all faces are flat,
all edges are tangent to the unit sphere, and the centroid of the polyhedron is
at the origin. As a consequence, all faces are convex.
There is no canonical form guaranteed to exist for general non-convex
polyhedra, however: in particular, there may be no position of the vertices
such that all the faces are flat or convex. The "Lozenge" operator in
:numref:`lozenge` creates concave faces when applied to a planar tiling.

Some operators can be applied to degenerate spherical polyhedra (dihedra and
hosohedra) with a result that is a convex polyhedron. Specifically, operators
with :math:`k > 1` may create a convex polyhedron from a dihedron, and
operators with :math:`\ell > 1` may create a convex polyhedron from a
hosohedron. (This is not guaranteed. For instance, try the lozenge operator on
a dihedron: the result won't even be 3-connected!)
For instance, a n-bipyramid is a kis n-dihedron, and (applying the
dual) an n-prism is a truncated n-hosohedron. Therefore the octahedron is a
kis 4-dihedron and the cube is a truncated 4-hosohedron.
This is interesting because the octahedron is also an ambo tetrahedron,
and the cube a join tetrahedron: if we admit degenerate polyhedra,
there are some polyhedra with two unequal reductions into operators and seeds.

Operators may also be applied to surfaces with boundary, although the behavior
of the operator at the boundary needs to be specified. In general,
this amounts to dropping incomplete faces or faces that cross over
the boundary, and dropping some related edges and vertices. We
lose the relationship with :math:`L_x` and :math:`M_x` by dropping those faces.

Extension - Operations that alter topology
------------------------------------------
In the topology of surfaces, the connected sum `A#B` of two surfaces `A` and `B`
can be thought of as removing a disk from A and B and stitching them together
along the created boundary.
If `B` has the topology of a sphere, then `A#B` has the topology of
`A`: a connected sum with a sphere does not change the topology. The
classification theorem of closed surfaces states that closed surfaces have the
topology of either a sphere or a connected sum of a number of toruses and/or
cross-caps.

In a topological sense, EROs and AEROs can be thought of as removing a disk from
a surface and replacing it with the chamber structure. In a more elaborate
sense, we can think of the operator chamber diagrams we've described so far
(even the alternating ones) as having the topology of a sphere: identify the two
edges on the left and the two edges on the right. Then, the operation can be
described as taking the connected sum of the operator chamber diagrams with each
face of the seed polyhedron. Thus assumption 2 in the list of assumptions at the
end of the "Operators on counts" section: taking the connected sum with a sphere
does not change the topology, so the operation does not change the Euler
characteristic.

.. _skeleton:
.. figure:: edge_chambers_skeleton.svg
   :align: right
   :figwidth: 25%

   Chambers of skeletonize operation.

However, operators that alter the topology can be described, introducing holes
or other features to a polyhedron. This may require us to think of the chamber
structure as having been extruded from a square into a square prism. One simple
operator of this kind makes nested or offset copies of the polyhedron:
obviously, this has :math:`M_x = n M_S = n I_3` where `n` is the number of
copies produced, and :math:`k = \ell = 1`. As expected, the Euler
characteristic of the result is the Euler characteristic of the seed times `n`.

Another operator is the skeletonize operator depicted in :numref:`skeleton`.
Edges and vertices are retained, but faces are removed. The red crosses
indicate that the base faces are not retained or replaced with vertices: they
are removed entirely. If `G` is the genus of the seed polyhedron, the genus of
the resulting "polyhedron" (inasmuch as an object with no faces can be
considered a polyhedron) is `G - f`. The :math:`M_x` form is obvious:

.. math::
   \begin{bmatrix}
   1 & 0 & 0 \\
   0 & 1 & 0 \\
   0 & 0 & 0 \end{bmatrix}

and :math:`k = \ell = 1`. (Technically :math:`\ell` could be any value, but it
makes sense to retain it as a measure of the hole created.)

Instead of annihilating the face completely, one can hollow out a space in its
center and leave behind a solid border. This can be done with the ``leonardo``
command in Antiprism, or the hollow/skeletonize/`h` operator in Polyhedronisme
(not to be confused with the skeletonize defined above, or the semi operator
from the last section). Although the operations differ in exactly how the new
faces are specified, topologically they both resemble a process like so:

* Duplicate the polyhedron as a slightly smaller polyhedron inside itself.
* For each face, remove the corresponding faces of the larger and smaller
  polyhedra. Take a torus and remove its outer half. Stitch the upper and lower
  boundary circles of this torus to the larger and smaller polyhedra where the
  faces were.

To represent this, we have to extrude the chamber structure out into a sort of
3d prism. The operator we'll describe here is essentially a process of replacing
each seed edge with a rectangular prism oriented with one edge along the seed
edge, somewhat like a 3d version of loft (`l`). (It is not the operation
performed by ``leonardo`` or Polyhedronisme, unfortunately; ``leonardo`` seems
to create overlapping faces.) In terms of invariants, :math:`k=\ell=1`,
:math:`b_4 = 2`, :math:`b'_4 = 4`, and :math:`M_x` is:

.. math::
   \begin{bmatrix}
   2 & 2 & 0 \\
   0 & 8 & 0 \\
   0 & 4 & 0 \end{bmatrix} .

If the seed polyhedron has Euler characteristic 2 (genus 0),
the result has Euler characteristic `4-2f`. The genus is `f-1`, not `f`,
because one torus is needed to connect the two copies of the sphere into
a (topologically) spherical surface.

One could also create operators that add arbitrary numbers of holes per edge.
(Operators that add cross-caps, e.g. based on a star polyhedron with Euler
characteristic 1 such as the tetrahemihexahedron, may be possible. Such
operators probably have more theoretical uses than aesthetic or practical ones,
and good luck getting the faces to be flat and not intersect awkwardly.)

Extensions - Multiple chambers
------------------------------
The concept of AEROs could be extended to k-partite graphs. :math:`k(k-1)/2`
interrelated chamber structures would have to be specified, which would get a
little unmanageable for large `k`. For example, if k=3, there would need to be
3 chambers: one on edges from set 1 to set 2, one from set 2 to set 3, and one
from set 1 to 3. By the four-color theorem, the largest `k` that is necessary
for a spherical tiling is 4, although larger `k` could be used.

Some EROs have forms where they are applied to only vertices or faces of
a certain order, such as :math:`t_3` to truncate vertices of order 3. These
could be described by a set of 3 chamber structures: on an edge between
order-3 vertices, on an edge from an order-3 vertex to a non-order-3 vertex
(or vice versa), and on an edge between non-order-3 vertices.

Neither of these schemes can be represented in the :math:`L_x` or :math:`M_x`
forms defined earlier.

Listing of operators and transformations
----------------------------------------
Where not specified, :math:`k` and :math:`\ell` are 1, and
:math:`b_i` and :math:`b'_i` are 0.

.. list-table:: EROs
   :header-rows: 1

   * - Operator `x`
     - Chiral?
     - Chambers of `x`
     - Matrix :math:`M_x`
     - :math:`k, \ell`, :math:`b_i`, :math:`b'_i`
     - Chambers of `dx`
     - Useful relations
   * - `S` (Seed, Identity)
     - N
     - .. image:: edge_chambers.svg
     - .. math::
        \begin{bmatrix}
        1 & 0 & 0 \\
        0 & 1 & 0 \\
        0 & 0 & 1 \end{bmatrix}
     -
     - .. image:: edge_chambers_dual.svg
     - `rr = S`, `dd = S`
   * - `j` (Join)
     - N
     - .. image:: edge_chambers_join.svg
     - .. math::
          \begin{bmatrix}
          1 & 0 & 1 \\
          0 & 2 & 0 \\
          0 & 1 & 0 \end{bmatrix}
     - :math:`b'_4=1`
     - .. image:: edge_chambers_ambo.svg
     - `j = jd = da = dad` (`jd=@j` and `ad=@a` if considering partitions)
   * - `k` (Kis)
     - N
     - .. image:: edge_chambers_kis.svg
     - .. math::
          \begin{bmatrix}
          1 & 0 & 1 \\
          0 & 3 & 0 \\
          0 & 2 & 0 \end{bmatrix}
     - :math:`k=2`, :math:`b'_3=2`
     - .. image:: edge_chambers_zip.svg
     - `k = nd = dz = dtd`
   * - `g` (Gyro)
     - Y
     - .. image:: edge_chambers_gyro.svg
     - .. math::
          \begin{bmatrix}
          1 & 2 & 1 \\
          0 & 5 & 0 \\
          0 & 2 & 0 \end{bmatrix}
     - :math:`b_3=2`, :math:`b'_5=2`
     - .. image:: edge_chambers_snub.svg
     - `g` = `rgdr` = `ds` = `rdsdr`
   * - `p` (Propeller)
     - Y
     - .. image:: edge_chambers_propeller-square_grid.svg
     - .. math::
          \begin{bmatrix}
          1 & 2 & 0 \\
          0 & 5 & 0 \\
          0 & 2 & 1 \end{bmatrix}
     - :math:`b_4=2`, :math:`b'_4=2`
     - .. image:: edge_chambers_dp.svg
     - `p = dpd`
   * - `c` (Chamfer)
     - N
     - .. image:: edge_chambers_chamfer.svg
     - .. math::
          \begin{bmatrix}
          1 & 2 & 0 \\
          0 & 4 & 0 \\
          0 & 1 & 1 \end{bmatrix}
     - :math:`b_3=2`, :math:`b'_6=1`
     - .. image:: edge_chambers_dc.svg
     - `c = dud`
   * - Lozenge
     - N
     - .. image:: edge_chambers_lozenge.svg
     - .. math::
          \begin{bmatrix}
          1 & 2 & 0 \\
          0 & 5 & 0 \\
          0 & 2 & 1 \end{bmatrix}
     - :math:`k=2`, :math:`\ell=2`, :math:`b_3=2`, :math:`b'_3=2`
     - .. image:: edge_chambers_dual_lozenge.svg
     - `x = dxd`
   * - `l` (Loft)
     - N
     - .. image:: edge_chambers_loft.svg
     - .. math::
          \begin{bmatrix}
          1 & 2 & 0 \\
          0 & 5 & 0 \\
          0 & 2 & 1 \end{bmatrix}
     - :math:`k=2`, :math:`b_3=2`, :math:`b'_4=2`
     - .. image:: edge_chambers_dual_loft.svg
     -
   * - `q` (Quinto)
     - N
     - .. image:: edge_chambers_quinto.svg
     - .. math::
          \begin{bmatrix}
          1 & 3 & 0 \\
          0 & 6 & 0 \\
          0 & 2 & 1 \end{bmatrix}
     - :math:`b_3=2`, :math:`b_4=1`, :math:`b'_5=2`
     - .. image:: edge_chambers_dual_quinto.svg
     -
   * - :math:`L_0` (Join-lace)
     - N
     - .. image:: edge_chambers_join-lace.svg
     - .. math::
          \begin{bmatrix}
          1 & 2 & 0 \\
          0 & 6 & 0 \\
          0 & 3 & 1 \end{bmatrix}
     - :math:`k=2`, :math:`b_4=2`, :math:`b'_3=2`, :math:`b'_4=1`
     - .. image:: edge_chambers_dual_lace0.svg
     -
   * - :math:`L` (Lace)
     - N
     - .. image:: edge_chambers_lace.svg
     - .. math::
          \begin{bmatrix}
          1 & 2 & 0 \\
          0 & 7 & 0 \\
          0 & 4 & 1 \end{bmatrix}
     - :math:`k=3`, :math:`b_4=2`, :math:`b'_3=4`
     - .. image:: edge_chambers_dual_lace.svg
     -
   * - :math:`K` (Stake)
     - N
     - .. image:: edge_chambers_stake.svg
     - .. math::
            \begin{bmatrix}
            1 & 2 & 1 \\
            0 & 7 & 0 \\
            0 & 4 & 0 \end{bmatrix}
     - :math:`k=3`, :math:`b_3=2`, :math:`b'_3=2`, :math:`b'_4=2`
     - .. image:: edge_chambers_dual_stake.svg
     -
   * - :math:`w` (Whirl)
     - Y
     - .. image:: edge_chambers_whirl.svg
     - .. math::
          \begin{bmatrix}
          1 & 4 & 0 \\
          0 & 7 & 0 \\
          0 & 2 & 1 \end{bmatrix}
     - :math:`b_3=4`, :math:`b'_6=2`
     - .. image:: edge_chambers_dual_whirl.svg
     -
   * - :math:`J=(kk)_0` (Join-kis-kis)
     - N
     - .. image:: edge_chambers_join-kis-kis.svg
     - .. math::
          \begin{bmatrix}
          1 & 2 & 1 \\
          0 & 8 & 0 \\
          0 & 5 & 0 \end{bmatrix}
     - :math:`k=3`, :math:`\ell=2`, :math:`b_3=2`, :math:`b'_3=4`, :math:`b'_4=1`
     - .. image:: edge_chambers_dual_kiskis0.svg
     -
   * - :math:`X` (Cross)
     - N
     - .. image:: edge_chambers_cross.svg
     - .. math::
          \begin{bmatrix}
          1 & 3 & 1 \\
          0 & 10 & 0 \\
          0 & 6 & 0 \end{bmatrix}
     - :math:`k=2`, :math:`b_4=2`, :math:`b_6=1`, :math:`b'_3=4`, :math:`b'_4=2`
     - .. image:: edge_chambers_dual_cross.svg
     -
   * - :math:`W` (Waffle) (New)
     - N
     - .. image:: edge_chambers_waffle.svg
     - .. math::
          \begin{bmatrix}
          1 & 4 & 1 \\
          0 & 9 & 0 \\
        0 & 4 & 0 \end{bmatrix}
     - :math:`b_3=2`, :math:`b_4=2`, :math:`b'_4=2`, :math:`b'_5=2`
     - .. image:: edge_chambers_dual_waffle.svg
     -
   * - :math:`B` (Bowtie) (New)
     - Y
     - .. image:: edge_chambers_bowtie.svg
     - .. math::
          \begin{bmatrix}
          1 & 5 & 1 \\
          0 & 10 & 0 \\
          0 & 4 & 0 \end{bmatrix}
     - :math:`b_3=4`, :math:`b_4=1`, :math:`b'_3=2`, :math:`b'_7=2`
     - .. image:: edge_chambers_dual_bowtie.svg
     - `rBr=Bd`

.. list-table:: ERO families
   :header-rows: 1

   * - Operator `x`
     - Chiral?
     - Matrix :math:`M_x`
     - :math:`k, \ell`, :math:`b_i`, :math:`b'_i`
     - Useful relations
   * - :math:`m_n` (Meta)
     - N
     - .. math::
          \begin{bmatrix}
          1 & n & 1 \\
          0 & 3n+3 & 0 \\
          0 & 2n+2 & 1 \end{bmatrix}
     - :math:`k=2`, :math:`\ell=n+1`, :math:`b_4=n`, :math:`b'_3=2n+2`
     - :math:`m_1 = m = kj`
   * - :math:`M_n` (Medial)
     - N
     - .. math::
          \begin{bmatrix}
          1 & n & 1 \\
          0 & 3n+1 & 0 \\
          0 & 2n & 1 \end{bmatrix}
     - :math:`\ell=n`, :math:`b_4=n`, :math:`b'_3=2n-2`, :math:`b'_4=2`
     - :math:`M_1 = o = jj`
   * - :math:`\Delta_{a,b}` if `T` divisible by 3
     - If :math:`a \ne b` and :math:`b \ne 0`
     - .. math::
          \begin{bmatrix}
          1 & T/3-1 & 1 \\
          0 & T & 0 \\
          0 & 2T/3 & 0 \end{bmatrix}
     - :math:`b_6=b`, :math:`b'_3=b'`
     - :math:`\Delta_{1,1} = n`,
       :math:`\Delta_{a,b}` :math:`= n \Delta_{(2a+b)/3, (b-a)/3}`
   * - :math:`\Delta_{a,b}` if `T` not divisible by 3
     - If :math:`a \ne b` and :math:`b \ne 0`
     - .. math::
          \begin{bmatrix}
          1 & (T-1)/3 & 0 \\
          0 & T & 0 \\
          0 & 2(T-1)/3 & 1 \end{bmatrix}
     - :math:`b_6=b`, :math:`b'_3=b'`
     - :math:`\Delta_{2,0} = u`, :math:`\Delta_{2,1} = dwd`
   * - :math:`\Box_{a,b}` if `T` even
     - If :math:`a \ne b` and :math:`b \ne 0`
     - .. math::
          \begin{bmatrix}
          1 & T/2-1 & 1 \\
          0 & T & 0 \\
          0 & T/2 & 0 \end{bmatrix}
     - :math:`b_4=b`, :math:`b'_4=b'`
     - :math:`\Box_{a,b} = \Box_{a,b}d`,
       :math:`\Box_{1,1} = j`, :math:`\Box_{2,0} = o = j^2`,
       :math:`\Box_{a,b}` :math:`= j\Box_{(a+b)/2,(b-a)/2}`,
       (:math:`\Box_{a,b}d = @\Box_{a,b}` if alternating vertices)
   * - :math:`\Box_{a,b}` if `T` odd
     - If :math:`a \ne b` and :math:`b \ne 0`
     - .. math::
          \begin{bmatrix}
          1 & (T-1)/2 & 0 \\
          0 & T & 0 \\
          0 & (T-1)/2 & 1 \end{bmatrix}
     - :math:`b_4` :math:`=b'_4` :math:`=b` :math:`=b'`
     - :math:`\Box_{a,b} = d\Box_{a,b}d`, :math:`\Box_{1,2} = p`

In the following two tables, when :math:`k^+=k^-`, both
are written as just :math:`k`.

.. list-table:: VAEROs
   :header-rows: 1

   * - Operator
     - Degree-2?
     - Chambers of `x`
     - Matrix
     - :math:`k_i, \ell_i`, :math:`b_i`, :math:`b'_i`
     - Chambers of `dx`
     - Useful relations
   * - Alternation, Hemi, Semi
     - Digons
     - .. image:: edge_chambers_alternating_semi.svg
     - .. math::
          \begin{bmatrix}
          1 & 0 & 0 & 0 \\
          0 & 0 & 1 & 0 \\
          0 & 1 & 0 & 1 \end{bmatrix}
     - :math:`k^+ = 2`,  :math:`\ell = 1/2`
     - .. image:: edge_chambers_alternating_dual_hemi.svg
     - `$xj = S`, `$dxj = d`
   * - Alternating Truncate (Pre-Chamfer)
     - N
     - .. image:: edge_chambers_alternating_truncate.svg
     - .. math::
          \begin{bmatrix}
          1 & 0 & 1 & 0 \\
          0 & 0 & 2 & 0 \\
          0 & 1 & 0 & 1 \end{bmatrix}
     - :math:`\ell = 3/2`, :math:`b_3=1`
     - .. image:: edge_chambers_alternating_dual_prechamfer.svg
     - `xj = c`, `dxjd = u`
   * - Pre-kis
     - Digons
     - .. image:: edge_chambers_alternating_bisect.svg
     - .. math::
          \begin{bmatrix}
          1 & 0 & 0 \\
          0 & 2 & 0 \\
          0 & 1 & 1 \end{bmatrix}
     - :math:`b'_3 = 1`, :math:`\ell = 1/2`, :math:`k^+ = 3`
     - .. image:: edge_chambers_alternating_dual_prekis.svg
     - `$xj = k`
   * - Pre-Join-Stake
     - N
     - .. image:: edge_chambers_alternating_prestake0.svg
     - .. math::
          \begin{bmatrix}
          1 & 1 & 0 \\
          0 & 3 & 0 \\
          0 & 1 & 1 \end{bmatrix}
     - :math:`k^+=2`, :math:`b_3=1`, :math:`b'_4=1`
     - .. image:: edge_chambers_alternating_dual_prestake0.svg
     - `xj = jk`
   * - Alternating Subdivide
     - N
     - .. image:: edge_chambers_alternating_subdivide.svg
     - .. math::
          \begin{bmatrix}
          1 & 1 & 0 \\
          0 & 3 & 0 \\
          0 & 1 & 1 \end{bmatrix}
     - :math:`\ell = 3/2`, :math:`b_4=1`, :math:`b'_3=1`
     - .. image:: edge_chambers_alternating_dual_subdivide.svg
     -
   * - Pre-Gyro
     - Degree-2 vertices
     - .. image:: edge_chambers_alternating_ortho.svg
     - .. math::
          \begin{bmatrix}
          1 & 1 & 1 \\
          0 & 3 & 0 \\
          0 & 1 & 0 \end{bmatrix}
     - :math:`\ell = 1/2`, :math:`b_3=1`, :math:`b'_6=1`
     - .. image:: edge_chambers_alternating_dual_pregyro.svg
     - `$xj = g`. Not the same as Pre-Join-Lace of dual.
   * - Pre-Join-Lace
     - N
     - .. image:: edge_chambers_alternating_prelace0.svg
     - .. math::
          \begin{bmatrix}
          1 & 0 & 1 & 0 \\
          0 & 0 & 3 & 0 \\
          0 & 1 & 1 & 1 \end{bmatrix}
     - :math:`k^+=2`, :math:`b_4=1`, :math:`b'_3=1`
     - .. image:: edge_chambers_alternating_dual_prejoinlace.svg
     - :math:`xj = L_0`. Not the same as pre-gyro of dual.
   * - Pre-Join-Kis-Kis
     - N
     - .. image:: edge_chambers_alternating_prekiskis0.svg
     - .. math::
          \begin{bmatrix}
          1 & 1 & 0 \\
          0 & 4 & 0 \\
          0 & 2 & 1 \end{bmatrix}
     - :math:`k^+=3`, :math:`k^-=2`, :math:`b_3=1`, :math:`b'_3=2`
     - .. image:: edge_chambers_alternating_dual_prekiskis0.svg
     - :math:`xj = (kk)_0`
   * - Pre-Cross
     - N
     - .. image:: edge_chambers_alternating_metaortho.svg
     - .. math::
          \begin{bmatrix}
          1 & 1 & 1 \\
          0 & 5 & 0 \\
          0 & 3 & 0 \end{bmatrix}
     - :math:`k^+=1`, :math:`k^-=2`, :math:`\ell = 3/2`,
       :math:`b_4=1`, :math:`b'_3=2`, :math:`b'_4=1`
     - .. image:: edge_chambers_alternating_dual_precross.svg
     - `xj = X`
   * - Alternating Meta/Join
     - N
     - .. image:: edge_chambers_alternating_metajoin.svg
     - .. math::
          \begin{bmatrix}
          1 & 1 & 1 \\
          0 & 5 & 0 \\
          0 & 3 & 0 \end{bmatrix}
     - :math:`k^+=1`, :math:`k^-=2`, :math:`\ell = 2`,
       :math:`b_3=1`, :math:`b'_3=3`
     - .. image:: edge_chambers_alternating_dual_mj.svg
     -
   * - Alternating Subdivide/Quinto
     - N
     - .. image:: edge_chambers_alternating_subdividequinto.svg
     - .. math::
          \begin{bmatrix}
          1 & 2 & 0 \\
          0 & 5 & 0 \\
          0 & 2 & 1 \end{bmatrix}
     - :math:`b_3=1`, :math:`b_5=1`, :math:`b'_4=2`
     - .. image:: edge_chambers_alternating_dual_uq.svg
     - `xj = jg`

Open questions
--------------
* Are there any irreducible EROs other than `j` that produce only
  quad faces?
* Are there any chiral EROs such that `rxr = dxd`? (They would have to be
  type 1 operators.)
* Are there other conditions that can be added to the invariants for
  :math:`L_x` to make the set of conditions sufficient as well as necessary?
* Is there an invariant related to the chirality of an operator?
* What other invariants need to be added to fully characterize EROs and AEROs?
