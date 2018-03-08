Notes on Conway operators
=========================

Conway originally defined a set of operations that could be performed on the
Platonic solids to obtain the Archimedean and Catalan solids. Although some of
the operations date back as early as Kepler, operations that obtain one
polyhedron using another one as a seed as now known as Conway operations or
Conway operators. Conway's original set of operations is denoted with the
letters `abdegjkmost`, and more are available depending on your software. For
instance, a truncated cube may be denoted `tC`, where `t` is the truncate
operation and `C` is a cube. Initially there was not much theory supporting
Conway operations, but [Brinkmann]_'s paper provides a framework. This text
is an attempt to use Brinkmann's work to find ways to better quantify and
analyze Conway operators.

Preliminaries
-------------
This assumes some familiarity with Conway operators. See [HartConway]_ for a
basic overview, or better still, spend some time playing with
[Polyhedronisme]_ (a web app) or ``conway`` in [Antiprism]_.
In general, this text uses the same terms as ``conway``.

Antiprism uses `S` (seed) for the identity operator. The dual operator is `d`,
and the mirror-image operator is `r`. Both `dd` and `rr` equal `S`. For every
operator `x`, there also exist operators `xd`, `dx`, and `dxd`, some of which
may be equal to `x`. Since the characteristics of the latter operators
are closely related to `x`, this text will generally pick one operator from
that collection and use that to represent the whole. (Usually the choice is
an operator that preserves the base vertices, for consistency's sake.)

An operator that cannot be expressed in terms of operators aside from `d` and
`r` is "primitive". For instance, `k` (Kis) and `j` (Join) are primitive,
but `m` (Meta) is not (it is equal to `kj`).

Chamber structure
-----------------
.. _facechambers:
.. figure:: Triangle_chambers.svg
   :align: right
   :figwidth: 25%

   Chambers of a triangular face.

.. _edgechambers:
.. figure:: edge_chambers.svg
   :align: right
   :figwidth: 25%

   Chambers adjacent to an edge.

[Brinkmann]_ et al. observed that Conway operators can be described in terms
of chambers. Each face may be divided into chambers by identifying the face
center and drawing lines from there to each vertex and edge midpoint, as in
:numref:`facechambers`. Similarly, each vertex of degree `n` is surrounded
by `n` white and `n` grey chambers. Each edge has a white and grey chamber on
each side of the edge, as shown in :numref:`edgechambers`.

The operator may then be specified by a structure of vertices and edges within
those chambers, possibly with edges crossing from one edge to another. If an
operator is achiral, the grey chamber is a reflection of the adjacent white
chamber. A list of operators with their chamber structures is listed at the
end of this text. In this text, we show the two chambers adjacent to each
side of a seed edge. Technically we only need the upper-left white chamber for
achiral operators or the upper white and grey chambers for chiral operators,
but showing both sides of the edge will make things easier later on.

For a given operator `x`, the chamber structure of `xd` is simply the
chamber structure of `x` rotated one quarter turn.

There is some freedom in where vertices are placed within the chambers.
This is more apparent with chiral operators. Often the operator is drawn
so that most of the vertices lie on the seed edge, but this is not necessary.
For instance, :numref:`propeller` is a chamber diagram for how George Hart
originally drew his propeller operator (see [HartPropeller]_),
but :numref:`propellersq` is topologically
equivalent and emphasizes the operator's relationship with a square grid.

.. _propeller:
.. figure:: edge_chambers_propeller.svg
   :align: center

   George Hart's drawing of the chambers for the operator `p` (Propeller)

.. _propellersq:
.. figure:: edge_chambers_propeller-square_grid.svg
   :align: center

   Grid drawing of the chambers for the operator `p` (Propeller)

Operators on counts
-------------------
In abstract algebraic terms, Conway operators form a monoid: a group without
an inverse, or a semigroup with an identity element. Let :math:`[v,e,f]` be
the count of vertices, edges, and faces of the seed, and :math:`v_i` and
:math:`f_i` be the count of vertices/faces of order `i` such that
:math:`\sum v_i = v` and :math:`\sum f_i = f`. There is a series of monoids and
homomorphisms between the monoids, as so:

* Conway operator `x` (acts on polyhedra)
* Infinite-dimensional linear operator :math:`L_x` (acts on :math:`v_i, e, f_i`)
* 3x3 matrix :math:`M_x` (acts on :math:`[v,e,f]`)
* Eigenvalues of :math:`M_x`

Each bullet will be handled in turn.

The action of the operator on the vertices of degree `i`, edges, and faces with
`i` sides can be described with an infinite linear operator :math:`L_x`. This
operator can be determined by counting elements off the chamber structure.
Step by step:

* Seed vertices are either retained or converted into faces centered on that
  vertex. (Other options are precluded by symmetry). Let `a = 1` if the
  seed vertices are retained, and 0 otherwise. Also, the degree of the vertex
  or face is either the same as the seed vertex, or a multiple of it; let `k`
  be that multiple.
* Seed face centers are either retained (possibly of in a smaller face) or
  converted into vertices. (Again, other options are precluded by symmetry).
  Let `c = 0` if the seed faces are retained, and 1 otherwise. Let :math:`\ell`
  serve a similar role as `k` above: the degree of the vertex or face
  corresponding to the seed face center is `k` times the degree of
* Except for the faces or vertices corresponding to the seed vertices and face
  centers, the added elements are in proportion to to the number of
  edges in the seed. `g` is the count of added edges (the edge multiplier or
  inflation rate from [Brinkmann]_ et al.),
  :math:`b_i` is the number of vertices of degree i added,
  and :math:`b'_i` is the number of faces of degree i added.

Count elements lying on or crossing the outer edge of the chamber structure as
half. It may help to draw an adjacent chamber, particularly when determining
the number of sides on a face. The result of the counting process can be
described in the following operator form;
variables in capital letters are the result of the operator.

.. math::
   E &= ge

   V_i &= a v_{i/k} + e b_i + c f_{i/\ell}

   F_i &= a' v_{i/k} + e b'_i + c' f_{i/\ell}

where `a`, :math:`a'`, `c`, and :math:`c'` are either 0 or 1, `g` is a
positive integer, all :math:`b_i` and :math:`b'_i` are nonnegative integers,
and `k` and :math:`\ell` are positive integers. The subscripted values like
:math:`v_{i/k}` should be interpreted as 0 if `i/k` is not an integer.

Under the constraint that the operator preserves the Euler characteristic,
it can be shown that :math:`a + a' = 1`, :math:`c + c' = 1`, and
:math:`g= b + b' + 1` where :math:`\sum b_i = b` and :math:`\sum b'_i = b'`.
Also, since :math:`b_i` and :math:`b'_i` are nonnegative integers, only a
finite number of their values can be non-zero. This makes the operator form
more manageable than the term "infinite linear operator" may suggest; in
reality, nearly all applications will only use a finite number of different
vertex and face degrees.

Applying the handshake lemma gives relations between the values:

.. math::
   2g &= 2ak + 2c\ell + \sum i b_i

   2g &= 2a'k + 2c'\ell + \sum i b'_i

If the polyhedron doesn't have degenerate features (e.g digons or degree-2
vertices), :math:`i \ge 3`. Together with characteristics from above, a
series of inequalities can be derived:

.. math::
   2k + 2\ell - 2 \le g + 1 \le 2a + 3b + 2c \le 2g

All these relations taken together  are necessary but not sufficient. The values
:math:`g=3`, :math:`a=1`, :math:`c=0`, :math:`k=1`, :math:`\ell=1`,
:math:`b_4=1`, :math:`b'_4=1` satisfy the relations, but do not appear
to correspond to any Conway operator. (However, see the "Extensions" section.)

The dual operator :math:`L_d` has the form :math:`E = e, V_i = f_i, F_i = v_i`.
With a little manipulation, it is easy to see that if :math:`L_x` has values
`a`, :math:`b_i`, `c`, etc, then applications of the dual operator have related
forms. :math:`L_x L_d`'s values exchange `a` with `c`, :math:`a'` with
:math:`c'`, and `k` with :math:`\ell`. :math:`L_d L_x`'s values exchange `a`
with :math:`a'`, `c` with :math:`c'`, and each :math:`b_i` with each
:math:`b'_i`. Finally, :math:`L_d L_x L_d`'s values exchange `a` with
:math:`c'`, and :math:`a'` with `c`, `k` with :math:`\ell`,
and each :math:`b_i` with each :math:`b'_i`.

The matrix form :math:`M_x` can be obtained from :math:`L_x` by summing
:math:`\sum v_i = v` and :math:`\sum f_i = f`, or from counting elements
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
   a & b & c \end{bmatrix}, \mathbf{M}_d \mathbf{M}_x \mathbf{M}_d = \begin{bmatrix}
   c' & b' & a' \\
   0 & g & 0 \\
   c & b & a \end{bmatrix}

The matrix :math:`M_x` has three eigenvalues: `1`, `g`, and `(a-c)`. Thus, its
determinant is `g(a-c)`. The first eigenvalue is constant and the second is the
edge multiplier defined earlier. The third is either equal to -1, 0, or 1.
The dual operator interchanges -1 and 1, which gives some motivation to using
operators with `a=1` as the representative operators over those with `a=-1`.
Operators can be thought of as having a parity based on `a` and `c`: if `a=c`,
the operator has even parity, otherwise it has odd parity. Like multiplication
of natural numbers, the composition of any operator with an even operator is
even, and the composition of two odd operators is odd.

For an operator `xy`, i.e. the composition of `x` and `y`, the expansion factor
`g` is the product of the `g` values for each operator, and the quantity `(a-c)`
is the product of each operator's `(a-c)`. For the matrix form, composition is
just the usual matrix multiplication: :math:`M_{xy} = M_x M_y`. Explicitly, let
:math:`g, a, b_i, b'_i, c, k, \ell` be the values for :math:`L_y`;
:math:`G, A, B_i, B'_i, C, K, L` for :math:`L_x`; and
:math:`\gamma, \alpha, \beta_i, \beta'_i, \sigma, \kappa, \lambda`
for :math:`L_{xy}`:

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

.. _waffle:
.. figure:: edge_chambers_waffle.svg
   :align: right
   :figwidth: 25%

   The waffle operator (W)

None of these homomorphisms are injections: there are certain
:math:`L_x` or :math:`M_x` that correspond to more than one Conway operator.
Examples for :math:`M_x` are easy to come by: where `n = kd`, :math:`M_k = M_n`.
For an example where the operators are not related by duality,
:math:`M_l = M_p`. For :math:`L_x`, :math:`L_{prp} = L_{pp}` but `prp` is not
the same as `pp` (one's chiral, one's not). For the operator depicted in
:numref:`waffle`, :math:`W \ne Wd`, but :math:`L_W = L_{Wd}`.
(This is a newly named operator, introduced in this text.)

Some further consequences of these representations:

* If `x=xd`, the operator is even. If `x=dxd`, the operator is odd.
* Operators where `g` is a prime number are primitive.
* There are no odd operators with `g=2`, so therefore odd operators
  with `g=2p`, where p is prime, are primitive.

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
* Two chiral operators can produce a chiral operator: `pp`, `pg`, `prg`

Further confusing things are operators where r and d interact.
The bowtie operator in :numref:`bowtie` is a simple example of this: `rB = Bd`.
(This is a newly named operator, introduced in this text.)

Relation to the Goldberg-Coxeter operation
------------------------------------------

The Goldberg-Coxeter operation can be fairly simply extended to a Conway
operator. In the master polygon, identify two vertices and the center: this is
the chamber structure of the operator.

Many of the named Conway operators are GC operations, or related by duality.
GC operators are also a good source of examples; in the 2-parameter families,
it's often easy to find an operator with a desired quality.
GC operators have an invariant `T`, the "trianglation number",
which is identical to the Conway operator edge factor `g`.

* :math:`\Box_{a,b}`: :math:`g = T = a^2 + b^2`
* :math:`\Delta_{a,b}`: :math:`g = T = a^2 + ab + b^2`

Extensions
----------
With some care, Conway operators can be applied to any polyhedron or tiling,
including those with holes. Chiral operators may only be applied to orientable
polyhedra. Planar tilings may be easier to analyze by taking a finite section
and treating it as a torus. There is no canonical form defined for non-spherical
polyhedra or tilings, however.

In [Coxeter8]_ (specifically section 8.6), Coxeter defines an alternation
operation `h` on polyhedra with only even-sided faces. Each face is replaced
with a face with half as many sides, and alternate vertices are either retained
as part of the faces or converted into vertices with number of sides equal to
the degree of the seed vertex. (He also defines a snub operation in section 8.4,
different from the `s` snub Conway defined, that is equivalent to `ht`.) The
alternation operation converts quadrilateral faces into digons. Usually the
digons are converted into edges, but for now, let digons be digons.


allow a, a', c, c' to be {0, 1/2, 1}, allow k and \ell in N/2

dealing with digons and order-2 vertices

Table of operators
------------------
Where not specified, :math:`k` and :math:`\ell` are 1, and
:math:`b_i` and :math:`b'_i` are 0.

.. list-table::

   * - Operator
     - Chiral?
     - Chambers
     - Matrix
     - :math:`k, \ell`, :math:`b_i`, :math:`b'_i`
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
     - `rr = S`
   * - `d` (Dual)
     - N
     - .. image:: edge_chambers_dual.svg
     - .. math::
          \begin{bmatrix}
          0 & 0 & 1 \\
          0 & 1 & 0 \\
          1 & 0 & 0 \end{bmatrix}
     -
     - `dd = S`
   * - `j` (Join)
     - N
     - .. image:: edge_chambers_join.svg
     - .. math::
          \begin{bmatrix}
          1 & 0 & 1 \\
          0 & 2 & 0 \\
          0 & 1 & 0 \end{bmatrix}
     - :math:`b'_4=1`
     - `j = jd = da = dad`
   * - `k` (Kis)
     - N
     - .. image:: edge_chambers_kis.svg
     - .. math::
          \begin{bmatrix}
          1 & 0 & 1 \\
          0 & 3 & 0 \\
          0 & 2 & 0 \end{bmatrix}
     - :math:`k=2`, :math:`b'_3=2`
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
     - `g = gd = ds = dsd`
   * - `p` (Propeller)
     - Y
     - .. image:: edge_chambers_propeller.svg
     - .. math::
          \begin{bmatrix}
          1 & 2 & 0 \\
          0 & 5 & 0 \\
          0 & 2 & 1 \end{bmatrix}
     - :math:`b_4=2`, :math:`b'_4=2`
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
     - `c = dud`
   * - `l` (Loft)
     - N
     - .. image:: edge_chambers_loft.svg
     - .. math::
          \begin{bmatrix}
          1 & 2 & 0 \\
          0 & 5 & 0 \\
          0 & 2 & 1 \end{bmatrix}
     - :math:`k=2`, :math:`b_3=2`, :math:`b'_4=2`
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
     -
   * - :math:`K_0` (Join-stake)
     - N
     - .. image:: edge_chambers_join-stake.svg
     - .. math::
          \begin{bmatrix}
          1 & 2 & 1 \\
          0 & 6 & 0 \\
          0 & 3 & 0 \end{bmatrix}
     - :math:`k=2`, :math:`b_3=2`, :math:`b'_4=3`
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
     - `rB=Bd`
   * - :math:`m_n` (Meta)
     - N
     -
     - .. math::
          \begin{bmatrix}
          1 & n & 1 \\
          0 & 3n+3 & 0 \\
          0 & 2n+2 & 1 \end{bmatrix}
     - :math:`k=2`, :math:`\ell=n+1`, :math:`b_4=n`, :math:`b'_3=2n+2`
     - :math:`m_1 = m = kj`
   * - :math:`M_n` (Medial)
     - N
     -
     - .. math::
          \begin{bmatrix}
          1 & n & 1 \\
          0 & 3n+1 & 0 \\
          0 & 2n & 1 \end{bmatrix}
     - :math:`\ell=n`, :math:`b_4=n`, :math:`b'_3=2n-2`, :math:`b'_4=2`
     - :math:`M_1 = o = jj`
   * - :math:`\Delta_{a,b}` if `T` divisible by 3
     - If :math:`a \ne b` and :math:`b \ne 0`
     -
     - .. math::
          \begin{bmatrix}
          1 & T/3-1 & 1 \\
          0 & T & 0 \\
          0 & 2T/3 & 0 \end{bmatrix}
     - :math:`b_6=b`, :math:`b'_3=b'`
     - :math:`\Delta_{2,0} = u`
   * - :math:`\Delta_{a,b}` if `T` not divisible by 3
     - If :math:`a \ne b` and :math:`b \ne 0`
     -
     - .. math::
          \begin{bmatrix}
          1 & (T-1)/3 & 0 \\
          0 & T & 0 \\
          0 & 2(T-1)/3 & 1 \end{bmatrix}
     - :math:`b_6=b`, :math:`b'_3=b'`
     - :math:`\Delta_{1,1} = n`, :math:`\Delta_{2,1} = dwd`
   * - :math:`\Box_{a,b}` if `T` even
     - If :math:`a \ne b` and :math:`b \ne 0`
     -
     - .. math::
          \begin{bmatrix}
          1 & T/2-1 & 1 \\
          0 & T & 0 \\
          0 & T/2 & 0 \end{bmatrix}
     - :math:`b_4=b`, :math:`b'_4=b'`
     - :math:`\Box_{a,b} = \Box_{a,b}d`,
       :math:`\Box_{1,1} = j`, :math:`\Box_{2,0} = o = j^2`
   * - :math:`\Box_{a,b}` if `T` odd
     - If :math:`a \ne b` and :math:`b \ne 0`
     -
     - .. math::
          \begin{bmatrix}
          1 & (T-1)/2 & 0 \\
          0 & T & 0 \\
          0 & (T-1)/2 & 1 \end{bmatrix}
     - :math:`b_4` :math:`=b'_4` :math:`=b` :math:`=b'`
     - :math:`\Box_{a,b} = d\Box_{a,b}d`, :math:`\Box_{1,2} = p`

Open questions
--------------
* Are there any operators such that `rx = dxd`? (They would have to be odd
  operators.)
* Is/are there an/other condition/s that can be added to the values for
  :math:`L_x` to make the set of conditions sufficient as well as necessary?
* Is there a good invariant related to the chirality of a Conway operator?
* What other invariants need to be added to fully characterize Conway operators?
