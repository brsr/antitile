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

Chamber structure
-----------------
.. figure:: Triangle_chambers.svg
   :align: right
   :figwidth: 25%

   Chambers of a triangular face.

.. figure:: edge_chambers.svg
   :align: right
   :figwidth: 25%

   Chambers adjacent to an edge.

[Brinkmann]_ observed that Conway operators can be described in terms of
chambers. Each face may be divided into chambers by identifying the face center
and drawing lines from there to each vertex and edge midpoint.
The operator may then be specified by a structure within those chambers. If an
operator is achiral, the grey chamber is a reflection of the adjacent white
chamber. Technically we only need the upper-left white chamber for achiral
operators or the upper white and grey chambers for chiral operators,
but showing both sides of the edge will make things easier later on.

There is some freedom in where vertices are placed within the chambers.
This is more apparent with chiral operators. Often the operator is drawn
so that most of the vertices lie on the seed edge, but this is not necessary.
For instance, the image on the left is a chamber diagram for how George Hart
originally drew his propeller operator (see [HartPropeller]_),
but the image on the right is topologically
equivalent and emphasizes the operator's relationship with a square grid.

.. image:: edge_chambers_propeller.svg

.. image:: edge_chambers_propeller-square_grid.svg

Operators on counts
-------------------
When `x` is the operator, :math:`[v,e,f]` are the vertices, edges, and faces of
the seed, and :math:`[v',e',f']` are the vertices, edges, and faces of the
result, then :math:`[v',e',f'] = \mathbf{M}_x [v,e,f]`.

.. math::
   \mathbf{M}_x = \begin{bmatrix}
   a & b & c \\
   0 & g & 0 \\
   a' & b' & c' \end{bmatrix}

where a + a' = 1, c + c' = 1, and g= b + b' + 1, and a, a', b, b', c, and c'
are all nonnegative integers. a, a', c, c' to be {0, 1}, and g is a positive
integer.


A more elaborate representation is as an infinite linear operator. Let `e` and
`e'` be the count of edges before and after like above, but now :math:`v_i` and
:math:`v'_i` are the count of vertexes of order `i` before and after, antipodal
:math:`f_i` and :math:`f'_i` are the count of faces with `i` sides.
:math:`\sum v_i = v`, and so on for the rest of these.

.. math::
   e' &= ge

   v'_i &= a v_{i/x} + e b_i + c f_{i/y}

   f'_i &= a' v_{i/x} + e b'_i + c' f_{i/y}

where :math:`\sum b_i = b`, :math:`\sum b'_i = b'`, all :math:`b_i` and
:math:`b'_i` are nonnegative integers, and x and y are positive integers. The
subscripts `i/x` should be interpreted as 0 if `i/x` is not an integer.

.. image:: edge_chambers_bowtie.svg
   :alt:   Chambers of the bowtie operator (B)

.. image:: edge_chambers_waffle.svg
   :alt:   Chambers of the waffle operator (W)

Relation to the Goldberg-Coxeter operation
------------------------------------------

Extensions
----------
allow a, a', c, c' to be {0, 1/2, 1}

.. math::
   e' &= ge

   v'_i &= a (v_{i/x_1} + v_{i/x_2})/2 + e b_i + c (f_{i/y_1} + f_{i/y_2})/2

   f'_i &= a' (v_{i/x_1} + v_{i/x_2})/2 + e b'_i + c'(f_{i/y_1} + f_{i/y_2})/2

dealing with digons and order-2 vertices

Table of values
---------------

======== ============ = =============== = = = ================================
Operator g            a b               c x y :math:`b_i` & :math:`b'_i`
======== ============ = =============== = = = ================================
S        1            1 0               0 1 1
d        1            0 0               1 1 1
j        2            1 0               1 1 1 :math:`b'_4=1`
k        3            1 0               1 2 1 :math:`b'_3=2`
g        5            1 2               1 1 1 :math:`b_3=2, b'_5=2`
p        5            1 2               0 1 1 :math:`b_4=2, b'_4=2`
c        4            1 2               0 1 1 :math:`b_3=2, b'_6=1`
l        5            1 2               0 2 1 :math:`b_3=2, b'_4=2`
q        6            1 3               0 1 1 :math:`b_3=2, b_4=1, b'_5=2`
K0       6            1 2               1 2 1 :math:`b_3=2, b'_4=3`
K        7            1 2               1 3 1 :math:`b_3=2, b'_3=2, b'_4=2`
L0       6            1 2               0 2 1 :math:`b_4=2, b'_3=2, b'_4=1`
L        7            1 2               0 3 1 :math:`b_4=2, b'_3=4`
w        7            1 4               0 1 1 :math:`b_3=4, b'_6=2`
J=(kk)0  8            1 2               1 3 2 :math:`b_3=2, b'_3=1, b'_4=4`
X        10           1 3               1 2 2 :math:`b_4=2, b_6=1, b'_3=4, b'_4=2`
m        :math:`3n+3` 1 `n`             1 2 n :math:`b_4=n, b'_3=2n+2`
M        :math:`3n+1` 1 `n`             1 1 n :math:`b_4=n, b'_3=2n-2, b'_4=2`
Δ even   `T`          1 :math:`T/3-1`   1 1 1 :math:`b_6=b, b'_3=b'`
Δ odd    `T`          1 :math:`(T-1)/3` 0 1 1 :math:`b_6=b, b'_3=b'`
☐ even   `T`          1 :math:`T/2-1`   1 1 1 :math:`b_4=b, b'_4=b'`
☐ odd    `T`          1 :math:`(T-1)/2` 0 1 1 :math:`b_4=b'_4=b`
Waffle   9            1 4               1 1 1 :math:`b_3=2, b_4=2, b'_4=2, b'_5=2`
Bowtie   10           1 5               1 1 1 :math:`b_3=4, b_4=1, b'_3=2, b'_7=2`
======== ============ = =============== = = = ================================
