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
   :figwidth: image

   Chambers of a triangular face.

.. figure:: edge_chambers.svg
   :align: right
   :figwidth: image

   Chamber structure of a Conway operator.

[Brinkmann]_ observed that Conway operators can be described in terms of
chambers. Each face may be divided into chambers by identifying the face center
and drawing lines from there to each vertex and edge midpoint. The operator
may then be specified by a structure within those chambers. Here
we'll draw the chamber structure of an operator as in the figure to the right:

* The two nodes on the left and right are vertices of the seed.
* The blue line from left to right is an edge of the seed.
* The two nodes on the top and bottom are face centers of the seed.
* The contents of the white and grey chambers are the same, rotated.

If an operator is achiral, the grey chamber is a reflection of the adjacent
white chamber. Technically we only need the upper-left white chamber for
achiral operators or the upper white and grey chambers for chiral operators,
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
subscripts `i/x` should be interpreted as 0 if i/x is not an integer.

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
