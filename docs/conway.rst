Notes on Conway operators
=========================

Chamber structure
-----------------
[Brinkmann]_ observed that Conway operators can be described in terms of
chambers

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
