Goldberg-Coxeter Operation on arbitrary tilings and polyhedra
=============================================================
This section will discuss tilings on arbitrary surfaces before we discuss the
spherical case, i.e. the behavior of the program with options
``gcopoly -n -p=flat``. Usually the GC operation is described on a usual
polyhedra of genus 0, but it can be applied to most surfaces. The main caveat
is that class III operators doesn't work on non-orientable surfaces: because
the operators are chiral, they depend on an underlying orientable surface.

Definition of GC operation
---------------------------
The GC operation comprises these steps:
* Choose a master polygon corresponding to :math:`\Box(a,b)` or
  :math:`\Delta(a,b)`.
* Replace each face of the base polyhedron with a stretched,
  transformed version of the master polygon that covers the base face.
* Stitch together any edges that might cross the edges of the base face.

There are some places in the literature where the GC operation is defined
as the dual of this operation. Particularly, in the fullerene literature,
most writers are interested in hexagonal faces, so they use the dual.

Coordinate form
---------------
To transform the master polygon to the target face polygon, this program
generally uses a standardized coordinate form.

Triangles can be are specified by barycentric coordinates
:math:`\beta_i` where :math:`\beta_1 + \beta_2 + \beta_3 = 1`, and then the
corresponding vertex is given by
:math:`\mathbf v = \beta_1\mathbf v_1+\beta_2\mathbf v_2+\beta_3\mathbf v_3`.
Given :math:`\mathbf v` and :math:`\mathbf v_i`, :math:`\beta_i` can be found
by using e.g. the Moore–Penrose inverse to solve the system given by
:math:`\beta_1 + \beta_2 + \beta_3 = 1` and
:math:`\mathbf v = \beta_1\mathbf v_1+\beta_2\mathbf v_2+\beta_3\mathbf v_3`.

Quadrilaterals are instead specified by what we'll call "xy coordinates"
where :math:`x` and :math:`y` are between 0 and 1. :math:`\mathbf v_1 +
(\mathbf v_2-\mathbf v_1) x + (\mathbf v_4-\mathbf v_1) y
+ (\mathbf v_1-\mathbf v_2+\mathbf v_3-\mathbf v_4)xy`. Unlike triangles,
quadrilaterals may have points that do not share a common plane: they may be
skew quadrilaterals. If the quadrilateral is a skew quadrilateral,
`x` and `y` smoothly parameterize a surface over that skew quadrilateral. In
this program, we don't need to find the xy coordinates for a point in an
arbitrary quadrilateral. There are some occasions where the xy coordinates are
determined for a square in the plane, in which case they can be converted
by rotation and scaling.

Linear index form
-----------------
Another way to think about the master polygon is to count lines in the plane
that cross it. This method can break down at the base vertices, but we already
know what those are supposed to be, so that's not a practical problem.

For :math:`\Delta(a,b)`, draw an equilateral triangle in the plane. Along one
edge, mark `a+b+1` evenly spaced points (including the base vertices). Going
clockwise, mark `a+1` points along the next edge, then `b+1` along the last.
Number the points from 0 to `a+b` along the first side, then along the
remaining two sides sequentially (counting the point shared between the other
two sides only once). Draw lines between points with the same number.
Then take those lines and make copies rotated 120° and 240° about the center.
The vertices of the master triangle correspond to the points where
3 lines intersect (not necessarily the marked points), and the linear index is
a 3-tuple of the number we gave to the original line, the 120° line, and the
240° line that intersect at that vertex.

For :math:`\Box(a,b)`, draw a square in the plane. Mark `b+1` vertices along
the left edge of the square (including the base vertices) and `a+1` vertices
along the bottom. Number the points from 0 to `a+b`, starting at the upper
left. Copy these points, rotate them by 180° about the center, and reverse
their numbering. Draw a line between each point with the same number. Then take
those lines and make a copy rotated by 90° about the center.
The vertices of the master square correspond to the points where
2 lines intersect (not necessarily the marked points), and the linear index is
a 2-tuple of the number we gave to the original line and the rotated line that
intersect that that vertex.

Given the Gaussian or Eisenstein coordinates of the vertices
in the master polygon, they can be directly converted to linear indexes.

Triangular: When the vertex is :math:`c + du`, in normal form,
and the master triangle is given by (a, b), then the linear indexes satisfy
:math:`\sum \ell_i = 2a+b`, with

.. math::
   \ell_1 = b + c + d
   \ell_2 = a - c
   \ell_3 = a - d

Quadrilateral: When the vertex is :math:`c + di`, in normal form,
and the master square is given by (a, b), then

.. math::
   \ell_1 = c + a
   \ell_2 = d

The primary use for linear index in the ``gcopoly`` program is "stitching"
neighboring polygons together at their edges. This operation is easy to
do by eye but a little complicated to automate. Linear indexes have the
advantage that they are integers, not floats, and therefore rounding error
isn't an issue and the stitching can be done exactly.

This program, when creating the master polygon, leaves in some vertices that
are outside the polygon but share edges with vertices inside it. This is
done in part to make stitching easier: the program determines which
vertices in adjacent faces match up, and collapses them into one.
In detail, the stitching procedure is as follows:

* Select two faces that share an edge. Make a copy of the linear indexes for
  each.
* Rotate the order of the linear index of both faces so that the vertices
  shared by the edge have the same lindex. (If a class I or II operation,
  can reverse the tuple order too to accomodate unoriented surfaces.)
* On one face, replace the linear index with :math:`\mathbf o - \mathbf \ell`,
  where :math:`\mathbf o = (a+b, a, 2a+b)` for :math:`\Delta(a,b)` and
  :math:`\mathbf o = (a+2b, b)` for :math:`\Box(a,b)`.
* Collapse vertices on different faces with the same linear index.

Surfaces with boundaries present a problem because there are tiles that don't
have a neighboring tile for each edge. In the language of [Deza2004]_, the
skeleton graph of this tiling has vertices with degree 2; that paper doesn't
even define their operation on such graphs.

As mentioned earlier, this program leaves in vertices that are outside the
master polygon but are connected to vertexes inside it. Correspondingly,
the program will retain those vertices when it transforms the master polygon.
This means that surfaces with boundaries may have faces outside the original
boundary if a Class II or III operator is applied.

(However, see the appendix gco-a-improper for information on the improper
*spherical* tilings that may have graphs with vertices of degree 2 or less.
These are not surfaces with boundaries, and can be handled differently.)

Composition: Topology vs Geometry
---------------------------------
On the plane, equality of operators holds in both a topological and geometric
sense: if :math:`\Box(a,b)\Box(c,d) = \Box(e,f)`, then the left-hand side
results in the same placement of vertices on the plane. This topological
equality carries over to tilings on arbitrary surfaces, but not
necessarily the geometric sense. Furthermore, commutativity and associativity
hold for topology, but not necessarily geometry.
So considering the exact position of vertices in space, where :math:`(a,b)`,
:math:`(c,d)`, and :math:`(e,f)` are complex numbers in the appropriate space:

* :math:`\Delta(a,b)\Delta(c,d)` may not equal :math:`\Delta(e,f)`,
  where :math:`(e,f)` is the normal form of the product of
  :math:`(a,b)` and :math:`(c,d)`
* :math:`\Delta(a,b)\Delta(c,d)` may not equal :math:`\Delta(c,d)\Delta(a,b)`
* :math:`\Delta(a,b)(\Delta(c,d)\Delta(e,f))` may not equal
  :math:`(\Delta(a,b)\Delta(c,d))\Delta(e,f)`

and the same for :math:`\Box`.

This fact actually turns out to be useful,
because it gives us some flexibility in vertex placement.