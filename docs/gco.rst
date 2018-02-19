Goldberg-Coxeter Operations on Polyhedra and Tilings
====================================================
The Goldberg-Coxeter (GC) operation can be used to subdivide the faces of a
polyhedron or tiling with triangular or square faces, by replacing the faces
with a portion of a lattice of triangular or square faces. Geodesic
subdivision is another name for the method on triangular faces; this method
can be attributed to [Goldberg]_, [Coxeter]_, Fuller, or [Caspar]_ and Klug,
depending on whether you're a mathematician or a scientist :) This method
often produces polyhedra with nice geometric qualities, for instance,
local symmetry preservation, minimal distortion, etc.

Strictly, the GC operation defined in e.g. [Deza2004]_ and [Deza2015]_ is
defined on a polyhedral graph, not a polyhedron. Embedding the graph into R3 produces
a polyhedron. In the embedding, it is necessary to consider geometric questions
like the exact placement of vertices that aren't an issue with graphs; those
are addressed here.

Notation
--------
This text uses :math:`\Delta(a,b)` for the triangular GC operator and
:math:`\Box(a,b)` for the quadrilateral GC operator. This contrasts with
[Deza2004]_, who use :math:`GC_{k,l}(G_0)` for both.

Contents
--------
.. toctree::

   gco-planar
   gco-arb
   gco-spherical
   gco-a-spherical
   gco-a-improper
   gco-a-mixed
