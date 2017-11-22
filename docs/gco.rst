Goldberg-Coxeter Operations on Polyhedra and Tilings
====================================================
The Goldberg-Coxeter (GC) operation can be used to subdivide the faces of a
polyhedron or tiling with triangular or square faces. Geodesic subdivision is
another name for the method on triangular faces; this method can be attributed
to Goldberg, Coxeter, Fuller, Caspar, and Klug. [Brinkmann]_ This method often
produces polyhedra with nice geometric qualities, for instance,
local symmetry preservation, minimal distortion, etc.

Strictly, the GC operation defined in e.g. [Deza]_ is defined on
a planar graph: however, planar graphs are closely related to polyhedra and
tilings. Where the literature describes the subdivision of vertices of a
planar graph, this program implements the subdivision of faces of a polyhedron
or tiling. One is (the skeleton of) the dual of the other; also, with 
polyhedra, it is necessary to consider geometric questions like the exact 
placement of vertices that aren't an issue with graphs.

Notation
--------
This text uses :math:`\Delta(a,b)` for the triangular GC operator and 
:math:`\Box(a,b)` for the quadrilateral GC operator. This program deviates 
from the parameterization of the triangular GC operation used in 
chemistry-oriented texts such as [Deza]_: it uses the traditional 
parameterization for geodesic domes and Goldberg polyhedra. 

Applications
------------
Applications and closely related topics include:

* Architecture: Geodesic domes
* Geodesy: Geodesic grids and the quad sphere
* Biology: The structure of viral capsids
* Chemistry: Buckminsterfullerene
* Computer graphics: Subdivision surfaces (Loop, Catmull–Clark)
* As a topic of study in itself: for example, the 2014 discovery that the
(icosahedral) Goldberg polyhedra form a fourth class of convex equilateral
polyhedra http://www.pnas.org/content/111/8/2920.abstract

.. toctree::

   gco-planar
   gco-arb
   gco-transform
   gco-spherical
   gco-a-dihedra
   gco-a-henahedron
   gco-a-mixed

References
----------
.. [Brinkmann] https://arxiv.org/abs/1705.02848
.. [Deza] Michel-Marie Deza, Mathieu Dutour Sikirić, Mikhail Ivanovitch 
   Shtogrin. (2015) Geometric Structure of Chemistry-Relevant Graphs: Zigzags 
   and Central Circuits. Springer. pp 131-148. 
   https://books.google.com/books?id=HLi4CQAAQBAJ&lpg=PA130&ots=ls1r5QkM51&dq=goldberg-coxeter&pg=PA130#v=onepage&q=goldberg-coxeter&f=false
