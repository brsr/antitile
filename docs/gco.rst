Goldberg-Coxeter Operations on Polyhedra and Tilings
====================================================
The Goldberg-Coxeter (GC) operation can be used to subdivide the faces of a
polyhedron or tiling with triangular or square faces. Geodesic subdivision is
another name for the method on triangular faces; this method can be attributed
to Goldberg, Coxeter, Fuller, Caspar, and Klug. [Goldberg]_ [Caspar]_ 
[Coxeter]_ [Brinkmann]_ This method often
produces polyhedra with nice geometric qualities, for instance,
local symmetry preservation, minimal distortion, etc.

Strictly, the GC operation defined in e.g. [Deza2015]_ is defined on
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
chemistry-oriented texts such as [Deza2015]_: it uses the traditional 
parameterization for geodesic domes and Goldberg polyhedra. 

Applications
------------
Applications and closely related topics include:

* Architecture: Geodesic domes
* Geodesy: Geodesic grids and the quad sphere
* Biology: The structure of viral capsids
* Chemistry: Buckminsterfullerene
* Computer graphics: Subdivision surfaces 
  (Loop [StamLoop]_, Catmull–Clark [StamCatmull]_)
* As a topic of study in itself: for example, the 2014 discovery that the
  (icosahedral) Goldberg polyhedra form a fourth class of convex equilateral
  polyhedra [Schein]_

Contents
--------
.. toctree::

   gco-planar
   gco-arb
   gco-spherical
   gco-a-spherical
   gco-a-improper
   gco-a-mixed

.. rubric:: References

.. [Brinkmann] https://arxiv.org/abs/1705.02848
.. [Deza2015] Michel-Marie Deza, Mathieu Dutour Sikirić, Mikhail Ivanovitch 
   Shtogrin. (2015) Geometric Structure of Chemistry-Relevant Graphs: Zigzags 
   and Central Circuits. Springer. pp 131-148. 
   https://books.google.com/books?id=HLi4CQAAQBAJ&lpg=PA130&ots=ls1r5QkM51&dq=goldberg-coxeter&pg=PA130#v=onepage&q=goldberg-coxeter&f=false   
.. [Altschuler] Altschuler, E. L. et al. 1997. Possible Global Minimum Lattice 
   Configurations for Thomson's Problem of Charges on a Sphere. Phys. Rev. 
   Lett. 78, 2681. [http://dx.doi.org/10.1103/PhysRevLett.78.2681 
   doi:10.1103/PhysRevLett.78.2681] 
   http://www.mcs.anl.gov/~zippy/publications/thomson/thomsonPRL.html
.. [Goldberg] Goldberg, M, 1937. A class of multi-symmetric polyhedra. Tohoku 
   Mathematical Journal.
.. [Hart1997] Hart, G, 1997. Calculating Canonical Polyhedra. Mathematica in 
   Education and Research, Vol 6 No. 3, Summer 1997, pp. 5-10. 
.. [Hart2012] Hart, G, 2012. Goldberg Polyhedra. In Senechal, Marjorie. 
   Shaping Space (2nd ed.). Springer. pp. 125–138. 
   [http://dx.doi.org/10.1007/978-0-387-92714-5_9 
   doi:10.1007/978-0-387-92714-5_9].
.. [Kenner] Kenner, H, 1976. Geodesic Math and How to Use It. University of 
   California Press.
.. [Schein] Schein & Gayed, 2014. Fourth class of convex equilateral 
   polyhedron with polyhedral symmetry related to fullerenes and viruses. 
   PNAS, Early Edition [http://dx.doi.org/10.1073/pnas.1310939111 
   doi:10.1073/pnas.1310939111]
.. [Oosterom] Van Oosterom, A & Strackee, J, 1983. The Solid Angle of a Plane 
   Triangle. IEEE Trans. Biom. Eng. BME-30 (2): 125–126. 
   [http://dx.doi.org/10.1109/TBME.1983.325207 doi:10.1109/TBME.1983.325207].   
.. [Caspar] D.L.D. Caspar and A. Klug. Physical principles in the construction 
   of regular viruses. In Cold Spring Harb Symp Quant Biol., volume 27, 
   pages 1–24, 1962.
.. [Coxeter] H.S.M. Coxeter. Virus macromolecules and geodesic domes. In J.C. 
   Butcher, editor, A spectrum of mathematics, pages 98–107. 
   Oxford University Press, 1971.
.. [Deza2004] M. Deza and M. Dutour. Goldberg-Coxeter constructions for 3- 
   and 4-valent plane graphs. The Electronic Journal of Combinatorics, 
   11, 2004. #R20.
.. [Tarnai] T. Tarnai, F. Kovacs, P.W. Fowler, and S.D. Guest. Wrapping the 
   cube and other polyhedra. Proceedings of the Royal Society A, 
   468:2652–2666, 2012.   
.. [StamLoop] Jos Stam: Evaluation of Loop Subdivision Surfaces, Computer 
   Graphics Proceedings ACM SIGGRAPH 1998,    
.. [StamCatmull] Stam, J. (1998). "Exact evaluation of Catmull-Clark 
   subdivision surfaces at arbitrary parameter values". Proceedings of the 
   25th annual conference on Computer graphics and interactive techniques - 
   SIGGRAPH '98. pp. 395–404. doi:10.1145/280814.280945. 
   ISBN 0-89791-999-8.   