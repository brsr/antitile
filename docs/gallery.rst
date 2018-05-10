Gallery
=======

.. list-table:: GC Operation

   * - .. image:: geodesic_8.png

       A standard (method 1) 8-frequency geodesic division of an icosahedron.
       ``gcopoly.py icosahedron.off -a 8 | off_color -f C | antiview -R 90,0,0``
     - .. image:: cube_4_3.png

       (4,3) GC operation on a cube, in canonical form. ``gcopoly.py cube.off
       -a 4 -b 3 | canonical | off_color -f C| antiview -R 90,0,0``
   * - .. image:: tet_5_3.png

       Standard (method 1) (5,3) subdivision of a tetrahedron. (Notice the
       vertices in the upper left and lower right.) ``gcopoly.py
       tetrahedron.off -a 5 -b 3 | off_color -f G | antiview -R 90,0,0``
     - .. image:: tet_5_3_opt.png

       (5,3) subdivision of a tetrahedron, but using `nslerp` projection and
       optimized to reduce variation in face area. ``gcopoly.py tetrahedron.off
       -a 5 -b 3 -p nslerp -k faces | off_color -f G | antiview -R 90,0,0``
   * - .. image:: cellular.gif

       Conway's Game of Life on the surface of a (5,3) subdivided cube.
       `Courtesy of Roger Kaufman
       <https://groups.google.com/d/msg/antiprism/hyKdjk-iwOI/tqi-Q1svAQAJ>`_.
     - .. image:: virus_1.png

       Repeatedly applying the (1,1) triangular subdivision using the flat
       projection and not normalizing gives this interesting organic shape.
       It resembles a viral capsid (moreso than geodesic polyhedra usually
       resemble a viral capsid), or a strange fruit. Produced by
       ``gcopoly.py icosahedron.off -a 9 -n -f |antiview``: ``-a 7`` and
       ``-a 5 -b 3`` are similar.
   * - .. image:: 3di.png

       Improper geodesic polyhedron on the 3-dihedron. Note dangling faces.
       ``gcopoly.py 4dihedron.off -a 5 -b 3 -p nslerp -k energy |
       off_color -f C | antiview -R 90,90,0``
     - .. image:: 4di.png

       GC operation on the 4-dihedron.
       ``gcopoly.py 4dihedron.off -a 5 -b 3 -p nslerp -k energy |
       off_color -f C | antiview -R 90,45,0``

.. list-table:: Other scripts

   * - .. image:: breakdown_4_4_nslerp2.png

       Quadrilateral `nslerp2` preserves diagonal lines across the
       face. ``breakdown.py 4 4 -q -p=nslerp2``
     - .. image:: breakdown_4_4_disk.png

       An attractive floral pattern produced by
       ``breakdown.py 4 4 -z 0 -q -p=disk``
   * - .. image:: balloon_10.svg

       A triangular balloon polyhedron. ``balloon.py 10 | off_color -f C |
       view_off.py -e 0 -d 5.5``
     - .. image:: balloon_10_qpl.svg

       A quadrilateral balloon polyhedron, which resembles the shape of a
       peeled coconut. ``balloon.py 10 -q -p -l | off_color -f C |
       view_off.py -d 5.6 -e 0``
