# antitile
Manipulation of polyhedra and tilings using Python. This package is designed to work with [Antiprism](https://github.com/antiprism/antiprism) but can be used on its own.

## Installation

    pip3 install git+git://github.com/brsr/antitile.git

## Usage
The package includes a number of scripts. These can be piped together with programs from Antiprism:
* `sgs.py`: Similar grid subdivision of tilings
* `balloon.py`: Balloon tiling of the sphere
* `cellular.py`: Colors polyhedra using cellular automata
* `sgsstats.py`: Statistics of the polyhedra/tiling, focused on the use of SGS tilings to model the sphere (see also `off_report` in Antiprism)
* `view_off.py`: A viewer for OFF files using matplotlib, allowing for export to SVG (see also `antiview` in Antiprism)

These are free-standing:
* `breakdown.py`: Visualize breakdown structures
* `factor.py`: Factors Gaussian, Eisenstein, Nietsnesie (Eisenstein based on the 6th root of unity instead of 3rd), and regular integers.

OFF files for the regular icosahedron, octahedron, tetrahedron, cube, and 3- and 4-edged dihedra are included in the `data` folder in the source.

## Examples
Statistics of a geodesic polyhedron (created using what geodesic dome people call Method 1):

    sgs.py -a 5 -b 3 icosahedron.off | sgsstats.py

Visualize a Goldberg polyhedron, with color:

    sgs.py -a 5 -b 3 icosahedron.off | off_color -v M -m group_color.txt | pol_recip | view_off.py

Create a quadrilateral-faced similar grid subdivision polyhedron, put it into canonical form (so the faces are all flat), and color it using Conway's Game of Life with random initial condition:

	sgs.py -a 5 -b 3 cube.off | canonical | off_color test.off -f n | cellular.py -v -b=3 -s=2,3
	view_off.py cellular100.off #or whatever if it reaches steady state early

`sgs.py` can subdivide (Class I and II) non-orientable surfaces too. Here, the base is a mobius strip-like surface with 12 faces:

	unitile2d 2 -s m -w 4 -l 1 | sgs.py -a 2 -b 2 -n | view_off.py	
	
A quadrilateral balloon polyhedra, which happens to resemble a peeled coconut:

    balloon.py 8 -pq | view_off.py

## For Contributors
This code makes heavy use of vectorized operations on NumPy multidimensional arrays, which are honestly pretty impenetrable until you get familiar with them. (And, uh, even after that.) I use the convention that the last axis of an array specifies the spatial coordinates:

    x, y, z = v[..., 0], v[..., 1], v[..., 2]
