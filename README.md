[![Documentation Status](https://readthedocs.org/projects/antitile/badge/?version=latest)](http://antitile.readthedocs.io/en/latest/?badge=latest)

# antitile
Manipulation of polyhedra and tilings using Python. This package is designed to work with [Antiprism](https://github.com/antiprism/antiprism) but can be used on its own.

## Installation

    pip3 install git+git://github.com/brsr/antitile.git

## Usage
The package includes a number of scripts.
* `gcopoly.py`: The Goldberg-Coxeter subdivision operation of tilings
* `balloon.py`: Balloon tiling of the sphere
* `cellular.py`: Colors polyhedra using cellular automata
* `gcostats.py`: Statistics of the polyhedra/tiling, focused on the use of GCO to model the sphere (see also `off_report` in Antiprism)
* `view_off.py`: A viewer for OFF files using matplotlib, allowing for export to SVG (see also `antiview` in Antiprism)
* `factor.py`: Factors Gaussian, Eisenstein, Steineisen (Eisenstein expressed with the 6th root of unity instead of 3rd), and regular integers.

Some Jupyter notebooks exploring various aspects of these programs are included in the `misc` folder in the source, as well as some OFF files for simple polyhedra including the 3- and 4-dihedra.

## Examples
Statistics of a geodesic polyhedron (created using what geodesic dome people call Method 1):

    gcopoly.py -a 5 -b 3 icosahedron.off | sgsstats.py

Visualize a Goldberg polyhedron, with color:

    gcopoly.py -a 5 -b 3 icosahedron.off | off_color -v M -m group_color.txt | pol_recip | view_off.py

Create a quadrilateral-faced similar grid subdivision polyhedron, put it into canonical form (so the faces are all flat), and color it using Conway's Game of Life with random initial condition:

	gcopoly.py -a 5 -b 3 cube.off | canonical | off_color test.off -f n | cellular.py -v -b=3 -s=2,3
	view_off.py cellular100.off
	# or whatever the last file is if it reaches steady state early

`gcopoly.py` can subdivide non-orientable surfaces too, at least for Class I and II subdivisions. Here, the base is a Möbius strip-like surface with 12 faces:

	unitile2d 2 -s m -w 4 -l 1 | gcopoly.py -a 2 -b 2 -n | view_off.py

A quadrilateral balloon polyhedra, which happens to resemble a peeled coconut:

    balloon.py 8 -pql | view_off.py

## For Contributors
This code makes heavy use of vectorized operations on NumPy multidimensional arrays, which are honestly pretty impenetrable until you get familiar with them. (And, uh, even after that.) I use the convention that the last axis of an array specifies the spatial coordinates:

    x, y, z = v[..., 0], v[..., 1], v[..., 2]
