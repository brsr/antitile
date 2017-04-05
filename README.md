# antitile
Manipulation of polyhedra and tilings using Python. This package is designed to work with [Antiprism](https://github.com/antiprism/antiprism) but can be used on its own.

## Installation

    pip3 install git+git://github.com/brsr/antitile.git

## Usage
There are currently 4 programs included with this package:
* `view_off.py`: A viewer for OFF files using matplotlib.
* `balloon.py`: Balloon tiling of the sphere
* `sgs.py`: Similar grid subdivision of tilings
* `sgsstats.py`: Statistics of SGS tilings
These can be chained together with programs from Antiprism.

OFF files for the regular icosahedron, octahedron, tetrahedron, cube, and 3- and 4-point dihedra are included (although where Python puts them may depend on your system).

## Examples
Statistics of a geodesic polyhedron (using what geodesic dome people call Method 1):

	sgs.py -a 5 -b 3 icosahedron.off | sgsstats.py

Visualize a Goldberg polyhedron, with color:

    sgs.py -a 5 -b 3 icosahedron.off | off_color -v M -m color_map.txt | pol_recip | view_off.py

Canonical form (no skew faces) of a quadrilateral-faced similar grid subdivision polyhedron:

    sgs.py -a 5 -b 3 cube.off | canonical | view_off.py
    
Quadrilateral balloon polyhedra kinda look like a peeled coconut:

    balloon.py 8 -pq | view_off.py

## For Contributors
This code makes heavy use of vectorized operations on NumPy multidimensional arrays, which are honestly pretty impenetrable until you get familiar with them. (And, uh, even after that.) I use the convention that the last axis of an array specifies the spatial coordinates:

    x, y, z = v[..., 0], v[..., 1], v[..., 2]
