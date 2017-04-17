# antitile
Manipulation of polyhedra and tilings using Python. This package is designed to work with [Antiprism](https://github.com/antiprism/antiprism) but can be used on its own.

## Installation

    pip3 install git+git://github.com/brsr/antitile.git

## Usage
These scripts are included with this package:
* `sgs.py`: Similar grid subdivision of tilings
* `balloon.py`: Balloon tiling of the sphere
* `cellular.py`: Colors polyhedra using cellular automata
* `sgsstats.py`: Statistics of the polyhedra/tiling, focused on the use of SGS tilings to model the sphere (see also `off_report` in Antiprism)
* `view_off.py`: A viewer for OFF files using matplotlib, allowing for export to SVG (see also `antiview` in Antiprism)

These can be piped together with programs from Antiprism.

OFF files for the regular icosahedron, octahedron, tetrahedron, cube, and 3- and 4-edged dihedra are included in the `data` folder in the source.

## Examples
Statistics of a geodesic polyhedron (created using what geodesic dome people call Method 1):

    sgs.py -a 5 -b 3 icosahedron.off | sgsstats.py

Visualize a Goldberg polyhedron, with color:

    sgs.py -a 5 -b 3 icosahedron.off | off_color -v M -m color_map.txt | pol_recip | view_off.py

Canonical form (no skew faces) of a quadrilateral-faced similar grid subdivision polyhedron:

    sgs.py -a 5 -b 3 cube.off | canonical | view_off.py
    
A quadrilateral balloon polyhedra, which happens to resemble a peeled coconut:

    balloon.py 8 -pq | view_off.py

## For Contributors
This code makes heavy use of vectorized operations on NumPy multidimensional arrays, which are honestly pretty impenetrable until you get familiar with them. (And, uh, even after that.) I use the convention that the last axis of an array specifies the spatial coordinates:

    x, y, z = v[..., 0], v[..., 1], v[..., 2]
