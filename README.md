# antitile
Manipulation of tilings using Python. This package is designed to work with [Antiprism](https://github.com/antiprism/antiprism) but can be used on its own.

Under construction.

## Installation

    pip3 install git+git://github.com/brsr/antitile.git

## Usage
There are currently 4 scripts included with this package:
* `view_off.py`: A viewer for OFF files using matplotlib.
* `balloon.py`: Balloon tiling of the sphere
* `sgs.py`: Similar grid subdivision of tilings
* `sgsstats.py`: Statistics of SGS tilings

## For Contributors
This code makes heavy use of vectorized operations on NumPy multidimensional arrays, which are honestly pretty impenetrable until you get familiar with them. (And, uh, even after that.) I use the convention that the last axis of an array specifies the spatial coordinates:

    x, y, z = v[..., 0], v[..., 1], v[..., 2]
