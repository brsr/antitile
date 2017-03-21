# antitile
Manipulation of tilings using Python and [Antiprism](https://github.com/antiprism/antiprism).

Under construction.

## Installation

    pip3 install git+git://github.com/brsr/antitile.git

## Usage
There are currently 3 scripts included with this package:
* `view_off.py` : A viewer for OFF files using matplotlib.
* `balloon.py` : Balloon tiling of the sphere
* `sgs.py` : Similar grid subdivision of tilings

## For Contributors
This code makes heavy use of vectorized operations on NumPy multidimensional arrays, which are honestly pretty impenetrable until you get familiar with them. (And, uh, even after that.) I use the convention that the last axis of an array specifies the spatial coordinates:

    x, y, z = v[..., 0], v[..., 1], v[..., 2]
