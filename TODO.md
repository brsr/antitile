# To Do
* Clean up the code
* Document
* Test cases
* Improve coloring
* Add 2d mode for `view_off.py`
* Port over the cellular automata on polyhedra faces script from the old `geogrid` module
* Refactor `figures_*.py`
* Replace the stitching logic in `antitile/sgs.py` so that for Class I and II subdivision, the faces do not need to be oriented. In theory this would allow SGS to operate on tilings on unorientable surfaces, although self-intersecting surfaces might have quirks to deal with.
* Allow "mixed" triangle & quad grids, at least for Class I. (Triangle class I and quad class I work together, but triangle class II works with particular quad class III subdivisions, not every triangle class III has a compatible quad class III, and quad class II doesn't work with triangle at all. Probably need to work out the theory of what's compatible with what.)
* Allow arbitrary face figures other than SGS (would probably need whole new program)
