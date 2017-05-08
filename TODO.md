# To Do
* Clean up/refactor
* Eliminate the use of `pandas`
* Document
* Test cases
* Make `setup.py` put the files in `data` somewhere accessible
* Replace the stitching logic in `antitile/sgs.py` so that for Class I and II subdivision, the faces do not need to be oriented. In theory this would allow SGS to operate on tilings on unorientable surfaces, although self-intersecting surfaces might have quirks to deal with.
* Work out a theory of mixed triangle & quad grids other than Class I. (Triangle (1,1) seems to work with square (2,1), for instance, while square (1,1) doesn't seem to work with any triangle breakdown.)
* Optimize (big frequencies)
* Allow arbitrary face figures other than SGS (would probably need whole new program)
