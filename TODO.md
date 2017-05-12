# To Do
* Clean up/refactor
* Eliminate the use of `pandas`
* Document
* Test cases
* Make `setup.py` put the files in `data` somewhere accessible
* Optimize (big frequencies)
* `antitile/sgs.py` sometimes creates duplicate faces with non-oriented faces for Class II divisions. This is probably due to overlapping faces with different orientations.
* Work out a theory of mixed triangle & quad grids other than Class I. (Triangle (1,1) seems to work with square (2,1), for instance, while square (1,1) doesn't seem to work with any triangle breakdown.)
* Allow arbitrary face figures other than SGS (would probably need whole new program)
