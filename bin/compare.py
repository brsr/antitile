#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create OFF files for each method and a range of frequencies
"""
import argparse
from sys import stdin
from antitile import sgs, off, tiling, projection, xmath

def all_SGS(infile, base, frequency, proj, tweak):
    poly = sgs.SGS(base, frequency, proj, tweak)
    if proj in projection.PARALLEL and frequency not in [(1,0),(1,1),(2,0)]:
        results = [('1', 1)]
        for measure in sgs.MEASURES:
            k = sgs.optimize_k(poly, base, sgs.MEASURES[measure],
                               not tweak, True)
            results.append((measure, k))
    else:
        results = [('0', 0)]
    pars = sgs.parallels(poly, base, exact=not tweak)
    for label, k in results:
        print(frequency, proj, tweak, label, k)
        pv = poly.vertices + k*pars
        pv = xmath.normalize(pv)
        result = off.write_off(pv, poly.faces)
        result += '#frequency = {}\n'.format(frequency)
        result += '#projection = {}\n'.format(proj)
        result += '#k = {} ({})\n'.format(k, label)
        result += '#tweak = {}\n'.format(tweak)
        result += '#normalized = True\n'
        outfile = '{:02d}-{:02d}'.format(*frequency)
        outfile += '-{}'.format(proj)
        outfile += '-t-' if tweak else '-n-'
        outfile += '{}-{}-'.format(label, k)
        outfile += infile
        with open(outfile, 'w') as f:
            f.write(result)

METHODS_3 = [('flat', None),
               ('other', None),
               ('gc', False),
               ('gc', True),
               ('nslerp', False),
               ('nslerp', True)]

METHODS_4 = [('flat', None),
               ('other', None),
               ('gc', False),
               ('gc', True),
               ('nslerp', False),
               ('nslerp', True),
               ('nslerp2', False),
               ('nslerp2', True)]

METHODS_D3 = [('nslerp', False),
               ('nslerp', True),
               ('disk', False),
               ('disk', True)]

METHODS_D4 = [('nslerp', False),
               ('nslerp', True),
               ('nslerp2', False),
               ('nslerp2', True),
               ('disk', False),
               ('disk', True)]


def main():
    parser = argparse.ArgumentParser(description="Create OFF files for each "
                                     "method and a range of frequencies")
    parser.add_argument("filename", help="filename")
    parser.add_argument("-a", default=9, help="max frequency", type=int)

    args = parser.parse_args()
    file = open(args.filename)
    with file as f:
        vertices, faces, _fc, _e, _ec, _v, _vc = off.load_off(f)
    base = tiling.Tiling(vertices, faces)
    methods = METHODS_3
    for a in range(1, args.a + 1):
        for b in range(a + 1):
            for proj, tweak in methods:
                all_SGS(args.filename.split('/')[-1], base,
                        (a, b), proj, tweak)


if __name__ == "__main__":
    main()
