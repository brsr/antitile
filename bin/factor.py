#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Factor integers in various Euclidean domains
"""
import argparse
from collections import Counter
from antitile import factor

DOMAINS = {"i": factor.Integer,
           "g": factor.Gaussian,
           "e": factor.Eisenstein,
           "n": factor.Nietsnesie}

DESC = """Factor various algebraic integers from Euclidean domains.
Input is of the form a + b*u, where u is the unit in whichever domain.
([g]aussian: the imaginary unit, [e]isenstein: third root of unity,
[n]ietsnesie: sixth root of unity, [i]nteger: zero)"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESC)
    parser.add_argument("a", type=int, help="First part of integer")
    parser.add_argument("b", nargs="?", default=0, type=int,
                        help="Second part of integer")
    parser.add_argument("-d", choices=DOMAINS, default="g",
                        help="Specify domain. Defaults to [g]aussian.")
    args = parser.parse_args()
    if args.a == 0 and args.b == 0:
        print('0')
        exit()
    constructor = DOMAINS[args.d]
    x = constructor(args.a, args.b)
    factors = x.factor()
    nfcount = Counter(factors)
    gen = ('{}'.format(x) if y == 1 else '{}^{}'.format(x, y)
           for x, y in sorted(nfcount.items()))
    print(' * '.join(gen))
