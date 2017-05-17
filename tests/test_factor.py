#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests for the factor module in antitile.
"""

import unittest
from operator import mul
from functools import reduce
from itertools import product
from antitile import factor

class TestFactor(unittest.TestCase):
    def setUp(self):
        self.domains = [factor.Integer,
                       factor.Gaussian,
                       factor.Eisenstein,
                       factor.Nietsnesie]

    def test_multiplication(self):
        for d in self.domains:
            for v1 in product(range(-5, 10), repeat=2):
                for v2 in product(range(-5, 10), repeat=2):
                    number1 = d(*v1)
                    number2 = d(*v2)
                    px = number1 * number2
                    anorm1 = number1.anorm()
                    anorm2 = number2.anorm()
                    self.assertEqual(px.anorm(), anorm1*anorm2)

    def test_factoring(self):
        for d in self.domains:
            for v in product(range(-5, 100), repeat=2):
                number = d(*v)
                fx = number.factor()
                backcalc_anorm = reduce(mul, (f.anorm() for f in fx), 1)
                self.assertEqual(number.anorm(), backcalc_anorm)
                unit = d(1)
                backcalc = reduce(mul, fx, unit)
                self.assertEqual(number, backcalc)

def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(factor))
    return tests


if __name__ == '__main__':
    unittest.main()
