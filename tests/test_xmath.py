#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests for the xmath module in antitile.
"""

import unittest
import doctest
from antitile import xmath

def load_tests(loader, tests, ignore):
    optionflags = (doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS |
                   doctest.IGNORE_EXCEPTION_DETAIL)
    tests.addTests(doctest.DocTestSuite(xmath, optionflags=optionflags))
    return tests

if __name__ == '__main__':
    unittest.main()
