# -*- coding: utf-8 -*-
"""
Objects to represent elements of Euclidean domains that can be expressed as
a + b u, where u is some unit and a and b are integers. Also implements
the Euclidean algorithm to factor those elements.
    Gaussian: u = j
    Eisenstein: u = exp(2pi/3)
    Steineisen: u = exp(pi/3) (isomorphic to Eisenstein but useful)
    (regular) Integers: b = 0, u = -1

See https://en.wikipedia.org/wiki/Euclidean_domain for more info.
"""

import abc
import math
from math import gcd

def smallest_prime_factor(n):
    """The smallest factor of n larger than 1. Will be n if n is prime.
    >>> smallest_prime_factor(7)
    7
    >>> smallest_prime_factor(15)
    3
    """
    if abs(n) <= 3:
        return n
    for i in range(2, n + 1):
        if n%i == 0:
            return i

class EuclideanInteger(metaclass=abc.ABCMeta):
    """Abstract base class for elements in Euclidean domains.

    Args:
        a: First component of the element.
        b: Second component (default: 0)

    Class-level attributes:
        symbol: Symbol to use for the unit (j, w, etc.)
        unit: Numeric form of the unit
        mod: Argument for modulus in Euclidean algorithm.
        unitmap: A map from units to a nice printable form for them."""
    symbol = NotImplemented
    unit = NotImplemented
    mod = NotImplemented
    unitmap = NotImplemented

    def __init__(self, a, b=0):
        self.a = int(a)
        self.b = int(b)

    @property
    def tuple(self):
        """The components of the element, in tuple form."""
        return (self.a, self.b)

    def __complex__(self):
        return self.a + self.unit * self.b

    def __eq__(self, other):
        return (self.a == other.a and self.b == other.b
                and self.unit == other.unit)

    def __lt__(self, other):
        return self.anorm() < other.anorm()

    def __bool__(self):
        return bool(self.a or self.b)

    def __hash__(self):
        return hash((self.tuple))

    def __sub__(self, other):
        return type(self)(self.a - other.a, self.b - other.b)

    def __neg__(self):
        return type(self)(-self.a, -self.b)

    def __pos__(self):
        return self

    def __abs__(self):
        return math.sqrt(self.anorm())

    @abc.abstractmethod
    def __mul__(self, other):
        return NotImplemented

    @abc.abstractmethod
    def __divmod__(self, other):
        return NotImplemented, NotImplemented

    @abc.abstractmethod
    def conjugate(self):
        """Conjugate of the element."""
        return NotImplemented

    @abc.abstractmethod
    def anorm(self):
        """Algebraic norm of the element."""
        return NotImplemented

    @staticmethod
    def _testfactor(p):
        """Test factor for use in factoring algorithm"""
        return NotImplemented

    def __floordiv__(self, other):
        return divmod(self, other)[0]

    def __mod__(self, other):
        return divmod(self, other)[1]

    def __str__(self):
        a, b = self.a, self.b
        if self.anorm() == 1:
            return self._print_unit()
        elif b == 0:
            return str(a)
        elif a == 0:
            return '{}'.format(b) + self.symbol
        elif b < 0:
            return '({} - {}'.format(a, abs(b)) + self.symbol + ')'
        else:
            return '({} + {}'.format(a, b) + self.symbol + ')'

    def _print_unit(self):
        return self.unitmap[self.tuple]

    def __repr__(self):
        return str(self)


    def gcd(self, other):
        """Greatest common denominator of this and another element."""
        #print(self, other)
        if not other:
            return self
        else:
            return other.gcd(self % other)

    def _factor(self):
        """Workhorse function for factoring. This is the Euclidean algorithm
        that Euclidean domains get their name from."""
        constructor = type(self)
        an = self.anorm()
        p = smallest_prime_factor(an)
        if an <= 1 or an == p:
            return [self]
        elif p%self.mod == self.mod-1:
            factor = constructor(p, 0)
            left = self//factor
        else:
            testfactor = self._testfactor(p)
            pfactor = constructor(p, 0)
            factor = testfactor.gcd(pfactor)
            if factor.anorm() == 1:
                raise Exception('failed to find factor')
            left, rem = divmod(self, factor)
            if rem:
                factor = factor.conjugate()
                left, rem = divmod(self, factor)
            if rem:
                raise Exception('failed to find divisor')
        return [factor] + left._factor()

    def factor(self):
        """Factor this element."""
        constructor = type(self)
        one = constructor(1)
        factors = self._factor()
        nf = [f.normal_form()[0] for f in factors if f.anorm() > 1]
        if nf:
            backcalc = nf[-1]
            for i in range(len(nf)-1):
                f = nf[i]
                if f.anorm() > 1:
                    backcalc *= f
        else:
            backcalc = one
        unit = self//backcalc
        if unit != one or not nf:
            nf.append(unit)
        return sorted(nf)

    def normal_form(self):
        """Returns normal form of this element and how many multiplications
        of the unit it takes to get from normal form to this form."""
        if (self.a > 0 and self.b >= 0) or (self.a == 0 and self.b == 0):
            return self, 0
        else:
            unit = type(self)(0, 1)
            nm = self//unit
            nf, n = nm.normal_form()
            return nf, 1 + n

class Integer(EuclideanInteger):
    """Class representing regular, run-of-the-mill integers.

    Args:
        a: The integer.
        b: Ignored, only present so the constructor has the same signature
            as the other classes that inherit from EuclideanInteger.

    See docs for EuclideanInteger for more information.
        """
    mod = 0
    unit = -1
    b = 0
    symbol = ''
    unitmap = {(1, 0): '1',
               (-1, 0): '-1'}

    def __init__(self, a, b=0):
        super().__init__(a, 0)

    def anorm(self):
        return abs(self.a)

    def conjugate(self):
        return self

    def __mul__(self, other):
        return Integer(self.a*other.a)

    def __truediv__(self, other):
        return Integer(self.a/other.a)

    def __divmod__(self, other):
        result = divmod(self.a, other.a)
        return Integer(result[0]), Integer(result[1])

    def gcd(self, other):
        return Integer(gcd(self.a, other.a))

    def factor(self):
        an = self.anorm()
        p = smallest_prime_factor(an)
        if an <= 1 or an == p:
            return [self]
        else:
            left = Integer(self.a/p)
            return [Integer(p)] + left.factor()

    def normal_form(self):
        if self.a < 0:
            return Integer(-self.a), 1
        else:
            return self, 0

class Gaussian(EuclideanInteger):
    """Class representing the Gaussian integers, a + bj where j is the
    imaginary unit sqrt(-1) and a and b are integers."""
    mod = 4
    unit = 1j
    symbol = 'j'
    unitmap = {(1, 0): '1',
               (0, 1): 'j',
               (-1, 0): '-1',
               (0, -1): '-j'}

    def anorm(self):
        return self.a**2 + self.b**2

    def conjugate(self):
        a, b = self.a, self.b
        return type(self)(a, -b)


    def __mul__(self, other):
        a, b = self.a, self.b
        c, d = other.a, other.b
        return type(self)(a*c - b*d, b*c + a*d)

    def __divmod__(self, other):
        a, b = self.a, self.b
        c, d = other.a, other.b
        an = other.anorm()
        if an == 0:
            raise ZeroDivisionError()
        x = (a*c+b*d)/an
        y = (b*c-a*d)/an
        x = math.floor(x) if x > 0 else math.ceil(x)
        y = math.floor(y) if y > 0 else math.ceil(y)
        div = type(self)(x, y)
        mod = self - div*other
        return div, mod

    @classmethod
    def _testfactor(cls, p):
        exponent = (p - 1)//2
        for n in range(1, p):
            ksq = n**exponent + 1
            if ksq % p == 0:
                k = n**(exponent//2)
                break
        else:
            msg = "couldn't find test factor for {} in the Gaussians".format(p)
            raise Exception(msg)
        return cls(k, 1)

class Eisenstein(EuclideanInteger):
    """Class representing the Eisenstein integers, a + bw where
    w = exp(j*2pi/3) and a and b are integers.

    See docs for EuclideanInteger for more information."""
    mod = 3
    unit = (-1 + 1j*math.sqrt(3))/2
    symbol = 'w'
    unitmap = {(1, 0): '1',
               (1, 1): '-w^2',
               (0, 1): 'w',
               (-1, 0): '-1',
               (-1, -1): 'w^2',
               (0, -1): '-w'}

    def anorm(self):
        a, b = self.a, self.b
        return a**2 - a*b + b**2

    def conjugate(self):
        a, b = self.a, self.b
        return Eisenstein(a - b, -b)

    def __mul__(self, other):
        a, b = self.a, self.b
        c, d = other.a, other.b
        return Eisenstein(a*c - b*d, b*c + a*d - b*d)

    def __divmod__(self, other):
        a, b = self.a, self.b
        c, d = other.a, other.b
        an = other.anorm()
        if an == 0:
            raise ZeroDivisionError()
        x = (a*c-a*d+b*d)/an
        y = (b*c-a*d)/an
        x = math.floor(x) if x > 0 else math.ceil(x)
        y = math.floor(y) if y > 0 else math.ceil(y)
        div = Eisenstein(x, y)
        mod = self - div*other
        return div, mod

    def normal_form(self):
        if ((self.a > 0 and self.b >= 0 and self.b <= self.a) or
                (self.a == 0 and self.b == 0)):
            return self, 0
        else:
            unit = Eisenstein(1, 1)
            nm = self//unit
            nf, n = nm.normal_form()
            return nf, 1 + n

    @classmethod
    def _testfactor(cls, p):
        for k in range(1, p**2):
            kf = k**2 - k + 1
            if kf % p == 0:
                break
        else:
            msg = "couldn't find test factor for {} in the Eisensteins".format(p)
            raise Exception(msg)
        return cls(k, 1)


class Steineisen(EuclideanInteger):
    """Class representing the Steineisen integers,
    a + bu where u = exp(j*pi/3) and a and b are integers. These are
    isomorphic to the Eisenstein integers but this parameterization
    is more convenient for use with gcopoly.py.

    See docs for EuclideanInteger for more information.
    """

    mod = 3
    unit = (1 + 1j*math.sqrt(3))/2
    symbol = 'u'
    unitmap = {(1, 0): '1',
               (0, 1): 'u',
               (-1, 1): 'u^2',
               (-1, 0): '-1',
               (0, -1): '-u',
               (1, -1): '-u^2'}

    def anorm(self):
        a, b = self.a, self.b
        return a**2 + a*b + b**2

    def conjugate(self):
        a, b = self.a, self.b
        return Steineisen(a + b, -b)

    def __mul__(self, other):
        a, b = self.a, self.b
        c, d = other.a, other.b
        return Steineisen(a*c - b*d, b*c + a*d + b*d)

    def __divmod__(self, other):
        a, b = self.a, self.b
        c, d = other.a, other.b
        an = other.anorm()
        if an == 0:
            raise ZeroDivisionError()
        x = (a*c + a*d + b*d)/an
        y = (b*c - a*d)/an
        x = math.floor(x) if x > 0 else math.ceil(x)
        y = math.floor(y) if y > 0 else math.ceil(y)
        div = Steineisen(x, y)
        mod = self - div*other
        return div, mod

    @classmethod
    def _testfactor(cls, p):
        for k in range(0, p**2):
            kf = k**2 + k + 1
            if kf % p == 0:
                break
        else:
            msg = "couldn't find test factor for {} in the Steineisens".format(p)
            raise Exception(msg)
        return cls(k, 1)
