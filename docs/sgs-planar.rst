Similar Grid Subdivision of planar tilings
==========================================

Square and triangular lattices
------------------------------

- triangle dual of hex
- stays similar even if skewed

Gaussian and Eisenstein (Nietsnesie) integers
---------------------------------------------
The vertices of square and triangular lattices match with the elements of
certain integral domains [#]_ in the complex plane. The Gaussian integers form
a square lattice, and the Eisenstein integers form a triangular lattice. The
Nietsnesie [#]_ integers also form a triangular lattice, and in fact are
isomorphic to the Eisenstein integers: they are introduced here because they
are a more convenient parameterization of the triangular lattice for the
purposes of Similar Grid Subdivision.

.. [#] A type of ring
.. [#] "Eisenstein" backwards; a nonce name to distinguish the two.

.. list-table::
   :header-rows: 1

   * - Type
     - Gaussian
     - Eisenstein
     - Nietsnesie
   * - Form
     - $a + b i$
     - $a + b w$
     - $a + b u$
   * - Adjoined element
     - $i =\sqrt{-1}$
     - $w = \frac{1}{2}(-1 + i\sqrt 3) = e^{2\pi i/3}$
     - $u = \frac{1}{2}(1 + i\sqrt 3) = e^{\pi i/3}$
   * - Units
       $(a, b)$
     - $(1, 0) = 1$,
       $(0, 1) = i$,
       $(-1, 0) = -1$,
       $(0, -1) = -i$
     - $(1, 0) = 1$
       $(1, 1) = -w^2 = 1 + w$,
       $(0, 1) = w$,
       $(-1, 0) = -1$,
       $(-1, -1) = w^2 = -1-w$,
       $(0, -1) = -w$
     - $(1, 0) = 1$
       $(0, 1) = u$,
       $(-1, 1) = u^2 = u-1$,
       $(-1, 0) = -1$,
       $(0, -1) = -u$,
       $(1, -1) = -u^2=1-u$
   * - Algebraic Norm (Euclidean Function)
       $|x|^2=x\overline x$
     - $a^2 + b^2$
     - $a^2 - ab + b^2$
     - $a^2 + ab + b^2$
   * - Multiplication
       $(a+bz) (c+dz)$
     - $(ac-bd) + (bc+ad)i$
     - $(ac-bd)+(bc+ad-bd)w$
     - $(ac-bd)+(bc+ad+bd)u$
   * - Division
       $(a+bz)/(c+dz)$
     - $\left({ac + bd \over c^2 + d^2}\right) +
       \left( {bc - ad \over c^2 + d^2} \right)i</math>$
     - $\left({ac-ad+bd \over c^2 - cd+ d^2}\right) +
       \left({bc-ad \over c^2 - cd+ d^2}\right)w$
     - $\left({ac+ad+bd \over c^2 + cd+ d^2}\right) +
       \left({bc-ad \over c^2 + cd+ d^2}\right)u$
   * - Divisor $t$ for prime test
     - 4
     - 3
     - 3
   * - Normal elements
     - $a > 0$, $b >= 0$
     - $a > 0$, $b >= 0$, $b <= a$
     - $a > 0$, $b >= 0$

Elements $x$ and $y$ of a domain such that $x = zy$, where $z$ is a unit, are
called associates. "Normal elements" are defined so that we have a "normal"
form for associates. All non-zero elements of the domain can be expressed as
a normal integer times a unit. All the units have $|z| = 1$. The normal
integers define a wedge of space in $\complex$, and the units define the
rotation of that space around zero.

All three of these are Euclidean domains. That means that all non-zero
elements have a unique factorization into prime elements and units, and
what's more, they can be factored using an extension of the Euclidean
division algorithm on the natural integers. A prime element can only be
divided by itself, an associate of itself, or a unit. More specifically,
non-zero elements of these domains can be factored as:

    $ x = zp_1^{n_1}p_2^{n_2} \cdots p_k^{n_k}$

where $z$ is a unit, $n_i$ is a natural integer, and $p_i$ is a normal prime
element of the domain.

An element of these domains is prime if:

- $x = pz$, where $z$ is a unit, $p$ is a natural prime, and 
  $p = -1\mod t$, or
- $|x|^2$ is a natural prime.

Subdivision of a triangle/square using a similar grid
-----------------------------------------------------
We'll use the operator $\Delta_{a,b}$ to denote triangle subdivision, and
$\Box_{a,b}$ for quadrilateral subdivision. Repeated application of an
operator will be denoted with a superscript, e.g. $\Delta^2_{a,b}$. Operators
apply from right to left, like Conway operators. We'll also denote the base
triangular and square lattices as $\Delta$ and $\Box$.

The points for the triangle are (0,0), (a,b), and (-b, a+b). The points for
the square are (0,0), (a,b), (a-b, a+b), and (-b,a).

In the terminology of geodesic domes,

- Class I operators have b=0.
- Class II operators have b=a.
- Class III operators are all others, and occur in chiral pairs.
  $\Delta_{a,b}$ is the chiral pair of $\Delta_{b,a}$.

$\Delta_{1,0}$ and $\Box_{1,0}$ are identity operators: where defined, they
cause no change to the subdivision.

Composition of subdivisions
---------------------------
Multiplication in the complex numbers can be interpreted geometrically as
rotation and scaling.
$\Delta_{a,b} \Delta_{c,d} = \Delta_{(a,b)*(c,d)}
$\Box_{a,b} \Box_{c,d} = \Box_{(a,b)*(c,d)}

Composition of these operators corresponds to multiplication in the
corresponding domain.

Associates and normal form
--------------------------
Some of the multiplications defined in the last section may result in a
non-normal value. You may have also noticed that there is more than one way
to draw these subdivisions on a lattice: for instance, the triangle defined
by (0,0), (-a,-b), and (b, -a-b) is just a rotated version of the triangle
defined by (0,0), (a,b), and (-b, a+b). Symbolically, if $x = (a,b)$ and
$y = xz$, where $z$ is a unit, $\Delta_x$ and $\Delta_y$ define the same
subdivision. To restore uniqueness, we can refine the definition of define
$\Delta_x$: $\Delta_x$ corresponds to the equivalence class of x and all its
associated elements in the Nietsnesie integers, and $\Box_x$ corresponds to
the equivalence class of x and all its associated elements in the Gaussian
integers. The "normal element" is the preferred label $x$ for that
equivalence class.

Chirality and conjugation
-------------------------
As mentioned earlier, $\Delta_{a,b}$ is the chiral pair of $\Delta_{b,a}$.
In the Gaussian and Nietsnesie domains, $(a,b)$ is an associate element to
the complex conjugate of $(b,a)$. Therefore, complex conjugation corresponds
to reflection of an operator.

Proof:
|    $x = a + bz$, $z$ a unit (so $|z| = 1$)
|    $xz$ is in the equivalence class of $z$
|    $z\overline x = z (a + b \overline z) = az + b z\overline z = b + az$

Interestingly, the notation used in Conway notation for the chiral pair of an
operator is an overbar, also used to denote complex conjugation.

A consequence of this is that the composed operator $\Delta_{a,b}\Delta_{b,a}$
will always be a Class I operator, since a complex number times its conjugate
is a real number. For every Class III (a,b), there will exist three operators
$\Delta_{a,b}\Delta_{b,a}$, $\Delta_{a,b}\Delta_{a,b}$, and
$\Delta_{b,a}\Delta_{b,a}$ with the same algebraic norm. The same holds for
$\Box_{a,b}$.

Operator factorization
----------------------
Since composition of operators corresponds to multiplication in the
corresponding domain, and these domains are Euclidean, we can "factor"
operators into smaller operators.

When (a,b) = a + bi is an element of the Gaussian integers, and

.. math::

   (a,b) = x = z p_1^{n_1}p_2^{n_2} \cdots p_k^{n_k}

then

.. math::

    \Box_{a,b} = \Box^{n_1}_{p_1}\Box^{n_2}_{p_2}\cdots\Box^{n_k}_{p_k}

Similarly, when (a,b) = a + bu is an element of the Nietsnesie integers, and

.. math::

    (a,b) = x = z p_1^{n_1}p_2^{n_2} \cdots p_k^{n_k}

then

.. math::

   \Delta_{a,b} = \Delta^{n_1}_{p_1}\Delta^{n_2}_{p_2}\cdots\Delta^{n_k}_{p_k}

When we move to non-planar tilings, these equalities will hold in a 
topological sense, but not necessarily in a geometric sense. That is, 
elements of topology like the number and connectivity of elements will be 
the same between factored and unfactored operators, but the exact position 
of vertices and length of edges may differ. Furthermore, the operators 
commute in topology but not necessarily geometry: for a non-planar tiling, 
$\Box_{a,b}\Box_{c,d}$ may not equal $\Box_{c,d}\Box_{a,b}$ in geometry. 
This actually turns out to be useful, because it allows us to tweak the 
geometry.