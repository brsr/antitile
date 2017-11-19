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
     - :math:`a + b i`
     - :math:`a + b w`
     - :math:`a + b u`
   * - Adjoined element
     - :math:`i =\sqrt{-1}`
     - :math:`w = \frac{1}{2}(-1 + i\sqrt 3) = e^{2\pi i/3}`
     - :math:`u = \frac{1}{2}(1 + i\sqrt 3) = e^{\pi i/3}`
   * - Units
       :math:`(a, b)`
     - :math:`(1, 0) = 1`,
       :math:`(0, 1) = i`,
       :math:`(-1, 0) = -1`,
       :math:`(0, -1) = -i`
     - :math:`(1, 0) = 1`
       :math:`(1, 1) = -w^2 = 1 + w`,
       :math:`(0, 1) = w`,
       :math:`(-1, 0) = -1`,
       :math:`(-1, -1) = w^2 = -1-w`,
       :math:`(0, -1) = -w`
     - :math:`(1, 0) = 1`
       :math:`(0, 1) = u`,
       :math:`(-1, 1) = u^2 = u-1`,
       :math:`(-1, 0) = -1`,
       :math:`(0, -1) = -u`,
       :math:`(1, -1) = -u^2=1-u`
   * - Algebraic Norm (Euclidean Function)
       :math:`|x|^2=x\overline x`
     - :math:`a^2 + b^2`
     - :math:`a^2 - ab + b^2`
     - :math:`a^2 + ab + b^2`
   * - Multiplication
       :math:`(a+bz) (c+dz)`
     - :math:`(ac-bd) + (bc+ad)i`
     - :math:`(ac-bd)+(bc+ad-bd)w`
     - :math:`(ac-bd)+(bc+ad+bd)u`
   * - Division
       :math:`(a+bz)/(c+dz)`
     - :math:`\left({ac + bd \over c^2 + d^2}\right) +
       \left( {bc - ad \over c^2 + d^2} \right)i`
     - :math:`\left({ac-ad+bd \over c^2 - cd+ d^2}\right) +
       \left({bc-ad \over c^2 - cd+ d^2}\right)w`
     - :math:`\left({ac+ad+bd \over c^2 + cd+ d^2}\right) +
       \left({bc-ad \over c^2 + cd+ d^2}\right)u`
   * - Divisor `t` for prime test
     - 4
     - 3
     - 3
   * - Normal elements
     - :math:`a > 0`, `b >= 0`
     - :math:`a > 0`, `b >= 0`, `b <= a`
     - :math:`a > 0`, `b >= 0`

Elements :math:`x` and :math:`y` of a domain such that :math:`x = zy`, where 
:math:`z` is a unit, are called associates. "Normal elements" are defined so 
that we have a "normal" form for associates. All non-zero elements of the 
domain can be expressed as a normal integer times a unit. All the units have 
:math:`|z| = 1`. The normal integers define a wedge of space in 
:math:`\mathbb{C}`, and the units define the rotation of that space around zero.

All three of these are Euclidean domains. That means that all non-zero
elements have a unique factorization into prime elements and units, and
what's more, they can be factored using an extension of the Euclidean
division algorithm on the natural integers. A prime element can only be
divided by itself, an associate of itself, or a unit. More specifically,
non-zero elements of these domains can be factored as:

.. math::
   x = zp_1^{n_1}p_2^{n_2} \cdots p_k^{n_k}

where :math:`z` is a unit, :math:`n_i` is a natural integer, and :math:`p_i` 
is a normal prime element of the domain.

An element of these domains is prime if:

- :math:`x = pz`, where :math:`z` is a unit, :math:`p` is a natural prime, 
  and :math:`p = -1\mod t`, or
- :math:`|x|^2` is a natural prime.

Subdivision of a triangle/square using a similar grid
-----------------------------------------------------
We'll use the operator :math:`\Delta_{a,b}` to denote triangle subdivision, and
:math:`\Box_{a,b}` for quadrilateral subdivision. Repeated application of an
operator will be denoted with a superscript, e.g. :math:`\Delta^2_{a,b}`. Operators
apply from right to left, like Conway operators. We'll also denote the base
triangular and square lattices as :math:`\Delta` and :math:`\Box`.

The points for the triangle are :math:`(0,0)`, :math:`(a,b)`, and 
:math:`(-b, a+b)`. The points for the square are :math:`(0,0)`, 
:math:`(a,b)`, :math:`(a-b, a+b)`, and :math:`(-b,a)`.

In the terminology of geodesic domes,

- Class I operators have :math:`b=0`.
- Class II operators have :math:`b=a`.
- Class III operators are all others, and occur in chiral pairs.
  :math:`\Delta_{a,b}` is the chiral pair of :math:`\Delta_{b,a}`.

:math:`\Delta_{1,0}` and :math:`\Box_{1,0}` are identity operators: 
where defined, they cause no change to the subdivision.

Composition of subdivisions
---------------------------
Multiplication in the complex numbers can be interpreted geometrically as
rotation and scaling.
:math:`\Delta_{a,b} \Delta_{c,d} = \Delta_{(a,b)*(c,d)}`
:math:`\Box_{a,b} \Box_{c,d} = \Box_{(a,b)*(c,d)}`

Composition of these operators corresponds to multiplication in the
corresponding domain.

Associates and normal form
--------------------------
Some of the multiplications defined in the last section may result in a
non-normal value. You may have also noticed that there is more than one way
to draw these subdivisions on a lattice: for instance, the triangle defined
by :math:`(0,0)`, :math:`(-a,-b)`, and :math:`(b, -a-b)` is just a rotated 
version of the triangle defined by :math:`(0,0)`, :math:`(a,b)`, and 
:math:`(-b, a+b)`. Symbolically, if :math:`x = (a,b)` and :math:`y = xz`, 
where :math:`z` is a unit, :math:`\Delta_x` and :math:`\Delta_y` define the 
same subdivision. To restore uniqueness, we can refine the definition of 
define :math:`\Delta_x`: :math:`\Delta_x` corresponds to the equivalence class 
of :math:`x` and all its associated elements in the Nietsnesie integers, and 
:math:`\Box_x` corresponds to the equivalence class of :math:`x` and all its 
associated elements in the Gaussian integers. The "normal element" is the 
preferred label :math:`x` for that equivalence class.

Chirality and conjugation
-------------------------
As mentioned earlier, :math:`\Delta_{a,b}` is the chiral pair of 
:math:`\Delta_{b,a}`. In the Gaussian and Nietsnesie domains, :math:`(a,b)` 
is an associate element to the complex conjugate of :math:`(b,a)`. 
Therefore, complex conjugation corresponds to reflection of an operator.

Proof:

* :math:`x = a + bz`, :math:`z` a unit (so :math:`|z| = 1`)
* :math:`xz` is in the equivalence class of :math:`z`
* :math:`z\overline x = z (a + b \overline z) = az + b z\overline z = b + az`

Interestingly, the notation used in Conway notation for the chiral pair of an
operator is an overbar, also used to denote complex conjugation.

A consequence of this is that the composed operator 
:math:`\Delta_{a,b}\Delta_{b,a}` will always be a Class I operator, since a 
complex number times its conjugate is a real number. For every Class III 
:math:`(a,b)`, there will exist three operators 
:math:`\Delta_{a,b}\Delta_{b,a}`, :math:`\Delta_{a,b}\Delta_{a,b}`, and
:math:`\Delta_{b,a}\Delta_{b,a}` with the same algebraic norm. 
The same holds for :math:`\Box_{a,b}`.

Operator factorization
----------------------
Since composition of operators corresponds to multiplication in the
corresponding domain, and these domains are Euclidean, we can "factor"
operators into smaller operators.

When :math:`(a,b) = a + bi` is an element of the Gaussian integers, and

.. math::

   (a,b) = x = z p_1^{n_1}p_2^{n_2} \cdots p_k^{n_k}

then

.. math::

    \Box_{a,b} = \Box^{n_1}_{p_1}\Box^{n_2}_{p_2}\cdots\Box^{n_k}_{p_k}

Similarly, when :math:`(a,b) = a + bu` is an element of the Nietsnesie 
integers, and

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
:math:`\Box_{a,b}\Box_{c,d}` may not equal :math:`\Box_{c,d}\Box_{a,b}` in 
geometry. This actually turns out to be useful, because it allows us to 
tweak the geometry.