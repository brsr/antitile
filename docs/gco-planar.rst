Goldberg-Coxeter Operation on planar tilings
==========================================

Square and triangular lattices
------------------------------

- triangle dual of hex
- stays similar even if skewed

Relation to Gaussian and Eisenstein integers
--------------------------------------------
The vertices of square and triangular lattices match with the elements of
certain integral domains [#]_ in the complex plane. The Gaussian integers form
a square lattice, and the Eisenstein integers form a triangular lattice. For
reasons that will be explained later, it is more convenient to parameterize the
Eisenstein integers with :math:`u = e^{\pi i/3}`
instead of the usual :math:`w = e^{2\pi i/3}`.
This defines exactly the same points in the complex plane, but for the sake of 
clarity, we'll call the version in terms of :math:`u` the Steineisen integers. 
The relationship between the two is simple: :math:`a + b w = a - b + bu`.

.. list-table::
   :header-rows: 1

   * - Type
     - Gaussian
     - Eisenstein
     - Steineisen
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
       :math:`(1, 1) = -w^2 = 1 + w = u`,
       :math:`(0, 1) = w`,
       :math:`(-1, 0) = -1`,
       :math:`(-1, -1) = w^2 = -1-w = -u`,
       :math:`(0, -1) = -w`
     - :math:`(1, 0) = 1`
       :math:`(0, 1) = u`,
       :math:`(-1, 1) = u^2 = u-1 = w`,
       :math:`(-1, 0) = -1`,
       :math:`(0, -1) = -u`,
       :math:`(1, -1) = -u^2=1-u = -w`
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
   * - Breakdown polygon 
     - 0, :math:`x`, :math:`x(1+i)`, :math:`xi`
     - 0, :math:`x`, :math:`x(1+w)`
     - 0, :math:`x`, :math:`xu`
   * - Breakdown polygon center
     - :math:`x(1+i)/2`
     - :math:`x(2+w)/3`
     - :math:`x(1+u)/3`
     
Elements :math:`x` and :math:`y` of a domain such that :math:`x = zy`, where 
:math:`z` is a unit, are called associates. "Normal elements" are defined so 
that we have a "normal" form for associates. All non-zero elements of the 
domain can be expressed as a normal element times a unit. All the units have 
:math:`|z| = 1`, so the associate elements can be thought of as the elements
rotated around the origin (by steps of 90 degrees for the Gaussians, 60 
degrees for the Steineisen/Eisenstein.) 

An equivalence class can be defined from the associate relationship (excluding 
the zero element): all associated elements belong to the same equivalence 
class. The "normal element" is the preferred label :math:`x` for that 
equivalence class. [#]_ 

All three of these domains are Euclidean domains. That means that all non-zero
elements have a unique factorization into prime elements and units, and
what's more, they can be factored using an extension of the Euclidean
division algorithm on the natural integers. A prime element can only be
divided by itself, an associate of itself, or a unit. More specifically,
non-zero elements of these domains can be factored as:

.. math::
   x = zp_1^{n_1}p_2^{n_2} \cdots p_k^{n_k}

where :math:`z` is a unit, :math:`n_i` is a natural integer, 
and :math:`p_i` is a normal prime element of the domain.

An element of these domains is prime if:

- :math:`x = pz`, where :math:`z` is a unit, :math:`p` is a natural prime, 
  and :math:`p = -1\mod t`, or
- :math:`|x|^2` is a natural prime.

Defining the operators
----------------------
We define :math:`\Delta(a,b)` for the triangular GC operator and 
:math:`\Box(a,b)` for the quadrilateral GC operator.
:math:`\Delta(a,b)` corresponds to the equivalence class of associates in the 
Steineisen integers, and :math:`\Box(a,b)` corresponds to the equivalence 
class of associates in the Gaussian integers. There are some very nice 
consequences of this:

* Composition of operators corresponds to multiplication of complex numbers.
  If :math:`\Delta(a,b)\Delta(c,d) = \Delta(e,f)`, then 
  :math:`(a + bu)(c + du) = z(e + fu)` for some unit z. A similar relation
  holds for :math:`\Box`.
* Since elements of these domains can be factored, operators can be "factored"
  into a sequence of operators. When :math:`x = a + bi` is an element of 
  the Gaussian integers, and :math:`x = z p_1^{n_1}p_2^{n_2} \cdots p_k^{n_k}`,
  then :math:`\Box(a,b) = 
  \Box^{n_1}_{p_1}\Box^{n_2}_{p_2}\cdots\Box^{n_k}_{p_k}`, and similarly for
  the Steineisen integers and :math:`\Delta(a,b)`.

:math:`\Delta_{1,0}` and :math:`\Box_{1,0}` are identity operators: 
where defined, they cause no change to the subdivision.

In the terminology of geodesic domes,

- Class I operators have :math:`b=0`.
- Class II operators have :math:`b=a`.
- Class III operators are all others, and occur in chiral pairs.
  :math:`\Delta(a,b)` is the chiral pair of :math:`\Delta(b,a)`. [#]_ 

Chiral pairs correspond to conjugation of the complex number :math:`a+bi` or 
:math:`a+bu`. That is, there exists some unit z such that 
:math:`z(a-bi) = b+ai`, or :math:`z(a+b\overbar{u}) = b+au`. 
Interestingly, the notation used in Conway notation for the chiral pair of an
operator is an overbar, also used to denote complex conjugation.
A consequence of this relation between chiral pairs is that the composed 
operator :math:`\Delta(a,b)\Delta_{b,a}` will always be a Class I operator, 
since a complex number times its conjugate is a real number. 
For every Class III :math:`(a,b)`, there will exist three operators 
:math:`\Delta(a,b)\Delta_{b,a}`, :math:`\Delta(a,b)\Delta(a,b)`, and
:math:`\Delta_{b,a}\Delta_{b,a}` with the same algebraic norm. 
The same holds for :math:`\Box(a,b)`.

On the plane, equality of operators holds in both a topological and geometric 
sense: if :math:`\Box(a,b)\Box(c,d) = \Box(e,f)`, then the left-hand side
results in the same placement of vertices on the plane. This topological
equality carries over to non-planar tilings and polyhedra, but not necessarily
the geometric sense. Furthermore, commutativity and associativity hold for 
topology, but not necessarily geometry. So considering the exact position of 
vertices in space:

* :math:`\Delta(a,b)\Delta(c,d)` may not equal :math:`\Delta(e,f)`, 
  where :math:`(e,f)` is the normal form of the product of 
  :math:`(a,b)` and :math:`(c,d)`
* :math:`\Delta(a,b)\Delta(c,d)` may not equal :math:`\Delta(c,d)\Delta(a,b)`
* :math:`\Delta(a,b)(\Delta(c,d)\Delta(e,f))` may not equal 
  :math:`(\Delta(a,b)\Delta(c,d))\Delta(e,f)`
  
This actually turns out to be useful, because it allows us to 
tweak the geometry.

Notes/References
----------------
.. [#] A type of ring.
.. [#] One benefit of using the Steineisen integers over the Eisenstein 
   integers is that the normal form for the Gaussian and Steineisen integers 
   is the same.
.. [#] Another benefit of using the Steineisen integers is that this 
   definition is the same for Gaussian and Steineisen integers, and 
   it reflects the symmetry of the operators. Were we to define this for 
   Eisenstein integers, class II would be :math:`b=2a`, and the chiral pairs
   would be :math:`(a, b)` and :math:`(a, a-b)`.)
