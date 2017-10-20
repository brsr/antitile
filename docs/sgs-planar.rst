Similar Grid Subdivision of planar tilings
==========================================

Square and triangular lattices
------------------------------

Gaussian and Eisenstein (Nietsnesie) integers
---------------------------------------------
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
   * - Algebraic Norm
       $|x|^2=x\overline x$
     - $a^2 + b^2$
     - $a^2 - ab+ b^2$
     - $a^2 + ab+ b^2$
   * - Multiplication
       $(a+bz) (c+dz)$
     - $(ac-bd) + (bc+ad)i$
	 - $(ac-bd)+(bc+ad-bd)w$
	 - $(ac-bd)+(bc+ad+bd)u$
   * - Division
       $(a+bz)/(c+dz)$
     - $\left({ac + bd \over c^2 + d^2}\right) + \left( {bc - ad \over c^2 + d^2} \right)i</math>$
     - $\left({ac-ad+bd \over c^2 - cd+ d^2}\right) +  \left({bc-ad \over c^2 - cd+ d^2}\right)w$
	 - $\left({ac+ad+bd \over c^2 + cd+ d^2}\right) +  \left({bc-ad \over c^2 + cd+ d^2}\right)u$
   * - Divisor $n$ for prime test 
     - 4
	 - 3
	 - 3

Subdivision of a triangle/square using a lattice
------------------------------------------------

Composition of subdivisions
---------------------------

Associates and normal form
--------------------------

Chirality and conjugation
-------------------------
$x = a + bz$, $|z| = 1$
$z\overline x = z (a + b \overline z) = az + b z\overline z = b + az$

As a Euclidean domain
---------------------
prime if:
* $x = pz$, where z is a unit, p is a natural prime, and $p = -1\mod n$, or
* $|x|^2$ is a natural prime
