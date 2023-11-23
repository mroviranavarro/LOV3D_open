     Wigner 3j-, 6j-, 9j-symbols
     ===========================
Support both numerical and symbolic computations
in a wide range of quantum numbers;
vector element-wise computing is supported.

The method is based on the Racah formula in a general case
and a set of known partial-case equations.

=======================================================================
     AUTHOR
Vladimir Borisovich SOVKOV
St. Petersburg State University (Russia),
Shanxi University (China)

Licensing provisions: BSD

=======================================================================
     USAGE

      Wigner3j(j1,j2,j3,m1,m2,m3)
wig = Wigner3j(j1,j2,j3,m1,m2,m3,ifs,ifcb)

      Wigner6j(j1,j2,j3,j4,j5,j6)
wig = Wigner6j(j1,j2,j3,j4,j5,j6,ifs,ifcb)

      Wigner9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)
wig = Wigner9j(j1,j2,j3,j4,j5,j6,j7,j8,j9,ifs,ifcb)


See details in the comments to the program codes.

=======================================================================
     DEFINITIONS:
(1)
R. N. Zare, Angular Momentum: Understanding Spatial Aspects in Chemistry and
Physics, John Wiley & Sons, New York–Chichester–Brisbane–Toronto–Sigapore, 1988. URL:
026http://eu.wiley.com/WileyCDA/WileyTitle/productCd-0471858927.ht
(2)
https://en.wikipedia.org/wiki/3-j_symbol
(3)
https://en.wikipedia.org/wiki/Clebsch–Gordan_coefficients
(4)
https://en.wikipedia.org/wiki/6-j_symbol
(5)
https://en.wikipedia.org/wiki/9-j_symbol

=======================================================================
     CONTENTS

readme.txt - this file

     MAIN PROGRAMS
Wigner3j.m - Wigner 3j-symbols or Clebsch-Gordan coefficients;
Wigner6j.m - Wigner 6j-symbols or corresponding recoupling coefficients;
Wigner9j.m - Wigner 9j-symbols or corresponding recoupling coefficients;

     AUXILIARY PROGRAMS
ArranA.m   - a simultaneous mutually adjusted sorting of several numerical vectors;
sper.m     - calculation of a permutation parity;
if3jc.m    - check of triangle rules and other necessary conditions on quantum numbers;
tc.m       - triangular coefficients.
