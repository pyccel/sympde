Linearity
*********

- must check that BilinearForm and LinearForm are really linear, by calling them and then applying atomize
- this verification must be done at the __new__ as long as the flag check_linearity is True. The later is an argument of the __new__ method for BilinearForm and LinearForm


Weak formulation
****************

- add the concept of weak formulation
- a weak formulation does not need to be linear/bilinear.
- add linearize function that returns either LinearForm or BilinearForm

Kron
****

the kronecker features must be upgraded as well as the tensorize function


