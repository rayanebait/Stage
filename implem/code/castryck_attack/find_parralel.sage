from time import time
from argparse import ArgumentParser


from sage.schemes.elliptic_curves.hom_frobenius import EllipticCurveHom_frobenius as frob
from richelot_aux import *

from utilities.supersingular import compute_point_order_D, torsion_basis
from utilities.find_prime import find_prime_with_torsion

import utilities.elliptic_curve as ec



parser = ArgumentParser()
parser.add_argument("-t", "--two_torsion", default = "128")
parser.add_argument("-T", "--three_torsion", default = "100")

args = parser.parse_args()

"""
data:
    -N is the rational torsion available
    -d is the isogeny degree
    -D is the torsion we take
"""
e2=Integer(args.two_torsion)
D = 2**(3*e2)
d = 3
N = d*D

e3 = ceil(log(2,3)*e2)

p, cofactor = find_prime_with_torsion(N)
while cofactor%3 == 0:
    p, cofactor = find_prime_with_torsion(N, cofactor+1)

print(f"Found prime {p}=2**(3*{e2})*3*({factor(cofactor)})-1\n\n")

Fp2d = GF(p**2, modulus=[1,0,1], name="i")
E0=EllipticCurve(Fp2d, [1,0])

P0, Q0 = torsion_basis(E0, 2**(3*e2))
PQ0 = P0+Q0

print(f"Now computing radical isogenies:\n\n")
t3 = time()

zeta3 = (-1 + Fp2d(-3).sqrt())/2

E1, xs = ec.chain_3radials(E0, [P0.xy()[0], Q0.xy()[0], PQ0.xy()[0]],zeta3, 2*e3)
t4 = time()

print(f"Computed new isogeny\
    {E0.j_invariant()}--->{E1.j_invariant()}\n")
print(f"In time: {t4-t3}")

# transform to a Montgomery curve
x4 = xs[0]
for _ in range((3*e2)-2):
    x4 = ec.x_onlyDoubling(E1, x4)
E1_, xs = ec.WeierstrassToMontgomery(E1, x4, xs, x_only=True)
Pd = E1_.lift_x(xs[0])
Qd = E1_.lift_x(xs[1])
if not (Pd + Qd).xy()[0] == xs[2]:
    Qd = -Qd


R = Fp2d['X']
X = R.gen()

a = E1_.a_invariants()

f = R([0, a[3], a[1], 1])
alpha = f.roots()[1][0]
alpha1 = alpha[1]
alpha0 = alpha[0]

assert alpha == Fp2d(alpha0 + Fp2d.gen()*alpha1)

f1 = R([-1, 2*alpha0/alpha1, 1])
f2 = R([-1,-2*alpha0/alpha1, 1])
f3 = R([-1,-2*(alpha0*(alpha0**2+alpha1**2-1))/(alpha1*(alpha0**2+alpha1**2+1)), 1])

print(f"f1={f1}\nf2={f2}, f3={f3}")

g = R(f1*f2*f3)
print(f"g={g}\n")

rs = [r[0] for r in g.roots()]

assert len(rs) == 6
isogeneous_prods = []

print(f"Searching for products (2,2)-isogeneous to the weil restriction of E_alpha\n")
for r1 in rs:
    for r2 in rs:
        if r1 != r2:
            for r3 in rs:
                if r3 not in [r1, r2]:
                    for r4 in rs:
                        if r4 not in [r1,r2,r3]:
                            g1 = R((X-r1)*(X-r2))
                            g2 = R((X-r3)*(X-r4))
                            g3 = R(g/(g1*g2))
                            try:
                                _, (E, E_) = FromJacToProd(g1, g2, g3)
                                isogeneous_prods.append(\
                                                        E_.j_invariant()\
                                                         )
                                isogeneous_prods.append(\
                                                        E.j_invariant(),\
                                                         )
                            except:
                                print(".", end ='')


s = list(set(isogeneous_prods))
assert E1_.j_invariant() in s

print(f"\n\nFound {list(set(isogeneous_prods))}")

