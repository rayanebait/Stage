from time import time
from argparse import ArgumentParser

from sage.schemes.elliptic_curves.weierstrass_morphism import *
from sage.schemes.elliptic_curves.hom_composite \
    import EllipticCurveHom_composite as hom_comp


from montgomery_isogenies.isogenies_x_only \
    import isogeny_from_scalar_x_only, evaluate_isogeny_x_only,\
    random_isogeny_x_only
from utilities.discrete_log import weil_pairing_pari, discrete_log_pari
from utilities.supersingular import compute_point_order_D, torsion_basis
from utilities.find_prime import find_prime_with_torsion

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
tt=Integer(args.two_torsion)
Tt=Integer(args.three_torsion)
d = 3**(Tt)
D = 2**(tt)
N = d*D


p, cofactor = find_prime_with_torsion(N)

print(f"Found prime {p}=2**{tt}*3**{Tt}*({factor(cofactor)})-1\n\n")

K.<J> = GF((p,2), modulus=[1,0,1])

E0=EllipticCurve(K, [1,0])
iota=WeierstrassIsomorphism(E0, [-J,0,0,0], E0)

P, Q = torsion_basis(E0, D)


t0 = time()
phi, E = random_isogeny_x_only(E0, d)
t1 = time()

"""
d is the isogeny degree, D is the point orders
"""

t2 = time()
P_, Q_ = evaluate_isogeny_x_only(phi, P, Q, D, d)
t3 = time()

print(f"Isogeny computation:\n {E0.j_invariant()}--->{E.j_invariant()}\
    \n\ntook {t1-t0} seconds.\n\n")

print(f"Isogeny evaluation at:\n\n {P},\n {Q}\n with image\n {P_}, {Q_}\
    \n\ntook {t3-t2} seconds.\n")

R,_ = torsion_basis(E0, d)
t4 = time()
phi =hom_comp(E0, R)
t5 = time()

t6 = time()
P_, Q_ = phi(P), phi(Q)
t7 = time()

print(f"Isogeny computation:\n {E0.j_invariant()}--->{E.j_invariant()}\
    \n\ntook {t5-t4} seconds.\n\n")

print(f"Isogeny evaluation at:\n\n {P},\n {Q}\n with image\n {P_}, {Q_}\
    \n\ntook {t7-t6} seconds.\n")
