from time import time
from argparse import ArgumentParser

from theta_structures.couple_point import CouplePoint
from theta_isogenies.product_isogeny_sqrt\
    import EllipticProductIsogenySqrt


from utilities.supersingular import compute_point_order_D, torsion_basis
from utilities.find_prime import find_prime_with_torsion
from utilities.parameter_generate import calcFields

import elliptic_curve as ec
import utilities.utilities_festa as utils
import utilities.endomorphism as End
import utilities.quaternion as quat 



parser = ArgumentParser()
parser.add_argument("-t", "--two_torsion", default = "128")
parser.add_argument("-T", "--three_torsion", default = "100")

args = parser.parse_args()

def Fp2ToFp2d(x, zeta2, Fp2d):
    return ZZ((x + x**p)/2) + ZZ((x - x**p)/(2*zeta2)) * Fp2d.gen()

def Fp2dToFp2(x, zeta2, Fp2d):
    i = Fp2d.gen()
    x_= ZZ((x + x**p)/2) + ZZ((x - x**p)/(2 * i)) * zeta2 
    return x_

def basis_change_ring(E0, basis, zeta2, Fp2d):
    basis =[E0(Fp2dToFp2(P[0], zeta2, Fp2d),\
               Fp2dToFp2(P[1], zeta2, Fp2d))\
                for P in basis]
    return basis

"""
TODO: use e2-2 instead of e2 for theta isogenies and EllipticProduct
Isogeny instead of EllipticProductIsogenySqrt
"""
def NonSmoothRandomIsog(e, N, basis2, action_matrices, strategy, use_theta=False):
    P, Q = basis2
    E0 = P.curve()
    p = E0.base_ring().characteristic()
    assert N % 2 == 1
    assert (p + 1) % 2**e == 0
    assert ((2**e)*P).is_zero() and ((2**e)*Q).is_zero()
    assert ((2**(e-1))*P).weil_pairing((2**(e-1))*Q, 2) == -1

    alpha = quat.FullRepresentInteger(N*(2**e - N), p)
    vP = End.action_by_matrices(alpha, [1, 0], action_matrices)
    alphaP = vP[0]*P + vP[1]*Q
    vQ = End.action_by_matrices(alpha, [0, 1], action_matrices)
    alphaQ = vQ[0]*P + vQ[1]*Q
 
    assert P.weil_pairing(Q, 2**e)**(N*(2**e - N)) == alphaP.weil_pairing(alphaQ, 2**e)
    if use_theta:
        X, Y, XY = d2isogeny.D2IsogenyImage(E0, E0, (2**e - N)*P, (2**e - N)*Q, alphaP, alphaQ, e, [(P, E0(0)), (Q, E0(0)), (P + Q, E0(0))], strategy, use_theta)
        if not (X[0] + Y[0] == XY[0] or X[0] + Y[0] == -XY[0]):
            Y[0] = -Y[0]
        Pd, Qd = X[0], Y[0]
        if not Pd.weil_pairing(Qd, 2**e) == P.weil_pairing(Q, 2**e)**N:
            if not (X[1] + Y[1] == XY[1] or X[1] + Y[1] == -XY[1]):
                Y[1] = -Y[1]
            Pd, Qd = X[1], Y[1]
        assert Pd.weil_pairing(Qd, 2**e) == P.weil_pairing(Q, 2**e)**N
    else:
        X, Y = d2isogeny.D2IsogenyImage(E0, E0, (2**e - N)*P, (2**e - N)*Q, alphaP, alphaQ, e, [(P, E0(0)), (Q, E0(0))], strategy, use_theta)
        Pd, Qd = X[0], Y[0]
        if not Pd.weil_pairing(Qd, 2**e) == P.weil_pairing(Q, 2**e)**N:
            Pd, Qd = X[1], Y[1]
        assert Pd.weil_pairing(Qd, 2**e) == P.weil_pairing(Q, 2**e)**N

    return Pd, Qd

"""
TODO: use e2-2 instead of e2 for theta isogenies and EllipticProduct
Isogeny instead of EllipticProductIsogenySqrt
"""
def get_matrices(E0, l, e, Fp2d, basis=None):
    p = E0.base_ring().characteristic()

    assert ((p+1)%(l**e) == 0)
    assert E0 == EllipticCurve(Fp2d, [1,0])

    if basis == None:
        basis = torsion_basis(E0, l**e)

    """
    Changing the field of definition to a subfield of Fp4 only for 
    action_matrices.
    """
    Fp4, Fp2, zeta2 = calcFields(p)
    E0=EllipticCurve(Fp2, [1,0])

    basis_bis = basis_change_ring(E0, basis, zeta2, Fp2d)

    Ms = End.action_matrices(basis_bis ,l**e ,zeta2, Fp4)

    return Ms

"""
Pushes an isogeny:
    E0--->EA
of degree q(2**e2-q) to a random 
    Es--->~Es
of the same degree through parallel 3**e3 isogenies.
I'm giving P0, Q0 because I'm not using the canonical 
2**e2-basis.
"""
"""
phi: E->E1 is a q-isogeny, psi:E1 -> E_ a 2^e-q-isogeny forming:
                  phi
                E----->E1
                |      |
                |      |
            psi'|      | psi
                |      |
                E2---->E_
                  phi'
with f = psi*phi 
    -We compute the isogeny with kernel <-q*P, psi*phi(P), P in E[2^e]>
        from ExE_ to E1xE2 of matrix (phi -hatpsi )
                                     (phi' hatpsi')
"""
def Kani_from_kernel(E, E_, P, Q, fP, fQ, q, e):
    kernel = [CouplePoint(q*P, fP), CouplePoint(q*Q, fQ)]
    Phi = EllipticProductIsogenySqrt(kernel, e)
    return Phi

"""Evaluates a non smooth isogeny of degree q in higher dimensional 
representation by Phi, on the l^e torsion"""
def eval_q_2a_q_from_kani(Phi, Pl, Ql, Pl_, Ql_, l, e, as_vector = True):
    O=Pl-Pl
    O_=Pl_-Pl_
    Rs = [[Pl, O_], [Ql, O_], [O, Pl_], [O, Ql_]]
    Rs = [CouplePoint(R[0], R[1]) for R in Rs]
    Rs = [Phi(R) for R in Rs]

    xP, yP = utils.BiDLP(Rs[0][0], Rs[2][0], Rs[3][0], l**e)
    xQ, yQ = utils.BiDLP(Rs[1][0], Rs[2][0], Rs[3][0], l**e)

    if as_vector:
        return [xP, yP], [xQ, yQ]
    else:
        return xP*Pl_+yP*Ql_, xQ*Pl_+yQ*Ql_

def eval_non_smooth_from_kani(Phi, Pl, Ql, Pl_, Ql_, q, e, l=2):
    O=Pl-Pl
    O_=Pl_-Pl_
    Rs = [[Pl, O_], [Ql, O_], [Pl+Ql, O_]]
    Rs = [CouplePoint(R[0], R[1]) for R in Rs]
    X, Y, XY = [Phi(R) for R in Rs]

    if not (X[0] + Y[0] == XY[0] or X[0] + Y[0] == -XY[0]):
        Y[0] = -Y[0]
    Pd, Qd = X[0], Y[0]
    if not Pd.weil_pairing(Qd, l**e) == Pl.weil_pairing(Ql, l**e)**q:
        if not (X[1] + Y[1] == XY[1] or X[1] + Y[1] == -XY[1]):
            Y[1] = -Y[1]
        Pd, Qd = X[1], Y[1]
    assert Pd.weil_pairing(Qd, l**e) == Pl.weil_pairing(Ql, l**e)**q

    return Pd, Qd

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

q = randint(0, 2**(2*e2))
while q*(2**(3*e2)-q)%2 == 0 or q*(2**(3*e2)-q)%3 == 0:
    q = randint(0, 2**(2*e2))


P0, Q0 = torsion_basis(E0, 2**(3*e2))
Ms = get_matrices(E0, 2, 3*e2, Fp2d, basis=[P0, Q0])

theta = quat.FullRepresentInteger(q*(2**(3*e2) - q), p)

print(f"Found endomorphism with coefficients: {theta}")
vthetaP0 = End.action_by_matrices(theta, [1,0], Ms)
vthetaQ0 = End.action_by_matrices(theta, [0,1], Ms)

thetaP0 = vthetaP0[0]*P0+vthetaP0[1]*Q0
thetaQ0 = vthetaQ0[0]*P0+vthetaQ0[1]*Q0

t0 = time()
Phi = Kani_from_kernel(E0, E0, P0, Q0, thetaP0, thetaQ0, q, 3*e2)
t1 = time()
P1, Q1 = eval_non_smooth_from_kani(Phi, P0, Q0, thetaP0,\
                                    thetaQ0, q, 3*e2)

E1 = P1.curve()

print(f"Computed isogeny\n\n\
        {E0.j_invariant()}--->{E1.j_invariant()}\n\n\
        In time {t1-t0}\n")

print(f"Now computing radical isogenies:\n\n")
sleep(1)
t3 = time()
PQ1 = P1+Q1

zeta3 = (-1 + Fp2d(-3).sqrt())/2

E1_, xs = ec.chain_3radials(E1, [P1.xy()[0], Q1.xy()[0], PQ1.xy()[0]],zeta3, 2*e3)
t4 = time()

print(f"Computed new isogeny {E1.j_invariant()}--->{E1_.j_invariant()}\n\
        ")
print(f"In time: {t4-t3}")
