

# This file was *autogenerated* from the file push_isog.sage
from sage.all_cmdline import *   # import sage library

_sage_const_2 = Integer(2); _sage_const_0 = Integer(0); _sage_const_1 = Integer(1); _sage_const_4 = Integer(4); _sage_const_3 = Integer(3)
from time import time
from argparse import ArgumentParser

from sage.schemes.elliptic_curves.weierstrass_morphism import *

from montgomery_isogenies.isogenies_x_only     import isogeny_from_scalar_x_only, evaluate_isogeny_x_only,    random_isogeny_x_only

from theta_structures.couple_point import CouplePoint
from theta_isogenies.product_isogeny_sqrt    import EllipticProductIsogenySqrt

from utilities.supersingular import compute_point_order_D, torsion_basis
from utilities.find_prime import find_prime_with_torsion
from utilities.parameter_generate import calcFields

import utilities.utilities_festa as utils
import utilities.endomorphism as End
import utilities.quaternion as quat 



parser = ArgumentParser()
parser.add_argument("-t", "--two_torsion", default = "128")
parser.add_argument("-T", "--three_torsion", default = "100")

args = parser.parse_args()

def Fp2ToFp2d(x, zeta2, Fp2d):
    return ZZ((x + x**p)/_sage_const_2 ) + ZZ((x - x**p)/(_sage_const_2 *zeta2)) * Fp2d.gen()

def Fp2dToFp2(x, zeta2, Fp2d):
    i = Fp2d.gen()
    x_= ZZ((x + x**p)/_sage_const_2 ) + ZZ((x - x**p)/(_sage_const_2  * i)) * zeta2 
    return x_

def basis_change_ring(E0, basis, zeta2, Fp2d):
    basis =[E0(Fp2dToFp2(P[_sage_const_0 ], zeta2, Fp2d),               Fp2dToFp2(P[_sage_const_1 ], zeta2, Fp2d))                for P in basis]
    return basis

"""
TODO: use e2-2 instead of e2 for theta isogenies and EllipticProduct
Isogeny instead of EllipticProductIsogenySqrt
"""
def isog_q_2a_q(E0, p, q, Fp2d, Fp2, Fp4, zeta2, verbose = False, verify=False):
    assert p % _sage_const_4  == _sage_const_3 

    """p+1 is really smooth so that will be fast"""
    factors = factor(p+_sage_const_1 )
    e2 = factors[_sage_const_0 ][_sage_const_1 ]
    e3 = factors[_sage_const_1 ][_sage_const_1 ]

    assert ((p+_sage_const_1 )%(_sage_const_2 **e2) == _sage_const_0 ) and ((p+_sage_const_1 )%(_sage_const_3 **e3) == _sage_const_0 )
    assert (q < _sage_const_2 **e2)
    assert (q*(_sage_const_2 **e2-q)%_sage_const_2  != _sage_const_0 ) and (q*(_sage_const_2 **e2-q)%_sage_const_3  != _sage_const_0 )

    """
    torsion_basis() needs E0.base_ring() to have generator a square 
    root of -1. But action_matrices needs the Fp2<Fp4 embedding with
    zeta2**2=-1.
    """
    assert E0 == EllipticCurve(Fp2d, [_sage_const_1 ,_sage_const_0 ])

    basis2 = torsion_basis(E0, _sage_const_2 **e2)
    basis3 = torsion_basis(E0, _sage_const_3 **e3)

    """
    Changing the field of definition to a subfield of Fp4 only for 
    action_matrices.
    """
    E0=EllipticCurve(Fp2, [_sage_const_1 ,_sage_const_0 ])

    basis2_bis = basis_change_ring(E0, basis2, zeta2, Fp2d)
    basis3_bis = basis_change_ring(E0, basis3, zeta2, Fp2d)

    Ms2 = End.action_matrices(basis2_bis ,_sage_const_2 **e2 ,zeta2, Fp4)
    Ms3 = End.action_matrices(basis3_bis ,_sage_const_3 **e3 ,zeta2, Fp4)


    """
    Go back to original field.
    """
    E0=EllipticCurve(Fp2d, [_sage_const_1 ,_sage_const_0 ])

    assert(q*(_sage_const_2 **e2-q)*(_sage_const_3 **e3)>p)
    theta = quat.FullRepresentInteger(q*(_sage_const_2 **e2-q)*(_sage_const_3 **e3), p)
    assert (quat.norm(theta,p)== q*(_sage_const_2 **e2-q)*(_sage_const_3 **e3))

    if verify:
        thetabar = quat.involution(theta)
        assert (quat.norm(thetabar,p)== q*(_sage_const_2 **e2-q)*(_sage_const_3 **e3))


    vector_thetaP =        End.action_by_matrices(theta, [_sage_const_1 , _sage_const_0 ], Ms2)
    vector_thetaQ =        End.action_by_matrices(theta, [_sage_const_0 , _sage_const_1 ], Ms2)
    vector_thetaR =        End.action_by_matrices(theta, [_sage_const_1 , _sage_const_0 ], Ms3)
    vector_thetaS =        End.action_by_matrices(theta, [_sage_const_0 , _sage_const_1 ], Ms3)

    if verify:
        vector_thetabarP =            End.action_by_matrices(thetabar, [_sage_const_1 , _sage_const_0 ], Ms2)
        vector_thetabarQ =            End.action_by_matrices(thetabar, [_sage_const_0 , _sage_const_1 ], Ms2)
        vector_thetabarR =            End.action_by_matrices(thetabar, [_sage_const_1 , _sage_const_0 ], Ms3)
        vector_thetabarS =            End.action_by_matrices(thetabar, [_sage_const_0 , _sage_const_1 ], Ms3)

        Mtheta2 = matrix([vector_thetaP, vector_thetaQ])
        Mtheta3 = matrix([vector_thetaR, vector_thetaS])

        Mthetabar2 = matrix([vector_thetabarP, vector_thetabarQ])
        Mthetabar3 = matrix([vector_thetabarR, vector_thetabarS])

        print(Mtheta3*Mthetabar3)

        assert (Mtheta2*Mthetabar2 == ((q*(_sage_const_2 **e2-q)*(_sage_const_3 **e3))%_sage_const_2 **e2)*matrix.identity(_sage_const_2 ))
        assert (Mtheta3*Mthetabar3 == _sage_const_0 *matrix.identity(_sage_const_2 ))
    """
    Compute image of the 2**e2 and 3**e3 torsion through an 
    endomorphism of degree q(2**e2-q)3**e3
    """

    P0, Q0 = basis2
    R0, S0 = basis3


    thetaP = vector_thetaP[_sage_const_0 ]*P0 + vector_thetaP[_sage_const_1 ]*Q0
    thetaQ = vector_thetaQ[_sage_const_0 ]*P0 + vector_thetaQ[_sage_const_1 ]*Q0

    thetaR = vector_thetaR[_sage_const_0 ]*R0 + vector_thetaR[_sage_const_1 ]*S0
    thetaS = vector_thetaS[_sage_const_0 ]*R0 + vector_thetaS[_sage_const_1 ]*S0

    print(factor(thetaR.order()), factor(thetaS.order()))
    print(factor(thetaP.order()), factor(thetaQ.order()))
    
    assert (thetaP.has_order(_sage_const_2 **e2))
    assert (thetaQ.has_order(_sage_const_2 **e2))
    assert (thetaR.has_order(_sage_const_3 **e3))
    assert (thetaS.has_order(_sage_const_3 **e3))

    """
    Find x,y such that x*theta(R0)+y*theta(S0)=O, so that 
    x*R0+y*S0 is the kernel of a 3**e3 isogeny. Alternatively,
    one of theta(R0) or theta(S0) is in the kernel of the 3**e3 
    isogeny ?
    """
    O=R0-R0
    x, y = utils.BiDLP(O, thetaR, thetaS, _sage_const_3 **e3)
    """Why is x,y always 0 0?"""
    print(x,y)

    if thetaR == O:
        psi, EA = isogeny_from_scalar_x_only(                        E0, _sage_const_3 **e3, vector_thetaS,                        basis=[R0, S0]                    )
    else:
        psi, EA = isogeny_from_scalar_x_only(                        E0, _sage_const_3 **e3, vector_thetaR,                        basis=[R0, S0]                    )
    """
    Need to multiply by inv3e3 for the kernel of the 
    (2**e2, 2**e2)-isogeny.
    """
    inv3e3 = ZZ(inverse_mod(_sage_const_3 **e3, _sage_const_2 **e2))

    """
    Should check sign. TODO
    """
    PA, QA = evaluate_isogeny_x_only(psi,                                       inv3e3*thetaP,                                       inv3e3*thetaQ,                                       _sage_const_2 **e2, _sage_const_3 **e3)
    
    return EA, PA, QA


"""
Pushes an isogeny:
    E0--->EA
of degree q(2**e2-q) to a random 
    Es--->~Es
of the same degree through parallel 3**e3 isogenies.
I'm giving P0, Q0 because I'm not using the canonical 
2**e2-basis.
"""
def push_isogeny(E0, EA, P0, Q0, PA, QA, q, _sage_const_2 , e):
    p = E0.base_ring().characteristic()
    assert p % _sage_const_4  == _sage_const_3 

    """p+1 is really smooth so that will be fast"""
    factors = factor(p+_sage_const_1 )
    e2 = factors[_sage_const_0 ][_sage_const_1 ]

    assert e<=e2

    e3 = factors[_sage_const_1 ][_sage_const_1 ]

    assert ((p+_sage_const_1 )%(_sage_const_2 **e2) == _sage_const_0 ) and ((p+_sage_const_1 )%(_sage_const_3 **e3) == _sage_const_0 )
    assert (q < _sage_const_2 **e2)
    assert (q*(_sage_const_2 **e2-q)%_sage_const_2  != _sage_const_0 ) and (q*(_sage_const_2 **e2-q)%_sage_const_2  != _sage_const_0 )

    m = randint(_sage_const_0 , _sage_const_3 **e3)
    R0, S0 = torsion_basis(E0, _sage_const_3 **e3)

    phi1, E1 = isogeny_from_scalar_x_only(                    E0, _sage_const_3 **e3, m,                    basis=[R0,S0]                )

    """
    Pushing the kernel of phi1 in the diagram through f=f1*f2 of 
    degrees q(2**e2-q)=deg(f):
                    f1
                E0------>E_A,1
                |          |
            f2' |          | f2
                |          |
                E_A,2--->E_A
    Once we have the imgs phi(ker(phi1), O) we get f1(R0+m*S0). We
    now want f2(f1(R0+m*S0)). We can take a basis of the 3**e3 torsion
    of EA, push it to E_A,1 and solve a BiDLP for f1(R0+m*S0).
    """

    kernel = [CouplePoint(q*P0, PA), CouplePoint(q*Q0, QA)]
    phi = EllipticProductIsogenySqrt(kernel, e)

    R3, S3 = torsion_basis(EA, _sage_const_3 **e3)
    Rs = [[R0+m*S0, R3-R3], [R0-R0, R3], [R0-R0, S3]]
    Rs = [CouplePoint(R[_sage_const_0 ], R[_sage_const_1 ]) for R in Rs]
    Rs = [phi(R) for R in Rs]

    x, y = utils.BiDLP(Rs[_sage_const_0 ][_sage_const_0 ], Rs[_sage_const_1 ][_sage_const_0 ], Rs[_sage_const_2 ][_sage_const_0 ], _sage_const_3 **e3)
    """
    Now we have x, y such that x*hatf2(R3)+y*hatf2(S3)=f1(R0+m*S0),
    so that (2**e2-q)^{-1}(x*R3+y*S3)=f(R0+m*S0), in particular, 
    x*R3+y*s3 and f(R0+m*S0) generate the same kernel isogeny!
    """

    #inv2e2 = inverse_mod(2**e2, 3**e3)
    #fK1 = inv2e2*(x*R3+y*S3)

    phi1_, E1_ = isogeny_from_scalar_x_only(                        EA, _sage_const_3 **e3,                        [x,y],                        basis=[R3, S3]                 )

    P1, Q1 = evaluate_isogeny_x_only(phi1,                                    P0,                                    Q0,                                    _sage_const_2 **e2, _sage_const_3 **e3)

    P2, Q2 = evaluate_isogeny_x_only(phi1_,                                    PA,                                    QA,                                    _sage_const_2 **e2, _sage_const_3 **e3)

    return E1, E1_, P1, Q1, P2, Q2

    

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
def eval_non_smooth_from_kani(Phi, Pl, Ql, Pl_, Ql_, l, e):
    O=Pl-Pl
    O_=Pl_-Pl_
    Rs = [[Pl, O_], [Ql, O_], [O, Pl_], [O, Ql_]]
    Rs = [CouplePoint(R[_sage_const_0 ], R[_sage_const_1 ]) for R in Rs]
    Rs = [Phi(R) for R in Rs]

    xP, yP = utils.BiDLP(Rs[_sage_const_0 ][_sage_const_0 ], Rs[_sage_const_2 ][_sage_const_0 ], Rs[_sage_const_3 ][_sage_const_0 ], l**e)
    xQ, yQ = utils.BiDLP(Rs[_sage_const_1 ][_sage_const_0 ], Rs[_sage_const_2 ][_sage_const_0 ], Rs[_sage_const_3 ][_sage_const_0 ], l**e)

    return [xP, yP], [xQ, yQ]

"""
data:
    -N is the rational torsion available
    -d is the isogeny degree
    -D is the torsion we take
"""
e2=Integer(args.two_torsion)
e3=Integer(args.three_torsion)
D = _sage_const_2 **(e2)
d = _sage_const_3 **(e3)
N = d*D

p, cofactor = find_prime_with_torsion(N)

print(f"Found prime {p}=2**{e2}*3**{e3}*({factor(cofactor)})-1\n\n")

factors = factor(p+_sage_const_1 )
e2 = factors[_sage_const_0 ][_sage_const_1 ]
e3 = factors[_sage_const_1 ][_sage_const_1 ]

Fp2d = GF(p**_sage_const_2 , modulus=[_sage_const_1 ,_sage_const_0 ,_sage_const_1 ], name="i")
E0=EllipticCurve(Fp2d, [_sage_const_1 ,_sage_const_0 ])

q = randint(_sage_const_0 , _sage_const_2 **e2)
while q*(_sage_const_2 **e2-q)%_sage_const_2  == _sage_const_0  or q*(_sage_const_2 **e2-q)%_sage_const_3  == _sage_const_0 :
    q = randint(_sage_const_0 , _sage_const_2 **e2)

Fp4, Fp2, zeta2 = calcFields(p)
EA, PA, QA = isog_q_2a_q(E0, p, q, Fp2d, Fp2, Fp4,zeta2,                                          verbose = True)

P0, Q0 = torsion_basis(E0, _sage_const_2 **e2)

t0 = time()
E1, E1_, P1, Q1, P1_, Q1_ = push_isogeny(E0, EA, P0, Q0, PA, QA, q, _sage_const_2 ,                                         e2)

#E2, E2_, P2, Q2, P2_, Q2_ = push_isogeny(E1, E1_, P1, Q1, P1_, Q1_, q, 2, e2)
t1 = time()
print(f"Computed starting isogeny\n\n\
        {E1.j_invariant()}--->{E1_.j_invariant()}\n\n\
        In time {t1-t0}\n")

P3, Q3 = torsion_basis(E1, _sage_const_3 **e3)
P3_, Q3_ = torsion_basis(E1_, _sage_const_3 **e3)

t3 = time()

Phi = Kani_from_kernel(E1, E1_, P1, Q1, P1_, Q1_, q, e2)
v1, v2 = eval_non_smooth_from_kani(Phi, P3, Q3, P3_, Q3_, _sage_const_3 , e3)

t4 = time()

print(f"Computed new isogeny {E1.j_invariant()}--->{E1_.j_invariant()}\n\
        ")
print(f"In time: {t4-t3} with image points {v1,v2}")



