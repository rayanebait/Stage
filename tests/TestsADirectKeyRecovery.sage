from sage.rings.finite_rings.integer_mod import square_root_mod_prime
from sage.schemes.elliptic_curves.weierstrass_morphism import *
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite as hom_comp
from sage.schemes.elliptic_curves.isogeny_small_degree import isogenies_2
from argparse import ArgumentParser
import time 
from pprint import pp

proof.arithmetic(False)

parser=ArgumentParser()

parser.add_argument('-n', '--parametersize', default='64')
parser.add_argument('-1', '--prime1', default='3')
parser.add_argument('-2', '--prime2', default='2')
parser.add_argument('-r', '--range', default='20')
args = parser.parse_args()

n=Integer(args.parametersize)

"""Assuming l_A>l_B"""
l_A=Integer(args.prime1)
l_B=Integer(args.prime2)
r=Integer(args.range)



def LBound(x,a=1/2,c=1): return exp(c*( (ln(x)**a) * (ln(ln(x))**(1-a)) ))

def find_smooth_f(l_A, l_B, n, fac_or_test_prime_or_test_prime_power=1):
    A=l_A**n
    ratio=float(log(l_A)/log(l_B))
    B=l_B**(n*(ceil(ratio)))
    print(f"log in base {l_B}:{n*(ceil(ratio))}")
    bound=float(LBound(B))
    print(B>A)

    G=cartesian_product([range(r), range(4)])
    if fac_or_test_prime_or_test_prime_power==0:
        for (i,j) in G:
            "Assume l_B==2"
            f=(B//(l_B**i))-(A//(l_A**j))

            factors=list(factor(f))
            print(factors)
            if max(factors)[0]<bound:
                return (f,(i,j))
    elif fac_or_test_prime_or_test_prime_power==1:
        for (i,j) in G:
            "Assume l_B==2"
            f=(B//(l_B**i))-(A//(l_A**j))
            if f.is_pseudoprime() and f>0:
                return (f,(i,j))
    else:
        for (i,j) in G:
            "Assume l_B==2"
            f=(B//(l_B**i))-(A//(l_A**j))
            if f.is_prime_power() and (not f.is_pseudoprime()):
                return (f,(i,j))
    return None

res=find_smooth_f(l_A, l_B, n, 1)

if res==None:
    print("Search failed")
else:
    f, (i,j)=res
    print(f"Found an f, {f}, with indices i={i}, j={j}")
