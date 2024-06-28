from time import time
from argparse import ArgumentParser
from pprint import pprint

from utilities.find_prime import find_prime_with_torsion
from utilities.quaternion import cornacchia_smith
from utilities.ecm_with_abort import ecm_with_abort

proof.arithmetic(False)

def two_squares_ecm(n, factors):
    if n < 0:
        return None
    if n == 0:
        return 0, 0
    if n == 1:
        return 1, 0
    a, b = 1, 0
    q_=1
    for q in factors:
        if q!=q_:
            e = 0
            q_ = q
        while n % q == 0:
            n //= q
            e += 1
        a *= q**(e//2)
        b *= q**(e//2)
        if e % 2 == 1:
            if q % 4 == 3:
                return None
            else:
                s, t = cornacchia_smith(q)
                at = a*s - b*t
                bt = a*t + b*s
                a, b = at, bt
    print(n)
    return a, b


parser = ArgumentParser()
parser.add_argument("-t", "--two_torsion", default = "128")

b=floor(float(log(2,3)*128))
parser.add_argument("-T", "--three_torsion", default = str(b))
parser.add_argument("-p", "--as_prime", action = "store_true")

args = parser.parse_args()

"""
data:
    -N is the rational torsion available
    -d is the isogeny degree
    -D is the torsion we take
    -as_prime, chooses wether we search fo 2**a-3**bd_A as a prime
        or factorize it with ecm_with_abort
"""
a=Integer(args.two_torsion)
b=Integer(args.three_torsion)

D = 2**(6*a)
d = 3**(2*b)
N = d*D
as_prime = args.as_prime


p, cofactor = find_prime_with_torsion(N)

print(f"Found prime {p}=2**(6*{a})*3**(2*{b})*({factor(cofactor)})-1\n\n")

"""2**a ~ p, 3**b = p**(1/4)"""
b = floor(log(2,3)*a/3)

""""""
r = ceil(log(p, 2))
dA = 3**b-r
n = 2**a - (3**(2*b))*dA
res = 1, 1

if as_prime:
    r = ceil(log(p))
    while dA < 3**b+r:
        print('.', end ='')
        if n.is_pseudoprime():
            print("Found pseudo prime {n}\n")
            if n in Primes() and n%4 == 1:
                print("Found candidate {n}\n")
                break
        n-=3**b
        dA+=1

    assert dA != 3**b+r
    x, y = cornacchia_smith(n)
else:
    fac = []
    while dA < 3**b+r:
        fac = ecm_with_abort(n, 0.5)
        print(fac)
        if fac != None:
            res = two_squares_ecm(n, fac)
            if res != None:
                break
        n-=3**b
        dA+=1

    assert dA != 3**b+r

    x, y = res

print(f"Found:\n\t x={x},\n\ty={y}")
print(f"With sum:\n\t x**2+y**2={x**2+y**2}")
print(f"For n={n}")
assert x**2+y**2 == n

print(f"And remainder/quotient r={r}, dA={dA}")




