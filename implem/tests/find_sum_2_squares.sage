from time import time
from argparse import ArgumentParser

from utilities.find_prime import find_prime_with_torsion

parser = ArgumentParser()
parser.add_argument("-t", "--two_torsion", default = "128")
parser.add_argument("-T", "--three_torsion", default = "1")

args = parser.parse_args()

"""
data:
    -N is the rational torsion available
    -d is the isogeny degree
    -D is the torsion we take
"""
tt=Integer(args.two_torsion)
Tt=Integer(args.three_torsion)

D = 2**(tt)
d = 3**(Tt)
N = d*D


p, cofactor = find_prime_with_torsion(N)

print(f"Found prime {p}=2**{tt}*3**{Tt}*({factor(cofactor)})-1\n\n")

"""2**a ~ p, 3**b = p**(1/3)"""
b = floor(log(2,3)*tt/3)

""""""
r = ceil(float(log(p)))
dA = 3**b-r
while dA < 3**b+r:
    fac = factor(2**a-(3**b)*dA)
    if [p[0]%4 for p in fac] == [1]*len(fac):
        break
    dA+=1

assert dA != 3**b+r

print(dA)
