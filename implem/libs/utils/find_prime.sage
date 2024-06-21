def find_prime_with_torsion(N, one_mod_three=False):
    cofactor = 1
    p = N*cofactor-1
    if one_mod_three:
        while not p.is_pseudoprime() and p&3==1:
            cofactor+=1
            p+=N
    else:
        while not p.is_pseudoprime():
            cofactor+=1
            p+=N

    assert(p in Primes()) 
    
    return p, cofactor
