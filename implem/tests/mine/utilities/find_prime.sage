def find_prime_with_torsion(N, cofactor=1):
    p = N*cofactor-1

    while (not p.is_pseudoprime()) or (p&3 == 1):
        cofactor+=1
        p+=N

    return p, cofactor

if __name__ == "__main__":
    p, cofactor = find_prime_with_torsion(102402393828)

    print(f"Is {p} prime? {p in Primes()}")

    p, cofactor= find_prime_with_torsion(102402393828)

    print(f"Is {p} prime? {p in Primes()}")
