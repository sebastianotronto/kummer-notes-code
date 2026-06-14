from sage.schemes.elliptic_curves.ell_generic import is_EllipticCurve

# Pre-tests on E and P
def suitable( E, P, ell ):
    if not is_EllipticCurve(E):
        print "E is not an elliptic curve"
        return False
    if E.base_field() != QQ:
        print "E is not defined over Q"
        return False
    if not P in E:
        print "P is not in E"
        return False
    if P.has_finite_order():
        print "P has finite order"
        return False
    if not is_prime(ell):
        print "ell is not prime"
        return False
    return True


# Stupid auxiliary function. Returns higest power of n dividing m.
def val( n, m ):
    ret = 0
    while m % (n^(ret+1)) == 0:
        ret += 1
    return ret

# Valuation of l-divisibility of the point P
# The parameters h and k describe the ell-part of E(F) (k>=h)
def divisibility( P, ell, k, h ):
    v = 0
    while v != k-1 and P.is_divisible_by(ell^(v+1)):
        v += 1
    return v

# label is the Cremona label of an elliptic curve over Q of rank >= 1.
# ell is a rational prime.
def test_label( label, ell ):
    E = EllipticCurve(label)
    if E.rank() < 1:
        print "E has rank 0"
        return
    P = E.gens()[0]
    test( E, P, ell )

# E is an elliptic curve over Q and P a point of infinite order on E.
# ell is a rational prime.
def test( E, P, ell ):
    if not suitable( E, P, ell ):
        return

    # Setting up for small primes cases
    print "Computations for small primes starting..."
    small_primes = []
    vl = []
    lpart = []
    for p in Primes():
        if p > 100:
            break
        if p == ell:
            print "Skipping", p, "because = ell"
            continue
        if E.discriminant() % p == 0:
            print "Skipping", p, "because E has bad reduction"
            continue
        if P.reduction(p).order() % ell != 0:
            print "Skipping", p, "because red of P is infinitely ell-divisible"
            continue
        small_primes.append(p)
        print "Working with prime p =", p

        # "torsion level"
        F = FiniteField(p)
        vlj = []
        lpj = []
        for i in range(3): ### I WOULD LIKE TO CHANGE THIS TO SOMETHING BIGGER
            print "Torsion level", i, "..."

            # We can build the division fields for increasing powers of l
            # incrementally. To get a division field, we first compute the
            # splitting field of the division polynomial. We may be off by
            # a degree 2 extension. If so, by finite field magic we know
            # exactly which degree 2 extension we need: the unique one!
            R.<x> = PolynomialRing(F)
            F.<a> = E.reduction(p).division_polynomial(ell^i).splitting_field()
            E_red = E.reduction(p).base_extend(F)
            # extend if necessary
            k = val(ell,E_red.gens()[0].order())
            h = val(ell,E_red.order()) - k
            if k < i or h < i:
                F.<a> = F.extension(2)
            E_red = E_red.base_extend(F)
            
            P_red = E_red.point(P.reduction(p))

            print "[Torsion field computed]"

            # l-part of the abelian group E(F_i)
            #gg = E_red.abelian_group().gens()
            #if len( gg ) == 1:
            #    lpj.append( ( val(ell,gg[0].order()), 0 ) )
            #else:
            #    lpj.append( (val(ell,gg[0].order()), val(ell,gg[1].order())) )
            k = val(ell,E_red.gens()[0].order())
            h = val(ell,E_red.order()) - k
           
            lpj.append( ( k, h ) )
            vlj.append(divisibility(P_red,ell,k,h))
            
        vl.append(vlj)
        lpart.append(lpj)

    # Output
    print "Done!"
    for i in range(3):
        print ""
        print "******************************************************"
        print "Torsion Level:", i
        print ""
        rows = [small_primes, [vl[j][i] for j in range(len(small_primes))],
                              [lpart[j][i] for j in range(len(small_primes))] ]
        print table(rows)
        print ""

# Wrapper for densities(E,P,ell)
def densities_label( label, ell ):
    E = EllipticCurve(label)
    if E.rank() < 1:
        print "E has rank 0"
        return
    P = E.gens()[0]
    densities( E, P, ell )

def ratio( h, k, d, ell ):
    x1 = max(k-d,0)
    y1 = max(h-d,0)
    x2 = max(k-d-1,0)
    y2 = max(h-d-1,0)
    up = ell^(x1+y1)-ell^(x2+y2)
    down = ell^(h+k)
    if d == 0:
        return RDF((up+1)/down)
    return RDF(up/down)


# Computes the densities dens[n] of primes such that P is ell^n-divisible mod p
def densities( E, P, ell ):
    if not suitable( E, P, ell ):
        return

    n_primes = 0
    divisible = []
    inf_divisible = 0

    # Array for counting divisibility with specified l-part
    div_part = [[[0 for i in range(10)] for j in range(10)] for k in range(10)]
    # total number of elements in the reductions that (do not) form the l-parts
    tot_l_part = 0
    tot_non_l_part = 0

    for p in Primes():
        if p > 10^2:
            break
        if p == ell or E.discriminant() % p == 0:
            continue

        n_primes +=1
        #print "Working with prime", p

        E_red = E.reduction(p)
        N = E_red.order()
        k = val(ell,E_red.gens()[0].order())
        h = val(ell,N) - k
        temp_non_l_part = (N / (ell^(h+k)))-1
        tot_non_l_part += temp_non_l_part
        temp_l_part = N - temp_non_l_part
        tot_l_part += temp_l_part
        # Debug
        # print "Group structure at", p, ":"
        # print E_red.abelian_group()
        # print "Our result:", N, (k,h), temp_l_part, temp_non_l_part

        if P.reduction(p).order() % ell != 0:
            inf_divisible += 1
            continue
        
        n = divisibility( P.reduction(p), ell, \
                val(ell,E.reduction(p).gens()[0].order()), 0 ) # Wrong parameters but ok
        while len(divisible) < n+1:
            divisible.append(0)
        divisible[n] += 1
        div_part[k][h][n] += 1

    N = len(divisible)
    divtotal = [0]*N
    divtotal[N-1] = divisible[N-1] + inf_divisible
    for i in range(2,N):
        divtotal[N-i] = divisible[N-i] + divtotal[N-i+1]

    print "Tested primes:", n_primes
    for i in range(10):
        for j in range(10):
            s = 0
            for l in range(10):
                s += div_part[i][j][l]
            if s != 0:
                print "l-part (%d,%d):"%(i,j)
                for l in range(10):
                    if div_part[i][j][l] != 0:
                        print "Exactly %d^%d-divisible: %d, expected %0.3f"%\
                                (ell,l,div_part[i][j][l],ratio(i,j,l,ell))
               
    for n in range(1,N):
        print "At least %d^%d-divisble: %0.3f density (%d times)"%(ell,n,\
                RDF(divtotal[n]/n_primes),divtotal[n]) 
    print "Infinitely-divisble: %0.3f density (%d times), expected %0.3f"%\
                (RDF(inf_divisible/n_primes),inf_divisible,\
                RDF(tot_non_l_part/(tot_l_part+tot_non_l_part)))



