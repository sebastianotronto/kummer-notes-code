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

def test_jump_kl2_label( label, ell ):
    E = EllipticCurve( label )
    if E.rank() < 1:
        print "E has rank 0"
        return
    P = E.gens()[0]
    test_jump_kl2( E, P, ell )

def test_jump_kl2( E, P, ell ):
    if not suitable( E, P, ell ):
        return

    flag = True
    for p in Primes():
        if p > 10^3:
            break
        if p == ell:                                                            
            # print "Skipping", p, "because = ell"                                
            continue                                                            
        if E.discriminant() % p == 0:                                           
            # print "Skipping", p, "because E has bad reduction"                  
            continue                                                            
        if P.reduction(p).order() % ell != 0:                                   
            # print "Skipping", p, "because red of P is infinitely ell-divisible" 
            continue

        # print "Working with prime p =", p

        F = FiniteField(p)

        for i in range(3):   
            # print "Torsion level", i, "..."                                     
                                                                                 
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

            flag = False
            if P_red.is_divisible_by(ell):
                # print "P ell-divisible in this torsion level"
                flag = True
                break

        if not flag:
            print "Point not divisible mod", p
            print "Stopping here"
            break
    if flag:
        print "*********************************"
        print "*** Candidate counterexample! ***"
        print "*********************************"

# Wrapper                                                
def test_jump_den_label( label, ell ):                                              
    E = EllipticCurve(label)                                                    
    if E.rank() < 1:                                                            
        print "E has rank 0"                                                    
        return                                                                  
    P = E.gens()[0]                                                             
    test_jump_den( E, P, ell ) 

def test_jump_den( E, P, ell ):
    if not suitable( E, P, ell ):
        return

    n_primes = 0
    inf_divisible = 0
    tot_l_part = 0
    tot_non_l_part = 0

    for p in Primes():
        if p > 5*10^3:
            break
        if p == ell or E.discriminant() % p == 0:
            continue

        n_primes += 1

        E_red = E.reduction(p)
        N = E_red.order()                                                       
        k = val(ell,E_red.gens()[0].order())                                    
        h = val(ell,N) - k                                                      
        temp_non_l_part = (N / (ell^(h+k)))-1                                   
        tot_non_l_part += temp_non_l_part                                       
        temp_l_part = N - temp_non_l_part
        tot_l_part += temp_l_part

        if P.reduction(p).order() % ell != 0:
            inf_divisible += 1

    found = RDF( inf_divisible / n_primes )
    expected = RDF( tot_non_l_part / (tot_l_part+tot_non_l_part) )
    print "Found: %0.3f, expected: %0.3f"%(found, expected)
    if abs(found-expected) > 0.05:
        print "Unexpected density! Checking divisibility in reductions..."
        test_jump_kl2( E, P, ell )

def test_from_file( filename ):
    attach(filename)
    for coord in data:
        label = EllipticCurve(coord).label()
        print "Trying curve", label, "with ell =", 3
        test_jump_den_label(label,3)
