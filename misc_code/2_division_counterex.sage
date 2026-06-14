
R1.<x> = PolynomialRing(QQ)
R2.<x,y> = PolynomialRing(QQ)

def extended_field( f, A, B, deg_mult ):
    # f: a polynomial whose roots are the x-coordinates of some points of an
    #    elliptic curve E: y^2 = x^3 + Ax + B.
    # return value: a field containing the x and y coordinates of those points
    # deg_mult: a positive integer known to be a multiple of the degree of the
    #    extended field.
    #
    # This function uses the properties of resultants (I can provide a pdf
    # explaining how it works).
    # 
    # It is useful to compute, e.g., the fields obtained by adjoining the
    # coordinates of the n-division points of a point (using the n-uplication
    # formulas to get the required polynomials).
    #
    # When used to compute the 2-division fields, it gives the same output as
    # E.division_field(2).

    
    g = y^2 - x^3 - A*x - B
    res = f.resultant(g,x)
    res = res.subs(y=x)

    #print "Splitting field of division pol:"
    #print R1(f).splitting_field('r')
    K.<b> = (R1(res*f)).splitting_field(degree_multiple=deg_mult)

    return K


A = -3
B = 17/4

E = EllipticCurve([0,0,0,A,B])
print E
print "of rank ", E.rank(), "with (some) points of infinite order:", E.gens()
print "CM:", E.has_cm()

rep = E.galois_representation()
print "mod 2 rep is surjective:", rep.is_surjective(2)

K_2.<a> = E.division_field(2)
print "2-division field:", K_2
P = [4,15/2]

# The following polynomial is derived from the duplication formula
# (Silverman, p.54) using the x-coordinate of the 2-division point as
# an indeterminate.x
f_P = R1(x^4 - 4*P[0]*x^3 - 2*A*x^2 - 4*(2*B + A*P[0])*x + A^2 - 4*B*P[0])

M = extended_field( f_P, A, B, 24 ) # The 2-division field of P
print "2-division field of P:", M

if M.degree() != 24:
    print "Stopping because 2-division of P is too small"
    exit()
