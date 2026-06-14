# Determines the proportion of real quadratic fields with odd class number.
def f(B):
    sqf = [ b for b in range(2,B) if b%4==1 and is_squarefree(b) ]
    total_odd = 0

    R.<x> = PolynomialRing(QQ)
    for a in sqf:
        if NumberField( x^2 - a, 'c' ).class_number() % 2 == 1:
            total_odd += 1

    print "=1 mod 4 squarefree integers up to", B, ": ", len(sqf)
    print "With odd class number:", total_odd
