# Prints a list of matrices, the elements of GL_2(Z/nZ)

def gcd( m, n ):
    if n == 0:
        return m
    return gcd( n, m%n )

def print_mat( a, b, c, d ):
    print ""
    print "|", a, " ", b, "|"
    print "|", c, " ", d, "|"

n = input("Choose n: ")

count = 0
for a in range(n):
    for b in range(n):
        for c in range(n):
            for d in range(n):
                if gcd( n, abs(a*d-b*c)) == 1:
                    print_mat(a,b,c,d)
                    count += 1

print "Found", count, "elements"
