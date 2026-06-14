E = EllipticCurve([0, 1, 0, -44, -84])
print "Rational torsion points:", E.torsion_points()
print "Generators of the free part:", E.gens()
P = E.gens()[0]

K = CyclotomicField(3)
EK = E.base_extend(K)
PK = EK.point(P)

print "2-division points over Q(zeta_3):", PK.division_points(2)
