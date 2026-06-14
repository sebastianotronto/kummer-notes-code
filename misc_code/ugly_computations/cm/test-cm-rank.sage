def test():
    attach("cm-rank+-all.sage")
    no = 0
    unk = 0
    for cf in data:
        E = EllipticCurve(cf)
        E_K = E.base_extend(NumberField(x^2-E.cm_discriminant(),'a'))
        try:
            if E_K.rank() != 2* E.rank():
                print "Non doubling:", cf, E.rank(), E_K.rank()
                no += 1
        except:
            print "Can't compute the rank for", cf
            unk += 1
    print "Total ok curves:", (len(data)-(no+unk)), "/", len(data), \
            "( unknown:", unk, ")"
