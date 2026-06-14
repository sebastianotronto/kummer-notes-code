read("400k-3.gp");
/* Testing all ECs with a 3-torsion point, conductor < 4*10^5 and rank > 0 */
/* Taken from LMFBD */

S = List();

{
for( i = 1, length( data ), 
  if( i % 100 == 1, print( "Tested ", i-1, " curves." ); );
  coeffs = data[i];
  E = ellinit( coeffs );
  f3 = elldivpol( E, 3 ); /* 3-division polynomial */
  eqn = y^2 + E[1]*x*y + E[3]*y - x^3 - E[2]*x^2 - E[4]*x - E[5]; /* eqn of E */
  res = substpol( polresultant( f3, eqn ), y, x ); /* Keep x for consistency */
  K3 = nfinit( polredbest( nfsplitting( res, 6 ) ) );
  E_ext = ellinit( coeffs, K3 );
  Egens = ellgenerators(E);
  for( j = 1, length( Egens ), 
    P = Egens[j];
    if( ellisdivisible( substpol(E_ext,x,y), P, 3, &Q ),
      cremona_label = ellidentify(E)[1][1];
      listput( S, [cremona_label, P, Q] );
      write( "pari-output-ext", cremona_label );
      write( "pari-output-ext", P );
      write( "pari-output-ext", Q );
      print( "Found EC ", cremona_label, " with point ", P );
      /* Check for 3-divisibility up to torsion */
      T = elltors(E)[3][1];
      if( ellisdivisible( E, elladd(E,P,T), 3 ),
        print( "P+T_1 is 3-divisible!" );
      );
      T = ellmul(E,T,2);    
      if( ellisdivisible( E, elladd(E,P,T), 3 ),
        print( "P+T_2 is 3-divisible!" );
      );  
    );
  );
);
}

