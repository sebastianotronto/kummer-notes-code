/* Candidate counterex mod3 not surj, unexpected density, rank >=1 */
label = [ "220a3", "220a4", "660c3", "660c4", "770f3", "770f4" ];

{
for( i = 1, length( data ), 
  print( "Testing curve ", i, ", label ", label[i] );
  E = ellinit( label[i] );
  f9 = elldivpol( E, 9 ); /* 3-division polynomial */
  eqn = y^2 + E[1]*x*y + E[3]*y - x^3 - E[2]*x^2 - E[4]*x - E[5]; /* eqn of E */
  res = substpol( polresultant( f9, eqn ), y, x ); /* Keep x for consistency */
  print( "Computing K9..." );
  K9 = nfinit( nfsplitting( res, 3^5*2 ) ); /* mo3 is at most 6 */
  print( "Computed!" );
  E_ext = ellinit( E[1..5], K9 );
  Egens = ellgenerators(E);
  for( j = 1, length( Egens ), 
    P = Egens[j];
    if( ellisdivisible( substpol(E_ext,x,y), P, 3, &Q ),
      print( "Found EC ", cremona_label, " with point ", P, 
             " and 3-division point ", Q, " over K9" );
      /* Check for 3-divisibility up to torsion */
      /* T = elltors(E)[3][1];
      if( ellisdivisible( E, elladd(E,P,T), 3 ),
        print( "P+T_1 is 3-divisible!" );
      );
      T = ellmul(E,T,2);    
      if( ellisdivisible( E, elladd(E,P,T), 3 ),
        print( "P+T_2 is 3-divisible!" );
      ); */ 
    );
  );
);
}

