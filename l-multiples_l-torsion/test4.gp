read("bunch_of_curves.gp");

{
for( i = 1, length( data ), 
  print( "Testing curve with coefficients ", data[i] );
  E = ellinit( data[i] );
  f4 = elldivpol( E, 4 ); /* 4-division polynomial */
  eqn = y^2 + E[1]*x*y + E[3]*y - x^3 - E[2]*x^2 - E[4]*x - E[5]; /* eqn of E */
  res = substpol( polresultant( f4, eqn ), y, x ); /* Keep x for consistency */
  print( "Computing K4..." );
  K4 = nfinit( nfsplitting( res, 96 ) );
  print( "Computed!" );
  E_ext = ellinit( E[1..5], K4 );
  Egens = ellgenerators(E);
  for( j = 1, length( Egens ), 
    P = Egens[j];
    if( ellisdivisible( substpol(E_ext,x,y), P, 2, &Q ),
      print( "Found EC ", cremona_label, " with point ", P, 
             " and 2-division point ", Q, " over K4" );
      if( ellisdivisible( substpol(E_ext,x,y), P, 4, &Q ),
        print( "Even 4-divisible! 4-division point: ", Q );
      );
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

