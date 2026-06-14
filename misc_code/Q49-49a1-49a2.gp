elltorsfield_deg( E, n, deg ) = {                                                        
  my( f, eqn, res );                                                            
  f = elldivpol( E, n );                                                        
  eqn = y^2 + E[1]*x*y + E[3]*y - x^3 - E[2]*x^2 - E[4]*x - E[5];               
  res = substpol( polresultant( f, eqn ), y, x );                               
  printf( "Computing %d-torsion field...\n", n );                               
  return( nfinit( nfsplitting( f*res, deg ) ) );                                     
} 

print("Computing 49-torsion field for 49a1, expecting maximal degree 42*49\n");
elltorsfield_deg( ellinit( "49a1" ), 49, 42*49 );

print("Computing 49-torsion field for 49a1, expecting maximal degree 42*49\n");
elltorsfield_deg( ellinit( "49a1" ), 49, 42*49 );
