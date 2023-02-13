\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\       Copyright (C) 2014 Denis Simon
\\
\\ Distributed under the terms of the GNU General Public License (GPL)
\\
\\    This code is distributed in the hope that it will be useful,
\\    but WITHOUT ANY WARRANTY; without even the implied warranty of
\\    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
\\    General Public License for more details.
\\
\\ The full text of the GPL is available at:
\\
\\                 http://www.gnu.org/licenses/
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

/*
  Auteur :
  Denis SIMON -> simon@math.unicaen.fr
  adresse du fichier :
  www.math.unicaen.fr/~simon/ellQ.gp

  *********************************************
  *          VERSION 13/01/2014               *
  *********************************************


  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  \\             English                      \\
  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  This package provides functions to compute the rank of elliptic
  curves over Q using 2-descent.
  This package requires the other package qfsolve.gp downloadable at
  www.math.unicaen.fr/~simon/qfsolve.gp
  It also requires the package ellcommon.gp  downloadable at
  www.math.unicaen.fr/~simon/ellcommon.gp

  They can be run under GP by the commands 
  gp > \r qfsolve.gp
  gp > \r ellcommon.gp
  gp > \r ellQ.gp

  The main function is ellQ_ellrank(), which takes as an argument
  any elliptic curve in the form [a1,a2,a3,a4,a6]
  the result is a vector [r,s,v], where
    r is a lower bound for the rank,
    s is the rank of the 2-Selmer group
    v is a set of independant points in E(Q)/2E(Q).

  Example:

  gp > ell = [1,2,3,4,5];
  gp > ellQ_ellrank(ell)
  %1 = [1, 1, [[1,2]]
  In this example, the rank is exactly 1, and [1,2] has infinite order.

  more details on the computations may be obtained by setting
  DEBUGLEVEL_ell = 1 (the higher value, the more details)

  Other functions: 

  ell2descent_complete, ell2descent_gen, ell2descent_viaisog,
  ellhalf, ellredgen, ellsort, elltors2, elltorseven,
  locallysoluble, ratpoint, redquartic,
  bnfpSelmer, reducemodsquares

  Quick information is obtained by typing
  gp > ?NameOfTheFunction

  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  \\    Comment utiliser ce programme ?       \\
  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  Ce module contient des fonctions pour calculer le rang des courbes
  elliptiques sur Q, en utilisant la 2-descente.
  langage : GP

  Ce module utilise les modules complementaires suivants :
    qfsolve.gp
    ellcommon.gp
  qui sont disponibles a l'adresse
    www.math.unicaen.fr/~simon/qfsolve.gp
    www.math.unicaen.fr/~simon/ellcommon.gp

  Pour l'utiliser, lancer gp, puis taper
  \r qfsolve.gp
  \r ellcommon.gp
  \r ellQ.gp

  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  \\  Description des principales fonctions   \\
  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  Explications succintes :
  La fonction ellQ_ellrank() accepte toutes les courbes sous la forme
  [a1,a2,a3,a4,a6]
  Les coefficients peuvent etre entiers ou non.
  L'algorithme utilise est celui de la 2-descente.
  La 2-torsion peut etre quelconque.
  Il suffit de taper : 

  gp > ell = [a1,a2,a3,a4,a6];
  gp > ellQ_ellrank(ell)

  Retourne un vecteur [r,s,v] ou
    r est le rang probable (c'est toujours une minoration du rang),
    s est le 2-rang du groupe de Selmer,
    v est une liste de points independants dans E(Q)/2E(Q).

  Exemple :

  gp > ell = [1,2,3,4,5];
  gp > ellQ_ellrank(ell)
  %1 = [1, 1, [[1,2]]
  Ici, le rang est exactement 1, et le point [1,2] est d'ordre infini.

  Courbes de la forme : k*y^2 = x^3+A*x^2+B*x+C
  sans 2-torsion, A,B,C entiers.
  gp > bnf = bnfinit(x^3+A*x^2+B*x+C);
  gp > ell = ellinit([0,A,0,B,C]);
  gp > rank = ell2descent_gen(ell,bnf,k);

  Courbes avec #E[2](Q) >= 2 :
  ell doit etre sous la forme 
  y^2 = x^3 + A*x^2 + B*x
  avec A et B entiers.
  gp > ell = [0,A,0,B,0]
  gp > ell2descent_viaisog(ell)
  = algorithme de la 2-descente par isogenies
  Attention A et B doivent etre entiers

  Courbes avec #E[2](Q) = 4 : y^2 = (x-e1)*(x-e2)*(x-e3)
  gp > ell2descent_complete(e1,e2,e3)
  = algorithme de la 2-descente complete
  Attention : les ei doivent etre entiers.

  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  \\     Autres fonctions utiles              \\
  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  elltors2(E)      : determine le groupe E[2](Q)
  elltorseven(E)   : determine le groupe E[2^*](Q)
  ellhalf(E,P)     : liste les points Q tels que 2Q = P
  ellredgen(E,v)   : reduction des points de v sur E

  locallysoluble(pol): teste si y^2=pol(x) est ELS
  ratpoint(pol,lim): cherche un point sur y^2=pol(x)
  redquartic(pol): reduction de la quartique pol


  Aide en ligne :
  ?NomDeLaFonction


  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  \\       Affichage des calculs              \\
  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  On peut avoir plus ou moins de details de calculs avec
  DEBUGLEVEL_ell = 0;
  DEBUGLEVEL_ell = 1; 2; 3;...

*/

\\
\\ Usual global variables
\\

global(DEBUGLEVEL_ell, LIM1, LIM3, LIMTRIV, ELLREDGENFLAG, COMPLETE):small;

  DEBUGLEVEL_ell = 0; \\ From 0 to 5: choose a higher value to have
                      \\ more details printed.
  LIM1 = 5;           \\ Limit for the search of trivial points on quartics
  LIM3 = 50;          \\ Limit for the search of points on ELS quartics
  LIMTRIV = 3;        \\ Limit for the search of trivial points on the elliptic curve
  ELLREDGENFLAG = 1;  \\ to reduce the generators at the end
  COMPLETE = 0;       \\ Use Complete 2-descent when full 2-torsion,
                      \\ otherwise 2-descent via isogenies.

\\
\\  Technical global variables
\\

global(MAXPROB, LIMBIGPRIME):small;

  MAXPROB = 20;
  LIMBIGPRIME = 30; \\ for primes larger than this limit: use a probabilistic test
                    \\ LIMBIGPRIME = 0 means: only deterministic tests

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\          SCRIPT                             \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\    MANIPULATION OF GLOBAL VARIABLES         \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{default_ellQ(
  DEBUGLEVEL_ell_val:small = 0,
  LIM1_val:small = 5,
  LIM3_val:small = 50,
  LIMTRIV_val:small = 3,
  ELLREDGENFLAG_val:small = 1,
  COMPLETE_val:small = 0,
  MAXPROB_val:small = 20,
  LIMBIGPRIME_val:small = 30
  ) = 

  DEBUGLEVEL_ell = DEBUGLEVEL_ell_val;
  print("  DEBUGLEVEL_ell = ",DEBUGLEVEL_ell);

  LIM1 = LIM1_val;
  print("  LIM1 = ",LIM1);

  LIM3 = LIM3_val;
  print("  LIM3 = ",LIM3);

  LIMTRIV = LIMTRIV_val;
  print("  LIMTRIV = ",LIMTRIV);

  ELLREDGENFLAG = ELLREDGENFLAG_val;
  print("  ELLREDGENFLAG = ",ELLREDGENFLAG);

  COMPLETE = COMPLETE_val;
  print("  COMPLETE = ",COMPLETE);

  MAXPROB = MAXPROB_val;
  print("  MAXPROB = ",MAXPROB);

  LIMBIGPRIME = LIMBIGPRIME_val;
  print("  LIMBIGPRIME = ",LIMBIGPRIME);
}

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\    FUNCTIONS FOR POLYNOMIALS                \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{ratpoint( pol, lim:small = 1, singlepoint = 1, infinity = 1) =
\\ Search for points on y^2=pol(x).
\\ The coeff of pol must be integers.
\\ If singlepoint >= 1, stop after a first point is found.

my(listpoints,point1,odd,deg4,pol16,tab16,pol9,tab9,pol5,tab5,pol0,vecz,vecx,lead,zz,xx,evpol,iz,factpol,deg,vz,epsmax);

if( DEBUGLEVEL_ell >= 4,
  print("    Starting ratpoint with pol = ",pol); 
  print("    lim = ",lim););

  if( !singlepoint, listpoints = []);
  point1 = [];

\\
\\          trivial solutions
\\

\\ the leading coeff is a square
  if( infinity && issquare(pollead(pol)),
    point1 = [ 1, sqrtrat(pollead(pol)), 0];
if( DEBUGLEVEL_ell >= 3, print("   trivial solution: lead(pol) is a square"));
    if( singlepoint,
if( DEBUGLEVEL_ell >= 4, print("    end of ratpoint"));
      return(point1));
    listpoints = concat(listpoints,[point1]));

\\ the constant coeff is a square
  if( issquare(polcoeff(pol,0)),
    point1 = [ 0, sqrtrat(polcoeff(pol,0)) ];
if( DEBUGLEVEL_ell >= 3, print("   trivial solution: pol(0) is a square"));
    if( singlepoint,
if( DEBUGLEVEL_ell >= 4, print("    end of ratpoint"));
      return(point1));
    listpoints = concat(listpoints,[point1]));

\\ roots of pol ?
  factpol = nfroots(,pol);
  if( #factpol, 
if( DEBUGLEVEL_ell >= 3, print("   trivial solution: roots of pol",factpol));
    if( singlepoint,
if( DEBUGLEVEL_ell >= 4, print("    end of ratpoint"));
      return([factpol[1],0]));
    listpoints = concat(listpoints,vector(#factpol,i,[factpol[i],0]))
  );

\\
\\       Sieve
\\

\\ initialisation of the sieve modulo 16, 9 and 5
\\ used only with even degree when lim is large

  deg = poldegree(pol);
  odd = deg%2;
  deg4 = ( !odd && lim > 20);
  if( deg4,

    pol16 = (Vec(pol)*Mod(1,16))~;
    tab16 = matrix(16,16);
    for(xx = 0, 16-1, 
      for(zz = 0, 16-1,
        tab16[xx+1,zz+1] = !issquare(vector(deg+1,i,xx^(deg+1-i)*zz^(i-1))*pol16)));
    pol9 = (Vec(pol)~)*Mod(1,9);
    tab9 = matrix(9,9);
    for(xx = 0, 9-1, 
      for(zz = 0, 9-1,
        tab9[xx+1,zz+1] = !issquare(vector(deg+1,i,xx^(deg+1-i)*zz^(i-1))*pol9)));
    pol5 = (Vec(pol)~)*Mod(1,5);
    tab5 = matrix(5,5);
    for(xx = 0, 5-1, 
      for(zz = 0, 5-1,
        tab5[xx+1,zz+1] = !issquare(vector(deg+1,i,xx^(deg+1-i)*zz^(i-1))*pol5)));
  );

\\ if the degree is odd, search only for square denominators
  if( odd, 
    vecz = vector(lim,i,i^2);
  ,
\\ if the degree is even, the leading coeff must be
\\ a square modulo zz.
    lead = pollead(pol);
    vecz = vector(lim);
    zz = 0;
    for( i = 1, lim,
      zz++; while( !issquare(Mod(lead,zz)),zz++); vecz[i] = zz
  ));

\\ the constant coeff must be a square modulo xx.
  pol0 = polcoeff(pol,0);
  vecx = vector(lim);
  xx = 0;
  for( i = 1, lim,
    xx++; while( !issquare(Mod(pol0,xx)),xx++); vecx[i] = xx);

if( DEBUGLEVEL_ell >= 4, print("    xmax = ",vecx[lim]));
if( DEBUGLEVEL_ell >= 4, print("    zmax = ",vecz[lim]));

if( DEBUGLEVEL_ell >= 5, print("    vecx = ",vecx));
if( DEBUGLEVEL_ell >= 5, print("    vecz = ",vecz));

  if( deg4,
    vz = vector(lim,i,Pol(
      vector(deg+1,j,polcoeff(pol,deg+1-j)*vecz[i]^(j-1))));
  );

\\ if deg is even and pol has no odd terms,
\\ it is enough to consider xx > 0
  if( !odd && subst(pol,x,-x)==pol, epsmax = 1, epsmax = 2);

\\ loop over x = xx/zz
\\ the loop on [xx,zz] is done diagonally
\\ to start with the smallest values of both xx and zz.
  for( somme = 2, 2*lim,
    for( ix = max(1,somme-lim), min(lim,somme-1),
      xx = vecx[ix]; iz = somme-ix; zz = vecz[iz];
      if( gcd(zz,xx) != 1, next);
      for( eps = 1, epsmax,
\\ when eps = 1, xx > 0; when eps = 2, xx < 0.
        if( deg4,
          if( tab16[xx%16+1,zz%16+1], xx=-xx;next);
          if( tab9[xx%9+1,zz%9+1],    xx=-xx;next);
          if( tab5[xx%5+1,zz%5+1],    xx=-xx;next);
          evpol = subst(vz[iz],'x,xx)
        ,
          evpol = subst(pol,variable(pol),xx/zz)
        );
        if( issquare(evpol),
          point1 = [xx/zz,sqrtrat(subst(pol,variable(pol),xx/zz))];
if( DEBUGLEVEL_ell >= 4, print("    point found by ratpoint = ",point1));
          if( singlepoint, break(3));
          listpoints = concat(listpoints,[point1])
        );
        xx = -xx
  )));

  if( point1 != [], 
if( DEBUGLEVEL_ell >= 3, print("   point found by ratpoint = ",point1));
if( DEBUGLEVEL_ell >= 4, print("    end of ratpoint "));
    if( singlepoint, return(point1), return(listpoints))
  );

return([]);
}
{redquartic(pol) =
\\ reduction of the quartic polynomial.
\\ (also works with other degrees)
my(local_prec,prec0,d,disc2,test,r,normderiv,disc2v,q,M);

if( DEBUGLEVEL_ell >= 4, print("    starting redquartic"));
if( DEBUGLEVEL_ell >= 3, print("   reduction of the quartic ",pol));

\\ choice of the real precision used in the computation
  local_prec = prec0 = default(realprecision);
  d = poldegree(pol);
  disc2 = poldisc(pol)^2;
  test = 0;
  while( test == 0,
if( DEBUGLEVEL_ell >= 4, print("    precision = ",local_prec));
    r = polroots(pol);
    normderiv = vector( d, i, norm(subst(pol',variable(pol),r[i])));
    disc2v = prod( i = 1, d, normderiv[i]) * pollead(pol)^(2*d-4);
    test = abs(disc2v-disc2) < 10^(-local_prec\2);
    if( !test, default(realprecision, local_prec *= 2))
  );

\\ former choice of the quadratic form
\\  q = Vec(sum( i = 1, d, norm(x-r[i])));
\\ Now, uses the quadratic form normalized as in Cremona-Stoll
  q = Vec(sum( i = 1, d, norm('x-r[i]) / normderiv[i]^(1/(d-2))));
  M = qfbreduce([q[1],q[2]/2;q[2]/2,q[3]]);
  pol = subst(pol,variable(pol),Pol(M[1,])/Pol(M[2,]))*Pol(M[2,])^poldegree(pol);

  if( local_prec != prec0, default(realprecision,prec0));

if( DEBUGLEVEL_ell >= 3, print("   reduced quartic = ",pol));
if( DEBUGLEVEL_ell >= 4, print("    end of redquartic"));

  return([pol,M]);
}
{listratpoint( pol, redflag = 0) =
my(list,i,K,ff,C,p,M,U,newpol,factpol,ll,listf,rr);

if( DEBUGLEVEL_ell >= 5, print("     Starting listratpoint with ",pol));
  list = [[pol,matid(2),1,1]];
  i = 1; 
  while( i <= #list,

    pol = list[i][1];

    K = abs(content(pol));
    if( K != 1,
      pol = (list[i][1] /= K);
      list[i][3] *= K
    );

    K = list[i][3];
    if( K == 1, i++; next);

    ff = factor(K);
    if( vecmax(ff[,2]) > 1,
      ff[,2] \= 2;
      C = factorback(ff);
      list[i][4] *= C;
      K = ( list[i][3] /= C^2);
      if( K == 1, i++; next);
      ff = factor(K);
    );

    p = ff[1,1];
    M = list[i][2];
    C = list[i][4];

    if( pollead(pol)%p == 0,
      U = M*[1,0;0,p];
      if( content(U) == 1,
        newpol = subst(pol,'x,'x/p)*p^(poldegree(pol)-1);
        list = concat(list, [[newpol,U,K/p,C*p]])
      )
    );

    factpol = centerlift(polrootsmod(pol,p));
    for( j = 1, #factpol,
      U = M*[p,factpol[j];0,1];
      if( content(U) == 1,
        newpol = subst(pol,'x,p*'x+factpol[j])/p;
        list = concat(list, [[newpol,U,K/p,C*p]])
    ));

    i++;
  );

  ll = sum( i = 1, #list, list[i][3] == 1);
  listf = vector(ll);
  i = 1;
  for( j = 1, #list,
    if( list[j][3] == 1,
      listf[i] = list[j]; i++));

  if( redflag,
    for( i = 1, #listf,
      rr = redquartic(listf[i][1]);
      listf[i][1] = rr[1];
      listf[i][2] = listf[i][2]*rr[2]
    )
  );

if( DEBUGLEVEL_ell >= 5, print("     Output of listratpoint = ",listf));

return(listf);
}
{ratpoint2( pol, lim:small = 1, singlepoint = 1, redflag = 0) =
my(listpoints,list,rr,y2,aux);

  listpoints = [];
  list = listratpoint(pol,redflag);
  for( i = 1, #list,
    rr = ratpoint(list[i][1],lim,singlepoint);
    if( singlepoint && #rr, rr=[rr]);
    for( j = 1, #rr,
      y2 = rr[j][2]*list[i][4];
      if( #rr[j] == 2,
        aux = [rr[j][1],1]~
      , aux = [rr[j][1],rr[j][3]]~
      );
      aux = list[i][2] * aux;
      if( aux[2] == 0,
        rr[j] = [aux[1],y2,0]
      , rr[j] = [aux[1]/aux[2],y2/aux[2]^(poldegree(pol)\2)]
      );
    );
    if( singlepoint && #rr, return(rr[1]));
    listpoints = concat(listpoints,rr);
  );
  listpoints = vecsort(listpoints,,2);
return(listpoints);
}
{polrealrootsisolate(pol) =
\\ pol is a squarefree polynomial in Z[x].
\\ Returns a list of vectors [a,b] with a and b rationals
\\ such that the intervals ]a,b] are disjoint and contain
\\ all the real roots of pol, and exactly one in each interval.
my(st,a,res,ind,b,c,stab,stac);

if( DEBUGLEVEL_ell >= 5, print("     starting polrealrootsisolate with pol = ",pol));
  st = polsturm(pol);
  if( !st, return([]));
  a = 1;
  while( polsturm(pol,-a,a) < st, a <<= 1);
  res = [[-a,a,st]];
  ind = 1;
  while( #res < st,
    while( res[ind][3] == 1, ind++);
    a = res[ind][1]; b = res[ind][2]; stab = res[ind][3];
    c = (a+b)/2;
    stac = polsturm(pol,a,c);
    if( stac == 0, res[ind][1] = c; next);
    if( stac == stab, res[ind][2] = c; next);
    res[ind] = [a,c,stac];
    res = concat(res,[[c,b,stab-stac]]);
  );
  res = vector(st,i,[res[i][1],res[i][2]]);
  res = vecsort(res,1);
if( DEBUGLEVEL_ell >= 5, print("     end of polrealrootsisolate with res = ",res));
  return(res);
}
{polrealrootsimprove(pol,v) =
\\ pol is a polynomial and v is a vector v=[a,b]
\\ such that pol contains exactly one root in the interval ]a,b].
\\ Returns another interval with the same property, but with half length.
\\ (dichotomy)
my(c,v2,vc);

  c = (v[1]+v[2])/2;
  v2 = subst(pol,variable(pol),v[2]);
  if( v2 == 0, return([c,v[2]]));
  vc = subst(pol,variable(pol),c);
  if( sign(v2)*sign(vc) >= 0, v[2] = c, v[1] = c);
  return(v);
}
{polrootsmodpn(pol,p) =
\\ Compute a list v. Each element of v is of the form
\\ [a,b], with maximal b <= valuation(poldisc(pol),p)
\\ and a is a root of pol modulo p^b.
my(vd,rac,i,pol2,r,newrac);

if( DEBUGLEVEL_ell >= 5, print("     starting polrootsmodpn ",p,":",pol));

  vd = valuation(poldisc(pol),p);
  rac = [[0,0]];
  i = 1;
  while (i <= #rac,
\\    if( rac[i][2] > vd, i++; next);
    if( rac[i][2] >= vd, i++; next);
    pol2 = subst(pol,'x,rac[i][1]+'x*p^rac[i][2]);
    pol2 /= content(pol2);
    r = lift(polrootsmod(pol2,p));
    if( #r == 0, i++; next);
    newrac = vector(#r,j,[rac[i][1] + p^rac[i][2]*r[j],rac[i][2]+1]);
    rac = concat(rac, vector(#r-1,j,newrac[j+1]));
    rac[i] = newrac[1];
  );
if( DEBUGLEVEL_ell >= 5, print("     end of polrootsmodpn ",rac));
  return(rac);
}

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\    FUNCTIONS FOR LOCAL COMPUTATIONS         \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{ppinit( nf, p) =
\\ a little more structure than idealprimedec()
my(pdec,pp);

  pdec = idealprimedec(nf,p);
  pp = vector(#pdec,i,
    [ pdec[i]
    , nfbasistoalg(nf,pdec[i][2])
    , if( p == 2, idealstar(nf,idealpow(nf,pdec[i],1+2*pdec[i].e)),0)
    , nfmodprinit(nf,pdec[i])
    ]);
  return(pp);
}
{nfpsquareodd( nf, a, pr) =
\\ pr is a prime ideal of nf as output by nfmodprinit
\\ a is an element of nf.
\\ Returns 1 if a is a square in the p-adics, 0 otherwise
\\ works only for p odd.
my(p,v,ap,den,norme);

if( DEBUGLEVEL_ell >= 5, print("     starting nfpsquareodd(",a,",",pr));
  if( a == 0,
if( DEBUGLEVEL_ell >= 5, print("     end of nfpsquareodd"));
    return(1));

  p = pr[3];
  v = idealval(nf,lift(a),p);
  if( v%2,
if( DEBUGLEVEL_ell >= 5, print("     end of nfpsquareodd"));
    return(0));
  ap = nfalgtobasis(nf,a/nfbasistoalg(nf,p[2])^v);
  den = valuation(denominator(content(ap)),p.p);
  if( den,
    den += den%2;
    ap = p.p^den*nfeltmul(nf,ap,nfeltpow(nf,p[2],-den*p.e))
  );

  norme = (p.p^p.f-1)/2;
  ap = nfeltpowmodpr(nf,ap,norme,pr);
  ap[1] -= 1;
  if( ap == 0,
if( DEBUGLEVEL_ell >= 5, print("     end of nfpsquareodd"));
    return(1));
  if( idealval(nf,ap,p) > 0,
if( DEBUGLEVEL_ell >= 5, print("     end of nfpsquareodd"));
    return(1));
if( DEBUGLEVEL_ell >= 5, print("     end of nfpsquareodd"));
  return(0);
}
{psquare( a, p) =
\\ a is an integer.
\\ p is a prime integer.
\\ Returns 1 if a is a square in the p-adics, 0 otherwise.
my(v,ap);

if( DEBUGLEVEL_ell >= 5, print("     starting psquare ",[a,p]));

  if( a == 0, 
if( DEBUGLEVEL_ell >= 5, print("     end of psquare 1"));
    return(1));

  v = valuation(a,p);
  if( v%2,
if( DEBUGLEVEL_ell >= 5, print("     end of psquare 0"));
    return(0));
  if( p == 2, 
    ap = (a>>v)%8-1,
    ap = kronecker(a/p^v,p)-1
  );
if( DEBUGLEVEL_ell >= 5, print("     end of psquare ", !ap));
  return(!ap);
}
{lemma6(pol, p, nu, xx) =
\\ technical lemma for local solubility of quartics
\\ Only for p <> 2.
my(gx,gpx,lambda,mu);

  gx = subst( pol, variable(pol), xx);
  if( psquare(gx,p), return(1));
  gpx = subst( pol', variable(pol), xx);
  lambda = valuation(gx,p); mu = valuation(gpx,p);

  if( gpx != 0 && lambda > 2*mu, return(1));
\\  if( (lambda >= mu+nu) && (nu > mu), return(1));
  if( (lambda >= 2*nu) && (mu >= nu), return(0));
  return(-1);
}
{lemma7( pol, nu, xx) =
\\ technical lemma for local solubility of quartics
\\ at p = 2.
my(gx,gpx,lambda,mu,q);

  gx = subst( pol, variable(pol), xx);
  if( psquare(gx,2), return(1));
  gpx = subst( pol', variable(pol), xx);
  lambda = valuation(gx,2); mu = valuation(gpx,2);
  if( gpx != 0 && lambda > 2*mu, return(1));
  if( nu > mu,
    if( lambda%2, return(-1));
    q = mu+nu-lambda;
    if( q == 1, return(1));
    if( q == 2 && (gx>>lambda)%4 == 1, return(1));
    return(-1));
  q = lambda-2*nu;
  if( q >= 0, return(0));
  if( q == -2 && (gx>>lambda)%4 == 1, return(0));
  return(-1);
}
{zp_soluble(pol, p, nu, pnu, x0, pnup) =
my(result,pol2,fact,x1);

if( DEBUGLEVEL_ell >= 5, print("     starting zp_soluble ",[pol,p,x0,nu]));
  if( p == 2,
    result = lemma7(pol,nu,x0),
    result = lemma6(pol,p,nu,x0));
  if( result == +1,
if( DEBUGLEVEL_ell >= 5, print("     end of zp_soluble 1 lemma"));
    return(1));
  if( result == -1,
if( DEBUGLEVEL_ell >= 5, print("     end of zp_soluble 0 lemma"));
    return(0));
  pnup = pnu*p;
  nu++;
  if( p < LIMBIGPRIME || !LIMBIGPRIME,
    for( i = 0, p-1,
      if( zp_soluble(pol,p,nu,pnup,x0+pnu*i),
if( DEBUGLEVEL_ell >= 5, print("     end of zp_soluble"));
        return(1)))
  ,
    pol2 = subst(pol,variable(pol),x0+pnu*variable(pol));
    pol2 /= content(pol2);
    pol2 = pol2*Mod(1,p);
    if( !poldegree(pol2), return(0));
    fact = factormod(pol2,p)[,1];
    for( i = 1, #fact,
      x1 = -centerlift(polcoeff(fact[i],0));
      if( zp_soluble(pol,p,nu,pnup,x0+pnu*x1),
if( DEBUGLEVEL_ell >= 5, print("     end of zp_soluble"));
        return(1)));
    for( i = 1, MAXPROB,
      x1 = random(p);
      if( zp_soluble(pol,p,nu,pnup,x0+pnu*x1),
if( DEBUGLEVEL_ell >= 5, print("     end of zp_soluble"));
        return(1)))
  );
if( DEBUGLEVEL_ell >= 2,
  if( p >= LIMBIGPRIME, 
    print("  ******* probabilistic test at p = ",p,"*******")));
if( DEBUGLEVEL_ell >= 5, print("     end of zp_soluble"));
  return(0);
}
{qp_soluble(pol, p) = 
if( DEBUGLEVEL_ell >= 5, 
  print("     starting qp_soluble ",p);
  print("     pol = ",pol));
  if( psquare(pollead(pol),p),
if( DEBUGLEVEL_ell >= 5, print("     end of qp_soluble 1"));
    return(1));
  if( psquare(polcoeff(pol,0),p),
if( DEBUGLEVEL_ell >= 5, print("     end of qp_soluble 1"));
    return(1));
  if( zp_soluble(pol,p,0,1,0),
if( DEBUGLEVEL_ell >= 5, print("     end of qp_soluble 1"));
    return(1));
  if( zp_soluble(polrecip(pol),p,1,p,0),
if( DEBUGLEVEL_ell >= 5, print("     end of qp_soluble 1"));
    return(1));
if( DEBUGLEVEL_ell >= 5, print("     end of qp_soluble 0"));
  return(0);
}
{locallysoluble(pol) =
\\ Determines if y^2 = pol(x,z) is everywhere locally soluble
my(c,disc0,plist,p,vc);

if( DEBUGLEVEL_ell >= 4, print("    starting locallysoluble: ",pol));

\\ real place
  if( !(poldegree(pol)%2) && sign(pollead(pol)) < 0 
         && sign(polcoeff(pol,0)) < 0 && polsturm(pol) == 0,
if( DEBUGLEVEL_ell >= 3, print("   not ELS at infinity"));
if( DEBUGLEVEL_ell >= 4, print("    end of locallysoluble"));
    return(0));

\\
\\ finite places
\\
  pol *= denominator(content(pol))^2;
  c = content(pol);

  disc0 = poldisc(pol);
  plist = factor (abs(2*disc0));
if( DEBUGLEVEL_ell >= 4, print("    list of bad primes = ",plist));
  for( i = 1, #plist[,1],
    p = plist[i,1];
if( DEBUGLEVEL_ell >= 4, print("    p = ",p));
    vc = valuation(c,p);
    if( vc >= 2, 
      pol /= p^(2*(vc\2));
      plist[i,2] -= 2*(vc\2)*(2*poldegree(pol)-2)
    );
    if( poldegree(pol) == 4 && p != 2 && plist[i,2] < 2, next);
    if( !qp_soluble(pol,p),
if( DEBUGLEVEL_ell >= 3, print("   not ELS at ",p));
if( DEBUGLEVEL_ell >= 4, print("    end of locallysoluble"));
      return(0)));

if( DEBUGLEVEL_ell >= 2, print("  quartic ELS: Y^2 = ",pol));
if( DEBUGLEVEL_ell >= 4, print("    end of locallysoluble"));
  return(1);
}
{LS2localimage(nf,gen,pp) =
my(p,LS2image,ph,ival,delta);

if( DEBUGLEVEL_ell >= 4, print("    starting LS2localimage"));

  p = pp[1][1].p;
  LS2image = matrix( if( p == 2, sum(i=1,#pp,1+#pp[i][3].cyc), 2*#pp), #gen);

  for( j = 1, #gen,
    ph = [];
    for( i = 1, #pp,
      ival = idealval(nf,gen[j],pp[i][1]);
      ph = concat(ph,[ival]); 
      delta = gen[j]/pp[i][2]^ival;
      if( p == 2,
        ph = concat(ph,ideallog(nf,delta,pp[i][3])~);
      , ph = concat(ph,[1-nfpsquareodd(nf,delta,pp[i][4])]);
      )
    );
    LS2image[,j] = ph~
  );
  LS2image *= Mod(1,2);

if( DEBUGLEVEL_ell >= 4, print("    LS2image = ",lift(LS2image)));
if( DEBUGLEVEL_ell >= 4, print("    end of LS2localimage"));
  return(LS2image);
}

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\    GENERIC FUNCTIONS FOR ELLIPTIC CURVES    \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{ellhalf(ell,P)=
\\ returns all the points Q on ell such that 2Q = P.
my(pol2,ratroots,half,x2,y2,P2);

  if(#ell < 13, ell=ellinit(ell));

  pol2 = Pol([4,ell.b2,2*ell.b4,ell.b6]); \\ 2-division polynomial

  if( P == [0], 
    ratroots = nfroots(,pol2);
    half = vector(#ratroots,i,[ratroots[i],-(ell.a1*ratroots[i]+ell.a3)/2]);
    half = concat( [[0]], half);
    return(half)
  );

  x2 = Pol([1,0,-ell.b4,-2*ell.b6,-ell.b8]); \\ x(2P) = x2/pol2 

  half = [];
  ratroots = nfroots(,x2-P[1]*pol2);  
  if( #ratroots == 0, return(half));
  for( i = 1, #ratroots,
    y2 = ellordinate(ell,ratroots[i]);
    for( j = 1, #y2,
      P2 = [ratroots[i],y2[j]];
      if( ellpow(ell,P2,2) == P, half = concat(half,[P2]))
    )
  );

  return(half);
}
{elltors2(ell)=
\\ Compute the 2-torsion subgroup of the elliptic curve ell.
my(tors2);

if( DEBUGLEVEL_ell >= 3, print("   computing the 2-torsion"));

  tors2 = ellhalf(ell,[0]);
  if( #tors2 == 1, 
    tors2 = [1, [], []],
  if( #tors2 == 2, 
    tors2 = [2, [2], [tors2[2]]]
  , tors2 = [4, [2,2], [tors2[2],tors2[3]]]
  ));
if( DEBUGLEVEL_ell >= 3, print("   E[2] = ",tors2));
  return(tors2);
}
{elltorseven(ell)=
\\ Compute the 2-Sylow subgroup of the torsion of the elliptic curve ell.
my(torseven,P2);

if( DEBUGLEVEL_ell >= 4, print("    computing the 2^n-torsion"));
  if(#ell < 13, ell=ellinit(ell));
  torseven = elltors2(ell);

  while( torseven[1] != 1,
    P2 = ellhalf(ell,torseven[3][1]);
    if( #P2 > 0,
       torseven[1] *= 2;
       torseven[2][1] *= 2;
       torseven[3][1] = P2[1];
       next
    );
    if( #torseven[3] == 1, break());

    P2 = ellhalf(ell,torseven[3][2]);
    if( #P2 > 0,
       torseven[1] *= 2;
       torseven[2][2] *= 2;
       torseven[3][2] = P2[1];
       next
    );
    P2 = ellhalf(ell,elladd(ell,torseven[3][1],torseven[3][2]));
    if( #P2 > 0,
       torseven[1] *= 2;
       torseven[2][1] *= 2;
       torseven[3][1] = P2[1];
       next
    );
    break()
  );
  
if( DEBUGLEVEL_ell >= 4, print("    E[2^n] = ",torseven)); 
  return(torseven);
}
{ellsort(listpts) =
\\ Sorting the points listpts on an elliptic curve
\\ using the naive height.
my(n,v,aux,ord);

  v = vector(n = #listpts);
  for( i = 1, n, 
    if( listpts[i] == [0], v[i] = [0,0,0]; next);
    aux = denominator(listpts[i][2])/denominator(listpts[i][1]);
    v[i] = vecsort(abs([listpts[i][1]*aux^2, listpts[i][2]*aux^3,aux]),,4);
  );
  ord = Vec(vecsort(v,,3)); \\ ord = vecsort(v,,3);
  return(vector(n,i,listpts[ord[i]]));
}
{ellremovetorsion(ell,listgen) =
\\ Extracting the points of infinite order from listgen
my(d,extra);

if( DEBUGLEVEL_ell >= 5, print("     removing torsion from ",listgen));
  d = #listgen;
  extra = 0;
  for( i = 1, d,
\\ points of order 1 or 2
    if( listgen[i] == [0]
     || listgen[i] == ellpow(ell,listgen[i],-1)
    , extra += 1<<(i-1);
      next
    );
\\ detection of infinite order points by looking at
\\ 8*9*5*7*P modulo the prime 1048583
    if( ell.disc%1048583 != 0
     && denominator(listgen[i])%1048583 != 0
     && ellpow(ell,listgen[i]*Mod(1,1048583),2520) != [0]
    , next
    );
\\ detection of torsion points by ellorder()
    if( ellorder(ell,listgen[i]),
      extra += 1<<(i-1)
    )
  );
  if( extra,
    listgen = vecextract(listgen,1<<#listgen-1-extra);
  );
if( DEBUGLEVEL_ell >= 5, print("     without torsion = ",listgen));
  return(listgen);
}
{ellredgen(ell0,listgen,K=1) =
\\ reduction of the generators of points in listgen
\\ on the elliptic curve ell = [a1,a2,a3,a4,a6]
\\ or K*y^2 = x^3 + a2*x^2 + a4*x + a6 (when a1 = a3 = 0);
\\ using the canonical height.
my(d,ell=ell0,sqrtK,urst,extra,M,U,listgen2,tors2,vt);

if( DEBUGLEVEL_ell >= 3, print("   Reduction of the generators ",listgen));
if( DEBUGLEVEL_ell >= 5, print("     ell=",ell));          
  d = #listgen;
  if( d == 0, return([]));

\\ removing torsion points from listgen
  listgen = ellremovetorsion(ell0,listgen);
  d = #listgen;
  if( d == 0, return([]));  

  if( #ell < 13, ell = ellinit(ell));

  if( K != 1,
    if( ell.a1 != 0 || ell.a3 != 0, error(" ellredgen: a1*a3 != 0"));
    ell[2] *= K; ell[4] *= K^2; ell[5] *= K^3;
    ell[6] *= K; ell[7] *= K^2; ell[8] *= K^3; ell[9] *= K^4;
    ell[10] *= K^2; ell[11] *= K^3; ell[12] *= K^6;
    sqrtK = sqrt(K);
    if( #ell == 19,
      ell[14] *= K;
      ell[15] /= sqrtK; ell[16] /= sqrtK;
      ell[17] *= sqrtK; ell[18] *= sqrtK;
      ell[19] /= K
    );

    for( i = 1, d,
      for( j = 1, #listgen[i],
        listgen[i][j] *= K^j))
  );

if( d == 1,
  urst = [1,0,0,0];
,
  if( #ell < 19, ell = ellinit(ell));
  ell = ellminimalmodel(ell,&urst);
  listgen = ellchangepoint(listgen,urst);
if( DEBUGLEVEL_ell >= 5, print("     ell = ",ell));
if( DEBUGLEVEL_ell >= 5, print("     listgen = ",listgen));

\\ Looking for relations between the points in listgen
\\ using LLL on the height matrix

  extra = 1;
  while( extra,
    M = ellheightmatrix(ell,listgen);
if( DEBUGLEVEL_ell >= 4, print("    height matrix = ",M));
    if( abs(matdet(M)) > 10^(-default(realprecision)+10),break);
    U = qflll(round(M*10^(default(realprecision)-10)),4);
    U = concat(U[1],U[2]);
if( DEBUGLEVEL_ell >= 4, print("    change of basis proposed by LLL = ",U));
\\ the columns of U that have very small coefficients
\\ are either exact relations or reductions (coeff <= 20)
\\ the other ones are irrelevant.
    extra = 0;
    for( i = 1, d,
      if( vecmax(abs(U[,i])) > 20, extra += 1<<(i-1))
    );
    U = vecextract(U,1<<d-1-extra);
    U = completebasis(U);
if( DEBUGLEVEL_ell >= 4, print("    change of basis 1 = ",U));

    listgen2 = vector(d);
    for( i = 1, d,
      listgen2[i] = [0];
      for( j = 1, d,
        listgen2[i] = elladd(ell,listgen2[i],ellpow(ell,listgen[j],U[j,i]))));
    listgen = listgen2;
  );

\\ Extracting the points of infinite order

\\ removing torsion points from listgen
    listgen = ellremovetorsion(ell,listgen);
    d = #listgen;
    if( d == 0, return([]));  
  );

if( DEBUGLEVEL_ell >= 3, print("   infinite order points = ",listgen));

\\ Now, the points should be of infinite order and independant
\\ Reducing the points of infinite order

  if( d > 1,
    M = ellheightmatrix(ell,listgen);
if( DEBUGLEVEL_ell >= 4, print("    height matrix = ",M));
    U = qflllgram(M);
if( DEBUGLEVEL_ell >= 4, print("    change of basis 2 = ",U));

    listgen2 = vector(d);
    for( i = 1, d,
      listgen2[i] = [0];
      for( j = 1, d,
        listgen2[i] = elladd(ell,listgen2[i],ellpow(ell,listgen[j],U[j,i]))));
    listgen = listgen2;
  );

if( DEBUGLEVEL_ell >= 3, print("   infinite order points = ",listgen));

  listgen = ellchangepoint(listgen,ellinverturst(urst));

\\ Reducing modulo the 2-torsion

  tors2 = elltorseven(ell0);
  if( tors2[1] > 1,
    vt = vector(tors2[2][1],j,ellpow(ell0,tors2[3][1],j-1));
    if( #tors2[2] == 2,
      vt = concat(vt,vector(#vt,j,elladd(ell0,vt[j],tors2[3][2])))
    );
    for( i = 1, d,
      listgen[i] = ellsort(vector(#vt,j,elladd(ell0,listgen[i],vt[j])))[1];
    );
  );

if( DEBUGLEVEL_ell >= 3, print("   infinite order points = ",listgen));
  
  if( K != 1,
    for( i = 1, d,
      for( j = 1, 2,
        listgen[i][j] /= K^j)));

\\ keep only the points (x,y) with y >= 0

  if( ell0.a1 == 0 && ell0.a3 == 0, 
    for( i = 1, d,
      if( #listgen[i] == 2,
        listgen[i][2] = abs(listgen[i][2]))));

if( DEBUGLEVEL_ell >= 2, print("  reduced generators = ",listgen));
  return(listgen);
}

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\    FUNCTIONS FOR NUMBER FIELDS              \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{reducemodsquares(delta,d) =
\\ Uses LLL to find z such that delta*z^2 has a small coefficient in x^d.
\\ delta must be a t_POLMOD
my(deg,xx,z,qd,Qd,reduc);

  deg = poldegree(component(delta,1));  \\ deg = poldegree(delta.mod);
  xx = Mod('x,component(delta,1)); \\ xx = Mod(x,delta.mod);
  z = subst(Pol(vector(deg,i,eval(Str("a"i)))),'x,xx);
  qd = polcoeff(lift(delta*z^2),d,'x);
  Qd = simplify(matrix(deg,deg,i,j,deriv(deriv(qd,eval(Str("a"i))),eval(Str("a"j)))/2));

  reduc = qflllgram_indef(Qd);
  if( #reduc == 2, reduc = reduc[2][,1]);

  return(delta*subst(Pol(reduc),'x,xx)^2);
}
{bnfpSelmer(bnf,S=1,p) =
\\ p is a prime integer and bnf a big number field.
\\ Compute the p-Selmer group of the number field bnf
\\ relative to the prime ideals dividing S.
\\ This group is denoted by LS2 in the sequel.
\\ Returns [gen,S'] where gen is a vector containing the generators
\\ of the p-Selmer group, represented has elements of bnf modulo p-powers,
\\ and S' is the support of gen.
my(S1,oddclass,multS,Slist,LS2gen,newprimes,newprimesval,kerval);

if( DEBUGLEVEL_ell >= 3, print("   Constructing the field Selmer group: L(S,",p,")"));
  S1 = idealhnf(bnf,S);

  oddclass = 0; multS = 1;
  while( !oddclass,
    if( multS != 1, S1 = idealmul(bnf,S1,multS));
    Slist = idealfactor(bnf,S1)[,1]~;
if( DEBUGLEVEL_ell >= 4, print("    constructing the S-units "));
if( DEBUGLEVEL_ell >= 4, print("    S1 = ",Slist));
    LS2gen = bnfsunit(bnf,Slist);

\\ If the class group is divisible by p,
\\ need to enlarge S1.
    oddclass = LS2gen[5].no % p;
    if( !oddclass,
if( DEBUGLEVEL_ell >= 4, print("    class group divisible by p = ",LS2gen[5].no));
      multS = idealmul(bnf,S,LS2gen[5].gen[1]);
    )
  );
  LS2gen = Mod(LS2gen[1],bnf.pol);

\\ The valuation of the generators must be divisible by p outside S.
  newprimes = [];
  for( i = 1, #Slist,
    if( idealadd(bnf,S,Slist[i]) == 1,
      newprimes = concat(newprimes,[Slist[i]])));
if( DEBUGLEVEL_ell >= 4, print("    newprimes = ",newprimes));
  newprimesval = matrix(#newprimes,#LS2gen,i,j,
    idealval(bnf,LS2gen[j],newprimes[i]));
if( DEBUGLEVEL_ell >= 4, print("    newprimesval = ",newprimesval));
  kerval = lift(matker(newprimesval*Mod(1,p)));
if( DEBUGLEVEL_ell >= 4, print("    kerval = ",kerval));
  LS2gen = vector(#kerval,i,
    prod( j = 1, #LS2gen,
      LS2gen[j]^kerval[j,i]));

\\ Add the units
  LS2gen = concat(bnf.fu,LS2gen); \\ LS2gen = concat(bnf.fu,LS2gen);
\\ Add also the torsion unit if its order is divisible by p.
  if( bnf[8][4][1]%p == 0, \\ if( bnf.tu[1]%p == 0,
    LS2gen = concat( [Mod(bnf[8][4][2],bnf.pol)], LS2gen)); \\ LS2gen = concat( [Mod(bnf.tu[2],bnf.pol)], LS2gen));
if( DEBUGLEVEL_ell >= 3, print("   #LS2gen = ",#LS2gen));
if( DEBUGLEVEL_ell >= 4, print("    LS2gen = ",LS2gen));
  return([LS2gen,Slist]);
}
{kersign(gen,rootapprox) =
\\ Determine the kernel of the sign map
\\ restricted to the subgroup generated by gen,
\\ and relative to the real embedding corresponding to
\\ the root of pol contained in the interval rootapprox.
my(signs,elt,elt2,d,st,kers,compt);

if( DEBUGLEVEL_ell >= 3, print("   Computing the kernel of the sign ",rootapprox));

\\ determination of the signs

  signs = vector(#gen);
  for( i = 1, #gen,
    elt = lift(gen[i]);
if( DEBUGLEVEL_ell >= 5, print("     Computing the sign of elt = ",elt));
    if( poldegree(elt) == 0, signs[i] = sign(simplify(elt)) < 0; next);
    d = poldisc(elt);
    if( poldegree(elt) == 2,
      if( d <= 0, signs[i] = sign(pollead(elt)) < 0; next));
    elt2 = if( d == 0, elt/gcd(elt,elt'), elt);
    st = 1;
    compt = 0;
    while( st,
      st = polsturm(elt2,rootapprox[1],rootapprox[2]);
      if( st,
        rootapprox = polrealrootsimprove(component(gen[i],1),rootapprox);
\\        rootapprox = polrealrootsimprove(gen[i].mod,rootapprox);
\\ if the sign of elt is too difficult to determine, 
\\ try a reduction modulo squares
        if( compt++ == 5, gen[i] = reducemodsquares(gen[i]); i--; next(2));
\\ if the sign of elt is still too difficult to determine, 
\\ try the sign of 1/elt.
        if( compt%5 == 0, gen[i] = 1/gen[i]; i--; next(2))
    ));
    signs[i] = sign(subst(elt,variable(elt),rootapprox[2])) < 0
  );
if( DEBUGLEVEL_ell >= 4, print("    signs = ",signs));

\\ construction of the kernel
  kers = matker(Mat(signs*Mod(1,2)))*Mod(1,2);
if( DEBUGLEVEL_ell >= 4, print("    kers = ",lift(kers)));
  return(kers);
}
{kernorm(gen,S,p) =
\\ gen is a generating set for a subgroup of K^* / K^p.
\\ Compute the kernel of the norm map.
\\ Uses the fact that all valuations are 0 mod p,
\\ except maybe at primes in S.
my(normgen,normmap,kern);
  
if( DEBUGLEVEL_ell >= 3, print("   Computing the kernel of the norm map"));

  if( p == 2, S = concat([-1],S));
  normgen = norm(gen);
if( DEBUGLEVEL_ell >= 4, print("    normgen = ",normgen));

\\ matrix of the norm map
  normmap = matrix(#S,#normgen,i,j,
    if( i == 1 && p == 2,
      sign(normgen[j]) < 0
    , valuation(normgen[j],S[i])));
if( DEBUGLEVEL_ell >= 4, print("    normmap = ",normmap));

\\ construction of the kernel
  kern = matker(normmap*Mod(1,p))*Mod(1,p);
if( DEBUGLEVEL_ell >= 4, print("    ker = ",lift(kern)));
  return(kern);
}

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\    FUNCTIONS FOR 2-DESCENT                  \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{elllocalimage( nf, pp, K = 1) =
\\ pol is the cubic polynomial defining the elliptic curve,
\\ nf is nfinit(pol),
\\ p is a prime integer, and pp = ppinit(p).
\\ Returns the image of the p-adic points 
\\ E(Qp)/2E(Qp) in Kp/Kp^2.
\\ The algorithm consists of choosing random p-adic points in E(Qp)
\\ until the number of images is equal to #E(Qp)[2] / |2|_p
my(X,p,prank,rac,pts,bound,essai,mrank,r,xx,delta,ph,delta2,local_prec,ival);

if( DEBUGLEVEL_ell >= 4, print("    starting elllocalimage"));

  X = Mod('x,nf.pol);
  p = pp[1][1][1];
  prank  = #pp - (p != 2);
if( DEBUGLEVEL_ell >= 4, print("    prank = ",prank));

  rac = polrootsmodpn(K*nf.pol,p);
if( DEBUGLEVEL_ell >= 5, print("     rac = ",rac));

  pts = matrix(0,0);
  bound = p+6;
  essai = 0;
  mrank = 0;
  while( mrank < prank,

    essai ++;
    if( essai%16 == 0,
      pts = matimage(pts);
      bound *= p;
    );

    r = random(#rac)+1; local_prec = random(rac[r][2]+3)-2;
    xx = rac[r][1]+p^local_prec*random(bound);
if( DEBUGLEVEL_ell >= 5, print("     xx = ",xx));
    delta = K*(xx-X);

\\ rem: K*pol(xx) = norm(delta) ( = y^2 for a point on the elliptic curve)
    if( !psquare(K*subst(nf.pol,'x,xx),p), next);
    ph = [];
    for( i = 1, #pp,
      ph = concat(ph,[ ival = idealval(nf,delta,pp[i][1])]);
      delta2 = delta/pp[i][2]^ival;
      if( p == 2,
        ph = concat(ph,ideallog(nf,delta2,pp[i][3])~);
      , ph = concat(ph,[1-nfpsquareodd(nf,delta2,pp[i][4])])
      )
    );
if( DEBUGLEVEL_ell >= 5, print("     ph = ",ph));

    pts = concat(pts,ph~*Mod(1,2));
    mrank = matrank(pts*Mod(1,2));
if( DEBUGLEVEL_ell >= 5, print("     pts = ",lift(pts)));
if( DEBUGLEVEL_ell >= 5, print("    matrank = ",mrank));
  );

  pts = matimage(pts);
if( DEBUGLEVEL_ell >= 5, print("     essai = ",essai));
if( DEBUGLEVEL_ell >= 4, print("    end of elllocalimage"));
  return(pts);
}
{ell2descent_gen(ell,bnf,K=1,help=[],redflag=0) =
\\ This algorithm performs 2-descent on the elliptic curve ell
\\ when ell has trivial 2-torsion.

\\ ell must be of the form K*y^2=x^3+A*x^2+B*x+C
\\ ie ell=[0,A,0,B,C], with K,A,B,C integers.
\\ bnf is bnfinit(x^3+A*x^2+B*x+C,1).

\\
\\ help is a list of known points (maybe empty) on ell.
\\ if redflag != 0, reduces the elements 
\\ of the field Selmer group modulo squares.

my(A,B,C,polrel,polprime,ttheta,badprimes,S,LS2,selmer,rootapprox,p,pp,locimage,LS2image,listpointstriv,listpoints,iwhile,expo,zc,liftzc,den,point,idealfactorzc,idealzc,baseidealzc,q2,sol,param,q1,pol,redq,q0,pointxx,point2,rang);

if( DEBUGLEVEL_ell >= 4, print("    starting ell2descent_gen"));

  if( #ell < 13, ell = ellinit(ell));

  if( ell.a1 != 0 || ell.a3 != 0,
    error(" ell2descent_gen: the curve is not of the form [0,a,0,b,c]"));
  if( denominator(ell.a2) > 1 || denominator(ell.a4) > 1 || denominator(ell.a6) >1,
    error(" ell2descent_gen: non integral coefficients"));

  A = ell.a2; if( DEBUGLEVEL_ell >= 2, print("  A = ",A));
  B = ell.a4; if( DEBUGLEVEL_ell >= 2, print("  B = ",B));
  C = ell.a6; if( DEBUGLEVEL_ell >= 2, print("  C = ",C));
  
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\      Construction of L(S,2)      \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

if( DEBUGLEVEL_ell >= 2, print(); print("  Computing L(S,2)"));

  polrel = Pol([1,A,B,C]);
  polprime = polrel';
  ttheta = Mod('x,polrel);

  if( !bnf, 
if( DEBUGLEVEL_ell >= 3, print("   bnfinit(",polrel,")"));
    bnf = bnfinit(polrel,1));

  badprimes = abs(K*idealadd(bnf,polprime,bnf.index));
if( DEBUGLEVEL_ell >= 5, print("     badprimes = ",badprimes[1,1]));
  S = bnfpSelmer(bnf,badprimes,2);
  LS2 = S[1]; S = S[2];

if( DEBUGLEVEL_ell >= 2, print("  L(S,2) = ",LS2));

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\   Construction of the Selmer group \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

if( DEBUGLEVEL_ell >= 2, print(); print("  Computing the Selmer group"));

\\ elements with square norm
  selmer = kernorm(LS2,vector(#S,i,S[i].p),2);
if( DEBUGLEVEL_ell >= 3, print("   selmer = ",lift(selmer)));

\\ the first real embedding must be > 0
\\ since the norm is a square, this is automatic
\\ if there is a single real embedding.
  if( bnf.r1 == 3,
    rootapprox = polrealrootsisolate(polrel)[1];
    selmer = matintersect(selmer,kersign(LS2,rootapprox))*Mod(1,2);
if( DEBUGLEVEL_ell >= 3, print("   selmer = ",lift(selmer)));
  );

\\ p-adic points
if( DEBUGLEVEL_ell >= 3, print("   p-adic points"));
  badprimes = factorint(badprimes[1,1]*2)[,1];
if( DEBUGLEVEL_ell >= 2, print("  badprimes = ",badprimes));
  for( i = 1, #badprimes,
    p = badprimes[i];
if( DEBUGLEVEL_ell >= 4, print("    p = ",p));
    pp = ppinit(bnf.nf,p);
    locimage = elllocalimage(bnf.nf,pp,K);
    LS2image = LS2localimage(bnf.nf,LS2,pp);
    locimage = matintersect(LS2image,locimage);
    selmer = matintersect(
      selmer,
      concat(
        matker(LS2image),
        matinverseimage(LS2image,locimage)*Mod(1,2)));
    selmer = matimage(selmer*Mod(1,2));
if( DEBUGLEVEL_ell >= 4, print("    selmer = ",lift(selmer)));
  if( !#selmer, break);
  );

if( DEBUGLEVEL_ell >= 2, print("  selmer = ",lift(selmer)));
if( DEBUGLEVEL_ell >= 2, print("  Selmer rank = ",#selmer));

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\   Search for trivial points      \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  if( #selmer,
if( DEBUGLEVEL_ell >= 2, print(); print("  Search for trivial points on the curve"));
    listpointstriv = ratpoint(K^3*subst(polrel,'x,'x/K),LIMTRIV,0);
    for( i = 1, #listpointstriv, 
      if( #listpointstriv[i] == 3,
        listpointstriv[i] = [0]
      , for( j = 1, 2, listpointstriv[i][j] /= K^j))
     );
    listpointstriv = concat(help,listpointstriv);
if( DEBUGLEVEL_ell >= 2, print("  Trivial points on the curve = ",listpointstriv));
  );

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\   Run through the Selmer group   \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

if( DEBUGLEVEL_ell >= 2, print(); print("  Run through the Selmer group"));

  listpoints = [];
  selmer = lift(selmer);
  iwhile = 1;
  while( iwhile < 1<<#selmer,
if( DEBUGLEVEL_ell >= 2, print());
if( DEBUGLEVEL_ell >= 4, print("   iwhile = ",iwhile));

\\ the next element zc as an algebraic number modulo squares

    expo = selmer*vectorv(#selmer,i,bittest(iwhile,i-1));
    zc = prod( i = 1, #LS2, LS2[i]^expo[i]);
if( DEBUGLEVEL_ell >= 2, print("  zc = ",zc));
    liftzc = lift(zc);

\\ Reduction modulo squares

    if( redflag, 
      zc = reducemodsquares(zc,2);
      liftzc = lift(zc);
      den = denominator(content(liftzc))^2;
      zc *= den; liftzc *= den;
if( DEBUGLEVEL_ell >= 2, print("  zc reduced = ",zc))
    );

\\ Does it come from a trivial point ?

    for( i = 1, #listpointstriv,
      point = listpointstriv[i];
      if( #point == 2,
        if( nfissquare(bnf.nf,K*(point[1]-'x)*liftzc),
if( DEBUGLEVEL_ell >= 2, print("  comes from the trivial point ",point));
          listpoints = concat(listpoints,[point]);
          iwhile = 1 << (degre(iwhile)+1);
          next(2)
    )));
    
if( DEBUGLEVEL_ell >= 2, print("  does not come from a trivial point"));

\\ Construction of the quadratic form q2
\\ Change the basis using the square factors of zc

    idealfactorzc = idealfactor(bnf,zc);
    idealfactorzc[,2] *= -1;
    idealfactorzc[,2] \= 2;
\\    idealzc = idealfactorback(bnf,idealfactorzc);
    idealzc = matid(3);
    for( i = 1, #idealfactorzc[,1],
      idealzc = idealmul(bnf,idealzc,idealpow(bnf,idealfactorzc[i,1],idealfactorzc[i,2]));
    );
    baseidealzc = vector(3,i,nfbasistoalg(bnf,idealzc[,i]));
    q2 = matrix(3,3,i,j,trace(zc*baseidealzc[i]*baseidealzc[j]/polprime));
if( DEBUGLEVEL_ell >= 4, print("    q2 = ",q2));
if( DEBUGLEVEL_ell >= 4, print("    q2/content(q2) = ",q2/content(q2)));

\\ Solution of the quadratic equation q2=0

    sol = qfsolve(q2/content(q2));
if( DEBUGLEVEL_ell >= 4,print("    sol = ",sol));
    if( type(sol) == "t_INT", 
      error(" ell2descent_gen: WRONG ELEMENT IN THE SELMER GROUP, please report"));


\\ Parametrizing the solutions of q2=0

    param = qfparam(q2,sol,3)*['x^2,'x,1]~;
    param /= content(param);

\\ Construction of the quartic

    q1 = -matrix(3,3,i,j,trace(zc*baseidealzc[i]*baseidealzc[j]*(ttheta+A)/polprime));
    pol = param~*q1*param;
if( DEBUGLEVEL_ell >= 2, print("  quartic: ",K,"*Y^2 = ",pol));
    redq = redquartic(pol);
    pol = redq[1];
    den = denominator(content(K*pol));
    pol *= den^2;
if( DEBUGLEVEL_ell >= 2, print("  reduced: ",K,"*Y^2 = ",pol));

\\ Search for points on the quartic

    point = ratpoint(K*pol,LIM1,1);
    if( point == [], point = ratpoint2(K*pol,LIM3,1));
    if( point == [], iwhile ++; next );

    if( #point == 2, point = concat(point,[1]));
if( DEBUGLEVEL_ell >= 2, print("  point on the reduced quartic = ",point));
    point = concat(redq[2]*[point[1],point[3]]~,[point[2]/den]~);
if( DEBUGLEVEL_ell >= 2, print("  point on the quartic = ",point));

\\ Construction of the point on the elliptic curve from the point on the quartic

    param = subst(param,'x,'x/'y)*'y^2;
    param = subst(subst(param,'x,point[1]),'y,point[2]);
    param *= K/point[3];
if( DEBUGLEVEL_ell >= 3, print("   reconstruction of the point on the curve"));
    q0 = matrix(3,3,i,j,trace(zc*baseidealzc[i]*baseidealzc[j]*(ttheta^2+A*ttheta+B)/polprime));
    pointxx = param~*q0*param/K;
    point2 = [ pointxx, sqrtrat(subst(polrel,'x,pointxx)/K)];
if( DEBUGLEVEL_ell >= 1, print(" point on the curve = ",point2));
    listpoints = concat(listpoints,[point2]);
    iwhile = 1 << (degre(iwhile)+1)
  );

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\      Conclusion report           \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  rang = #listpoints;

if( DEBUGLEVEL_ell >= 2,
  print();
  print("  rank of found points     = ",#listpoints);
  print("  rank of the Selmer group = ",#selmer));

if( DEBUGLEVEL_ell >= 1, afficheselmer(rang,#selmer));
  
  if(  (#selmer - rang)%2,
    rang ++;
if( DEBUGLEVEL_ell >= 1, 
      print(" III should be a square, hence ");
      afficheselmer(rang,#selmer));
  );

\\ Verification

if( DEBUGLEVEL_ell >= 1, print("listpoints = ",listpoints));
  for( i = 1, #listpoints, 
    if( subst(polrel,'x,listpoints[i][1])-K*listpoints[i][2]^2,
      error(" ell2descent_gen: WRONG POINT = ",listpoints[i]," please report")));

\\ Reduction of the points

  listpoints = vecsort(listpoints,,2);
  if( #listpoints >= 2 && ELLREDGENFLAG,
    listpoints = ellredgen(ell,listpoints,K));

if( DEBUGLEVEL_ell >= 4, print("    end of ell2descent_gen"));
  return([rang,#selmer,listpoints]);
}
{afficheselmer(m1,m2,tors2) =

    print("#E(Q)[2]      = ",1<<tors2);
    print("#S(E/Q)[2]    = ",1<<m2);
  if( m1+tors2 == m2,
    print("#E(Q)/2E(Q)   = ",1<<(m1+tors2));
    print("#III(E/Q)[2]  = 1");
    print("rank(E/Q)     = ",m1);
  ,
    print("#E(Q)/2E(Q)  >= ",1<<(m1+tors2));
    print("#III(E/Q)[2] <= ",1<<(m2-m1-tors2));
    print("rank(E/Q)    >= ",m1)
  );
}
{ellQ_ellrank(ell,help=[]) =
\\ Algorithm of 2-descent on the elliptic curve ell.
\\ help is a list of known points on ell.
my(urst,urst1,den,eqell,tors2,bnf,rang,time1);

if( DEBUGLEVEL_ell >= 3, print("   starting ellQ_ellrank"));
  if( #ell < 13, ell = ellinit(ell));

\\ kill the coefficients a1 and a3
  urst = [1,0,0,0];
  if( ell.a1 != 0 || ell.a3 != 0,
    urst1 = [1,0,-ell.a1/2,-ell.a3/2];
    ell = ellchangecurve(ell,urst1);
    urst = ellcomposeurst(urst,urst1)
  );

\\ kill denominators
  while( (den = denominator([ell.a2,ell.a4,ell.a6])) > 1,
    den = factor(den); den[,2] = vectorv(#den[,2],i,1);
    den = factorback(den);
    urst1 = [1/den,0,0,0];
    ell = ellchangecurve(ell,urst1);
    urst = ellcomposeurst(urst,urst1)
  );

  help = ellchangepoint(help,urst);
  eqell = Pol([1,ell.a2,ell.a4,ell.a6]);
if( DEBUGLEVEL_ell >= 1, print(" Elliptic curve: Y^2 = ",eqell));

\\ Choice of the algorithm depending on the 2-torsion structure

  tors2 = ellhalf(ell,[0]);
if( DEBUGLEVEL_ell >= 1, print(" E[2] = ",tors2));

  if( #tors2 == 1,                              \\ case 1: 2-torsion trivial
if( DEBUGLEVEL_ell >= 3, print1("   bnfinit(",eqell,")"));
if( DEBUGLEVEL_ell >= 4, gettime());
    bnf = bnfinit(eqell,1);
if( DEBUGLEVEL_ell >= 4, time1 = gettime());
if( DEBUGLEVEL_ell >= 3, print(" done"));
    rang = ell2descent_gen(ell,bnf,1,help);
if( DEBUGLEVEL_ell >= 4, print("    time for bnfinit  = ",time1));
if( DEBUGLEVEL_ell >= 4, print("    time for the rest = ",gettime()));
  ,
  if( #tors2 == 2 || !COMPLETE,                 \\ case 2: 2-torsion >= Z/2Z
    if( ell.a6 != 0,
      urst1 = [1,tors2[2][1],0,0];
      ell = ellchangecurve(ell,urst1);
      urst = ellcomposeurst(urst,urst1)
    );
    eqell = Pol([1,ell.a2,ell.a4,ell.a6]);
if( DEBUGLEVEL_ell >= 1, print(" Elliptic curve: Y^2 = ",eqell));

    rang = ell2descent_viaisog(ell,help)
  ,                                             \\ case 3: 2-torsion = Z/2Z*Z/2Z
    rang = ell2descent_complete(tors2[2][1],tors2[3][1],tors2[4][1])
  ));

  rang[3] = ellchangepoint(rang[3],ellinverturst(urst));
if( DEBUGLEVEL_ell >= 3, print("   end of ellQ_ellrank"));

  return(rang);
}
{ell2descent_complete(e1,e2,e3,help) =
\\ Compute the rank of the elliptic curve
\\ E: Y^2 = (x-e1)*(x-e2)*(x-e3)
\\ using the complete 2-descent algorithm (see J.Silverman).
\\ Returns [r,s,v] with
\\   r is a lower bound for the rank of E(Q)
\\   s is the rank of Selmer[2]
\\   v is a system of independant points on E(Q)/2E(Q)

\\ e1, e2, e3 must be integers.
\\ help is a list of known points on E.

my(ee,d32,d31,d21,G1,G2,G3,vect1,vect2,vect3,selmer,rang,listpoints,b1,b2,q1,sol1,param1,param1x,quart,point,z1,solx,soly,strange,ell);

if( DEBUGLEVEL_ell >= 2, print("  Algorithm of complete 2-descent"));

\\ sort the integers e1, e2, e3 in increasing order

  ee = vecsort([e1,e2,e3]);
  e1 = ee[1]; e2 = ee[2]; e3 = ee[3];

\\ Computation of the groups G1 and G2

  d32 = e3-e2; d31 = e3-e1; d21 = e2-e1;
  G1 = factor(d31*d21)[,1];   \\ (G1 > 0)
  G2 = factor(-d32*d21)[,1];  \\ (G2 < 0)
  G3 = d31*d32;

if( DEBUGLEVEL_ell >= 3, print("   G1 = ",G1));
if( DEBUGLEVEL_ell >= 3, print("   G2 = ",G2));

\\ Run through G1*G2

  vect1 = vector(#G1,i,[0,1]);
  vect2 = vector(#G2,i,[0,1]);
  selmer = 0;
  rang = 0;
  listpoints = [];

  forvec( X = vect1,
    b1 = prod( i = 1, #G1, G1[i]^X[i]);

\\ b1*b2*b3 must be a square, where b3 is a divisor of d32*d31
    vect3 = vect2;
    for( i = 2, #G2, 
      if( G3%G2[i] !=0,
        vect3[i] = [1,1]*valuation(b1,G2[i])));

    forvec( Y = vect3,
      b2 = prod( i = 1, #G2, G2[i]^Y[i]);

if( DEBUGLEVEL_ell >= 3, print("   [b1,b2] = ",lift([b1,b2])));

\\ Trivial points coming from the 2-torsion

      if( b1==1 && b2==1,
if( DEBUGLEVEL_ell >= 4, print("    trivial point [0]"));
        selmer++; rang++; next);
      if( issquare(-d21*b2) && issquare(d31*d21*b1),
if( DEBUGLEVEL_ell >= 3, print("   trivial point [e1,0]"));
        selmer++; rang++; listpoints = concat(listpoints,[[e1,0]]); next);
      if( issquare(d21*b1) && issquare(-d32*d21*b2),
if( DEBUGLEVEL_ell >= 3, print("   trivial point [e2,0]"));
        selmer++; rang++; listpoints = concat(listpoints,[[e2,0]]); next);
      if( issquare(d31*b1) && issquare(d32*b2),
if( DEBUGLEVEL_ell >= 3, print("   trivial point [e3,0]"));
        selmer++; rang++; listpoints = concat(listpoints,[[e3,0]]); next);

\\ Trivial points coming from help

      for( i = 1, #help,
        if( #help[i] != 2 || help[i][2] == 0, next);
        if( issquare(b1*(help[i][1]-e1)) && issquare(b2*(help[i][1]-e2)),
if( DEBUGLEVEL_ell >= 3, print("   trivial point from help ",help[i]));
          selmer++; rang++; 
          listpoints = concat(listpoints,[help[i]]); next(2));
      );

\\ If one can solve 2 quadratic equations
\\ (1) q1: b1*z1^2-b2*z2^2 = e2-e1
\\ (2) q2: b1*z1^2-b1*b2*z3^2 = e3-e1
\\ then (x,y) = (b1*z1^2+e1,b1*b2*z1*z2*z3) is a point on E
\\ we also have
\\ (3) q3 = q1-q2: b1*b2*z3^2-b2*z2^2=e2-e3

\\ Solution of the q1

      q1 = matdiagonal([b1,-b2,-d21]);
if( DEBUGLEVEL_ell >= 3, print("   q1 = ",q1));
      sol1 = qfsolve(q1);
      if( type(sol1) == "t_INT",
if( DEBUGLEVEL_ell >= 3, print("   q1 not ELS at ",sol1));
        next);
if( DEBUGLEVEL_ell >= 3, print("   solution of q1 = ",sol1));
      param1 = qfparam(q1,sol1,1);
      param1 /= content(param1);
if( DEBUGLEVEL_ell >= 3, print("   parametrization of q1 = ",param1));
      param1x = param1*['x^2,'x,1]~;

\\ Solution of the q2
\\ only useful to detect local non solubility

\\ my(q2,sol2);
\\      q2 = matdiagonal([b1,-b1*b2,-d31]);
\\if( DEBUGLEVEL_ell >= 3, print("   q2 = ",q2));
\\      sol2 = qfsolve(q2);
\\      if( type(sol2) == "t_INT",
\\if( DEBUGLEVEL_ell >= 3, print("   q2 not ELS at ",sol2));
\\        next);

\\ Construction of the quartic

      quart = b1*b2/gcd(b1,b2)^2*(b1*param1x[1]^2-d31*param1x[3]^2);
if( DEBUGLEVEL_ell >= 3, print("   quart = ",quart));

\\ Search for trivial points on the quartic

      point = [];
\\      point = ratpoint(quart,LIM1,1);

\\ Local solubility of the quartic
   
      if( point == [] && !locallysoluble(quart),
if( DEBUGLEVEL_ell >= 3, print("   quartic not ELS "));
        next);
if( DEBUGLEVEL_ell >= 2, print("  y^2 = ",quart));
      selmer++;

\\ Search for points on the quartic

      if( point == [], point = ratpoint2(quart,LIM3,1));
      if( point != [], 
if( DEBUGLEVEL_ell >= 2, print("  point found on the quartic !!"));
if( DEBUGLEVEL_ell >= 3, print("   ",point));
        if( #point == 2,
          z1 = subst(param1x[1],'x,point[1])/subst(param1x[3],'x,point[1])
        , z1 = param1[1,1]/param1[3,1]);
        solx = b1*z1^2+e1;
        soly = sqrtrat((solx-e1)*(solx-e2)*(solx-e3));
        listpoints = concat(listpoints,[[solx,soly]]);
if( DEBUGLEVEL_ell >= 1, print(" point on the elliptic curve = ",[solx,soly]));
        rang++
      ,
if( DEBUGLEVEL_ell >= 2, print("  no point found on the quartic"))
      )
    )
  );

\\ end

if( DEBUGLEVEL_ell >= 1, 
  print("#S^(2)      = ",selmer));
  if( rang > selmer/2, rang = selmer);
if( DEBUGLEVEL_ell >= 1,
  strange = rang != selmer;
  if( strange,
  print("#E[K]/2E[K]>= ",rang)
, print("#E[K]/2E[K] = ",rang));
  print("#E[2]       = 4");
);
  rang = ceil(log(rang)/log(2))-2;
  selmer = valuation(selmer,2);
if( DEBUGLEVEL_ell >= 1,
  if( strange,
  print(selmer-2," >= rank  >= ",rang)
, print("rank        = ",rang));
  if( rang, print("points = ",listpoints));
);
  ell = ellinit([0,-(e1+e2+e3),0,e1*e2+e2*e3+e3*e1,-e1*e2*e3],1);
  listpoints = vecsort(listpoints,,2);
  if( ELLREDGENFLAG, listpoints = ellredgen(ell,listpoints));
  listpoints = concat(ellsort(elltorseven(ell)[3]),listpoints);

  return([rang,selmer,listpoints]);
}
{ellcount( c, d, KS2gen, listpointstriv=[]) =
my(found,listgen,listpointscount,m1,m2,lastloc,mask,i,d1,iaux,j,triv,pol,point,qf,solqf,para,point1,v);

if( DEBUGLEVEL_ell >= 4, print("    starting ellcount ",[c,d]));
if( DEBUGLEVEL_ell >= 4, print("    KS2gen = ",KS2gen));
if( DEBUGLEVEL_ell >= 4, print("    listpointstriv = ",listpointstriv));

  found = 0;
  listgen = KS2gen;
  listpointscount = [];

  m1 = m2 = 0; lastloc = -1;

  mask = 1 << #KS2gen;
  i = 1;
  while( i < mask,
    d1 = 1; iaux = i; j = 1;
    while( iaux, 
      if( iaux%2, d1 *= listgen[j]);
      iaux >>= 1; j++);
if( DEBUGLEVEL_ell >= 3, print("   d1 = ",d1));
    triv = 0;
    for( j = 1, #listpointstriv,
      if( listpointstriv[j][1] && issquare(d1*listpointstriv[j][1]),
        listpointscount = concat(listpointscount,[listpointstriv[j]]);
if( DEBUGLEVEL_ell >= 2, print("  trivial point"));
        triv = 1; m1++;
        if( degre(i) > lastloc, m2++);
        found = 1; lastloc = -1; break));
    if( !triv,
    pol = Pol([d1,0,c,0,d/d1]);
if( DEBUGLEVEL_ell >= 3, print("   quartic = y^2 = ",pol));
    point = ratpoint(pol,LIM1,1);
    if( point != [],
if( DEBUGLEVEL_ell >= 2, print("  point on the quartic"));
if( DEBUGLEVEL_ell >= 3, print(   point));
      m1++;
      listpointscount = concat(listpointscount,[d1*point[1]*point]);
      if( degre(i) > lastloc, m2++);
      found = 1; lastloc = -1
    ,
      if( locallysoluble(pol),
        if( degre(i) > lastloc, m2++; lastloc = degre(i));
\\        point = ratpoint2(pol,LIM3,1);
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Instead of solving directly y^2 = d1*x^4+c*x^2+d/d1,
\\ we solve first y^2 = d1*X^2+c*X+d/d1, then solve the quartic X = x^2
\\ which gives a new quartic
        qf = 2*[d1,c/2,0;c/2,d/d1,0;0,0,-1];
        solqf = qfsolve(qf);
        para = qfparam(qf,solqf,2)*['x^2,'x,1]~;
if( DEBUGLEVEL_ell >= 3, 
  print("   the conic y^2 = ",Pol([d1,c,d/d1]));
  print("   is parametrized by [x,y] = "subst([para[1]/para[2],para[3]/para[2]],'x,'t)));
        point1 = ratpoint2(para[1]*para[2],LIM3,1);
        if( point1 != [],
          if(#point1 == 2,
            para = subst(para,'x,point1[1])
          , point1 = [1,point1[2]/point1[1]^2,0];
            para = vector(3,ii,polcoeff(para[ii],2))
          );
          point = [point1[2]/para[2],para[3]/para[2]];
        , point = [];
        );
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        if( point != [],
if( DEBUGLEVEL_ell >= 2, print("  point on the quartic"));
if( DEBUGLEVEL_ell >= 3, print(   point));
          m1++;
          listpointscount = concat(listpointscount,[d1*point[1]*point]);
          if( degre(i) > lastloc, m2++);
          found = 1; lastloc = -1
        ,
if( DEBUGLEVEL_ell >= 2, print("  no point found on the quartic"))
          ))));
    if( found,
      found = 0;
      v = 0; iaux = (i>>1);
      while( iaux, iaux >>= 1; v++);
      mask >>= 1;
      listgen = vecextract(listgen,(1<<#listgen)-(1<<v)-1);
      i = (1<<v)
    , i++)
  );

  for( i = 1, #listpointscount,
   if( #listpointscount[i] > 1,
    if( subst('x^3+c*'x^2+d*'x,'x,listpointscount[i][1])-listpointscount[i][2]^2 != 0,
      error(" ellcount: WRONG POINT, please report ",i))));
if( DEBUGLEVEL_ell >= 4, print("    end of ellcount"));

  return([listpointscount,[m1,m2]]);
}
{ell2descent_viaisog(ell,help=[]) =
\\ Computation of the rank of the elliptic curve ell
\\ having rational 2-torsion, using the algorithm via 2-isogenies.
\\
\\ ell must be on the form
\\ y^2=x^3+ax^2+bx -> ell = [0,a,0,b,0]
\\ with a and b integers.

my(P,Pfact,tors,listpointstriv,KS2prod,KS2gen,listpoints,pointgen,n1,n2,certain,apinit,bpinit,np1,np2,listpoints2,aux1,aux2,certainp,rang,strange);

if( DEBUGLEVEL_ell >= 2, print("  Algorithm of 2-descent via isogenies"));
  if( #ell < 13, ell = ellinit(ell));

  if( ell.disc == 0,
    error(" ell2descent_viaisog: singular curve !!"));
  if( ell.a1 != 0 || ell.a3 != 0 || ell.a6 != 0, 
    error(" ell2descent_viaisog: the curve is not on the form [0,a,0,b,0]"));
  if( denominator(ell.a2) > 1 || denominator(ell.a4) > 1,
    error(" ell2descent_viaisog: non-integral coefficients"));

\\
\\ Working with the initial curve
\\

\\ Construction of trivial points: torsion

  P = Pol([1,ell.a2,ell.a4]);
  Pfact = factor(P)[,1];
  tors = #Pfact;
  listpointstriv = concat(help,elltorseven(ell)[3]);

\\ Construction of trivial points: small naive height

if( DEBUGLEVEL_ell >= 3, print("   Search for trivial points on the curve"));
  P *= 'x;
if( DEBUGLEVEL_ell >= 3, print("   Y^2 = ",P));
  listpointstriv = concat( listpointstriv, ratpoint(P,LIMTRIV,0));
if( DEBUGLEVEL_ell >= 1, print(" trivial points on E(Q) = ",listpointstriv); print());

  KS2prod = -abs(ell.a4);
  if( ell.a2^2-4*ell.a4 < 0, KS2prod *= -1);
  KS2gen = factor(KS2prod)[,1];

if( DEBUGLEVEL_ell >= 2,
  print("  #K(b,2)gen          = ",#KS2gen);
  print("  K(b,2)gen = ",KS2gen));

  listpoints = ellcount(ell.a2,ell.a4,KS2gen,listpointstriv);
  pointgen = listpoints[1];
if( DEBUGLEVEL_ell >= 1, print(" points on E(Q) = ",pointgen); print());
  n1 = listpoints[2][1]; n2 = listpoints[2][2];

  certain = (n1 == n2);
if( DEBUGLEVEL_ell >= 1,
  if( certain, 
    print("[E(Q):phi'(E'(Q))]  = ",1<<n1);
    print("#S^(phi')(E'/Q)     = ",1<<n2);
    print("#III(E'/Q)[phi']    = 1"); print()
  ,
    print("[E(Q):phi'(E'(Q))] >= ",1<<n1);
    print("#S^(phi')(E'/Q)     = ",1<<n2);
    print("#III(E'/Q)[phi']   <= ",1<<(n2-n1)); print())
);

\\
\\ Working with the isogeneous curve
\\

  apinit = -2*ell.a2; bpinit = ell.a2^2-4*ell.a4;
  KS2prod = -abs(bpinit);
  if( ell.a4 < 0, KS2prod *= -1);
  KS2gen = factor(KS2prod)[,1];

if( DEBUGLEVEL_ell >= 2,
  print("  #K(a^2-4b,2)gen     = ",#KS2gen);
  print("  K(a^2-4b,2)gen     = ",KS2gen));

\\ Construction of trivial points: torsion

  P = Pol([1,apinit,bpinit]);
  listpointstriv = elltorseven([0,apinit,0,bpinit,0])[3];
   
\\ Construction of trivial points: small naive height

if( DEBUGLEVEL_ell >= 3, print("   Search for trivial points on the curve"));
  P *= 'x;
if( DEBUGLEVEL_ell >= 3, print(" Y^2 = ",P));
  listpointstriv = concat( listpointstriv, ratpoint(P,LIMTRIV,0));
if( DEBUGLEVEL_ell >= 1, print(" trivial points on E'(Q) = ",listpointstriv); print());

  listpoints = ellcount(apinit,bpinit,KS2gen,listpointstriv);

if( DEBUGLEVEL_ell >= 1, print(" points on E'(Q) = ",listpoints[1]));
  np1 = listpoints[2][1]; np2 = listpoints[2][2];
  listpoints2 = vector(#listpoints[1],i,0);
  for( i = 1, #listpoints[1],
    listpoints2[i] = [0,0];
    aux1 = listpoints[1][i][1]^2;
    if( aux1 != 0,
      aux2 = listpoints[1][i][2];
      listpoints2[i][1] = aux2^2/aux1/4;
      listpoints2[i][2] = aux2*(bpinit-aux1)/aux1/8
    , listpoints2[i] = listpoints[1][i]));
if( DEBUGLEVEL_ell >= 1, print(" points on E(Q) = ",listpoints2); print());
  pointgen = concat(pointgen,listpoints2);

  certainp = (np1 == np2);
if( DEBUGLEVEL_ell >= 1,
  if( certainp,
    print("[E'(Q):phi(E(Q))]   = ",1<<np1);
    print("#S^(phi)(E/Q)       = ",1<<np2);
    print("#III(E/Q)[phi]      = 1"); print()
  ,
    print("[E'(Q):phi(E(Q))]  >= ",1<<np1);
    print("#S^(phi)(E/Q)       = ",1<<np2);
    print("#III(E/Q)[phi]     <= ",1<<(np2-np1)); print());

  if( !certain && (np2 > np1), print1(1<<(np2-np1)," <= "));
  print1("#III(E/Q)[2]       ");
  if( certain && certainp, print1(" "), print1("<"));
  print("= ",1<<(n2+np2-n1-np1));

  print("#E(Q)[2]            = ",1<<tors);
);
  rang = n1+np1-2;
if( DEBUGLEVEL_ell >= 1,
  if( certain && certainp,
    print("#E(Q)/2E(Q)         = ",(1<<(rang+tors)));
    print("rank                = ",rang); print()
  ,
    print("#E(Q)/2E(Q)        >= ",(1<<(rang+tors))); print();
    print("",rang," <= rank          <= ",n2+np2-2); print()
  ));

  strange = (n2+np2-n1-np1)%2;
  if( strange, 
if( DEBUGLEVEL_ell >= 1,
      print(" !!! III should be a square !!!"); print("hence"));
    if( certain, 
      np1++;
      certainp = (np1 == np2);
if( DEBUGLEVEL_ell >= 1,
        if( certainp,
          print("[E'(Q):phi(E(Q))]   = ",1<<np1);
          print("#S^(phi)(E/Q)       = ",1<<np2);
          print("#III(E/Q)[phi]      = 1"); print()
        ,
          print("[E'(Q):phi(E(Q))]  >= ",1<<np1);
          print("#S^(phi)(E/Q)       = ",1<<np2);
          print("#III(E/Q)[phi]     <= ",1<<(np2-np1)); print())
      )
    ,
    if( certainp,
      n1++;
      certain = (n1 == n2);
if( DEBUGLEVEL_ell >= 1, 
        if( certain,
          print("[E(Q):phi'(E'(Q))]   = ",1<<n1);
          print("#S^(phi')(E'/Q)       = ",1<<n2);
          print("#III(E'/Q)[phi']      = 1"); print()
        ,
          print("[E(Q):phi'(E'(Q))]  >= ",1<<n1);
          print("#S^(phi')(E'/Q)      = ",1<<n2);
          print("#III(E'/Q)[phi']    <= ",1<<(n2-n1)); print())
      )
    , n1++)
  );

if( DEBUGLEVEL_ell >= 1,
  print("#S^(2)(E/Q)           = ",1<<(n2+np2-1));
  if( !certain && (np2 > np1), print1(" ",1<<(np2-np1)," <= "));
  print1("#III(E/Q)[2]       ");
  if( certain && certainp, print1(" "), print1("<"));
  print("= ",1<<(n2+np2-n1-np1));
  print("#E(Q)[2]            = ",1<<tors);
);
  rang = n1+np1-2;
if( DEBUGLEVEL_ell >= 1,
  if( certain && certainp, 
    print("#E(Q)/2E(Q)         = ",(1<<(rang+tors))); print();
    print("rank                = ",rang); print()
  ,
    print("#E(Q)/2E(Q)        >= ",(1<<(rang+tors))); print();
    print(rang," <= rank          <= ",n2+np2-2); print())
  ));

\\ end of strange 
  
  pointgen = vecsort(pointgen,,2);
  if( ELLREDGENFLAG, pointgen = ellredgen(ell,pointgen));
  pointgen = concat(ellsort(elltorseven(ell)[3]),pointgen);
if( DEBUGLEVEL_ell >= 1, print("points = ",pointgen));
if( DEBUGLEVEL_ell >= 3, print("   end of ell2descent_viaisog"));

  return([rang,n2+np2-2+tors,pointgen]);
}

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\            HELP MESSAGES                \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
\\                  functions for elliptic curves
  addhelp(ell2descent_complete,
    "ell2descent_complete(e1,e2,e3): Performs a complete 2-descent on the elliptic curve y^2 = (x-e1)*(x-e2)*(x-e3). See ?ellQ_ellrank for the format of the output.");
  addhelp(ell2descent_gen,
    "ell2descent_gen((E,bnf,k=1,help=[]): E is a vector of the form [0,A,0,B,C], (or the result of ellinit of such a vector) A,B,C integers such that x^3+A*x^2+B*x+C; bnf is the corresponding bnfinit(,1); Performs 2-descent on the elliptic curve Ek: k*y^2=x^3+A*x^2+B*x+C. See ?ellQ_ellrank for the format of the output.");
  addhelp(ell2descent_viaisog,
    "ell2descent_viaisog(E,help=[]): E is an elliptic curve of the form [0,a,0,b,0], with a, b integers. Performs a 2-descent via isogeny on E. See ?ellQ_ellrank for the format of the output.");
  addhelp(ellQ_ellrank,
    "ellQ_ellrank(E,help=[]): E is any elliptic curve defined over Q. Returns a vector [r,s,v], where r is a lower bound for the rank of E, s is the rank of its 2-Selmer group and v is a list of independant points in E(Q)/2E(Q). If help is a vector of nontrivial points on E, the result might be faster. This function might be used in conjunction with elltors2(E). See also ?default_ellQ");
  addhelp(ellhalf,
    "ellhalf(E,P): returns the vector of all points Q on the elliptic curve E such that 2Q = P");
  addhelp(ellredgen,
    "ellredgen(E,v): returns a vector of smallest possible points on the elliptic curve E generating the same subgroup as v, up to torsion.");
  addhelp(ellsort,
    "ellsort(v): v being a vector of points on some elliptic curve, returns the vector v sorted according to the naive height.");
  addhelp(elltors2,
    "elltors2(E): for an elliptic curve E, returns the group E(K)[2], where K is the field of definition of the coefficients of E (Q, R, Qp or Fp)."); 
  addhelp(elltorseven,
    "elltorseven(E): for an elliptic curve E, returns 2-Sylow subgroup of E(K)_tors, where K is the field of definition of the coefficients of E: (Q, R, Qp or Fp).");


\\                  functions for polynomials
  addhelp(locallysoluble,
    "locallysoluble(pol): returns 1 if y^2=pol(x) is everywhere locally soluble, 0 otherwise.");
  addhelp(ratpoint,
    "ratpoint(pol,lim=1,singlepoint=1): search for rational points on y^2=pol(x), for about within the bounds given by lim. The coefficients of pol must be integral. If singlepoint=1, returns at most one point, otherwise as many as possible.");
  addhelp(redquartic,
    "redquartic(pol): reduction of the quartic pol using Cremona-Stoll algorithm. Returns [p,M], where p is the reduced quartic and M is the GL2(Z) transformation. Also works with other degree polynomials.");


\\                  functions for number fields
  addhelp(bnfpSelmer,
    "bnfpSelmer(K,S,p): K being a number field given by bnfinit, S an ideal of K, and p a prime number, computes a set of generators of the group K(S,p) = { x in K^* / K^*^p, v_P(x) = 0 (mod p) for all P coprime to S}");
  addhelp(reducemodsquares,
    "reducemodsquares(delta,d): delta being a t_POLMOD, returns another delta'=delta*z^2, such that delta' has a small coefficient in x^d.");


\\                  others
  addhelp(default_ellQ,
    "default_ellQ(DEBUGLEVEL_ell, LIM1, LIM3, LIMTRIV, ELLREDGENFLAG, COMPLETE, MAXPROB, LIMBIGPRIME): set the value of the global variables used for ellQ_ellrank() and other related functions. DEBUGLEVEL_ell: 0-5: choose the quantity of information printed during the computation (default=0: print nothing); LIM1 (resp LIM3): search limit for easy (resp hard) points on quartics; LIMTRIV: search limit for trivial points on elliptic curves; ELLREDGENFLAG: if != 0, try to reduce the generators at the end; COMPLETE: if != 0 and full 2-torsion, use complete 2-descent, otherwise via 2-isogeny; MAXPROB, LIMBIGPRIME: technical.");
/*  addhelp(DEBUGLEVEL_ell,
    "DEBUGLEVEL_ell: Choose a higher value of this global variable to have more details of the computations printed during the 2-descent algorithm. 0 = don't print anything; 1 = (default) just print the result; 2 = print more details including the Selmer group and the nontrivial quartics.");
*/
}

