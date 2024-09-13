\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\       Copyright (C) 2011 Denis Simon
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
  www.math.unicaen.fr/~simon/ell.gp

  *********************************************
  *          VERSION 06/04/2011               *
  *********************************************

  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  \\    Comment utiliser ce programme ?       \\
  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  Programme de calcul du rang des courbes elliptiques
  dans les corps de nombres.
  langage: GP
  pour l'utiliser, lancer gp, puis taper
  \r ell.gp

  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  \\  Description des principales fonctions   \\
  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  Explications succintes :
  definition du corps :
  bnf=bnfinit(y^2+1);
  (il est indispensable que la variable soit y).
  on peut ensuite poser : 
  X = Mod(y,bnf.pol);

  La fonction bnfellrank() accepte toutes les courbes sous la forme
  [a1,a2,a3,a4,a6]
  Les coefficients peuvent etre entiers ou non.
  L'algorithme utilise est celui de la 2-descente.
  La 2-torsion peut etre quelconque.
  Il suffit de taper : 

  gp > ell = [a1,a2,a3,a4,a6];
  gp > bnfellrank(bnf,ell)

  Retourne un vecteur [r,s,vec]
  ou r est le rang probable (c'est toujours une minoration du rang),
  s est le 2-rang du groupe de Selmer,
  vec est une liste de points dans E(K)/2E(K).

  Courbes avec #E[2](K) >= 2:
  ell doit etre sous la forme 
  y^2 = x^3 + A*^2 + B*x
  avec A et B entiers algebriques
  gp > ell = [0,A,0,B,0]
  gp > bnfell2descent_viaisog(ell)
  = algorithme de la 2-descente par isogenies
  Attention A et B doivent etre entiers

  Courbes avec #E[2](K) = 4: y^2 = (x-e1)*(x-e2)*(x-e3)
  -> bnfell2descent_complete(bnf,e1,e2,e3);
  = algorithme de la 2-descente complete
  Attention: les ei doivent etre entiers algebriques.


  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  \\       Affichage des calculs              \\
  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  On peut avoir plus ou moins de details de calculs avec
  DEBUGLEVEL_ell = 0;
  DEBUGLEVEL_ell = 1; 2; 3;...

*/


nf_scalar_or_multable_to_alg(nf, z) = {
  if (type(z) == "t_MAT", nfbasistoalg(nf, z[,1]), z);
}

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\          SCRIPT                             \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\    MANIPULATION OF GLOBAL VARIABLES         \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{default_ell(
  DEBUGLEVEL_ell_val:small = 0,
  LIM1_val:small = 2,
  LIM3_val:small = 4,
  LIMTRIV_val:small = 2,
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

  MAXPROB = MAXPROB_val;
  print("  MAXPROB = ",MAXPROB);

  LIMBIGPRIME = LIMBIGPRIME_val;
  print("  LIMBIGPRIME = ",LIMBIGPRIME);
}

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\    COMMON FUNCTIONS TO ell.gp AND ellQ.gp   \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{ellcomposeurst(urst1,urst2) =
local(u1 = urst1[1], r1 = urst1[2], s1 = urst1[3], t1 = urst1[4],
      u2 = urst2[1], r2 = urst2[2], s2 = urst2[3], t2 = urst2[4]);
  [u1*u2,u1^2*r2+r1,u1*s2+s1,u1^3*t2+s1*u1^2*r2+t1];
}
{ellinverturst(urst) =
local(u = urst[1], r = urst[2], s = urst[3], t = urst[4]);
  [1/u,-r/u^2,-s/u,(r*s-t)/u^3];
}
{mysubst(polsu,subsx) =
  if( type(lift(polsu)) == "t_POL",
    return(simplify(subst(lift(polsu),variable(lift(polsu)),subsx)))
  , return(simplify(lift(polsu))));
}
{degre(idegre) =
local(ideg = idegre, jdeg = 0);

  while( ideg >>= 1, jdeg++);
  return(jdeg);
}
{nfissquare(nf, a) = #nfsqrt(nf,a) > 0;
}
{nfsqrt( nf, a) =
\\ if a is a square in the number field nf returns [sqrt(a)], otherwise [].
local(alift,ta,py,pfact);

if( DEBUGLEVEL_ell >= 5, print("     starting nfsqrt ",a));
  if( a==0 || a==1, 
if( DEBUGLEVEL_ell >= 5, print("     end of nfsqrt ",a));
    return([a]));

  alift = lift(a);
  ta = type(a);
  if( !poldegree(alift), alift = polcoeff(alift,0));

  if( type(alift) != "t_POL",
    if( issquare(alift), 
if( DEBUGLEVEL_ell >= 5, print("     end of nfsqrt ",sqrtrat(alift)));
      return([sqrtrat(alift)])));

  if( poldegree(nf.pol) <= 1,
if( DEBUGLEVEL_ell >= 5, print("     end of nfsqrt ",[]));
    return([]));
  if( ta == "t_POL", a = Mod(a,nf.pol));

\\ the norm should be a square

  if( !issquare(norm(a)),
if( DEBUGLEVEL_ell >= 5, print("     end of nfsqrt ",[]));
      return([]));

\\ the real embeddings must all be >0

  for( i = 1, nf.r1,
    if( nfrealsign(nf,a,i) < 0,
if( DEBUGLEVEL_ell >= 5, print("     end of nfsqrt ",[]));
      return([])));

\\ factorization over nf of the polynomial X^2-a

  if( variable(nf.pol) == 'x,
    py = subst(nf.pol,'x,'y);
    pfact = lift(factornf('x^2-mysubst(alift,Mod('y,py)),py)[1,1])
  ,
    pfact = lift(factornf('x^2-a,nf.pol)[1,1]));
  if( poldegree(pfact) == 2,
if( DEBUGLEVEL_ell >= 5, print("     end of nfsqrt ",[]));
    return([]));
if( DEBUGLEVEL_ell >= 5, print("     end of nfsqrt ",pfact));
  return([subst(polcoeff(pfact,0),'y,Mod(variable(nf.pol),nf.pol))]);
}
{nfrealsign(nf,a,i) =
\\ return the sign of the algebraic number a in the i-th real embedding.
local(nf_roots,ay,prec0);

  if( a == 0, return(0));

  a = lift(a);
  if( type(a) != "t_POL",
    return(sign(a)));

  nf_roots = nf.roots;
  prec0 = default(realprecision);

  ay = 0;
  while( ay == 0 || precision(ay) < 38,

    ay = subst(a,variable(a),nf_roots[i]);

    if( ay == 0 || precision(ay) < 38,
if( DEBUGLEVEL_ell >= 3, 
  print(" **** Warning: doubling the real precision in nfrealsign **** ",
        2*default(realprecision)));
      default(realprecision,2*default(realprecision));
      nf_roots = real(polroots(nf.pol))
    )
  );
  default(realprecision,prec0);

  return(sign(ay));
}
{sqrtrat(a) =
  sqrtint(numerator(a))/sqrtint(denominator(a));
}

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\      FUNCTIONS SPECIFIC TO ell.gp          \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{nfpolratroots(nf,pol) =
local(f,ans);
  f = nffactor(nf,lift(pol))[,1];
  ans = [];
  for( j = 1, #f, 
    if( poldegree(f[j]) == 1,
      ans = concat(ans,[-polcoeff(f[j],0)/polcoeff(f[j],1)])));
  return(ans);
}
{mynfhilbert2(nf,a,b,p) =
\\ p is a prime output by idealprimedec() above 2
local(v,res);

if( DEBUGLEVEL_ell >= 5, print("     starting mynfhilbert2"));
  v = idealval(nf,a,p)\2;
  if( v, a *= nfbasistoalg(nf,p[2])^(-2*v));

  v = idealval(nf,b,p)\2;
  if( v, b *= nfbasistoalg(nf,p[2])^(-2*v));

  if( nfqp_soluble(nf,a*'x^2+b,initp(nf,p)), res = 1, res = -1);
if( DEBUGLEVEL_ell >= 5, print("     end of mynfhilbert2",res));
  return(res);
}
{mynfhilbertp(nf,a,b,p) =
\\ calcule le symbole de Hilbert quadratique local (a,b)_p
\\ * en l'ideal premier p du corps nf,
\\ * a et b sont des elements non nuls de nf, sous la forme
\\ * de polmods ou de polynomes, et p renvoye par idealprimedec().

if( DEBUGLEVEL_ell >= 5, print("     starting mynfhilbertp at ",p));
  if( a == 0 || b == 0, error("mynfhilbertp: argument = 0"));
  if( p.p == 2,
if( DEBUGLEVEL_ell >= 5, print("     end of mynfhilbertp"));
    return(mynfhilbert2(nf,a,b,p)));
  if( type(a) != "t_POLMOD", a = Mod(a,nf.pol));
  if( type(b) != "t_POLMOD", b = Mod(b,nf.pol));
return(nfhilbert(nf,a,b,p));
}
{ideallistfactor(nf,listfact) =
local(Slist,S1,test,k);

if( DEBUGLEVEL_ell >= 5, print("     starting ideallistfactor"));
  Slist = []; test = 1;
  for( i = 1, #listfact,
    if( listfact[i] == 0, next);
    S1 = idealfactor(nf,listfact[i])[,1];
    for( j = 1, #S1, k = #Slist;
      for( k = 1, #Slist,
        if( Slist[k] == S1[j], test = 0; break));
      if( test, Slist = concat(Slist,[S1[j]]), test = 1);
     ));
if( DEBUGLEVEL_ell >= 5, print("     end of ideallistfactor"));
  return(Slist);
}
{mynfhilbert(nf,a,b) =
\\ calcule le symbole de Hilbert quadratique global (a,b):
\\ =1 si l'equation X^2-aY^2-bZ^2=0 a une solution non triviale,
\\ =-1 sinon,
\\ a et b doivent etre non nuls.
local(al,bl,S);

if( DEBUGLEVEL_ell >= 4, print("    starting mynfhilbert ",[a,b]));
  if( a == 0 || b == 0, error("mynfhilbert : argument = 0"));
  al = lift(a); bl = lift(b);

\\ solutions locales aux places reelles 

  for( i = 1, nf.r1,
    if( nfrealsign(nf,al,i) < 0 && nfrealsign(nf,bl,i) < 0,
if( DEBUGLEVEL_ell >= 3, print("   mynfhilbert: no solution at infinity"));
if( DEBUGLEVEL_ell >= 4, print("    end of mynfhilbert"));
      return(-1))
  );

  if( type(a) != "t_POLMOD", a = Mod(a,nf.pol));
  if( type(b) != "t_POLMOD", b = Mod(b,nf.pol));

\\  solutions locales aux places finies (celles qui divisent 2ab)

  S = ideallistfactor(nf,[2,a,b]);
  forstep ( i = #S, 2, -1,
\\ d'apres la formule du produit on peut eviter un premier 
    if( mynfhilbertp(nf,a,b, S[i]) == -1,
if( DEBUGLEVEL_ell >= 3, print("   mynfhilbert: no solution at: ",S[i]));
if( DEBUGLEVEL_ell >= 4, print("    end of mynfhilbert"));
      return(-1)));
if( DEBUGLEVEL_ell >= 4, print("    end of mynfhilbert"));
  return(1);
}
{initp( nf, p) =
\\   pp[1] est l'ideal sous forme reduite 
\\   pp[2] est un entier de Zk avec une valuation 1 en p
\\   pp[3] est la valuation de 2 en p
\\   pp[4] sert a detecter les carres dans Qp 
\\     si p|2 il faut la structure de Zk/p^(1+2v) d'apres Hensel
\\     sinon il suffit de calculer x^(N(p)-1)/2 
\\   pp[5] est un systeme de representants de Zk/p 
\\     c'est donc un ensemble de cardinal p^f .
local(idval,pp);

if( DEBUGLEVEL_ell >= 5, print("     starting initp for p = ",p));
  idval = idealval(nf,2,p);
  pp=[ p, nfbasistoalg(nf,p[2]), idval, 0, repres(nf,p) ];
  if( idval,
    pp[4] = idealstar(nf,idealpow(nf,p,1+2*idval)),
    pp[4] = p.p^p.f\2 );
if( DEBUGLEVEL_ell >= 5, print("     end of initp"));
  return(pp);
}
{deno(num) = 
\\ calcule un denominateur du polynome num

  if( num == 0, return(1));
  if( type(num) == "t_POL",
     return(denominator(content(num))));
  return(denominator(num));
}
{nfratpoint(nf,pol,lim,singlepoint=1) =
\\ Si singlepoint == 1, cherche un seul point, sinon plusieurs.
local(compt1,compt2,deg,n,AA,point,listpoints,vectx,evpol,sq,xpol);

if( DEBUGLEVEL_ell >= 4,
  print("    starting nfratpoint with pol = ",pol);
  print("    lim = ",lim));

  compt1 = 0; compt2 = 0;
  deg = poldegree(pol); n = poldegree(nf.pol);
  AA = lim<<1;
  if( !singlepoint, listpoints = []);

\\ cas triviaux 
  sq = nfsqrt(nf,polcoeff(pol,0));
  if( sq!= [],
    point = [ 0, sq[1], 1];
    if( singlepoint,
if( DEBUGLEVEL_ell >= 4, print("    end of nfratpoint"));
      return(point));
    listpoints = concat(listpoints,[point])
  );
  sq = nfsqrt(nf,pollead(pol));
  if( sq != [],
    point = [ 1, sq[1], 0];
    if( singlepoint,
if( DEBUGLEVEL_ell >= 4, print("    end of nfratpoint"));
      return(point));
    listpoints = concat(listpoints,[point])
  );
  
\\ boucle generale
  point = [];
  vectx = vector(n,i,[-lim,lim]);
  for( denoz = 1, lim,
    if( poldegree(pol)%2 == 0 &&
       !issquare(Mod(norm(pollead(pol)),denoz)), next);
    forvec( xx = vectx,
      if( denoz == 1 || gcd(content(xx),denoz) == 1,
        xpol = nfbasistoalg(nf,xx~);
        evpol = subst(pol,'x,xpol/denoz);
        sq = nfsqrt(nf,evpol);
        if( sq != [],
          point = [xpol/denoz, sq[1], 1];
          if( singlepoint, break(2));
          listpoints = concat(listpoints,[point])));
  ));

  if( singlepoint, listpoints = point);
if( DEBUGLEVEL_ell >= 3, print("   points found by nfratpoint = ",listpoints));
if( DEBUGLEVEL_ell >= 4, print("    end of nfratpoint"));
  return(Vec(listpoints));
}
{repres(nf,p) =
\\ calcule un systeme de representants Zk/p
local(fond,mat,f,rep,pp,ppi,pp2,jppi,gjf);

if( DEBUGLEVEL_ell >= 5, print("     starting repres"));
  fond = [];
  mat = idealhnf(nf,p);
  for( i = 1, #mat,
    if( mat[i,i] != 1, fond = concat(fond, nf.zk[i])));
  f = #fond;
  pp = p.p;
  rep = vector(pp^f,i,0);
  rep[1] = 0;
  ppi = 1;
  pp2 = pp\2;
  for( i = 1, f,
    for( j = 1, pp-1,
      if( j <= pp2, gjf = j*fond[i], gjf = (j-pp)*fond[i]);
      jppi = j*ppi;
      for( k = 0, ppi-1, rep[jppi+k+1] = rep[k+1]+gjf ));
    ppi *= pp);
if( DEBUGLEVEL_ell >= 5, print("     end of repres"));
  return(Mod(rep,nf.pol));
}
{val(nf,num,p) =
  if( num == 0, 32000, idealval(nf,lift(num),p));
}
{nfissquaremodpodd( nf, a, p) =
\\ Return 1 if a is a p-adic square, 0 otherwise.
\\ Only for a prime ideal p coprime to 2.
\\ a = t_POLMOD
local(v,ap,norme,den);

if( DEBUGLEVEL_ell >= 5, print("     starting nfissquaremodpodd"));
  if( a == 0,
if( DEBUGLEVEL_ell >= 5, print("     end of nfissquaremodpodd"));
    return(1));
  v = idealval(nf,lift(a),p);
  if( v%2,
if( DEBUGLEVEL_ell >= 5, print("     end of nfissquaremodpodd"));
    return(0));
  ap = a/nfbasistoalg(nf,p[2])^v;

  norme = (p.p^p.f-1)/2;
  den = denominator(content(lift(ap)))%p.p;
  if( sign(den), ap *= Mod(1,p.p));
  ap = ap^norme-1;
  if( ap == 0,
if( DEBUGLEVEL_ell >= 5, print("     end of nfissquaremodpodd"));
    return(1));
  ap = lift(lift(ap));
  if( idealval(nf,ap,p) > 0,
if( DEBUGLEVEL_ell >= 5, print("     end of nfissquaremodpodd"));
    return(1));
if( DEBUGLEVEL_ell >= 5, print("     end of nfissquaremodpodd"));
  return(0);
}
{nfissquaremodp( nf, a, p, zinit) =
\\ a is an algebraic integer of nf
\\ returns 1 if a is a square modulo the prime ideal p, 0 otherwise
\\ a = t_POLMOD
local(valap,zlog);

if( DEBUGLEVEL_ell >= 5, print("     starting nfissquaremodp ",[a,p,zinit]));
  if( a == 0,
if( DEBUGLEVEL_ell >= 5, print("     end of nfissquaremodp"));
    return(1));

  if( p.p != 2,
if( DEBUGLEVEL_ell >= 5, print("     end of nfissquaremodp"));
    return(nfissquaremodpodd(nf,a,p)));

  valap = idealval(nf,a,p);
  if( valap%2,
if( DEBUGLEVEL_ell >= 5, print("     end of nfissquaremodp"));
    return(0));
  if( valap,
    zlog = ideallog(nf,a*(nf_scalar_or_multable_to_alg(nf,p[5])/p.p)^valap,zinit)
  ,
    zlog = ideallog(nf,a,zinit));
  for( i = 1, #zinit[2][2],
    if( !(zinit[2][2][i]%2) && (zlog[i]%2),
if( DEBUGLEVEL_ell >= 5, print("     end of nfissquaremodp"));
      return(0)));
if( DEBUGLEVEL_ell >= 5, print("     end of nfissquaremodp"));
  return(1);
}
{nfissquaremodpq( nf, a, p, q) =
\\ cette fonction renvoie 1 si a est un carre
\\ ?inversible? modulo P^q et 0 sinon.
\\ P divise 2, et ?(a,p)=1?.
local(vala,zinit,zlog);

if( DEBUGLEVEL_ell >= 5, print("     starting nfissquaremodpq ",[a,p,q]));
  if( a == 0,
if( DEBUGLEVEL_ell >= 5, print("     end of nfissquaremodpq"));
    return(1));
  vala = idealval(nf,a,p);
  if( vala >= q,
if( DEBUGLEVEL_ell >= 5, print("     end of nfissquaremodpq"));
    return(1));
  if( vala%2,
if( DEBUGLEVEL_ell >= 5, print("     end of nfissquaremodpq"));
    return(0));
  zinit = idealstar(nf,idealpow(nf,p,q-vala),2);
  zlog = ideallog(nf,a*nf_scalar_or_multable_to_alg(nf,p[5]/2)^vala,zinit);
  for( i = 1, #zinit[2][2],
    if( !(zinit[2][2][i]%2) && (zlog[i]%2),
if( DEBUGLEVEL_ell >= 5, print("     end of nfissquaremodpq"));
      return(0)));
if( DEBUGLEVEL_ell >= 5, print("     end of nfissquaremodpq"));
  return(1);
}
{nfsqrtmodpq(nf,a,p,q) =
\\ suppose que a est un carre modulo p^q
\\ et renvoie x tel que  x^2 = a mod p^q.
local(p_hnf,f,aaa,qq,e,xx,yy,r,aux,b,m,vp,inv2x,zinit,zlog,expo,p_ini,non_sq,id);

if( DEBUGLEVEL_ell >= 5, print("    starting nfsqrtmodpq ",a,p,q));
  if( a == 0 || a == 1,
if( DEBUGLEVEL_ell >= 5, print("    end of nfsqrtmodpq"));
    return(a));
  f = idealval(nf,a,p);
  if( f >= q,
if( DEBUGLEVEL_ell >= 5, print("     end of nfsqrtmodpq"));
    return(0));
if( f%2, error("nfsqrtmodpq: a is not a square, odd valuation"));
  a = nfalgtobasis(nf,a);
  if( f, aaa = nfeltpow(nf,nfeltdiv(nf,a,nf_scalar_or_multable_to_alg(nf,p[5]/p.p)),f), aaa = a);
  p_hnf = idealhnf(nf,p);
  p_ini = nfmodprinit(nf,p);
if( DEBUGLEVEL_ell >= 5, print("     p_hnf = ",p_hnf));
if( DEBUGLEVEL_ell >= 5, print("     p_ini = ",p_ini));

  if( p.p != 2,
\\ first case : p is odd

\\ Shanks sqrt algorithm
if( DEBUGLEVEL_ell >= 5, print("     Shanks sqrt algorithm"));
    non_sq = nfrandintmodid(nf,p_hnf);
    while( nfissquaremodpodd(nf,non_sq,p), non_sq = nfrandintmodid(nf,p_hnf));
    non_sq = nfalgtobasis(nf,non_sq);
    qq = ( p.p^p.f -1) \ 2;
    e = 1; while( !(qq%2), e++; qq \= 2);
    non_sq = nfeltpowmodpr(nf,non_sq,qq,p_ini);
    yy = non_sq; r = e;
    xx  = nfeltpowmodpr(nf,aaa,qq\2,p_ini);
    aux = nfeltmulmodpr(nf,aaa,xx,p_ini);
    b   = nfeltmulmodpr(nf,aux,xx,p_ini); 
    xx = aux;
    aux = b; m = 0;
    while( !val(nf,nfbasistoalg(nf,aux)-1,p),
       m++;
       aux = nfeltpowmodpr(nf,aux,2,p_ini)
    );
    while( m,
      if( m == r, error("nfsqrtmodpq: m = r"));
      aux = nfeltpowmodpr(nf,yy,1<<(r-m-1),p_ini);
      yy  = nfeltpowmodpr(nf,aux,2,p_ini);
      r = m;
      xx = nfeltmulmodpr(nf,xx,aux,p_ini);
      b  = nfeltmulmodpr(nf,b,yy,p_ini);
      aux = b;m = 0;
      while( !val(nf,nfbasistoalg(nf,aux)-1,p),
        m++;
        aux = nfeltpowmodpr(nf,aux,2,p_ini)
      )
    );

\\ lift de Hensel
\\

    xx = nfbasistoalg(nf,xx);
    aaa= nfbasistoalg(nf,aaa);
    if( q > 1,
if( DEBUGLEVEL_ell >= 5, print("     Hensel lifting"));
      vp = idealval(nf,xx^2-aaa,p);
      if( vp < q-f,
        yy = 2*xx;
        inv2x = nfbasistoalg(nf,idealaddtoone(nf,yy,p)[1])/yy;
        while( vp < q, vp++; xx -= (xx^2-aaa)*inv2x);
      );
      if( f, xx *= nfbasistoalg(nf,p[2])^(f\2));
    );
    xx = mynfeltreduce(nf,xx,idealpow(nf,p,q))
  ,
\\ cas ou p divise 2 */
    if( q-f > 1, id = idealpow(nf,p,q-f), id = p_hnf);
    zinit = idealstar(nf,id,2);
    zlog = ideallog(nf,aaa,zinit);
    xx = 1;
    for( i = 1, #zlog,
      expo = zlog[i];
      if( expo,
        if( !expo%2,
          expo = expo>>1
        , aux = zinit[2][i];
          expo = expo*((aux+1)>>1)%aux
        );
        xx *= nfbasistoalg(nf,zinit[2][3][i])^expo
      )
    );
    if( f,
      xx *= nfbasistoalg(nf,p[2])^(f>>1);
      id = idealpow(nf,p,q));
    xx = mynfeltreduce(nf,xx,id);
  );
if( DEBUGLEVEL_ell >= 5, print("    end of nfsqrtmodpq ",xx));
  return(xx);
}
{nflemma6( nf, pol, p, nu, xx) =
local(gx,gpx,lambda,mu);

if( DEBUGLEVEL_ell >= 5, print("     starting nflemma6"));
  gx = subst( pol, 'x, xx);
  if( nfissquaremodpodd(nf,gx,p),
if( DEBUGLEVEL_ell >= 5, print("     end of nflemma6"));
    return(1));
  gpx = subst( pol', 'x, xx);
  lambda = val(nf,gx,p);mu = val(nf,gpx,p);

  if( lambda>2*mu,
if( DEBUGLEVEL_ell >= 5, print("     end of nflemma6"));
    return(1));
  if( (lambda >= 2*nu)  && (mu >= nu),
if( DEBUGLEVEL_ell >= 5, print("     end of nflemma6"));
    return(0));
if( DEBUGLEVEL_ell >= 5, print("     end of nflemma6"));
  return(-1);
}
{nflemma7( nf, pol, p, nu, xx, zinit) =
local(gx,gpx,v,lambda,mu,q);

if( DEBUGLEVEL_ell >= 5, print("entree dans nflemma7 ",[xx,nu]));
  gx = subst( pol, 'x, xx);
  if( nfissquaremodp(nf,gx,p,zinit),
if( DEBUGLEVEL_ell >= 5, print("fin de nflemma7"));
    return(1));
  gpx = subst( pol', 'x, xx);
  v = p[3];
  lambda = val(nf,gx,p);mu = val(nf,gpx,p);
  if( lambda>2*mu,
if( DEBUGLEVEL_ell >= 5, print("fin de nflemma7"));
    return(1));
  if( nu > mu,
    if( lambda%2,
if( DEBUGLEVEL_ell >= 5, print("fin de nflemma7"));
      return(-1));
    q = mu+nu-lambda;
    if( q > 2*v,
if( DEBUGLEVEL_ell >= 5, print("fin de nflemma7"));
      return(-1));
    if( nfissquaremodpq(nf,gx*nf_scalar_or_multable_to_alg(nf,p[5]/2)^lambda,p,q),
if( DEBUGLEVEL_ell >= 5, print("fin de nflemma7"));
      return(1))
  ,
    if( lambda >= 2*nu,
if( DEBUGLEVEL_ell >= 5, print("fin de nflemma7"));
      return(0));
    if( lambda%2,
if( DEBUGLEVEL_ell >= 5, print("fin de nflemma7"));
      return(-1));
    q = 2*nu-lambda;
    if( q > 2*v,
if( DEBUGLEVEL_ell >= 5, print("fin de nflemma7"));
      return(-1));
    if( nfissquaremodpq(nf,gx*nf_scalar_or_multable_to_alg(nf,p[5]/2)^lambda,p,q),
if( DEBUGLEVEL_ell >= 5, print("fin de nflemma7"));
      return(0))
  );
if( DEBUGLEVEL_ell >= 5, print("fin de nflemma7"));
  return(-1);
}
{nfzp_soluble( nf, pol, p, nu, pnu, x0) =
local(result,pnup,lrep);

if( DEBUGLEVEL_ell >= 5, print("entree dans nfzp_soluble ",[lift(x0),nu]));
  if( p[3] == 0,
    result = nflemma6(nf,pol,p[1],nu,x0),
    result = nflemma7(nf,pol,p[1],nu,x0,p[4]));
  if( result == +1,
if( DEBUGLEVEL_ell >= 5, print("fin de nfzp_soluble"));
    return(1));
  if( result == -1,
if( DEBUGLEVEL_ell >= 5, print("fin de nfzp_soluble"));
    return(0));
  pnup = pnu*p[2];
  lrep = #p[5];
  nu++;
  for( i = 1, lrep,
    if( nfzp_soluble(nf,pol,p,nu,pnup,x0+pnu*p[5][i]),
if( DEBUGLEVEL_ell >= 5, print("fin de nfzp_soluble"));
      return(1)));
if( DEBUGLEVEL_ell >= 5, print("fin de nfzp_soluble"));
  return(0);
}
{mynfeltmod(nf,a,b) =
local(qred);

  qred = round(nfalgtobasis(nf,a/b));
  qred = a-b*nfbasistoalg(nf,qred);
  return(qred);
}
{mynfeltreduce(nf,a,id) =
  nfbasistoalg(nf,nfeltreduce(nf,nfalgtobasis(nf,a),id));
}
{nfrandintmodid( nf, id) =
local(res);

if( DEBUGLEVEL_ell >= 5, print("entree dans nfrandintmodid"));
  res = 0;
  while( !res,
    res = nfrandint(nf,0);
    res = mynfeltreduce(nf,res,id));
if( DEBUGLEVEL_ell >= 5, print("fin de nfrandintmodid"));
  return(res);
}
{nfrandint( nf, borne) =
local(d,res);

if( DEBUGLEVEL_ell >= 5, print("     starting nfrandint"));
  d = poldegree(nf.pol);
  res = vectorv(d,i,if( borne, random(borne<<1)-borne, random()));
  res = nfbasistoalg(nf,res);
if( DEBUGLEVEL_ell >= 5, print("     end of nfrandint"));
  return(res);
}
{nfqp_solublebig( nf, pol, p,ap=0,b=1) =
local(deg,xx,z,Px,cont,pi,pol2,Roots);

if( DEBUGLEVEL_ell >= 4, print("    starting nfqp_solublebig avec ",p.p));
  deg = poldegree(pol);

  if( nfissquaremodpodd(nf,polcoeff(pol,0),p),
if( DEBUGLEVEL_ell >= 4, print("    end of nfqp_solublebig"));
    return(1));
  if( nfissquaremodpodd(nf,pollead(pol),p),
if( DEBUGLEVEL_ell >= 4, print("    end of nfqp_solublebig"));
    return(1));

\\ on tient compte du contenu de pol
  cont = idealval(nf,polcoeff(pol,0),p);
  for( i = 1, deg,
    if( cont, cont = min(cont,idealval(nf,polcoeff(pol,i),p))));
  if( cont, pi = nf_scalar_or_multable_to_alg(nf,p[5]/p.p));
  if( cont > 1, pol *= pi^(2*(cont\2)));

\\ On essaye des valeurs de x au hasard
  if( cont%2,
    pol2 = pol*pi
  , pol2 = pol;
    for( i = 1, MAXPROB,
      xx = nfrandint(nf,0);
      z = 0; while( !z, z = random());
      xx = -ap*z+b*xx;
      Px=polcoeff(pol,deg);
      forstep (j=deg-1,0,-1,Px=Px*xx+polcoeff(pol,j));
      Px *= z^deg;
      if( nfissquaremodpodd(nf,Px,p),
if( DEBUGLEVEL_ell >= 4, print("    end of nfqp_solublebig"));
        return(1));
    )
  );

\\ On essaye les racines de pol
  Roots = nfpolrootsmod(nf,pol2,p);
  pi = nfbasistoalg(nf,p[2]);
  for( i = 1, #Roots,
    if( nfqp_solublebig(nf,subst(pol,'x,pi*'x+Roots[i]),p),
if( DEBUGLEVEL_ell >= 4, print("fin de nfqp_solublebig"));
      return(1)));

if( DEBUGLEVEL_ell >= 4, print("    end of nfqp_solublebig"));
  return(0);
}
{nfpolrootsmod(nf,pol,p) =
\\ calcule les racines modulo l'ideal p du polynome pol.
\\ p est un ideal premier de nf, sous la forme idealprimedec
local(factlist,sol);

  factlist = nffactormod(nf,pol,p)[,1];
  sol = [];
  for( i = 1, #factlist,
    if( poldegree(factlist[i]) == 1,
      sol = concat(sol, [-polcoeff(factlist[i],0)/polcoeff(factlist[i],1)])));
  return(sol);
}
{nfqp_soluble( nf, pol, p) =
\\ p is a prime output by initp()
\\ return 1 if pol(x) = y^2 has a p-adic rational solution
\\ (possibly oo)
\\ coeff of nfqfsoluble must be integers.

if( DEBUGLEVEL_ell >= 5, print("     starting nfqp_soluble ",p));
if( DEBUGLEVEL_ell >= 5, print("     pol = ",pol));
  if( nfissquaremodp(nf,pollead(pol),p[1],p[4]),
if( DEBUGLEVEL_ell >= 5, print("     end of nfqp_soluble"));
    return(1));
  if( nfissquaremodp(nf,polcoeff(pol,0),p[1],p[4]),
if( DEBUGLEVEL_ell >= 5, print("     end of nfqp_soluble"));
    return(1));
  if( nfzp_soluble(nf,pol,p,0,1,0),
if( DEBUGLEVEL_ell >= 5, print("     end of nfqp_soluble"));
    return(1));
  if( nfzp_soluble(nf,polrecip(pol),p,1, p[2],0),
if( DEBUGLEVEL_ell >= 5, print("     end of nfqp_soluble"));
    return(1));
if( DEBUGLEVEL_ell >= 5, print("     end of nfqp_soluble"));
  return(0);
}
{nflocallysoluble( nf, pol, r=0,a=1,b=1) =
\\ Test whether Y^2 = pol is Everywhere Locally Soluble
local(pol0,plist,add,ff,p,Delta,c,s,t);

if( DEBUGLEVEL_ell >= 4, print("    starting nflocallysoluble ",[pol,r,a,b]));
  pol0 = pol;

\\ places finies

  pol *= deno(content(lift(pol)))^2;
  for( ii = 1, 3,
    if( ii == 1, plist = idealprimedec(nf,2));
    if( ii == 2 && r, plist = idealfactor(nf,poldisc(pol0/pollead(pol0))/pollead(pol0)^6/2^12)[,1]);
    if( ii == 2 && !r, plist = idealfactor(nf,poldisc(pol0))[,1]);
    if( ii == 3,
      add = idealadd(nf,a,b);
      ff = factor(idealnorm(nf,add))[,1];
      addprimes(ff);
if( DEBUGLEVEL_ell >= 4, print("    list of primes = ",ff));
      plist = idealfactor(nf,add)[,1]);
      for( i = 1, #plist,
        p =  plist[i];  
if( DEBUGLEVEL_ell >= 3, print("   p = ",p));
        if( p.p < LIMBIGPRIME || !LIMBIGPRIME,
          if( !nfqp_soluble(nf,pol,initp(nf,p)),
if( DEBUGLEVEL_ell >= 2, print("  not ELS at ",p));
if( DEBUGLEVEL_ell >= 4, print("    end of nflocallysoluble"));
            return(0)),
          if( !nfqp_solublebig(nf,pol,p,r/a,b),
if( DEBUGLEVEL_ell >= 2, print("  not ELS at ",p.p," ( = big prime )"));
if( DEBUGLEVEL_ell >= 4, print("    end of nflocallysoluble"));
            return(0))));
  );

\\ places reelles 

/*
  Peter Bruin, August 2014: code adapted to fix precision problems.

  We want to count the number of real roots of a polynomial of the
  form x^4 + s*x^2 + a*x + b.  The action of complex conjugation on
  the complex roots determines the sign of the discriminant Delta:

    cycle type () or (1 2)(3 4), i.e. 4 or 0 real roots: Delta > 0
    cycle type (1 2), i.e. exactly 2 real roots: Delta < 0

  We may assume Delta > 0 and have to decide whether there are 0 or 4
  real roots.  We apply Sturm's theorem to the following sequence
  (with t = 2*s^3 - 8*b*s + 9*a^2):

    p0 = x^4 + s*x^2 + a*x + b
    p1 = 4*x^3 + 2*s*x + a
    p2 = x*p1 - 4*p0
       = -2*s*x^2 - 3*a*x - 4*b
    p3 = (3*a - 2*s*x)*p2 - s^2*p1
       = -t*x - a*(s^2 + 12*b)
    p4 = Delta

  Note that the case of 4 real roots can only occur if both s and t
  are nonzero, otherwise the Sturm sequence is shorter than this one.
  By Sturm's theorem, the number of roots equals

    (sign changes in [1, -4, -2*s, t, Delta]) -
    (sign changes in [1, 4, -2*s, -t, Delta]).

  Hence the only way to have 4 real roots is if s < 0 and t < 0.
*/
  if( nf.r1,
    c = pollead(pol);
    pol /= c;
    pol = subst(pol, 'x, 'x - polcoeff(pol, 3)/4);
    s = polcoeff(pol, 2);
    t = 2*s^3 - 8*s*polcoeff(pol, 0) + 9*polcoeff(pol, 1)^2;
    Delta = poldisc(pol);
    for( i = 1, nf.r1,
      if( nfrealsign(nf,c,i) > 0, next);
      if( nfrealsign(nf,Delta,i) < 0, next);
      if( nfrealsign(nf,s,i) < 0 && nfrealsign(nf,t,i) < 0, next);
if( DEBUGLEVEL_ell >= 2, print("  not ELS at infinity"));
if( DEBUGLEVEL_ell >= 4, print("    end of nflocallysoluble"));
      return(0);
  ));
if( DEBUGLEVEL_ell >= 2, print("  quartic ELS "));
if( DEBUGLEVEL_ell >= 4, print("    end of nflocallysoluble"));
  return(1);
}
{nfellcount( nf, c, d, KS2gen, pointstriv) =
local(found,listgen,listpointscount,m1,m2,lastloc,mask,i,d1,iaux,j,triv,pol,point,deuxpoints,aux,v);

if( DEBUGLEVEL_ell >= 4, print("    starting nfellcount ",[c,d]));
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
if( DEBUGLEVEL_ell >= 2, print("  d1 = ",d1));
    triv = 0;
    for( j = 1, #pointstriv,
      if( pointstriv[j][3]*pointstriv[j][1]
        && nfissquare(nf,d1*pointstriv[j][1]*pointstriv[j][3]),
        listpointscount = concat(listpointscount,[pointstriv[j]]);
if( DEBUGLEVEL_ell >= 2, print("  comes from a trivial point"));
        triv = 1; m1++;
        if( degre(i) > lastloc, m2++);
        found = 1; lastloc = -1; break
    ));
    if( !triv,
    pol = Pol([d1,0,c,0,d/d1]);
if( DEBUGLEVEL_ell >= 3, print("   quartic: y^2 = ",pol));
    point = nfratpoint(nf,pol,LIM1,1);
    if( point != [],
if( DEBUGLEVEL_ell >= 2, print("  point on the quartic"));
if( DEBUGLEVEL_ell >= 3, print("   ",point));
      m1++;
      if( point[3] != 0,
        aux = d1*point[1]/point[3]^2;
        deuxpoints = [ aux*point[1], aux*point[2]/point[3] ]
      ,
        deuxpoints = [0]);
      listpointscount = concat(listpointscount,[deuxpoints]);
      if( degre(i) > lastloc, m2++);
      found = 1; lastloc = -1
    ,
      if( nflocallysoluble(nf,pol),
        if( degre(i) > lastloc, m2++; lastloc = degre(i));
        point = nfratpoint(nf,pol,LIM3,1);
        if( point != [],
if( DEBUGLEVEL_ell >= 2, print("  point on the quartic"));
if( DEBUGLEVEL_ell >= 3, print("   ",point));
          m1++;
          aux = d1*point[1]/point[3]^2;
          deuxpoints = [ aux*point[1], aux*point[2]/point[3] ];
          listpointscount = concat(listpointscount,[deuxpoints]);
          if( degre(i) > lastloc, m2++);
          found = 1; lastloc = -1
        ,
if( DEBUGLEVEL_ell >= 2, print("  no point found on the quartic"));
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
      error("nfellcount: WRONG POINT = ",listpointscount[i]))));
if( DEBUGLEVEL_ell >= 4, print("    end of nfellcount"));
  return([listpointscount,[m1,m2]]);
}
{bnfell2descent_viaisog( bnf, ell) =
\\ Calcul du rang des courbes elliptiques avec 2-torsion
\\ dans le corps de nombres bnf
\\ par la methode des 2-isogenies.
\\
\\ ell = [a1,a2,a3,a4,a6]
\\ y^2+a1xy+a3y=x^3+a2x^2+a4x+a6
\\
\\ ell doit etre sous la forme
\\ y^2=x^3+ax^2+bx -> ell = [0,a,0,b,0]
\\ avec a et b entiers.
local(P,Pfact,tors,pointstriv,apinit,bpinit,plist,KS2prod,oddclass,KS2gen,listpoints,pointgen,n1,n2,certain,np1,np2,listpoints2,aux1,aux2,certainp,rang,strange);

if( DEBUGLEVEL_ell >= 2, print("  Algorithm of 2-descent via isogenies"));
if( DEBUGLEVEL_ell >= 3, print("   starting bnfell2descent_viaisog"));
  if( variable(bnf.pol) != 'y,
    error("bnfell2descent_viaisog: the variable of the number field must be y"));
  ell = ellinit(Mod(lift(ell),bnf.pol));

  if( ell.disc == 0,
    error("bnfell2descent_viaisog: singular curve !!"));
  if( ell.a1 != 0 || ell.a3 != 0 || ell.a6 != 0,
    error("bnfell2descent_viaisog: curve not of the form [0,a,0,b,0]"));
  if( denominator(nfalgtobasis(bnf,ell.a2)) > 1 || denominator(nfalgtobasis(bnf,ell.a4)) > 1,
    error("bnfell2descent_viaisog: non integral coefficients"));

  P = Pol([1,ell.a2,ell.a4])*Mod(1,bnf.pol);
  Pfact = factornf(P,bnf.pol)[,1];
  tors = #Pfact;
  if( #Pfact > 1,
    pointstriv = [[0,0,1],[-polcoeff(Pfact[1],0),0,1],[-polcoeff(Pfact[2],0),0,1]]
  , pointstriv = [[0,0,1]]);

  apinit = -2*ell.a2; bpinit = ell.a2^2-4*ell.a4;

\\ calcul des ideaux premiers de plist
\\ et de quelques renseignements associes
  plist = idealfactor(bnf,6*ell.disc)[,1];

if( DEBUGLEVEL_ell >= 3, print("   Search for trivial points on the curve"));
  P *= 'x;
if( DEBUGLEVEL_ell >= 3, print("   Y^2 = ",P));
  pointstriv = concat( pointstriv, nfratpoint(bnf.nf,P,LIMTRIV,0));
if( DEBUGLEVEL_ell >= 1, print(" trivial points on E(K) = ");
  print(lift(pointstriv)); print());

  KS2prod = ell.a4;
  oddclass = 0;
  while( !oddclass,
    KS2gen = bnfsunit(bnf,idealfactor(bnf,KS2prod)[,1]~);
    oddclass = KS2gen[5][1]%2;
    if( !oddclass,
      KS2prod = idealmul(bnf,KS2prod,(KS2gen[5][3][1])));
 );
  KS2gen = KS2gen[1];
  for( i = 1, #KS2gen,
    KS2gen[i] = nfbasistoalg(bnf, KS2gen[i]));
  KS2gen = concat(Mod(lift(concat(bnf.tu[2], bnf.fu)),bnf.pol),KS2gen);
if( DEBUGLEVEL_ell >= 2,
  print("  #K(b,2)gen          = ",#KS2gen);
  print("  K(b,2)gen = ",KS2gen));

  listpoints = nfellcount(bnf.nf,ell.a2,ell.a4,KS2gen,pointstriv);
  pointgen = listpoints[1];
if( DEBUGLEVEL_ell >= 1,
  print(" points on E(K) = ",lift(pointgen));
  print());
  n1 = listpoints[2][1]; n2 = listpoints[2][2];

  certain = (n1 == n2);
if( DEBUGLEVEL_ell >= 1,
  if( certain,
    print("[E(K):phi'(E'(K))]  = ",1<<n1);
    print("#S^(phi')(E'/K)     = ",1<<n2);
    print("#III(E'/K)[phi']    = 1"); print()
  ,
    print("[E(K):phi'(E'(K))] >= ",1<<n1);
    print("#S^(phi')(E'/K)     = ",1<<n2);
    print("#III(E'/K)[phi']   <= ",1<<(n2-n1)); print())
);

  KS2prod = bpinit;
  oddclass = 0;
  while( !oddclass,
    KS2gen = bnfsunit(bnf,idealfactor(bnf,KS2prod)[,1]~);
    oddclass = (KS2gen[5][1]%2);
    if( !oddclass,
      KS2prod = idealmul(bnf,KS2prod,(KS2gen[5][3][1]))));
  KS2gen = KS2gen[1];
  for( i = 1, #KS2gen,
    KS2gen[i] = nfbasistoalg(bnf, KS2gen[i]));
  KS2gen = concat(Mod(lift(concat(bnf.tu[2], bnf.fu)),bnf.pol),KS2gen);
if( DEBUGLEVEL_ell >= 2,
  print("  #K(a^2-4b,2)gen     = ",#KS2gen);
  print("  K(a^2-4b,2)gen     = ",KS2gen));

  P = Pol([1,apinit,bpinit])*Mod(1,bnf.pol);
  Pfact= factornf(P,bnf.pol)[,1];
  if( #Pfact > 1,
    pointstriv=[[0,0,1],[-polcoeff(Pfact[1],0),0,1],[-polcoeff(Pfact[2],0),0,1]]
  , pointstriv = [[0,0,1]]);

if( DEBUGLEVEL_ell >= 3, print("   Search for trivial points on the curve"));
  P *= 'x;
if( DEBUGLEVEL_ell >= 3, print("   Y^2 = ",P));
  pointstriv = concat( pointstriv, nfratpoint(bnf.nf,P,LIMTRIV,0));
if( DEBUGLEVEL_ell >= 1, print(" trivial points on E'(K) = ");
  print(lift(pointstriv)); print());

  listpoints = nfellcount(bnf.nf,apinit,bpinit,KS2gen,pointstriv);
if( DEBUGLEVEL_ell >= 1, print(" points on E'(K) = ",lift(listpoints[1])));
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
if( DEBUGLEVEL_ell >= 1, print(" points on E(K) = ",lift(listpoints2)); print());
  pointgen = concat(pointgen,listpoints2);

  certainp = (np1 == np2);
if( DEBUGLEVEL_ell >= 1,
  if( certainp,
    print("[E'(K):phi(E(K))]   = ",1<<np1);
    print("#S^(phi)(E/K)       = ",1<<np2);
    print("#III(E/K)[phi]      = 1"); print()
  ,
    print("[E'(K):phi(E(K))]  >= ",1<<np1);
    print("#S^(phi)(E/K)       = ",1<<np2);
    print("#III(E/K)[phi]     <= ",1<<(np2-np1)); print());

  if( !certain && (np2>np1), print1(1<<(np2-np1)," <= "));
  print1("#III(E/K)[2]       ");
  if( certain && certainp, print1(" "), print1("<"));
  print("= ",1<<(n2+np2-n1-np1));

  print("#E(K)[2]            = ",1<<tors);
);
  rang = n1+np1-2;
if( DEBUGLEVEL_ell >= 1,
  if( certain && certainp,
    print("#E(K)/2E(K)         = ",(1<<(rang+tors)));
    print("rank                = ",rang); print()
  ,
    print("#E(K)/2E(K)        >= ",(1<<(rang+tors))); print();
    print(rang," <= rank          <= ",n2+np2-2); print()
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
          print("[E'(K):phi(E(K))]   = ",1<<np1);
          print("#S^(phi)(E/K)       = ",1<<np2);
          print("#III(E/K)[phi]      = 1"); print()
        ,
          print("[E'(K):phi(E(K))]  >= ",1<<np1);
          print("#S^(phi)(E/K)       = ",1<<np2);
          print("#III(E/K)[phi]     <= ",1<<(np2-np1)); print())
      )
    ,
    if( certainp,
      n1++;
      certain = ( n1 == n2);
if( DEBUGLEVEL_ell >= 1,
        if( certain,
          print("[E(K):phi'(E'(K))]   = ",1<<n1);
          print("#S^(phi')(E'/K)       = ",1<<n2);
          print("#III(E'/K)[phi']      = 1"); print()
        ,
          print("[E(K):phi'(E'(K))]  >= ",1<<n1);
          print("#S^(phi')(E'/K)      = ",1<<n2);
          print("#III(E'/K)[phi']    <= ",1<<(n2-n1)); print())
      )
    , n1++)
  );

if( DEBUGLEVEL_ell >= 1,
  if( !certain && (np2>np1), print1(1<<(np2-np1)," <= "));
  print1("#III(E/K)[2]       ");
  if( certain && certainp, print1(" "), print1("<"));
  print("= ",1<<(n2+np2-n1-np1));
  print("#E(K)[2]            = ",1<<tors);
);
  rang = n1+np1-2;
if( DEBUGLEVEL_ell >= 1,
  if( certain && certainp,
    print("#E(K)/2E(K)         = ",(1<<(rang+tors))); print();
    print("rank                = ",rang); print()
  ,
    print("#E(K)/2E(K)        >= ",(1<<(rang+tors))); print();
    print(rang," <= rank          <= ",n2+np2-2); print())
  ));

\\ end of strange

if( DEBUGLEVEL_ell >= 1, print(" points = ",pointgen));
if( DEBUGLEVEL_ell >= 3, print("   end of bnfell2descent_viaisog"));
  return([rang,n2+np2-2+tors,pointgen]);
}
{nfchinese( nf, b, fact) =
\\ Chinese Remainder Theorem
local(l,fact2);

if( DEBUGLEVEL_ell >= 4, print("    starting nfchinese"));
  l = #fact[,1];
  fact2 = vector(l,i,idealdiv(nf,b,idealpow(nf,fact[i,1],fact[i,2])));
  fact2 = idealaddtoone(nf,fact2);
  for( i = 1, l,
    fact2[i] = nfbasistoalg(nf,fact2[i]));
if( DEBUGLEVEL_ell >= 4, print("    end of nfchinese"));
  return(fact2);
}
{bnfqfsolve2(bnf, aleg, bleg, aut=['y]) =
\\ Solves Legendre Equation x^2-aleg*Y^2=bleg*Z^2
\\ Using quadratic norm equations
\\ aut contains the Galois automorphisms of bnf ( as polynomials in y)
\\ with aut[1] = y.
local(aux,solvepolrel,auxsolve,solvepolabs,exprxy,rrrnf,bbbnf,SL0,SL1,SL,sunL,fondsunL,normfondsunL,SK,sunK,fondsunK,vecbleg,matnorm,matnormmod,expsolution,solution,reste,carre,verif,x0,x1);

if( DEBUGLEVEL_ell >= 3, print("   starting bnfqfsolve2"));
  solvepolrel = 'x^2-aleg;
if( DEBUGLEVEL_ell >= 4, print("    aleg = ",aleg));
if( DEBUGLEVEL_ell >= 4, print("    bleg = ",bleg));

  if( nfissquare(bnf,aleg),
if( DEBUGLEVEL_ell >= 4, print("    aleg is a square !!"));
    aleg = nfsqrt(bnf,aleg)[1];
    return([aleg,1,0]~)
  );

  if( #aut > 1,
if( DEBUGLEVEL_ell >= 4, print("    factorization of the discriminant using the automorphisms of bnf"));
    for( i = 2, #aut,
      aux = abs(polresultant(lift(aleg)-subst(lift(aleg),'y,aut[i]),bnf.pol));
      if( aux, addprimes(factor(aux)[,1]))));

  auxsolve = rnfequation(bnf,solvepolrel,1);
  solvepolabs = auxsolve[1];
  exprxy = auxsolve[2];
  if( auxsolve[3],
if( DEBUGLEVEL_ell >= 5, print("     case with auxsolve[3] != 0")));
if( DEBUGLEVEL_ell >= 4, print("    bbbnfinit ",solvepolabs));
  rrrnf = rnfinit(bnf,solvepolrel);
  bbbnf = bnfinit(solvepolabs,1);
if( DEBUGLEVEL_ell >= 4, print("    done"));
  SL0 = 1;
if( DEBUGLEVEL_ell >= 4, print("    bbbnf.clgp = ",bbbnf.clgp));
  for( i = 1, #bbbnf.clgp[2],
    if( bbbnf.clgp[2][i]%2 == 0,
      SL0 = idealmul(bbbnf,SL0,bbbnf.clgp[3][i][1,1])));
  SL1 = idealmul(bbbnf,SL0,rnfeltup(rrrnf,bleg));
  SL = idealfactor(bbbnf,SL1)[,1]~;
  sunL = bnfsunit(bbbnf,SL);
  fondsunL = concat(concat(bbbnf.fu, bbbnf.tu[2]),vector(#sunL[1],i,nfbasistoalg(bbbnf,sunL[1][i])));
  normfondsunL = vector(#fondsunL, i, norm(rnfeltabstorel(rrrnf,fondsunL[i])));
  SK = idealfactor(bnf,idealnorm(bbbnf,SL1))[,1]~;
  sunK = bnfsunit(bnf,SK);
  fondsunK = concat(concat(bnf.fu, bnf.tu[2]),vector(#sunK[1],i,nfbasistoalg(bnf,sunK[1][i])));
  vecbleg = bnfissunit(bnf,sunK,bleg);
  matnorm = matrix(#fondsunK,#normfondsunL,i,j,0);
  for( i = 1, #normfondsunL,
    matnorm[,i] = lift(bnfissunit( bnf,sunK,normfondsunL[i] )));
  matnormmod = matnorm*Mod(1,2);
  expsolution = lift(matinverseimage( matnormmod, vecbleg*Mod(1,2)));
  if( !length(expsolution), error("bnfqfsolve2 : NO SOLUTION !! "));
  solution = prod( i = 1, #expsolution, fondsunL[i]^expsolution[i]);
  solution = rnfeltabstorel(rrrnf,solution);
  reste = (lift(vecbleg) - matnorm*expsolution)/2;
  carre = prod( i = 1, #vecbleg, fondsunK[i]^reste[i]);
  solution *= carre;
  x1 = polcoeff(lift(solution),1,'x);
  x0 = polcoeff(lift(solution),0,'x);
  verif = x0^2 - aleg*x1^2-bleg;
if( verif, error("bnfqfsolve2: WRONG POINT"));
if( DEBUGLEVEL_ell >= 3, print("   end of bnfqfsolve2"));
  return([x0,x1,1]~);
}
{bnfqfsolve(bnf, aleg, bleg, flag3, aut=['y]) =
\\ cette fonction resout l'equation X^2-aleg*Y^2=bleg*Z^2
\\ dans le corps de nombres nf.
\\ la solution est [X,Y,Z],
\\ [0,0,0] sinon.
local(nf,aa,bb,na,nb,maxnb,mat,resl,t,sq,pol,vecrat,alpha,xx,yy,borne,test,sun,fact,suni,f,l,aux,alpha2,maxnbiter,idbb,rem,nbiter,mask,oldnb,newnb,bor,testici,de,xxp,yyp,rap,verif);

if( DEBUGLEVEL_ell >= 5, print("     starting bnfqfsolve"));
if( DEBUGLEVEL_ell >= 3, print("   (a,b) = (",aleg,",",bleg,")"));
  nf = bnf.nf;
  aleg = Mod(lift(aleg),nf.pol); aa = aleg;
  bleg = Mod(lift(bleg),nf.pol); bb = bleg;

  if( aa == 0, 
if( DEBUGLEVEL_ell >= 5, print("     end of bnfqfsolve"));
    return([0,1,0]~));
  if( bb == 0,
if( DEBUGLEVEL_ell >= 5, print("     end of bnfqfsolve"));
    return([0,0,1]~));

  na = abs(norm(aa)); nb = abs(norm(bb));
  if( na > nb, maxnb = na, maxnb = nb);
  maxnb <<= 20;
  mat = Mod(matid(3),nf.pol); borne = 1;
  test = 0; nbiter = 0;

  while( 1,
    if( flag3 && bnf.clgp[1]>1, resl = bnfqfsolve2(bnf,aa,bb,aut); break);
if( DEBUGLEVEL_ell >= 4, print("    (na,nb,a,b) = ",lift([na,nb,aa,bb,norm(aa),norm(bb)])));
if( DEBUGLEVEL_ell >= 5, print("     ***",nb,"*** "));
    if( nb >= maxnb, 
      mat = Mod(matid(3),nf.pol);
      aa = aleg; bb = bleg; na = abs(norm(aleg)); nb = abs(norm(bleg)));
    if( aa == 1, resl = [1,1,0]~; break);
    if( bb == 1, resl = [1,0,1]~; break);
    if( aa+bb == 1, resl = [1,1,1]~; break);
    if( aa+bb == 0, resl = [0,1,1]~; break);
    if( aa == bb  && aa != 1,
      t = aa*mat[,1];
      mat[,1] = mat[,3]; mat[,3] = t;
      aa = -1; na = 1);
    if( issquare(na),
      sq = nfsqrt(nf,aa);
      if( sq != [], resl = [sq[1],1,0]~; break));
    if( issquare(nb),
      sq = nfsqrt(nf,bb);
      if( sq != [], resl = [sq[1],0,1]~; break));
    if( na > nb,
      t = aa; aa = bb; bb = t;
      t = na; na = nb; nb = t;
      t = mat[,3]; mat[,3] = mat[,2]; mat[,2] = t);
    if( nb == 1,
if( DEBUGLEVEL_ell >= 4, print("    (a,b) = ",lift([aa,bb])));
if( DEBUGLEVEL_ell >= 4, print("    (na,nb) = ",lift([na,nb])));
      if( aleg == aa && bleg == bb, mat = Mod(matid(3),nf.pol));
      if( flag3, resl = bnfqfsolve2(bnf,aa,bb,aut); break);
      pol = aa*'x^2+bb;  
      vecrat = nfratpoint(nf,pol,borne++,1);
      if( vecrat != 0, resl=[vecrat[2],vecrat[1],vecrat[3]]~; break);

      alpha = 0;
if( DEBUGLEVEL_ell >= 4, print("    bound = ",borne));
      while( alpha==0,
        xx = nfrandint(nf,borne); yy = nfrandint(nf,borne);
        borne++; 
        alpha = xx^2-aa*yy^2 );
      bb *= alpha; nb *= abs(norm(alpha));
      t = xx*mat[,1]+yy*mat[,2];
      mat[,2] = xx*mat[,2]+aa*yy*mat[,1];
      mat[,1] = t;
      mat[,3] *= alpha
    ,
      test = 1;
if( DEBUGLEVEL_ell >= 4, print("on factorise bb = ",bb));
      sun = bnfsunit(bnf,idealfactor(bnf,bb)[,1]~);
      fact = lift(bnfissunit(bnf,sun,bb));
if( DEBUGLEVEL_ell >= 4, print("fact = ",fact));
      suni = concat(concat(bnf.fu, bnf.tu[2]),vector(#sun[1],i,nfbasistoalg(bnf,sun[1][i])));
      for( i = 1, #suni,
        if( (f = fact[i]>>1), 
          test =0;
          for( k = 1, 3, mat[k,3] /= suni[i]^f);
          nb /= abs(norm(suni[i]))^(2*f);
          bb /= suni[i]^(2*f)));
if( DEBUGLEVEL_ell >= 4, print("    factorization of bb = ",bb));
      fact = idealfactor(nf,bb);
if( DEBUGLEVEL_ell >= 4, print("    fact = ",fact));
      l = #fact[,1];

      if( test,
        aux = 1;
        for( i = 1, l,
          if( (f = fact[i,2]>>1) &&
               !(fact[i,1][1]%2) && !nfissquaremodpodd(nf,aa,fact[i,1]),
            aux=idealmul(nf,aux,idealpow(nf,fact[i,1],f))));
        if( aux != 1,
          test = 0;
          alpha = nfbasistoalg(nf,idealappr(nf,idealinv(nf,aux)));
          alpha2 = alpha^2;
          bb *= alpha2; nb *= abs(norm(alpha2));
          mat[,3] *= alpha));
      if( test,
        maxnbiter = 1<<l;
        sq = vector(l,i,nfsqrtmodpq(nf,aa,fact[i,1],fact[i,2]));
        l = #sq;
if( DEBUGLEVEL_ell >= 4,
  print("    sq = ",sq);
  print("    fact = ",fact);
  print("    l = ",l));
        if( l > 1, 
          idbb = idealhnf(nf,bb);
          rem = nfchinese(nf,idbb,fact));
        test = 1; nbiter = 1;
        while( test && nbiter <= maxnbiter,
          if( l > 1,
            mask = nbiter; xx = 0;
            for( i = 1, l,
              if( mask%2, xx += rem[i]*sq[i], xx -= rem[i]*sq[i] ); mask >>= 1)
          , 
            test = 0; xx = sq[1]);
          xx = mynfeltmod(nf,xx,bb);
          alpha = xx^2-aa;
          if( alpha == 0, resl=[xx,1,0]~; break(2));
          t = alpha/bb;
if( DEBUGLEVEL_ell >= 4, print("    [alpha,bb] = ",[alpha,bb]));
          oldnb = nb;
          newnb = abs(norm(t));
if( DEBUGLEVEL_ell >= 4, print("    [oldnb,newnb,oldnb/newnb] = ",[oldnb,newnb,oldnb/newnb+0.]));
          while( nb > newnb,
            mat[,3] *= t;
            bb = t; nb = newnb;
            t = xx*mat[,1]+mat[,2];
            mat[,2] = aa*mat[,1] + xx*mat[,2];
            mat[,1] = t;
            xx = mynfeltmod(nf,-xx,bb);
            alpha = xx^2-aa;
            t = alpha/bb;
            newnb = abs(norm(t));
          );
        if( nb == oldnb, nbiter++, test = 0);
        );
        if( nb == oldnb,
          if( flag3, resl = bnfqfsolve2(bnf,aa,bb,aut); break);
          pol = aa*'x^2+bb;
          vecrat =nfratpoint(nf,pol,borne++<<1,1);
          if( vecrat != 0, resl = [vecrat[2],vecrat[1],vecrat[3]]~; break);

          bor = 1000; yy = 1; testici = 1;
          for( i = 1, 10000, de = nfbasistoalg(nf,vectorv(poldegree(nf.pol),j,random(bor)));
            if( idealadd(bnf,de,bb) != matid(poldegree(bnf.pol)),next);
            xxp = mynfeltmod(bnf,de*xx,bb); yyp = mynfeltmod(bnf,de*yy,bb);
            rap = (norm(xxp^2-aa*yyp^2)/nb^2+0.);
            if( abs(rap) < 1,
if( DEBUGLEVEL_ell >= 4, print("    ********** \n \n MIRACLE ",rap," \n \n ***"));
              t = (xxp^2-aa*yyp^2)/bb;
              mat[,3] *= t;
              bb = t; nb = abs(norm(bb));
if( DEBUGLEVEL_ell >= 4, print("    newnb = ",nb));
              t = xxp*mat[,1]+yyp*mat[,2];
              mat[,2] = aa*yyp*mat[,1] + xxp*mat[,2];
              mat[,1] = t;
              xx = xxp; yy = -yyp; testici = 0;
              ));

          if( testici,  
            alpha = 0;
            while( alpha == 0,
              xx = nfrandint(nf,4*borne); yy = nfrandint(nf,4*borne);
              borne++;
              alpha = xx^2-aa*yy^2);
            bb *= alpha; nb *= abs(norm(alpha));
            t = xx*mat[,1] + yy*mat[,2];
            mat[,2] = xx*mat[,2]+aa*yy*mat[,1];
            mat[,1] = t;
            mat[,3] *= alpha;)))));
  resl = lift(mat*resl);
if( DEBUGLEVEL_ell >= 5, print("     resl1 = ",resl));
if( DEBUGLEVEL_ell >= 5, print("     content = ",content(resl)));
  resl /= content(resl);
  resl = Mod(lift(resl),nf.pol);
if( DEBUGLEVEL_ell >=5, print("     resl3 = ",resl));
  fact = idealadd(nf,idealadd(nf,resl[1],resl[2]),resl[3]);
  fact = bnfisprincipal(bnf,fact,3);
  resl *=1/nfbasistoalg(nf,fact[2]);
if( DEBUGLEVEL_ell >= 5, print("     resl4 = ",resl));
if( DEBUGLEVEL_ell >= 3, print("   resl = ",resl));
  verif = (resl[1]^2-aleg*resl[2]^2-bleg*resl[3]^2 == 0);
  if( !verif, error("bnfqfsolve: WRONG POINT"));
if( DEBUGLEVEL_ell >= 3, print("   end of bnfqfsolve"));
  return(resl);
}
{bnfredquartique( bnf, pol, r,a,b) =
\\ reduction d'une quartique issue de la 2-descente
\\ en connaissant les valeurs de r, a et b. 
local(gcc,princ,rp,pol2);

if( DEBUGLEVEL_ell >= 4, print("    starting bnfredquartique"));
if( DEBUGLEVEL_ell >= 4, print("    ",[r,a,b]));
if( DEBUGLEVEL_ell >= 3, print("   reduction of the quartic ",pol));

  if( a == 0,
    rp = 0
  ,
    gcc = idealadd(bnf,b,a);
    if( gcc == 1,
      rp = nfbasistoalg(bnf,idealaddtoone(bnf.nf,a,b)[1])/a;
      rp = mynfeltmod(bnf,r*rp,b)
    ,
      princ = bnfisprincipal(bnf,gcc,3);
      if( princ[1] == 0, gcc = nfbasistoalg(bnf,princ[2])
      ,
if( DEBUGLEVEL_ell >= 3, print("   quartic not reduced"));
if( DEBUGLEVEL_ell >= 4, print("    end of bnfredquartique"));
        return([pol,0,1]));
      rp = nfbasistoalg(bnf,idealaddtoone(bnf.nf,a/gcc,b/gcc)[1])/(a/gcc);
      rp = mynfeltmod(bnf,r*rp,b)/gcc;
      b /= gcc;
    )
  );
  pol2 = subst(pol/b,'x,rp+b*'x)/b^3;
if( DEBUGLEVEL_ell >= 3, print(" quartic reduced: ",pol2));
if( DEBUGLEVEL_ell >= 4, print("    end of bnfredquartique"));
  return([pol2,rp,b]);
}
{bnfell2descent_gen( bnf, ell, ext, help=[], bigflag=1, flag3=1, aut=['y]) =
\\ bnf a un polynome en y.
\\ si ell= y^2=P(x), alors ext est
\\ ext[1] est une equation relative du corps (=P(x)),
\\ ext[2] est le resultat rnfequation(bnf,P,1);
\\ ext[3] est le bnfinit (sur Q) de l'extension.
\\ dans la suite ext est note L = K(theta).
\\ help est une liste de points deja connus sur ell.
\\ si bigflag !=0 alors on applique bnfredquartique.
\\ si flag3 ==1 alors on utilise bnfqfsolve2 (equation aux normes) pour resoudre Legendre
\\ aut est une liste d'automorphismes connus de bnf
\\ (ca peut aider a factoriser certains discriminants).
\\ ell est de la forme y^2=x^3+A*x^2+B*x+C
\\ ie ell=[0,A,0,B,C], avec A,B et C entiers.
\\
local(nf,unnf,ellnf,A,B,C,S,plist,Lrnf,SLprod,LS2,LS2gen,polrel,alpha,ttheta,KS2gen,normLS2,normcoord,LS2coordtilda,LS2tilda,aux,listgen,listpointstriv,listpoints,m1,m2,loc,lastloc,maskwhile,iwhile,zc,iaux,liftzc,ispointtriv,point,c,b,a,found,alphac,r,denc,dena,cp,alphacp,beta,mattr,vec,z1,cont,d,e,polorig,pol,redq,transl,multip,UVW,pointxx,point2,rang,listELS,listnotELS);

if( DEBUGLEVEL_ell >= 4, print("    starting bnfell2descent_gen"));

  nf = bnf.nf;
  unnf = Mod(1,nf.pol);
  ellnf = ell*unnf;
  if( #ellnf <= 5, ellnf = ellinit(ellnf));

  A = ellnf.a2; if( DEBUGLEVEL_ell >= 2, print("  A = ",A));
  B = ellnf.a4; if( DEBUGLEVEL_ell >= 2, print("  B = ",B));
  C = ellnf.a6; if( DEBUGLEVEL_ell >= 2, print("  C = ",C));

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\      Construction of L(S,2)         \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

if( DEBUGLEVEL_ell >= 2, print(); print("  Computing L(S,2)"));

  polrel = ext[1];
  alpha = Mod(Mod('y,nf.pol),polrel); \\ alpha est l'element primitif de K
  ttheta = Mod('x,polrel);            \\ ttheta est la racine de P(x)

  S = 6*lift(ellnf.disc);
  plist = idealfactor(nf,S)[,1];
  Lrnf = ext[3];
  SLprod = subst(lift(polrel'),'y,lift(ext[2][2]));
if( DEBUGLEVEL_ell >= 3, print("   ",ext[2]));

  while( 1,
\\ Construction des S-unites
    LS2gen = bnfsunit(Lrnf, idealfactor(Lrnf,SLprod)[,1]~);
if( DEBUGLEVEL_ell >= 4, print("    LS2gen = ",LS2gen));
\\ si le groupe de classes est impair, on a fini.
    if( LS2gen[5][1]%2, break);
if( DEBUGLEVEL_ell >= 3, print("   2-class group ",LS2gen[5][3][1][1,1]));
    S *= LS2gen[5][3][1][1,1];
    SLprod = idealmul(Lrnf,SLprod,(LS2gen[5][3][1]));
  );

  KS2gen = bnfsunit(bnf,idealfactor(nf,S)[,1]~);
 
if( DEBUGLEVEL_ell >= 3, print("   #KS2gen = ",#KS2gen[1]));
if( DEBUGLEVEL_ell >= 3, print("    KS2gen = ",KS2gen[1]));

  LS2gen = LS2gen[1];
  LS2 = vector(#LS2gen,i,lift(nfbasistoalg(Lrnf,LS2gen[i])));
  LS2 = concat(lift(concat(Lrnf.fu, Lrnf.tu[2])),LS2);

  LS2 = subst(LS2,'x,ttheta);
  LS2 = LS2*Mod(1,polrel);
if( DEBUGLEVEL_ell >= 3, print("   #LS2 = ",#LS2));
if( DEBUGLEVEL_ell >= 3, print("    LS2 = ",LS2));

if( DEBUGLEVEL_ell >= 2, print("  L(S,2) = ",LS2));

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\   Construction of the Selmer group \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

if( DEBUGLEVEL_ell >= 2, print(); print("  Computing the Selmer group"));

\\ dans LS2gen, on ne garde que ceux dont la norme est un carre.

  normLS2 = norm(LS2);
if( DEBUGLEVEL_ell >= 4, print("    normLS2 = ",normLS2));

\\ matrice de l'application norme

  normcoord = matrix(#KS2gen[1]+#bnf[8][5]+1,#normLS2,i,j,0);
  for( i = 1, #normLS2,
    normcoord[,i] = bnfissunit(bnf,KS2gen,normLS2[i]));
if( DEBUGLEVEL_ell >= 4, print("    normcoord = ",normcoord));

\\ construction du noyau de la norme

  LS2coordtilda = lift(matker(normcoord*Mod(1,2)));
if( DEBUGLEVEL_ell >= 4, print("    LS2coordtilda = ",LS2coordtilda));
  LS2tilda = vector(#LS2coordtilda[1,],i,0);
  for( i = 1, #LS2coordtilda[1,],
    aux = 1;
    for( j = 1, #LS2coordtilda[,i],
      if( sign(LS2coordtilda[j,i]),
        aux *= LS2[j]));
    LS2tilda[i] = aux;
  );

if( DEBUGLEVEL_ell >= 3, print("   LS2tilda = ",LS2tilda));
if( DEBUGLEVEL_ell >= 3, print("   norm(LS2tilda) = ",norm(LS2tilda)));

\\ Fin de la construction de L(S,2)

  listgen = LS2tilda;
if( DEBUGLEVEL_ell >= 2, print("  #LS2gen = ",#listgen));
if( DEBUGLEVEL_ell >= 2, print("   LS2gen = ",listgen));

if( DEBUGLEVEL_ell >= 3, print("   (A,B,C) = ",[A,B,C]));

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\   Recherche de points triviaux   \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

if( DEBUGLEVEL_ell >= 2, print("  Search for trivial points on the curve"));
  listpointstriv = nfratpoint(nf,'x^3+A*'x^2+B*'x+C,LIMTRIV,0);
  listpointstriv = concat(help,listpointstriv);
if( DEBUGLEVEL_ell >= 1, print(" Trivial points on the curve = ",listpointstriv));

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\          parcours de L(S,2)         \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  listpoints = [];
  m1 = 0; m2 = 0; lastloc = -1;
  maskwhile = 1<<#listgen;
  listELS = [0]; listnotELS = [];
  iwhile = 1;

  while( iwhile < maskwhile,
if( DEBUGLEVEL_ell >= 4,
  print("    iwhile = ",iwhile);
  print("    listgen = ",listgen)
);

\\ utilise la structure de groupe pour detecter une eventuelle solubilite locale.
    loc = 0;
    for( i = 1, #listELS,
      for( j = 1, #listnotELS,
        if( bitxor(listELS[i],listnotELS[j]) == iwhile,
if( DEBUGLEVEL_ell >= 3, print("   Not ELS from group structure"));
          listnotELS = concat(listnotELS,[iwhile]);
          iwhile++;
          next(3)
    )));

    for( i = 1, #listELS,
      for( j = i+1, #listELS,
        if( bitxor(listELS[i],listELS[j]) == iwhile,
if( DEBUGLEVEL_ell >= 3, print("   ELS from group structure"));
          listELS = concat(listELS,[iwhile]);
          loc = 1;
          break(2)
    )));

    iaux = vectorv(#listgen,i,bittest(iwhile,i-1));
    iaux = (LS2coordtilda*iaux)%2;
    zc = unnf*prod( i = 1, #LS2, LS2[i]^iaux[i]);

if( DEBUGLEVEL_ell >= 2, print("  zc = ",zc));
    liftzc = lift(zc);

\\ Est-ce un point trivial ?
    found = 0;
    ispointtriv = 0;
    for( i = 1, #listpointstriv,
      point = listpointstriv[i];
      if( #point == 2 || point[3] != 0,
        if( nfissquare(Lrnf.nf,subst((lift(point[1])-'x)*lift(liftzc),'y,lift(ext[2][2]))),
if( DEBUGLEVEL_ell >= 2, print("  comes from the trivial point ",point));
          listpoints = concat(listpoints,[point]);
          found = 1; ispointtriv = 1; break
    )));

\\ \\\\\\\\\\\\\
\\ On cherche a ecrire zc sous la forme a-b*theta
\\ \\\\\\\\\\\\\

    if( !found, \\ si ce n'est pas un point trivial
      a = polcoeff(liftzc,0);
      b =-polcoeff(liftzc,1);
      c = polcoeff(liftzc,2);

      if( c != 0,
        alphac = (A*b+B*c-a)*c+b^2;
if( DEBUGLEVEL_ell >= 3, print("   alphac = ",alphac));
        r = nfsqrt(nf,norm(zc))[1];
        if( alphac == 0,
\\ cas particulier
if( DEBUGLEVEL_ell >= 3, print("   continuing with 1/zc"));
          zc = norm(zc)*(1/zc);
if( DEBUGLEVEL_ell >= 2, print("  zc = ",zc))
        ,
\\ Il faut resoudre une forme quadratique
\\ Existence (locale = globale) d'une solution :
          denc = deno(lift(c));
          if( denc != 1, cp = c*denc^2, cp = c);
          dena = deno(lift(alphac));
          if( dena != 1, alphacp = alphac*dena^2, alphacp = alphac);
if( DEBUGLEVEL_ell >= 2, print("  Hilbert symbol (",alphacp,",",cp,") = "));
          if( !loc && mynfhilbert(nf, alphacp,cp) < 0,
if( DEBUGLEVEL_ell >= 3, print("   no local solution"));
            listnotELS = concat(listnotELS,[iwhile]);
            iwhile++;
            next
          );
\\ Ici on a l'existence locale
          beta = A*(A*b*c+B*c^2+b^2)-C*c^2+a*b;
          mattr = matid(3);
          mattr[1,1] = c ;mattr[2,2] = alphac ;
          mattr[3,3] = r ;mattr[2,3] = -beta;
          mattr[1,2] = -(b +A*c) ;mattr[1,3] = a-B*c+A*(A*c+b);
if( DEBUGLEVEL_ell >= 2, print1("  sol of quadratic equation = "));
          vec = bnfqfsolve(bnf,alphacp,cp,flag3,aut);
if( DEBUGLEVEL_ell >= 2, print(lift(vec)));
          aux = vec[2]*dena;
          vec[2] = vec[1];vec[1] = aux;
          vec[3] = vec[3]*denc;
          vec = (mattr^(-1))*vec;
          vec /= content(lift(vec));
          z1 = (vec[3]*ttheta+vec[2])*ttheta+vec[1];
if( DEBUGLEVEL_ell >= 3, print("   z1 = ",z1));
          zc *= z1^2;
if( DEBUGLEVEL_ell >= 2, print("  zc*z1^2 = ",zc));
         )
      )
    );

\\ \\\\\\\\\\
\\ Maintenant zc est de la forme a-b*theta
\\ \\\\\\\\\\

    if( !found,

if( DEBUGLEVEL_ell >= 3, print("   zc = ",zc));
      liftzc = lift(zc);
      a = polcoeff(liftzc,0);
      b =-polcoeff(liftzc,1);
      c = polcoeff(liftzc,2);
if( c, error("bnfell2descent_gen : c <> 0"));

\\ remove denominators and try to simplify zc
      cont = idealadd(nf,a,b);
      if( cont != 1,
        cont = idealfactor(nf,cont);
        cont[,2] \= 2;
        aux = 1;
        for( i = 1, #cont[,1],
\\ *****************************************************************
\\          aux = idealmul(nf,aux,idealpow(nf,cont[i,1],cont[i,2]))
          aux = idealmul(nf,aux,idealpow(nf,idealhnf(nf,cont[i,1]),cont[i,2]))
\\ *****************************************************************
        );
        cont = nfbasistoalg(nf,bnfisprincipal(bnf,aux)[2])^2;
        a /= cont;
        b /= cont;
        zc/= cont;
        liftzc = lift(zc);
if( DEBUGLEVEL_ell >= 3, print("   new zc = ",zc));
      );

      if( nfissquare(nf,b),
if( DEBUGLEVEL_ell >= 3, print("   b is a square"));
        point = [a/b,nfsqrt(nf,(a/b)^3+A*(a/b)^2+B*(a/b)+C)[1]];
if( DEBUGLEVEL_ell >= 2, print("  point found = ",point));
        listpoints = concat(listpoints,[point]);
        found = 1; ispointtriv = 1
      )
    );

\\ \\\\\\\\\\\
\\ Construction de la quartique 
\\ \\\\\\\\\\\

    if( !found, \\ si ce n'est pas un point trivial
      r = nfsqrt(nf,norm(zc))[1];
if( DEBUGLEVEL_ell >= 4, print("    r = ",r));
      c = -2*(A*b+3*a);
if( DEBUGLEVEL_ell >= 4, print("    c = ",c));
      d = 8*r;
if( DEBUGLEVEL_ell >= 4, print("    d = ",d));
      e = (A^2*b^2 - 2*A*a*b-4*B*b^2-3*a^2);
if( DEBUGLEVEL_ell >= 4, print("    e = ",e));
      polorig = b*('x^4+c*'x^2+d*'x+e)*unnf;
if( DEBUGLEVEL_ell >= 2, print("  quartic: (",lift(b),")*Y^2 = ",lift(polorig/b)));
      pol = polorig;
      if( bigflag,
        redq = bnfredquartique(bnf,pol,r,a,b);
if( DEBUGLEVEL_ell >= 2, print("  reduced: Y^2 = ",lift(redq[1])));
        pol = redq[1]; transl = redq[2]; multip = redq[3]
      );

\\ Search for a point on the quartic
      point = nfratpoint(nf,pol,LIM1,1);
      found = point != [];
      if( found,
        loc = 1
      );
\\ If the quartic is not known to be ELS, check if it is
      if( !loc,
        if( bigflag,
          loc = nflocallysoluble(nf,pol,r,a,b)
        , loc = nflocallysoluble(nf,pol,0,1,1)
      ));
      if( !loc,
        listnotELS = concat(listnotELS,[iwhile]);
        iwhile++;
        next
      )
    );

\\ If no point is found, search harder
    if( !found,
if( DEBUGLEVEL_ell >= 2, print("quartic is ELS"));
      point = nfratpoint(nf,pol,LIM3,1);
      found = point != []
    );

    if( found && !ispointtriv,
      if( bigflag,
        point[1] = point[1]*multip+transl;
        point[2] = nfsqrt(nf,subst(polorig,'x,point[1]/point[3]))[1]);
      mattr = matid(3);
      mattr[1,1] = -2*b^2; mattr[1,2] = (A*b+a)*b;
      mattr[1,3] = a^2+(2*B-A^2)*b^2; mattr[2,2] = -b;
      mattr[2,3] = a+A*b; mattr[3,3] =r;
      UVW = [point[1]^2,point[3]^2,point[1]*point[3]]~;
      vec = (mattr^(-1))*UVW;
      z1 = (vec[3]*ttheta+vec[2])*ttheta+vec[1];
      zc *= z1^2;
      zc /= -polcoeff(lift(zc),1);
if( DEBUGLEVEL_ell >= 3, print("   zc*z1^2 = ",zc));
      pointxx = polcoeff(lift(zc),0);
      point2 = [ pointxx, nfsqrt(nf,subst('x^3+A*'x^2+B*'x+C,'x,pointxx))[1]];
if( DEBUGLEVEL_ell >= 1, print("  point found = ",point2));
      listpoints = concat(listpoints,[point2]);
    );

    listELS = concat(listELS,[iwhile]);
    if( degre(iwhile) > lastloc, m2++; lastloc = degre(iwhile));

    if( found,
      m1++;
      found = 0; lastloc = -1;
      iwhile = 1<<degre(iwhile);
      maskwhile >>= 1;
      LS2coordtilda = vecextract(LS2coordtilda,1<<#listgen-iwhile-1);
      listgen = vecextract(listgen,1<<#listgen-iwhile-1);
      while( listELS[#listELS] >= iwhile,
        listELS = vecextract(listELS,1<<(#listELS-1)-1));
      while( #listnotELS && listnotELS[#listnotELS] >= iwhile,
        listnotELS = vecextract(listnotELS,1<<(#listnotELS-1)-1))
    , iwhile ++
    )
  );

if( DEBUGLEVEL_ell >= 2,
  print("  m1 = ",m1);
  print("  m2 = ",m2));
if( DEBUGLEVEL_ell >= 1,
  print("#S(E/K)[2]    = ",1<<m2));
  if( m1 == m2,
if( DEBUGLEVEL_ell >= 1,
    print("#E(K)/2E(K)   = ",1<<m1);
    print("#III(E/K)[2]  = 1");
    print("rank(E/K)     = ",m1));
    rang = m1
  ,
if( DEBUGLEVEL_ell >= 1,
    print("#E(K)/2E(K)  >= ",1<<m1);
    print("#III(E/K)[2] <= ",1<<(m2-m1));
    print("rank(E/K)    >= ",m1));
    rang = m1;
    if( (m2-m1)%2,
if( DEBUGLEVEL_ell >= 1,
      print(" III should be a square, hence ");
      if( m2-m1 > 1,
        print("#E(K)/2E(K)  >= ",1<<(m1+1));
        print("#III(E/K)[2] <= ",1<<(m2-m1-1));
        print("rank(E/K)    >= ",m1+1)
      ,
        print("#E(K)/2E(K)  = ",1<<(m1+1));
        print("#III(E/K)[2] = 1");
        print("rank(E/K)    = ",m1+1)));
      rang = m1+1)
  );

if( DEBUGLEVEL_ell >= 1, print(" listpoints = ",listpoints));
  for( i = 1, #listpoints,
    if( #listpoints[i] == 3, 
      listpoints[i] = vecextract(listpoints[i],3));
    if( !ellisoncurve(ellnf,listpoints[i]), 
      error("bnfell2descent: WRONG POINT ")));
if( DEBUGLEVEL_ell >= 4, print("    end of bnfell2descent_gen"));
  return([rang,m2,listpoints]);
}
{bnfellrank(bnf,ell,help=[],bigflag=1,flag3=1) =
\\ Algorithme de la 2-descente sur la courbe elliptique ell.
\\ help est une liste de points connus sur ell.

\\ attention bnf a un polynome en y.
\\ si bigflag !=0, on reduit les quartiques
\\ si flag3 != 0, on utilise bnfqfsolve2
local(urst,urst1,den,factden,eqtheta,rnfeq,bbnf,ext,rang,f);

if( DEBUGLEVEL_ell >= 3, print("   starting bnfellrank"));
  if( #ell < 5, ell = ellinit(ell));
  ell = vector(5, i, ell[i]);

\\ removes the coefficients a1 and a3
  urst = [1,0,0,0];
  if( ell.a1 != 0 || ell.a3 != 0,
    urst1 = [1,0,-ell.a1/2,-ell.a3/2];
    ell = ellchangecurve(ell,urst1);
    urst = ellcomposeurst(urst,urst1)
  );

\\ removes denominators
  while( (den = idealinv(bnf,idealadd(bnf,idealadd(bnf,1,ell.a2),idealadd(bnf,ell.a4,ell.a6))))[1,1] > 1,
    factden = idealfactor(bnf,den)[,1];
    den = 1;
    for( i = 1, #factden,
      den = idealmul(bnf,den,factden[i]));
    den = den[1,1];
    urst1 = [1/den,0,0,0];
    ell = ellchangecurve(ell,urst1);
    urst = ellcomposeurst(urst,urst1);
  );

  help = ellchangepoint(help,urst);

\\ choix de l'algorithme suivant la 2-torsion
  ell *= Mod(1,bnf.pol);
  eqtheta = Pol([1,ell.a2,ell.a4,ell.a6]);
if( DEBUGLEVEL_ell >= 1, print(" elliptic curve: Y^2 = ",eqtheta));
  f = nfpolratroots(bnf,eqtheta);

  if( #f == 0,                                  \\ cas 1: 2-torsion triviale
    rnfeq = rnfequation(bnf,eqtheta,1);
    urst1 = [1,-rnfeq[3]*Mod('y,bnf.pol),0,0];
    if( rnfeq[3] != 0, 
      ell = ellchangecurve(ell,urst1);
      urst = ellcomposeurst(urst,urst1);
      eqtheta = subst(eqtheta,'x,'x-rnfeq[3]*Mod('y,bnf.pol));
      rnfeq = rnfequation(bnf,eqtheta,1);
if( DEBUGLEVEL_ell >= 2, print("  translation: working with Y^2 = ",eqtheta));
    );
if( DEBUGLEVEL_ell >= 3, print1("   bbnfinit "));
    bbnf = bnfinit(rnfeq[1],1);
if( DEBUGLEVEL_ell >= 3, print("   done"));
    ext = [eqtheta, rnfeq, bbnf];
    rang = bnfell2descent_gen(bnf,ell,ext,help,bigflag,flag3)
  ,
  if( #f == 1,                                  \\ cas 2: 2-torsion = Z/2Z
    if( f[1] != 0, 
      urst1 = [1,f[1],0,0];
      ell = ellchangecurve(ell,urst1);
      urst = ellcomposeurst(urst,urst1)
    );
    rang = bnfell2descent_viaisog(bnf,ell)
  ,                                             \\ cas 3: 2-torsion = Z/2Z*Z/2Z
    rang = bnfell2descent_complete(bnf,f[1],f[2],f[3],flag3)
  ));
 
  rang[3] = ellchangepoint(rang[3],ellinverturst(urst));
if( DEBUGLEVEL_ell >= 3, print("   end of bnfellrank"));

return(rang);
}
{bnfell2descent_complete(bnf,e1,e2,e3,flag3=1,aut=['y]) =
\\ calcul du rang d'une courbe elliptique
\\ par la methode de 2-descente complete.
\\ Y^2 = (x-e1)*(x-e2)*(x-e3);
\\ en suivant la methode decrite par J.Silverman
\\ si flag3 ==1 alors on utilise bnfqfsolve2 (equation aux normes)
\\ pour resoudre Legendre 
\\ on pourra alors utiliser aut=nfgaloisconj(bnf.pol)

\\ e1, e2 et e3 sont des entiers algebriques de bnf.

local(KS2prod,oddclass,KS2gen,vect,selmer,rang,b1,b2,vec,z1,z2,d31,quart0,quart,cont,fa,point,solx,soly,listepoints,strange);

if( DEBUGLEVEL_ell >= 2, print("  Algorithm of complete 2-descent"));

\\ calcul de K(S,2)

  KS2prod = (e1-e2)*(e2-e3)*(e3-e1)*2;
  oddclass = 0;
  while( !oddclass,
    KS2gen = bnfsunit(bnf,idealfactor(bnf,KS2prod)[,1]~);
    oddclass = (KS2gen[5][1]%2);
    if( !oddclass,
      KS2prod = idealmul(bnf,KS2prod,(KS2gen[5][3][1])));
  );
  KS2gen = KS2gen[1];
  for( i = 1, #KS2gen,
    KS2gen[i] = nfbasistoalg(bnf, KS2gen[i]));
  KS2gen = concat(Mod(lift(concat(bnf.tu[2], bnf.fu)),bnf.pol),KS2gen);
if( DEBUGLEVEL_ell >= 2,
  print("  #K(S,2)gen = ",#KS2gen);
  print("   K(S,2)gen = ",KS2gen)
);

\\ parcours de K(S,2)*K(S,2)

  vect = vector(#KS2gen,i,[0,1]);
  selmer = 0;
  rang = 0;
  listepoints = [];

  forvec( X = vect,
    b1 = prod( i = 1, #KS2gen, KS2gen[i]^X[i]);
    forvec( Y = vect,
      b2 = prod( i = 1, #KS2gen, KS2gen[i]^Y[i]);

if( DEBUGLEVEL_ell >= 3, print("   [b1,b2] = ",lift([b1,b2])));

\\ points triviaux provenant de la 2-torsion

      if( b1==1 && b2==1, 
if( DEBUGLEVEL_ell >= 2, print("  trivial point [0]"));
        selmer++; rang++; next);
      if( nfissquare(bnf.nf,(e2-e1)*b1)
        && nfissquare(bnf.nf,(e2-e3)*(e2-e1)*b2),
if( DEBUGLEVEL_ell >= 2, print("  trivial point [e2,0]"));
        selmer++; rang++; next);
      if( nfissquare(bnf.nf,(e1-e2)*b2)
        && nfissquare(bnf.nf,(e1-e3)*(e1-e2)*b1),
if( DEBUGLEVEL_ell >= 2, print("  trivial point [e1,0]"));
        selmer++; rang++; next);
      if( nfissquare(bnf.nf,(e3-e1)*b1)
        && nfissquare(bnf.nf,(e3-e2)*b2),
if( DEBUGLEVEL_ell >= 2, print("  trivial point [e3,0]"));
        selmer++; rang++; next);

\\ premier critere local : sur les formes quadratiques

      if( mynfhilbert(bnf.nf,b1*b2,b1*(e2-e1)) < 0
        || mynfhilbert(bnf.nf,b2,b1*(e3-e1)) < 0
        || mynfhilbert(bnf.nf,b1,b2*(e3-e2)) < 0
        ,
if( DEBUGLEVEL_ell >= 3, print("   not ELS"));
        next);

if( DEBUGLEVEL_ell >= 2, print("  [b1,b2] = ",lift([b1,b2])));
if( DEBUGLEVEL_ell >= 2, print("  quadratic forms locally soluble"));

\\ solution de la premiere forme quadratique

      if( b1 != b2,
        vec = bnfqfsolve(bnf,b1*b2,b1*(e2-e1),flag3);
if( DEBUGLEVEL_ell >= 3, print("   sol part = ",vec));
        if( vec[3] == 0, error("bnfell2descent_complete: BUG !!! : vec[3]=0 "));
        z1 = vec[1]/vec[3]/b1; z2 = vec[2]/vec[3]
      ,
        z1 = (1+(e2-e1)/b1)/2; z2 = z1-1
      );
      d31 = e3-e1;
      quart0 = b2^2*(z1^2*b1 - d31)*'x^4 - 4*z1*b2^2*z2*b1*'x^3
           + 2*b1*b2*(z1^2*b1 + 2*b2*z2^2 + d31)*'x^2 - 4*z1*b2*z2*b1^2*'x
           + b1^2*(z1^2*b1 - d31);
      quart = quart0*b1*b2;
if( DEBUGLEVEL_ell >= 4, print("    quart = ",quart));
      quart *= denominator(simplify(content(quart)))^2;
      cont = simplify(content(lift(quart)));
      fa = factor(cont);
      for( i = 1, #fa[,1], quart /= fa[i,1]^(2*(fa[i,2]\2)));
if( DEBUGLEVEL_ell >= 3, print("   quartic reduced = ",quart));

\\ la quartique est-elle localement soluble ?
   
      if( !nflocallysoluble(bnf.nf,quart),
if( DEBUGLEVEL_ell >= 2, print("  quartic not ELS"));
        next);
      selmer++;

\\ recherche de points sur la quartique.

      point = nfratpoint(bnf.nf,quart,LIM3,1);
      if( point != [], 
if( DEBUGLEVEL_ell >= 2, print("  point found on the quartic !!"));
if( DEBUGLEVEL_ell >= 3, print("   ",point));
        if( point[3], 
          point /= point[3];
          z1 = (2*b2*point[1]*z2-z1*(b1+b2*point[1]^2))/(b1-b2*point[1]^2);
          solx = b1*z1^2+e1;
          soly = nfsqrt(bnf.nf,(solx-e1)*(solx-e2)*(solx-e3))[1];
          listepoints = concat(listepoints,[[solx,soly]]);
if( DEBUGLEVEL_ell >= 1, print(" point on the elliptic curve =",[solx,soly]))
        );
        rang++
      ,
if( DEBUGLEVEL_ell >= 2, print("  no point found on the quartic"))
      )
    )
  );

\\ fin

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
  if( rang, print("points = ",listepoints));
);

return([rang,selmer,listepoints]);
}

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\            HELP MESSAGES                \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
  addhelp(bnfellrank,
    "bnfellrank(K,E,help=[]): E is a 5-component vector defining an elliptic curve defined over the number field K (output by bnfinit()). Returns a vector [r,s,v], where r is a lower bound for the rank of E(K), s is the rank of its 2-Selmer group and v is a list of independant points in E(K)/2E(K). If help is a vector of nontrivial points on E(K), the result might be faster. See also ?default_ell");
\\                  others
  addhelp(default_ell,
    "default_ell(DEBUGLEVEL_ell, LIM1, LIM3, LIMTRIV, MAXPROB, LIMBIGPRIME): output/set the value of the global variables used for ellrank() and other related functions. DEBUGLEVEL_ell: 0-5 : choose the quantity of information printed during the computation (default=0: print nothing); LIM1 (resp LIM3): search limit for easy (resp hard) points on quartics; LIMTRIV: search limit for trivial points on elliptic curves; MAXPROB, LIMBIGPRIME: technical.");
}
