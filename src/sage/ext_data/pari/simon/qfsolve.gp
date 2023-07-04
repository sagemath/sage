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
  Author:
  Denis SIMON -> simon@math.unicaen.fr
  address of the file:
  www.math.unicaen.fr/~simon/qfsolve.gp

  *********************************************
  *          VERSION 09/01/2014               *
  *********************************************

  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  \\                English help                \\
  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  This package provides functions to solve quadratic equations over Q.
  language: GP
  It can be run under GP by the command
  gp > \r qfsolve.gp

  This package contains 4 main functions:

  - qfsolve(G,{factD}): Solve over Q the quadratic equation X~*G*X = 0.
  G must be a symmetric matrix n*n, with coefficients in Z.
  The solution might be a single vector (vectorv)
  or a matrix (whose columns generate a totally isotropic subspace).
  If no solution exists, the output is a prime number
  indicating that there is no solution in the local field Q_p
  (-1 for the reals, p for Q_p).
  factD is an optional parameter. If present, it must be equal to
  the factorization of -abs(2*matdet(G)). This saves a lot of time.

  Example:
  gp > G = [1,0,0;0,1,0;0,0,-34];
  gp > qfsolve(G)
  %1 = [-3, -5, 1]~

  - qfparam(G,sol,fl): Coefficients of quadratic forms that parametrize the
  solutions of the ternary quadratic form G, using the particular
  solution sol.
  fl is optional and can be 1, 2, or 3, in which case the 'fl'th form is
  reduced. The default is fl=3.

  Example:
  gp > qfparam(G,[-3,-5,1]~)
  %2 = 
  [ 3 -10 -3]

  [-5  -6  5]

  [ 1   0  1]
  Indeed, the solutions can be parametrized as
  [3*x^2 - 10*y*x - 3*y^2, -5*x^2 - 6*y*x + 5*y^2, x^2 + y^2]~
  
  - qflllgram_indef(G,c): Solve or reduce the quadratic form G with
  integral coefficients. G might be definite or indefinite.
  This is an lll-type algorithm with a constant 1/4<c<=1.
  c is optional and the default is c=1.

  - quadclass2(d,factd): Compute the 2-Sylow of the (narrow) class group
  of discriminant d. d must be a fondamental discriminant.
  factD is an optional parameter. If present, it must be equal to
  the factorization of abs(2*d). In this case, the 
  algorithm runs in polynomial time.

  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  \\         Description des fonctions          \\
  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  Programme de resolution des equations quadratiques
  langage: GP
  pour l'utiliser, lancer gp, puis taper
  \r qfsolve.gp

  Ce fichier contient 4 principales fonctions:

  - qfsolve(G,{factD}): pour resoudre l'equation quadratique X^t*G*X = 0
  G doit etre une matrice symetrique n*n, a coefficients dans Z.
  S'il n'existe pas de solution, la reponse est un entier
  indiquant un corps local dans lequel aucune solution n'existe
  (-1 pour les reels, p pour Q_p).
  factD est optionnel. Il doit etre egal a -abs(2*matdet(G)),
  Cela entraine un gain de temps.

  Exemple:
  gp > G = [1,0,0;0,1,0;0,0,-34];
  gp > qfsolve(G)
  %1 = [-3, -5, 1]~

  - qfparam(G,sol,fl): pour parametrer les solutions de la forme
  quadratique ternaire G, en utilisant la solution particuliere sol.
  fl est optionnel et peut prendre les valeurs 1,2 ou 3.
  Il indique que la 'fl'eme forme quadratique est reduite.
  Par defaut, fl=3

  Exemple:
  gp > qfparam(G,[-3,-5,1]~)
  %2 = 
  [ 3 -10 -3]

  [-5  -6  5]

  [ 1   0  1]
  Ici, les solutions sont parametrees par
  [3*x^2 - 10*y*x - 3*y^2, -5*x^2 - 6*y*x + 5*y^2, x^2 + y^2]~

  - qflllgram_indef(G,c): pour resoudre ou reduire la forme quadratique
  G a coefficients entiers. Il s'agit d'un algorithme type LLL, avec la
  constante 1/4<c<=1.
  c est optionnel et par defaut c=1.

  - quadclass2(d,factd): determine le 2-Sylow du (narrow) groupe de classes de
  discriminant d, ou d est un discriminant fondamental.
  factd est optionnel. Il doit etre egal a -abs(2*d),
  et dans ce cas le reste de l'algorithme est polynomial.

*/

\\
\\ Usual global variables
\\ 

global(DEBUGLEVEL_qfsolve):small;

  DEBUGLEVEL_qfsolve = 0;

\\ use the variable DEBUGLEVEL_qfsolve :
\\ From 0 to 5 : choose a higher value to have
\\ more details printed.


\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\          SCRIPT                             \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{default_qfsolve(DEBUGLEVEL_qfsolve_val:small = 0) =

  DEBUGLEVEL_qfsolve = DEBUGLEVEL_qfsolve_val;
  print("  DEBUGLEVEL_qfsolve = ",DEBUGLEVEL_qfsolve);
}

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\          TYPE CONVERSIONS                   \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ THIS FUNCTION WILL BE REPLACED BY matconcat(matdiagonal(v))
\\ IN VERSION >= 2.6.0
\\ build the matrix whose diagonal blocks are listed in the vector v.
{matdiagonalblock(v) =
my(M);
my(lv,lt=0);

  if( type(v) != "t_VEC" && type(v) != "t_COL",
    error("wrong type in matdiagonalblock()"));

  lv = length(v);
  for( i = 1, lv, lt += length(v[i]));
  M = matrix(lt,lt);
  lt = 0;
  for( i = 1, lv,
    my( lvi = length(v[i]));
    for( j = 1, lvi,
      for( k = 1, lvi,
        M[lt+j,lt+k] = v[i][j,k]));
    lt += lvi
  );
return(M);
}

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\          LINEAR ALGEBRA                     \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ Gives a unimodular matrix with the last column equal to v.
\\ redflag = 0 or 1. If redflag = 1, then the n-#v first columns are reduced.
{completebasis(v,redflag=0) =
my(Mv,U,re);
my(n);

  if( type(v) != "t_COL" && type(v) != "t_MAT",
    error("wrong type in completebasis"));

  if( type(v) == "t_COL", Mv = Mat(v), Mv = v);
  n = length(Mv[,1]);
  if( n == length(Mv), return(Mv));
  U = (mathnf(Mv~,1)[2]~)^-1;
  if( n==1 || !redflag, return(U));
\\ extract the n-#v columns and LLL-reduce them
  re = qflll(vecextract(U,1<<n-1,1<<(n-#Mv)-1));
\\ apply the reduction
  re = U*matdiagonalblock([re,matid(#Mv)]);

return(re);
}

\\ Compute the kernel of M mod p.
\\ returns [d,U], where
\\ d = dim (ker M mod p)
\\ U in GLn(Z), and its first d columns span the kernel.
{kermodp(M,p) =
my(U);
my(n,d);

  n = length(M);
  U = centerlift(matker(M*Mod(1,p)));
  d = length(U);
  U = completebasis(U,0);
  U = matrix(n,n,i,j,U[i,n+1-j]);

return([d,U]);
}

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\          INVARIANTS COMPUTATIONS            \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ Compute the Hilbert symbol at p
\\ where p = -1 means real place and not p = 0 as in gp
{myhilbert(a,b,p) =
  if( sign(p) < 0, return(hilbert(a,b,0)));
  return(hilbert(a,b,p));
}

\\ Given a symmetric matrix G over Z, compute the local invariant
\\ (=Witt invariant) of G at the prime p (at real place if p = -1)
\\ Assume that none of the determinant G[1..i,1..i] is 0.
{qflocalinvariant(G,p) =
my(vdet,diag);
my(n,c);

  n = length(G);
\\ Diagonalize G first.
  vdet  = vector( n+1, i, matdet(matrix(i-1,i-1,j,k,G[j,k])));
  diag = vector( n, i, vdet[i+1]*vdet[i]);

\\ Then compute the product of the Hilbert symbols
\\ (diag[i],diag[j])_p for i < j
  c = 1;
  for( i = 1, n,
    for( j = i+1, n,
          c *= myhilbert( diag[i], diag[j], p)));
return(c);
}

\\ G is a quadratic form, or a symmetrix matrix,
\\ or a list of quadratic forms with the same discriminant.
\\ If given, fa must be equal to factor(-abs(2*matdet(G)))[,1].
{qflocalinvariants(G,fa=[]) =
my(vG,sol,vdet);
my(lG);

if( DEBUGLEVEL_qfsolve >= 4, print("    starting qflocalinvariants ",G));

\\ convert G into a vector of symmetric matrices
  if( type(G) == "t_VEC", vG = G , vG = [G]);
  lG = length(vG);
  for( j = 1, lG,
    if( type(vG[j]) == "t_QFI" || type(vG[j]) == "t_QFR",
      vG[j] = Mat(vG[j])));

\\ compute the factorization of -2*abs(det(G))
  if( !length(fa), 
    fa = (factor(-abs(2*matdet(vG[1])))[,1]));

\\ in dimension 2, each invariant is a single Hilbert symbol.
  if( length(vG[1]) == 2,
    sol = matrix(length(fa),lG,i,j,
      myhilbert(vG[j][1,1],-matdet(vG[1]),fa[i]) < 0);
if( DEBUGLEVEL_qfsolve >= 4, print("    end of qflocalinvariants"));
    return([fa,sol])
  );

  sol = matrix(length(fa),lG);
  for( j = 1, lG,
    my( n = length(vG[j]));
\\ in dimension n, we need to compute a product of n Hilbert symbols.
    vdet = vector(n+1, i, matdet(matrix(i-1,i-1,k,m,vG[j][k,m])));
    for( i = 1, length(fa),
      my( p = fa[i]);
      my( h = 1);
      for( k = 1, n-1, h *= myhilbert(-vdet[k],vdet[k+1],p));
      h *= myhilbert(vdet[n],vdet[n+1],p);
      sol[i,j] = h < 0;
    )
  );
if( DEBUGLEVEL_qfsolve >= 4, print("    end of qflocalinvariants"));
return([fa,sol]);
}

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\        QUADRATIC FORMS REDUCTION            \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ M = [a,b;b;c] has integral coefficients.
\\ Gauss reduction of the binary quadratic form
\\   qf = (a,b,c)=a*X^2+2*b*X*Y+c*Y^2
\\ Returns the reduction matrix with det = +1.
{qfbreduce(M) =
my(H,di,aux);
my(a,b,c,q,r,nexta,nextb,nextc);
my(test);

if( DEBUGLEVEL_qfsolve >= 5, print("     starting qfbreduce with ",M));

  a = M[1,1]; b = M[1,2]; c = M[2,2];

  H = matid(2); test = 1;
  while( test && a,
    di = divrem(b,a); q = di[1]; r = di[2];
    if( 2*r > abs(a), r -= abs(a); q += sign(a));
    H[,2] -= q*H[,1];
    nextc = a; nextb = -r; nexta= (nextb-b)*q+c;

    test = abs(nexta) < abs(a);
    if( test,
      c = nextc; b = nextb; a = nexta;
      aux = H[,1]; H[,1] = -H[,2]; H[,2] = aux
    )
  );

if( DEBUGLEVEL_qfsolve >= 5, print("     end of qfbreduce with ",H));
return(H);
}

\\ Performs first a LLL reduction on a positive definite
\\ quadratic form QD bounding the indefinite G.
\\ Then finishes the reduction with qfsolvetriv().
{qflllgram_indef(G,c=1,base=0) =
my(M,QD,M1,S,red);
my(n);

if( DEBUGLEVEL_qfsolve >= 4, print("    qflllgram_indef with G = ",G));
  n = length(G);
  M = matid(n);
  QD = G;
  for( i = 1, n-1,
    if( !QD[i,i],
      return(qfsolvetriv(G,base))
    );
    M1 = matid(n);
    M1[i,] = -QD[i,]/QD[i,i];
    M1[i,i] = 1;
    M = (M*M1);
    QD = (M1~*QD*M1)
  );
  M = (M^(-1));
  QD = (M~*abs(QD)*M);
  S = qflllgram(QD/content(QD));
  if( #S < n, S = completebasis(S));
  red = qfsolvetriv(S~*G*S,base);
  if( type(red) == "t_COL",
    red = (S*red);
    return(red));
  red[2] = S*red[2];
  if( length(red) == 3,
    red[3] = S*red[3]);
return(red);
}

\\ LLL reduction of the quadratic form G (Gram matrix)
\\ where we go on, even if an isotropic vector is found.
{qflllgram_indefgoon(G,c=1) =
my(red,U1,G2,U2,G3,U3,G4,U,V,B,U4,G5,U5,G6);
my(n);

  red = qflllgram_indef(G,c,1);
\\ If no isotropic vector is found, nothing to do.
  if( length(red) == 2, return(red));
\\ otherwise a solution is found:
  U1 = red[2];
  G2 = red[1]; \\ On a G2[1,1] = 0
  U2 = mathnf(Mat(G2[1,]),4)[2];
  G3 = U2~*G2*U2;
\\ The first line of the matrix G3 only contains 0,
\\ except some 'g' on the right, where g^2| det G.
  n = length(G);
  U3 = matid(n); U3[1,n] = round(-G3[n,n]/G3[1,n]/2);
  G4 = U3~*G3*U3;
\\ The coeff G4[n,n] is reduced modulo 2g
  U = vecextract(G4,[1,n],[1,n]);
  if( n == 2,
    V = matrix(2,0)
  , V = vecextract(G4,[1,n],1<<(n-1)-2));
  B = round(-U^-1*V);
  U4 = matid(n);
  for( j = 2, n-1,
    U4[1,j] = B[1,j-1];
    U4[n,j] = B[2,j-1]
  );
  G5 = U4~*G4*U4;
\\ The last column of G5 is reduced
  if( n < 4, return([G5,U1*U2*U3*U4]));

  red = qflllgram_indefgoon(matrix(n-2,n-2,i,j,G5[i+1,j+1]),c);
  U5 = matdiagonalblock([Mat(1),red[2],Mat(1)]);
  G6 = U5~*G5*U5;
return([G6,U1*U2*U3*U4*U5]);
}

\\ LLL reduction of the quadratic form G (Gram matrix)
\\ in dim 3 only, with detG = -1 and sign(G) = [2,1];
{qflllgram_indefgoon2(G,c=1) =
my(red,U1,G2,bez,U2,G3,U3);
my(cc);

  red = qflllgram_indef(G,c,1);
\\ We always find an isotropic vector.
  U1 = [0,0,1;0,1,0;1,0,0];
  G2 = U1~*red[1]*U1;
\\ G2 has a 0 at the bottom right corner.
  bez = bezout(G2[3,1],G2[3,2]);
  U2 = [bez[1],G2[3,2]/bez[3],0;bez[2],-G2[3,1]/bez[3],0;0,0,-1];
  G3 = U2~*G2*U2;
\\ G3 has 0 under the co-diagonal.
  cc = (G3[1,1])%2;
  U3 = [1,0,0;  cc,1,0;
        round(-(G3[1,1]+cc*(2*G3[1,2]+G3[2,2]*cc))/2/G3[1,3]),
        round(-(G3[1,2]+cc*G3[2,2])/G3[1,3]),1];

return([U3~*G3*U3,red[2]*U1*U2*U3]);
}

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\          CLASS GROUP COMPUTATIONS           \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ Compute the square root of the quadratic form qfb.
\\ This function is not fully implemented.
\\ For the moment it only works for detqfb squarefree
\\ (except at 2, where the valuation is 2 or 3).
\\ factdet must be given and equal to factor(2*abs(det qfb))
{qfbsqrt(qfb,factdet) =
my(G,m,n);
my(a,b,c,d,p);
my(aux,Q1,M);

if( DEBUGLEVEL_qfsolve >=3, print("   starting qfbsqrt with ",qfb,factdet));
  G = Vec(qfb);
  a = G[1]; b = (G[2]/2); c = G[3]; d = a*c-b^2;

\\ 1st step: solve m^2 = a (d), m*n = -b (d), n^2 = c (d)
  m = n = Mod(1,1);
  factdet[1,2] -= 3;
  for( i = 1, length(factdet[,1]),
    if( !factdet[i,2], next);
    p = factdet[i,1];
    if( gcd(a,p) == 1,
      aux = sqrt(Mod(a,p));
      m = chinese(m,aux);
      n = chinese(n,-b/aux)
    ,
      aux = sqrt(Mod(c,p));
      n = chinese(n,aux);
      m = chinese(m,-b/aux)
    )
  );
  m = centerlift(m);
  n = centerlift(n);
if( DEBUGLEVEL_qfsolve >=4, print("    m = ",m); print("    n = ",n));

\\ 2nd step: build Q1, with det=-1 such that Q1(x,y,0) = G(x,y)
  Q1 = [(n^2-c)/d, (m*n+b)/d, n ;
        (m*n+b)/d, (m^2-a)/d, m ;
        n,         m,         d ];
  Q1 = -matadjoint(Q1);

\\ 3rd step: reduce Q1 to [0,0,-1;0,1,0;-1,0,0]
  M = qflllgram_indefgoon2(Q1)[2][3,];
  if( M[1] < 0, M = -M);
if( DEBUGLEVEL_qfsolve >=3, print("   end of qfbsqrt"));
  if( M[1]%2,
    return(Qfb(M[1],2*M[2],2*M[3]))
  , return(Qfb(M[3],-2*M[2],2*M[1])));
}

\\ Implementation of Shanks/Bosma-Stevenhagen algorithm
\\ to compute the 2-Sylow of the class group of discriminant D.
\\ Only works for D = fundamental discriminant.
\\ When D = 1(4), work with 4D.
\\ If given, factdetG must be equal to factor(abs(2*D)).
\\ Apart from this factorization, the algorithm is polynomial time.
\\ If Winvariants is given, the algorithm stops as soon as
\\ an element having these W-invariants is found.
{quadclass2(D,factdetG,Winvariants,U2) =
my(factD,listgen,E,invgen,im,Ker,Kerim,listgen2,G2,clgp2,red,auxg);
my(p,aux);
my(n,rang,m,vD,vp);

if( DEBUGLEVEL_qfsolve >= 1, print(" Construction of the 2-class group of discriminant ",D));
  if( D%4 == 2 || D%4 == 3, error("quadclass2: Discriminant not congruent to 0,1 mod 4"));

  if( D==-4, return([[1],[Qfb(1,0,1)]]));

  if( !factdetG, factdetG = factor(2*abs(D)));
  factD = concat([-1],factdetG[,1]~);
  if( D%4 == 1, D *= 4; factdetG[1,2] += 2);

  n = length(factD); rang = n-3;
  if(D>0, m = rang+1, m = rang);
  if(m<0, m=0);
if( DEBUGLEVEL_qfsolve >= 3, print("   factD = ",factD));
  listgen = vector(m);

  vD = valuation(D,2);
  if( vD,
    E = Qfb(1,0,-D/4)
  , E = Qfb(1,1,(1-D)/4)
  );
if( DEBUGLEVEL_qfsolve >= 3, print("   E = ",E));

  if( type(Winvariants) == "t_COL"
    && (Winvariants == 0
        || length(matinverseimage(U2*Mod(1,2),Winvariants))>0)
    , return([[1],[E]]));

  for( i = 1, m, \\ no need to look at factD[1]=-1, nor factD[2]=2
    p = factD[i+2];
    vp = valuation(D,p);
    aux = (p^vp);
    if( vD,
      listgen[i] = Qfb(aux,0,-D/4/aux)
    , listgen[i] = Qfb(aux,aux,(aux-D/aux)/4))
  );
  if( vD == 2 && D%16 != 4,
    m++; rang++; listgen = concat(listgen,[Qfb(2,2,(4-D)/8)]));
  if( vD == 3,
    m++; rang++; listgen = concat(listgen,[Qfb(2^(vD-2),0,-D/2^vD)]));

if( DEBUGLEVEL_qfsolve >= 3, print("   listgen = ",listgen));
if( DEBUGLEVEL_qfsolve >= 2, print("  rank = ",rang));

  if( !rang, return([[1],[E]]));

 invgen = qflocalinvariants(listgen,factD)[2]*Mod(1,2);
if( DEBUGLEVEL_qfsolve >= 3, print("   invgen = ",lift(invgen)));

  clgp2 = vector(m,i,2);
  im = lift(matinverseimage(invgen,matimage(invgen)));
  while( (length(im) < rang) 
  || (type(Winvariants) == "t_COL" && length(matinverseimage(concat(invgen,U2),Winvariants) == 0)),
    Ker = lift(matker(invgen));
    Kerim = concat(Ker,im);
    listgen2 = vector(m);
    for( i = 1, m,
      listgen2[i] = E;
      for( j = 1, m,
        if( Kerim[j,i],
          listgen2[i] = qfbcompraw(listgen2[i],listgen[j])));
      if( norml2(Kerim[,i]) > 1,
        red = qfbreduce(auxg=Mat(listgen2[i]));
        auxg = red~*auxg*red;
        listgen2[i] = Qfb(auxg[1,1],2*auxg[1,2],auxg[2,2]))
    );
    listgen = listgen2;
    invgen = invgen*Kerim;

if( DEBUGLEVEL_qfsolve >= 4, print("    listgen = ",listgen));
if( DEBUGLEVEL_qfsolve >= 4, print("    invgen = ",lift(invgen)));

    for( i = 1, length(Ker),
      G2 = qfbsqrt(listgen[i],factdetG);
      clgp2[i] <<= 1;
      listgen[i] = G2;
      invgen[,i] = qflocalinvariants(G2,factD)[2][,1]*Mod(1,2)
    );

if( DEBUGLEVEL_qfsolve >= 3, print("   listgen = ",listgen));
if( DEBUGLEVEL_qfsolve >= 3, print("   invgen = ",lift(invgen)));
if( DEBUGLEVEL_qfsolve >= 3, print("   clgp2 = ",clgp2));

    im = lift(matinverseimage(invgen,matimage(invgen)))
  );

  listgen2 = vector(rang);
  for( i = 1, rang,
    listgen2[i] = E;
    for( j = 1, m,
      if( im[j,i],
        listgen2[i] = qfbcompraw(listgen2[i],listgen[j])));
    if( norml2(im[,i]) > 1,
      red = qfbreduce(auxg=Mat(listgen2[i]));
      auxg = red~*auxg*red;
      listgen2[i] = Qfb(auxg[1,1],2*auxg[1,2],auxg[2,2]))
  );
  listgen = listgen2;
\\ listgen = vector(rang,i,listgen[m-rang+i]);
  clgp2 = vector(rang,i,clgp2[m-rang+i]);

if( DEBUGLEVEL_qfsolve >= 2, print("  listgen = ",listgen));
if( DEBUGLEVEL_qfsolve >= 2, print("  clgp2 = ",clgp2));

return([clgp2,listgen]);
}

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\          QUADRATIC EQUATIONS                \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ Try to solve G = 0 with small coefficients
\\ This is proved to work if
\\ -  det(G) = 1, dim <= 6 and G is LLL reduced
\\
\\ Returns [G,matid] if no solution is found.
\\ Exit with a norm 0 vector if one such is found.
\\ If base == 1 and norm 0 is obtained, returns [H~*G*H,H,sol] where
\\   sol is a norm 0 vector and is the 1st column of H.
{qfsolvetriv(G,base=0) =
my(n);
my(H,sol,GG);

  n = length(G);
  H = matid(n);

\\ case 1: A basis vector is isotropic
  for( i = 1, n,
    if( G[i,i] == 0,
      sol = H[,i];
      if( base == 0, return(sol));
      H[,i] = H[,1]; H[,1] = sol;
      return([H~*G*H,H,sol])
    )
  );

\\ case 2: G has a block +- [1,0;0,-1] on the diagonal
  for( i = 2, n,
    if( G[i-1,i] == 0 && G[i-1,i-1]*G[i,i] == -1,
      H[i-1,i] = -1; sol = H[,i];
      if( base == 0, return(sol));
      H[,i] = H[,1]; H[,1] = sol;
      return([H~*G*H,H,sol])
    )
  );

\\ case 3: a principal minor is 0
  for( i = 1, n,
    GG = vecextract(G,1<<i-1,1<<i-1);
    if( matdet(GG) != 0, next);
    sol = matker(GG)[,1];
    sol /= content(sol);
    sol = concat(sol,vectorv(n-i));
    if( base == 0, return(sol));
    H = completebasis(sol);
    H[,n] = -H[,1]; H[,1] = sol;
    return([H~*G*H,H,sol])
  );

return([G,H]);
}

\\ p a prime number. 
\\ finds a solution mod p for the quadratic form G
\\ such that det(G) !=0 mod p and dim G = n>=3;
{qfsolvemodp(G,p) =
my(vdet,G2,sol,x1,x2,x3,N1,N2,N3,s);
my(r);
my(n);

  n = length(G);
  vdet = [0,0,0];
  for( i = 1, 3,
    G2 = vecextract(G,1<<i-1,1<<i-1)*Mod(1,p);
    vdet[i] = matdet(G2);
    if( !vdet[i],
      sol = kermodp(lift(G2),p)[2][,1];
      sol = vectorv(n, j, if( j <= i, sol[j], 0));
      return(sol)
    )
  );

\\ now, solve in dimension 3...
\\ reduction to the diagonal case:

  x1 = [1,0,0]~;
  x2 = [-G2[1,2],G2[1,1],0]~;
  x3 = [G2[2,2]*G2[1,3]-G2[2,3]*G2[1,2],G2[1,1]*G2[2,3]-G2[1,3]*G2[1,2],G2[1,2]^2-G2[1,1]*G2[2,2]]~;
  while(1,
    if( issquare( N1 = -vdet[2]),                 s = sqrt(N1); sol = s*x1+x2; break);
    if( issquare( N2 = -vdet[3]/vdet[1]),         s = sqrt(N2); sol = s*x2+x3; break);
    if( issquare( N3 = -vdet[2]*vdet[3]/vdet[1]), s = sqrt(N3); sol = s*x1+x3; break);
    r = 1;
    while( !issquare( s = (1-N1*r^2)/N3), r = random(p));
    s = sqrt(s); sol = x1+r*x2+s*x3; break
  );
  sol = vectorv(n, j, if( j <= 3, sol[j]));
return(sol);
}

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\            HELP MESSAGES                \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
addhelp(default_qfsolve,
  "default_qfsolve(DEBUGLEVEL_qfsolve):
  output or set the value of the global variable DEBUGLEVEL_qfsolve.
  The higher the value, the more information you get about intermediate
  results concerning functions related to qfsolve.
  default is 0: print nothing.");
addhelp(qflllgram_indef,
  "qflllgram_indef(G,{c}): Solve or reduce the quadratic form G with integral coefficients. G might be definite or indefinite. This is an lll-type algorithm with a constant 1/4<c<=1.
  c is optional and the default is c=1.
  The output is either a vectorv which is a solution of v~*G*v=0, or a 2-component vector [H,U], where U is a unimodular matrix such that H = U~*G*U is LLL-reduced.
  Example:
  gp > G=[1637490518557, -9118398255553, -17114399686722; -9118398255553, -40039266946520, 44178901566187; -17114399686722, 44178901566187, 150094052078168];
  gp > qflllgram_indef(G)
  %1 = [-24749181067550, 1904107022307, -3382470700136]~
");
}
