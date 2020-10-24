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
  www.math.unicaen.fr/~simon/ell.gp

  *********************************************
  *          VERSION 13/01/2014               *
  *********************************************

*/

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\          SCRIPT                             \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\    COMMON FUNCTIONS TO ell.gp AND ellQ.gp   \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{ellcomposeurst(urst1,urst2) =
my(u1 = urst1[1], r1 = urst1[2], s1 = urst1[3], t1 = urst1[4],
      u2 = urst2[1], r2 = urst2[2], s2 = urst2[3], t2 = urst2[4]);
  [u1*u2,u1^2*r2+r1,u1*s2+s1,u1^3*t2+s1*u1^2*r2+t1];
}
{ellinverturst(urst) =
my(u = urst[1], r = urst[2], s = urst[3], t = urst[4]);
  [1/u,-r/u^2,-s/u,(r*s-t)/u^3];
}
{mysubst(polsu,subsx) = 
  if( type(lift(polsu)) == "t_POL",
    return(simplify(subst(lift(polsu),variable(lift(polsu)),subsx)))
  , return(simplify(lift(polsu))));
}
{degre(idegre) =
my(ideg = idegre, jdeg = 0);

  while( ideg >>= 1, jdeg++);
  return(jdeg);
}
{nfrealsign(nf,a,i) =
\\ return the sign of the algebraic number a in the i-th real embedding.
my(nf_roots,ay,prec0);

  if( a == 0, return(0));

  a = lift(a);
  if( type(a) != "t_POL",
    return(sign(a)));

  nf_roots = nf.roots;
  prec0 = default(realprecision);

  ay = 0;
  while( ay == 0 || precision(ay) < 10,

    ay = subst(a,variable(a),nf_roots[i]);

    if( ay == 0 || precision(ay) < 10,
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
{nfsqrt( nf, a) =
\\ if a is a square in the number field nf returns [sqrt(a)], otherwise [].
my(alift,ta,minpola,py,pfact);

  if( a==0 || a==1, return([a]));

  alift = lift(a);
  ta = type(a);
  if( !poldegree(alift), alift = polcoeff(alift,0));

  if( type(alift) != "t_POL",
    if( issquare(alift), return([sqrtrat(alift)])));

  if( poldegree(nf.pol) <= 1, return([]));
  if( ta == "t_POL", a = Mod(a,nf.pol));

\\ the norm should be a square

  if( !issquare(norm(a)), return([]));

\\ the real embeddings must all be >0

  minpola = minpoly(a);
  if( polsturm(minpola,,0), return([]));

\\ factorization over nf of the polynomial X^2-a

  if( variable(nf.pol) == 'x,
    py = subst(nf.pol,'x,'y);
    pfact = lift(factornf('x^2-mysubst(alift,Mod('y,py)),py)[1,1])
  ,
    pfact = lift(factornf('x^2-a,nf.pol)[1,1]));
  if( poldegree(pfact) == 2, return([]));
  return([subst(polcoeff(pfact,0),'y,Mod(variable(nf.pol),nf.pol))]);
}
{nfissquare(nf, a) = #nfsqrt(nf,a) > 0;
}
{sqrtrat(a) = 
  sqrtint(numerator(a))/sqrtint(denominator(a));
}

