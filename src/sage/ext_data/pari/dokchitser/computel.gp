/************* ComputeL v1.3.8, 2001-2018, (c) Tim Dokchitser ************/
/**************** computing special values of L-functions ****************/
/* arXiv.org/abs/math.NT/0207280, Exper. Math. 13 (2004), no. 2, 137-150 */
/****** Questions/comments welcome! -> tim.dokchitser@bristol.ac.uk ******/

\\ ACKNOWLEDGEMENTS: I'd like to thank Mark Watkins, Steve Donnelly,
\\ William Stein, Anton Mellit, Almasa Odzak, Karim Belabas, Myoungil Kim,
\\ Chris King, F. Patrick Rabarison, Neil Dummigan, Maciej Radziejewski,
\\ François Brunault and Alex Best, for examples, bug fixes and suggestions

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Distributed under the terms of the GNU General Public License (GPL)
\\    This code is distributed in the hope that it will be useful,
\\    but WITHOUT ANY WARRANTY; without even the implied warranty of
\\    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
\\    GNU General Public License for more details.
\\ The full text of the GPL is available at:
\\                 http://www.gnu.org/licenses/
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ USAGE:   Take an L-function L(s) = sum of a(n)/n^s over complex numbers
\\          e.g. Riemann zeta-function, Dedekind zeta-function,
\\               Dirichlet L-function of a character, L-function
\\               of a curve over a number field, L-function of a modular form,
\\               any ``motivic'' L-function, Shintani's zeta-function etc.
\\          assuming L(s) satisfies a functional equation of a standard form,
\\          this package computes L(s) or its k-th derivative for some k
\\          for a given complex s to required precision
\\          - a short usage guide is provided below
\\          - or (better) just look at the example files ex-*
\\            they are hopefully self-explanatory
\\
\\ ASSUMED: L^*(s) = Gamma-factor * L(s) satisfies functional equation
\\          L^*(s) = sgn * L^*(weight-s),
\\            [ more generally L^*(s) = sgn * Ldual^*(weight-s) ]
\\          Gamma-factor = A^s * product of Gamma((s+gammaV[k])/2)
\\            where A = sqrt(conductor/Pi^d)
\\
\\      gammaV    = list of Gamma-factor parameters,
\\                  e.g. [0] for Riemann zeta, [0,1] for ell.curves
\\      conductor = exponential factor (real>0, usually integer),
\\                  e.g. 1 for Riemann zeta and modular forms under SL_2(Z)
\\                  e.g. |discriminant| for number fields
\\                  e.g. conductor for H^1 of curves/Q
\\      weight    = real > 0      (usually integer, =1 by default)
\\                  e.g. 1 for Riemann zeta, 2 for H^1 of curves/Q
\\      sgn       = complex number   (=1 by default)
\\
\\ 1. Read the package (\rcomputel)
\\ 2. Set the required working precision (say \p28)
\\
\\ 3. DEFINE gammaV, conductor, weight, sgn,
\\           Lpoles = vector of points where L^*(s) has (simple) poles
\\             Only poles with Re(s)>weight/2 are to be included
\\           Lresidues = vector of residues of L^*(s) in those poles
\\             or set Lresidues = automatic (default value; see ex-nf)
\\           if necessary, re-define coefgrow(), MaxImaginaryPart (see below)
\\
\\ [4.] CALL cflength()   determine how many coefficients a(n) are necessary
\\      [optional]        to perform L-function computations
\\
\\ 5. CALL initLdata(cfstr) where cfstr (e.g. "(-1)^k") is a string which
\\         evaluates to k-th coefficient a(k) in L-series, e.g.
\\      N    = cflength();                      \\ say returns N=10
\\      Avec = [1,-1,0,1,-1,0,1,-1,0,1,-1,0];   \\ must be at least 10 long
\\      initLdata("Avec[k]");
\\    If Ldual(s)<>L(s), in other words, if the functional equation involves
\\    another L-function, its coefficients are passed as a 3rd parameter,
\\      initLdata("Avec[k]",,"conj(Avec[k])"); see ex-chgen as an example
\\
\\ [7.] CALL checkfeq()     verify how well numerically the functional
\\      [optional]          equation is satisfied
\\                          also determines the residues if Lpoles!=[]
\\                          and Lresidues=automatic
\\    More specifically: for T>1 (default 1.2), checkfeq(T) should ideally
\\    return 0 (with current precision, e.g. 3.2341E-29 for \p28 is good)
\\      * if what checkfeq() returns does not look like 0 at all,
\\        probably functional equation is wrong
\\        (i.e. some of the parameters gammaV, conductor etc., or the coeffs)
\\      * if checkfeq(T) is to be used, more coefficients have to be
\\        generated (approximately T times more), e.g. call
\\           cflength(1.3), initLdata("a(k)",1.3), checkfeq(1.3)
\\      * T=1 always (!) returns 0, so T has to be away from 1
\\      * default value T=1.2 seems to give a reasonable balance
\\      * if you don't have to verify the functional equation or the L-values,
\\           call cflength(1) and initLdata("a(k)",1),
\\           you need slightly less coefficients then
\\
\\ 8. CALL L(s0)    to determine the value of L-function L(s) in s=s0
\\    CALL L(s0,c)  with c>1 to get the same value with a different cutoff
\\                  point (c close to 1); should return the same answer,
\\                  good to check if everything works with right precision
\\                  (if it doesn't, email me!)
\\                  needs generally more coefficients for larger ex
\\                  if L(s0,ex)-L(s0) is large, either the functional eq.
\\                  is wrong or loss of precision (should get a warning)
\\    CALL L(s0,,k) to determine k-th derivative of L(s) in s=s0
\\                  see ex-bsw for example
\\    CALL Lseries(s,,k) to get first k terms of Taylor series expansion
\\                        L(s)+L'(s)S+L''(s)*S^2/2!+...
\\                  faster than k calls to L(s)
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\      Default values for the L-function parameters                      \\
\\      All may be (and conductor and gammaV must be) re-defined          \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ MUST be re-defined, gives error if unchanged
conductor = automatic;
gammaV    = automatic;

\\ MAY be re-defined
weight    = 1;         \\ by default L(s)<->L(1-s)
sgn       = 1;         \\            with sign=1 in functional equation
Lpoles    = [];        \\            and L*(s) has no poles
Lresidues = automatic; \\ if this is not changed to [r1,r2,...] by hand,
                       \\ checkfeq() tries to determine residues automatically
                       \\ see ex-nf for instance

{
coefgrow(n) = if(length(Lpoles),    \\ default bound for coeffs. a(n)
   1.5*n^(vecmax(real(Lpoles))-1),  \\ you may redefine coefgrow() by hand
   sqrt(4*n)^(weight-1));           \\ if your a(n) have different growth
}                                   \\ see ex-delta for example

\\ - For s with large imaginary part there is a lot of cancellation when
\\ computing L(s), so a precision loss occurs. You then get a warning message
\\ - If you want to compute L(s), say, for s=1/2+100*I,
\\ set MaxImaginaryPart=100 before calling initLdata()
\\ - global variable PrecisionLoss holds the number of digits lost in
\\ the last calculation (independently of the MaxImaginaryPart setting)

MaxImaginaryPart = 0;    \\ re-define this if you want to compute L(s)
                         \\ for large imaginary s (see ex-zeta2 for example)

MaxAsympCoeffs  = 40;    \\ At most this number of terms is generated
                         \\ in asymptotic series for phi(t) and G(s,t)
                         \\ default value of 40 seems to work generally well


/******************* IMPLEMENTATION OF THE PACKAGE ************************/


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\         Some helfpul functions                                         \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ Extraction operations on vectors
vecleft(v,n)  = vecextract(v,concat("1..",Str(n)));
vecright(v,n) = vecextract(v,concat(Str(length(v)-n+1),concat("..",Str(length(v)))));

\\ Tabulate a string to n characters, e.g. StrTab(3,2)="3 ";
StrTab(x,n) = x=Str(x);while(length(x)<n,x=concat(x," "));x

\\ Concatenate up to 4 strings
concatstr(s1="",s2="",s3="",s4="")=concat(Str(s1),concat(Str(s2),concat(Str(s3),Str(s4))))

\\ Print a ``small error'', e.g. 0.00000013 as "1E-7"
{
errprint(x)=if(type(x)=="t_COMPLEX",x=abs(x));
   if(x==0,concatstr("1E-",default(realprecision)+1),
   concatstr(truncate(x/10^floor(log(abs(x))/log(10))),"E",floor(log(abs(x))/log(10))));
}


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ gammaseries(z0,terms)                                                   \\
\\ Taylor series expansion of Gamma(z0+x) around 0, z0 arbitrary complex   \\
\\ - up to O(x^(terms+1))                                                  \\
\\ - uses current real precision                                           \\
\\ See Luke "Mathematical functions and their approximations", section 1.4  \
\\ note a misprint there in the recursion formulas [(z-n) term in c3 below] \
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
gammaseries(z0,terms)=
  local(Avec,Bvec,Qvec,n,z,err,res,c0,c1,c2,c3,sinser,reflect,digts,srprec,negint);

  srprec=default(seriesprecision);
  if (z0==real(round(z0)),z0=real(round(z0)));    \\ you don't want to know
  negint=type(z0)=="t_INT" && z0<=0;              \\ z0 is a pole
  default(seriesprecision,terms+1+negint);
  if (terms==0 && !negint,res=gamma(z0)+O(x),     \\ faster to use
  if (terms==1 && !imag(z0) && !negint,           \\   built-in functions
      res=gamma(z0)*(1+psi(z0)*x+O(x^2)),         \\   in special cases
  if (z0==0, res=gamma(1+x)/x,
  if (z0==1, res=gamma(1+x),
  if (z0==2, res=gamma(1+x)*(1+x),
  \\ otherwise use Luke's rational approximations for psi(x)
  digts=default(realprecision);      \\ save working precision
  default(realprecision,digts+3);    \\   and work with 3 digits more
  reflect=real(z0)<0.5;              \\ left of 1/2 use reflection formula
  if (reflect,z0=1-z0);
  z=subst(Ser(precision(1.*z0,digts+3)+X),X,x);
    \\ work with z0+x as a variable gives power series in X as an answer
  Avec=[1,(z+6)/2,(z^2+82*z+96)/6,(z^3+387*z^2+2906*z+1920)/12];
  Bvec=[1,4,8*z+28,14*z^2+204*z+310];
  Qvec=[0,0,0,Avec[4]/Bvec[4]];
  n=4;
  until(err<0.1^(digts+1.5),         \\ Luke's recursions for psi(x)
    c1=(2*n-1)*(3*(n-1)*z+7*n^2-9*n-6);
    c2=-(2*n-3)*(z-n-1)*(3*(n-1)*z-7*n^2+19*n-4);
    c3=(2*n-1)*(n-3)*(z-n)*(z-n-1)*(z+n-4);
    c0=(2*n-3)*(n+1);
    Avec=concat(Avec,[(c1*Avec[n]+c2*Avec[n-1]+c3*Avec[n-2])/c0]);
    Bvec=concat(Bvec,[(c1*Bvec[n]+c2*Bvec[n-1]+c3*Bvec[n-2])/c0]);
    Qvec=concat(Qvec,Avec[n+1]/Bvec[n+1]);
    err=vecmax(abs(Vec(Qvec[n+1]-Qvec[n])));
    n++;
  );
  res=gamma(z0)*exp(intformal( psi(1)+2*(z-1)/z*Qvec[n] )); \\ psi->gamma
  if (reflect,                        \\ reflect if necessary
    sinser=Vec(sin(Pi*z));
    if (negint,sinser[1]=0);          \\ taking slight care at integers<0
    res=subst(Pi/res/Ser(sinser),x,-x);
  );
  default(realprecision,digts);
  )))));
  default(seriesprecision,srprec);
  res;
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ fullgamma(ss) - the full gamma factor (at s=ss)                        \\
\\   vA^s*Gamma((s+gammaV[1])/2)*...*Gamma((s+gammaV[d])/2)               \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
fullgamma(ss) =
  if(ss!=lastFGs,lastFGs=ss;
    lastFGval=prod(j=1,length(gammaV),gamma((ss+gammaV[j])/2),vA^ss)
  );
  lastFGval;
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ fullgammaseries(ss,extraterms) - Laurent series for the gamma factor   \\
\\                                  without the exponential factor, i.e.  \\
\\ Gamma((s+gammaV[1])/2)*...*Gamma((s+gammaV[d])/2)                      \\
\\ around s=ss with a given number of extra terms. The series variable is S.
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
fullgammaseries(ss,extraterms)=
  local(digts,GSD);

  digts=lfundigits;
  if (lastFGSs!=ss || lastFGSterms!=extraterms,
    GSD=sum(j=1,numpoles,(abs((ss+poles[j])/2-round(real((ss+poles[j])/2)))<10^(2-digts)) * PoleOrders[j] )+extraterms;
    lastFGSs=ss;
    lastFGSterms=extraterms;
    lastFGSval=subst(prod(j=1,length(gammaV),gammaseries((ss+gammaV[j])/2,GSD)),x,S/2);
  );
 lastFGSval;
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\         RecursionsAtInfinity(gammaV)                                   \\
\\   Recursions for the asymptotic expansion coefficients                 \\
\\   for phi(x) and G(s,x) at infinity.                                   \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
RecursionsAtInfinity(gammaV)=
  local(d,p,j,k,symvec,modsymvec,deltapol,recF,recG);

  \\ d = number of Gamma-factors in question
  \\ gammaV[k] = Gamma-factors
  \\ symvec = vector of elementary symmetric functions
  \\   1, gammaV[1]+...+gammaV[d], ... , gammaV[1]*...*gammaV[d], 0
  \\ modsymvec = symmetric expressions used in the formula

  d      = length(gammaV);
  symvec = concat(Vec(prod(k=1,d,(x+gammaV[k]))),[0]);

  modsymvec = vector(d+2,k,0);
  for (j=0,d,
  for (k=0,j,
    modsymvec[j+1]+=(-symvec[2])^k*d^(j-1-k)*binomial(k+d-j,k)*symvec[j-k+1];
  ));

  \\ Delta polynomials
  OldSeriesPrecision = default(seriesprecision);
  default(seriesprecision,2*d+2);
  deltapol=subst(Vec( (sinh(x)/x)^tvar ),tvar,x);
  default(seriesprecision,OldSeriesPrecision);

  \\ recursion coefficients for phi at infinity
  recF=vector(d+1,p,
    -1/2^p/d^(p-1)/n*sum(m=0,p,modsymvec[m+1]*prod(j=m,p-1,d-j)*
    sum(k=0,floor((p-m)/2),(2*n-p+1)^(p-m-2*k)/(p-m-2*k)!*subst(deltapol[2*k+1],x,d-p))));

  \\ recursion coefficients for G at infinity
  recG=vector(d,p,recF[p+1]-(symvec[2]+d*(s-1)-2*(n-p)-1)/2/d*recF[p]);

  [vector(d-1,p,recF[p+1]),recG]  \\ return them both
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\                 SeriesToContFrac(vec)                                  \\
\\    Convert a power series vec[1]+vec[2]*x+vec[3]*x^2+...               \\
\\    into a continued fraction expansion.                                \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
SeriesToContFrac(vec)=
  local(res=[],ind);

  vec=1.*vec;
  while (1,
    res=concat(res,[vec[1]]);
    ind=0;
    until(ind==length(vec) || abs(vec[ind+1])>10^-asympdigits,ind++;vec[ind]=0);
    if(ind>=length(vec),break);
    res=concat(res,[ind]);
    vec=Vec(x^ind/Ser(vec));
  );
  res;
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\                 EvaluateContFrac(cfvec,terms,t)                        \\
\\ Evaluate a continued fraction at x=t, using a given number of terms    \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
EvaluateContFrac(cfvec,terms,t)=
  local(res);

  if (terms<0 || terms>length(cfvec)\2,terms=length(cfvec)\2);
  res=cfvec[2*terms+1];
  while(terms>0,res=if(res,cfvec[2*terms-1]+t^cfvec[2*terms]/res,10^asympdigits);terms--);
  res;
}


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ cflength( cutoff=1.2 )                                                \\
\\                                                                        \\
\\ number of coefficients necessary to work with L-series with            \\
\\ current Gamma-factor gammaV, conductor, weight and working precision   \\
\\ - CUTOFF specifies largest t used as a cutoff point in checkfeq(t)     \\
\\   and L(...,t,...). Default is 1.2. Set it to 1 if checkfeq()          \\
\\   is not to be used.                                                   \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
cflength(cutoff=1.2)=
  local(vA,d,expdifff,asympconstf,err,t,t1=0,t2=2,tt,res,lfundigits);

  vA = sqrt(conductor/Pi^length(gammaV));
  d  = length(gammaV);
  lfundigits = default(realprecision) +
     max(ceil(-log(abs(fullgamma(0.7*weight+MaxImaginaryPart*I)))/log(10)),0);
  expdifff = (sum(k=1,d,gammaV[k])+1)/d-1;
  asympconstf = 2*prod(k=1,d,gamma(k/d));
  err = 10^(-lfundigits-0.5);
  until (t2-t1<=1,
    t  = if(t1,(t1+t2)\2,t2);
    tt = t/cutoff/vA;
    res = coefgrow(t) * asympconstf*exp(-d*tt^(2/d))*tt^expdifff;
    if (t1,if(abs(res)>err,t1=t,t2=t),if(abs(res)<err,t1=t2/2,t2*=2));
  );
  ceil(t2)
}


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ initLdata(cfstr, cutoff=1.2 [,cfdualstr])                             \\
\\                                                                        \\
\\ - to be called before the L-values computations.                       \\
\\ - gammaV, conductor and weight must be initialized by now.             \\
\\   also coefgrow(), MaxImaginaryPart and MaxAsympCoeffs are used here   \\
\\ - CFSTR must be a string which evaluates to a function of k            \\
\\   which gives k-th coefficient of the L-series, e.g. "(-1)^k"          \\
\\ - CUTOFF specifies largest t used as a cutoff point in checkfeq(t)     \\
\\   and L(...,t,...). Default is 1.2. Set it to 1 if checkfeq()          \\
\\   is not to be used.                                                   \\
\\ - if cutoff<0, force the number of coefficients to be -cutoff          \\
\\ - CFDUALSTR must evaluate (like cfstr) to the k-th coefficient of      \\
\\   the dual L-function if it is different from L(s),                    \\
\\   for instance initLdata("a(k)",,"conj(a(k))")   (see e.g. ex-chgen)   \\
\\ - uses current real precision to determine the desired precision       \\
\\   for L-values                                                         \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
initLdata(vstr,cutoff=1.2,vdualstr="")=
  local(len,d,pordtmp,recF,recG,terms);

  if (type(gammaV)!="t_VEC",
    error("Gamma-factor gammaV has to be defined before calling initLdata()"));
  if (type(conductor)!="t_INT" && type(conductor)!="t_REAL",
    error("conductor has to be defined before calling initLdata()"));

  len=if(cutoff<0,-cutoff,cflength(cutoff));
  cfvec=vector(len,k,eval(vstr));
  if (vdualstr=="",cfdualvec=cfvec,cfdualvec=vector(len,k,eval(vdualstr)));
  if (cutoff<0,len=cflength());

  lastFGs  = -1E99;         \\ globals to track what was calculated last
  lastFGSs = -1E99;         \\ to avoid re-calculating same values each time
  lastGs   = [-1E99,-1E99];
  lastLSs  = [-1E99,-1E99];

  d       = length(gammaV);                   \\ d = number of gamma factors

  \\ Calculate the necessary amount of extra digits

  answerdigits = default(realprecision);
  vA=sqrt(conductor/Pi^d);
  lfundigits = answerdigits +
     max(ceil(-log(abs(fullgamma(0.7*weight+MaxImaginaryPart*I)))/log(10)),0);
  termdigits   = lfundigits + floor(weight-1);
  taylordigits = 2*termdigits;
  asympdigits  = termdigits;

  \\ Exponential factor defined to maximal precision

  default(realprecision,taylordigits);
  vA     = sqrt(precision(conductor,taylordigits)/Pi^d); \\ exp. factor
  lastt  = len/vA;

  pordtmp = vector(d,k,1);                    \\ locate poles and their orders
  for(j=1,d,for(k=1,d,if(j!=k,\
    if(type(diff=gammaV[j]-gammaV[k])=="t_INT" & (diff%2==0) & (diff<=0),\
      pordtmp[j]+=pordtmp[k];pordtmp[k]=0))));
  poles      = [];
  PoleOrders = [];
  for(j=1,d,if(pordtmp[j]!=0,\
    poles=concat(poles,gammaV[j]);PoleOrders=concat(PoleOrders,pordtmp[j])));
  numpoles   = length(poles);

  \\ Initialize the asymptotic coefficients at infinity

  default(realprecision,asympdigits);

  recFG=RecursionsAtInfinity(gammaV);
  recF=recFG[1];
  recG=recFG[2];
  kill(recFG);

  \\ Maximum number of terms in the asymptotic expansion

  ncoeff=MaxAsympCoeffs;

  \\ Asymptotic behaviour at infinity for phi and G

  expdifff    = (sum(k=1,d,gammaV[k])+1)/d-1;
  asympconstf = 2*prod(k=1,d,gamma(k/d));
  expdiffg    = (sum(k=1,d,gammaV[k])+1)/d-1-2/d;
  asympconstg = prod(k=1,d,gamma(k/d));

  \\ Coefficients for the asymptotic expansion of phi(t) and G(t)

  Fvec=vector(d+ncoeff,X,0);
  Fvec[d]=1.;
  for(y=1,ncoeff,Fvec[d+y]=1.*sum(j=1,d-1,subst(recF[j],n,y)*Fvec[d+y-j]));

  Gvec=vector(d+ncoeff,X,0);
  Gvec[d]=1.;
  for(y=1,ncoeff,Gvec[d+y]=1.*sum(j=1,d,subst(recG[j],n,y)*Gvec[d+y-j]));

  \\ Convert the Fvec (Taylor asymptotic) coefficients into fcf (contfrac coeffs)

  fcf=SeriesToContFrac(vector(ncoeff+1,k,Fvec[d+k-1]));
  fncf=length(fcf)\2;     \\ at most ncoeff+1, less if terminates

  \\ Taylor series coefficients of phi(t) around t=infinity

  if (lastt<35,termstep=1,termstep=floor(lastt^(1/3)));
  phiinfterms=vector(round(lastt/termstep)+1,k,-1);

  terms=fncf;
  PhiCaseBound=0;
  for (k=1,length(phiinfterms),
    t1=(k-1)*termstep;
    while ((k>1)&&(terms>1)&&
      (abs(phiinf(t1,terms-1)-phiinf(t1,terms)))<10^(-termdigits-1),terms-=1);
    if (sum(j=1,terms,fcf[2*j])<ncoeff,phiinfterms[k]=terms,PhiCaseBound=k*termstep);
  );

  \\ Recursions for phi(t) and G(t,s) at the origin

  default(realprecision,taylordigits);

  \\ Initial values of the gamma factors for recursions

  InitV = vector(numpoles,j,prod(k=1,d,
          subst(gammaseries((-poles[j]+gammaV[k])/2,PoleOrders[j]-1),x,X/2)));

  \\ Taylor series coefficients of phi(t) around t=0 -> phiVser

  phiV    = [];
  phiVnn  = 0;
  phiVser = InitV;
  until((phiVnn>3)&&(vecmax(abs(phiV[phiVnn]))*((PhiCaseBound+1)*vA)^(2*phiVnn)<10^(-termdigits-1)),
    RecursephiV());

  default(realprecision,answerdigits);
}


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ phi(t) - inverse Mellin transform of the product of Gamma-factors      \\
\\          computed either with Taylor in 0 or from asymptotics          \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

phi(t)=if(t<PhiCaseBound,phi0(t),phiinf(t,phiinfterms[min(1+floor(abs(t)/termstep),length(phiinfterms))]));

{
RecursephiV()=     \\ compute one more term for the recursions at the origin
  phiVnn++;
  phiV=concat(phiV,[matrix(numpoles,vecmax(PoleOrders),j,k,polcoeff(phiVser[j],-k))]);
  for (j=1,numpoles,for(k=1,length(gammaV),
    phiVser[j]/=(X/2-poles[j]/2-phiVnn+gammaV[k]/2)));
}

{
phi0(t)=                              \\ phi(t) using series expansion at t=0
  local(t2,LogTTerm,TPower,res=0,nn=0,totalold);

  default(realprecision,taylordigits);
  t        = precision(t,taylordigits);
  t2       = t^2;
  LogTTerm = vector(vecmax(PoleOrders),k,(-log(t))^(k-1)/(k-1)!)~;
  TPower   = 1.0*vector(numpoles,j,t^poles[j]);
  until (abs(res-totalold)<10^-(termdigits+1)&&(nn>3),
    totalold=res;
    nn++;
    if(nn>phiVnn,RecursephiV());
    res+=TPower*phiV[nn]*LogTTerm;
    TPower*=t2;
  );
  default(realprecision,termdigits);
  res;
}

{                          \\ phi(t) using asymptotic expansion at infinity
phiinf(t,ncf=fncf)=
  local(res,d,td2);

  default(realprecision,asympdigits);
  t=precision(t,asympdigits);
  d=length(gammaV);
  td2=t^(-2/d);
  res=EvaluateContFrac(fcf,ncf,td2);
  res=res*asympconstf*exp(-d/td2)*t^expdifff;
  default(realprecision,termdigits);
  res;
}


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ G(t,s)   - incomplete Mellin transform of phi(t) divided by x^s        \\
\\            computed either with Taylor in 0 or from asymptotics        \\
\\ G(t,s,k) - its k-th derivative (default 0)                             \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
G(t,ss,der=0)=
  local(nn);

  if(lastGs!=[ss,der] || type(Ginfterms)!="t_VEC",
    initGinf(ss,der);lastGs=[ss,der]);
  if(t<GCaseBound,G0(t,ss,der),
    nn=min(1+floor(abs(t)/termstep),length(Ginfterms));
    Ginf(t,ss,der,Ginfterms[nn]));
}


LogInt(i,j,logt,der)=\
  if (abs(i)<10^(2-termdigits),0,sum(k=0,j-1,binomial((k-j),der)*der!*(-i*logt)^k/k!)/i^(j+der));

{
MakeLogSum(ss,der)=
  local(nn,V,logsum);

  if(length(LogSum)==phiVnn,                      \\ more phiV's necessary
    for (j=1,floor(phiVnn/10)+1,RecursephiV()));
  for (nn=length(LogSum)+1,phiVnn,                \\ generate logsums
    V=phiV[nn];
    logsum=vector(numpoles,j,sum(k=1,PoleOrders[j],V[j,k]*LogInt(poles[j]+2*(nn-1)+ss,k,lt,der)))~;
    LogSum=concat(LogSum,[logsum]);
  );
  lastLSs=[ss,der];
}

{
G0(t,ss,der)=                 \\ G(t,s,der) computed using Taylor series at 0
  local(t2,LT,TPower,res,nn,term,gmser,gmcf,dgts);

  default(realprecision,taylordigits);
  ss     = precision(ss,taylordigits);
  if ([ss,der]!=lastLSs,LogSum=[]);
  t      = precision(t,taylordigits);
  t2     = t^2;       \\ time
  LT     = log(t);    \\ = money
  TPower = vector(numpoles,j,t^poles[j]);
  res    = 0;
  nn     = 0;
  term   = 1;
  until ((nn>3) && abs(term)<10^-(termdigits+1),
    nn++;
    if(nn>length(LogSum),MakeLogSum(ss,der));
    term=TPower*subst(LogSum[nn],lt,LT);
    res+=term;
    TPower*=t2;
  );
  gmser=fullgammaseries(ss,der)/t^(S+ss);
  gmcf=polcoeff(gmser,der,S)*der!;
  res=(gmcf-res);
  default(realprecision,termdigits);
  res;
}


{
Ginf(t,ss,der,ncf=-1)=     \\ G(t,s,der) computed using asymptotic expansion
  local(res,d,tt);         \\ at infinity and associated continued fraction

  default(realprecision,asympdigits);
  ss=precision(ss,asympdigits);
  t=precision(t,asympdigits);
  if (ncf==-1,ncf=gncf);
  d=length(gammaV);
  tt=t^(-2/d);
  res=EvaluateContFrac(gcf,ncf,tt);
  res=asympconstg*exp(-d/tt)*t^expdiffg*tt^der*res;
  default(realprecision,termdigits);
  res;
}


{
initGinf(ss,der)=          \\ pre-compute asymptotic expansions for a given s

  local(d,gvec,gncf,terms,t1);
  default(realprecision,asympdigits);
  ss=precision(ss,asympdigits);
  d=length(gammaV);
  gvec=Gvec;
  for (k=1,der,gvec=deriv(gvec,s);gvec=concat(vecright(gvec,length(gvec)-1),1));
  gcf=SeriesToContFrac(vector(ncoeff+1,k,subst(gvec[d+k-1],s,ss)));
  gncf=length(gcf)\2;
  Ginfterms=vector(round(lastt/termstep)+1,k,-1);
  terms=gncf;
  GCaseBound=0;
  for (k=1,length(Ginfterms),
    t1=(k-1)*termstep;
    while ((k>1)&&(terms>1)&&
      (abs(Ginf(t1,ss,der,terms-1)-Ginf(t1,ss,der,terms)))<10^(-termdigits-2),terms-=1);
    if (sum(j=1,terms,gcf[2*j])<ncoeff,Ginfterms[k]=terms,GCaseBound=k*termstep);
  );
  default(realprecision,termdigits);
}


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ ltheta(t) = sum a(k)*phi(k*t/vA)                                        \\
\\ satisfies ltheta(1/t)=sgn*t^weight*ldualtheta(t) + residue contribution \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

ltheta(t) = t=precision(t/vA,taylordigits);\
  sum(k=1,min(floor(lastt/t+1),length(cfvec)),cfvec[k]*if(cfvec[k],phi(k*t),0));
ldualtheta(t) = t=precision(t/vA,taylordigits);\
  sum(k=1,min(floor(lastt/t+1),length(cfvec)),cfdualvec[k]*if(cfdualvec[k],phi(k*t),0));

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ checkfeq(t=1.2) = verify the functional equation by evaluating LHS-RHS  \\
\\                   for func.eq for ltheta(t), should return approx. 0    \\
\\ - also determines residues if Lresidues is set to automatic             \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
checkfeq(t=6/5)=                \\ determine residues if they are not yet set
  local(nlp,lpx,lpv,lpm,res);   \\ .. and check the functional equation

  if (Lresidues==automatic && Lpoles!=[],
    nlp=length(Lpoles);
    lpx=vector(nlp,k,1.15+if(k==1,0,k-1)/10);
    lpv=vector(nlp,k,tq=lpx[k];sgn*tq^weight*ldualtheta(tq)-ltheta(1/tq));
    lpm=matrix(nlp,nlp,k,j,tq=lpx[k];ss=Lpoles[j];tq^ss-sgn*tq^(weight-ss));
    Lresidues=matsolve(lpm,lpv~)~;
  );
  res=ltheta(t)-sgn*t^(-weight)*ldualtheta(1/t)+
    sum(k=1,length(Lpoles),Lresidues[k]*(t^-Lpoles[k]-sgn*t^(-weight+Lpoles[k])));
  default(realprecision,answerdigits);
  res;
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ L(ss,cutoff,k) = k-th derivative of L(s) at s=ss                       \\
\\                  cutoff = 1 by default (cutoff point), >=1 in general  \\
\\                  must be equal to 1 if k>0                             \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

L(ss,cutoff=1,der=0)=polcoeff(Lseries(ss,cutoff,der),der)*der!


{       \\ Lseries(s,,der) = L(s)+L'(s)S+L''(s)*S^2/2!+... ,first der terms
Lseries(ss,cutoff=1,der=0)=

  local(FGSeries,LGSeries,res);

  default(realprecision,lfundigits);
  FGSeries = fullgammaseries(ss,der)*vA^(ss+S);
  if (length(Lpoles) && (vecmin(abs(vector(length(Lpoles),k,Lpoles[k]-ss)))<10^(-answerdigits) ||
    vecmin(abs(vector(length(Lpoles),k,weight-Lpoles[k]-ss)))<10^(-answerdigits)),
    error("L*(s) has a pole at s=",ss));
\\  if (vecmin(abs(vector(length(Lpoles),k,Lpoles[k]-ss)))<10^(-answerdigits) ||
\\    vecmin(abs(vector(length(Lpoles),k,weight-Lpoles[k]-ss)))<10^(-answerdigits),
\\    error("L*(s) has a pole at s=",ss));
  PrecisionLoss = ceil(-log(vecmax(abs(Vec(FGSeries))))/log(10))
    -(lfundigits-answerdigits);
  if (PrecisionLoss>1,
    print("Warning: Loss of ",PrecisionLoss," digits due to cancellation"));
  LSSeries = sum(k=0,der,Lstar(ss,cutoff,k)*S^k/k!)+O(S^(der+1));
  res=LSSeries/FGSeries;
  default(realprecision,answerdigits);
  res;
}


{
Lstar(ss,cutoff=1,der=0)=                \\ Lstar(s) = L(s) * Gamma factor
  local(res,ncf1,ncf2);

  if (der & (cutoff!=1),error("L(s,cutoff,k>0) is only implemented for cutoff=1"));
  ss     = precision(ss,taylordigits);
  cutoff = precision(cutoff,taylordigits);
  ncf1   = min(round(lastt*vA*cutoff),length(cfvec));
  ncf2   = min(round(lastt*vA/cutoff),length(cfvec));
  default(realprecision,termdigits);
  res=(-sum(k=1,length(Lpoles),(-1)^der*der!*Lresidues[k]/(ss-Lpoles[k])^(der+1)*cutoff^(-Lpoles[k]))
       -sgn*sum(k=1,length(Lpoles),Lresidues[k]*der!/(weight-Lpoles[k]-ss)^(der+1)*cutoff^(-weight+Lpoles[k]))
       +sgn*sum(k=1,ncf1,if(cfdualvec[k],cfdualvec[k]*(-1)^der*G(k/vA/cutoff,weight-ss,der),0)/cutoff^weight)
       +sum(k=1,ncf2,if(cfvec[k],cfvec[k]*G(k*cutoff/vA,ss,der),0))
  )*cutoff^ss;
  res;
}
