\\ DimensionSk.g, V2.2: program to compute dimensions of spaces of cusp
\\ forms with non-trivial character.
\\
\\ By Kevin Buzzard (buzzard@ic.ac.uk)
\\
\\ Please report any problems to buzzard@ic.ac.uk .
\\
\\ Needs: genusn.g function S0(n,k)
\\
\\ This started off as my port of William's DimensionCuspForms magma script,
\\ and has grown since then.
\\
\\
\\ Version history:
\\
\\ V2.2: Fixed a bug in localevaluate() pointed out by Fernando Gouvea:
\\       s was used as a local variable but not declared as such.
\\       (this did actually cause problems with another script, Tpprog.g,
\\       which happened to use s for something else), 14/2/2002.
\\
\\ V2.1: minor tidying up and added a few \\ lines, 10/8/2001.
\\
\\ V2.0: severe tidying up, e.g. got rid of all the \s at end
\\ of lines and used the { trick instead; used local() for local
\\ variables instead of the old version, and so on, 9/8/2001.
\\
\\ V1.0, 26/7/00, just found it in the mess.
\\
\\ TODO: Optimise CohenOesterle2 and CohenOesterle3, currently they
\\       are ridiculous. See the magma version.
\\
\\ Remarks: has a superugly way of storing Dirichlet characters.
\\
\\ Syntax for characters:
\\
\\ eps is [N, i x 3 matrix], where eps[2][,1] is the primes dividing
\\ N, eps[2][,2] is the powers of these primes that divide N, and eps[2][,3]
\\ is the following: for p odd, p^n||N, it's t such that znprimroot(p^n)
\\ gets sent to exp(2*pi*i/phi(p^n))^t. And for p=2, it's
\\ 0 for 2^1, it's 0 (trivial) or -1 (non-trivial) for 2^2, and for p^n>=8
\\ it's either t>=0 for the even char sending 5 to exp(2*pi*i/p^(n-2))^t,
\\ or t<=-1 for the odd char sending 5 to exp(2*pi*i/p^(n-2))^(-1-t).
\\ (so either 0<=t<2^(n-2) or -1>=t>-1-2^(n-2) )

\\ Examples of creation:
\\
\\ TrivialCharacter(100) creates the trivial character of level 100.
\\ DirichletCharacter(45,[1,2]) creates the character of level 45
\\ which is the product of the character of level 9 sending
\\ znprimroot(9)=2 to exp(2*Pi*I/6)^1 and the character of level 5
\\ sending znprimroot(5)=2 to exp(2*Pi*I/4)^2=-1.

print("DimensionCuspForms(eps,k) computes the dimension of the space");
print("of cusp forms of character eps and weight k (note that eps knows");
print("its level!) The syntax for characters is rather complicated, see");
print("the detailed remarks at the beginning of the source code for");
print("more information.");


nearestinteger(z)=
{
  if(abs(z-round(z))>10^(-10),
    error("argument ",z," not near to an integer in nearestinteger()")
  ,
    real(round(z))
  )
}

DimensionCuspForms(eps,k)=
{
  local(N=0);
  if(k<=0,
    0
  ,
    if(k==1,
      error("Can't compute dimension of a space of cusp forms of weight one")
    );

    N=eps[1];
    if(istrivial(eps),
      S0(N,k)
    ,
      if(nearestinteger(evaluate(eps,N-1))!=(-1)^k,
        0
      ,
        nearestinteger(idxG0(N)*(k-1)/12 \
         +CohenOesterle1(eps,k)+CohenOesterle2(eps,k)+CohenOesterle3(eps,k))
      )
    )
  )
}

DimensionModularForms(eps,k)=
{
  local(N=0);
  if(k<=0,
    0
  ,
    if(k==1,
      error("Can't compute dimension of a space of modular
       forms of weight one")
    );
    \\ k>=2

    N=eps[1];
    if(nearestinteger(evaluate(eps,N-1))!=(-1)^k,
      0
    ,
      if(istrivial(eps),
        S0(N,k)+c0(N)-(k==2)
      ,
        nearestinteger(idxG0(N)*(k-1)/12 \
         -CohenOesterle1(eps,k)+CohenOesterle2(eps,k)+CohenOesterle3(eps,k) \
         -(k==2))
      )
    )
  )
}

CohenOesterle1(eps,k)=
{
  local(facN,f,facf);
  facN=matrix(matsize(eps[2])[1],2,i,j,eps[2][i,j]);
  f=charconductor(eps);
  facf=facN;facf[,2]=vectorv(matsize(facN)[1],i,valuation(f,facN[i,1]));

  (-1/2)*prod(i=1,matsize(facN)[1],lambda(facN[i,2],facf[i,2],facN[i,1]))
}

CohenOesterle2(eps,k)=
{
  local(N,gamma_k=0);
  N=eps[1];
  if(k%4==2,gamma_k=-1/4,if(k%4==0,gamma_k=1/4));
  if(gamma_k==0,
    0
  ,
    gamma_k*sum(i=0,N-1,if((i^2+1)%N==0,evaluate(eps,i),0))
  )
}

CohenOesterle3(eps,k)=
{
  local(N,mu_k=0);
  N=eps[1];
  if(k%3==2,mu_k=-1/3,if(k%3==0,mu_k=1/3));
  if(mu_k==0,
    0
  ,
    mu_k*sum(i=0,N-1,if((i^2+i+1)%N==0,evaluate(eps,i),0))
  )
}

lambda(r,s,p)=
{
  if(2*s<=r,
    if(r%2==0,
      p^(r\2)+p^((r\2)-1)
    ,
      2*p^((r-1)\2)
    )
  ,
   2*p^(r-s)
  )
}

evaluate(eps,n)=
{
  if(gcd(n,eps[1])>1,
    0
  ,
    if(istrivial(eps),
      1
    ,
      \\ main case here.
      prod(i=1,matsize(eps[2])[1],localevaluate(eps[2][i,],n))
    )
  )
}

localevaluate(v,n)=
{
  local(p,e,s,t);
  p=v[1];e=v[2];t=v[3];
  if(p==2, \\ deal with 2 separately
    if(e==1,
      1
    ,
      s=if(n%4==1,1,-1);
      if(t>=0,
        exp(2*Pi*I/2^(e-2)*t*z2log(Mod(n*s,2^e)))
      ,
        s*exp(2*Pi*I/2^(e-2)*(-1-t)*z2log(Mod(n*s,2^e)))
      )
    )
  ,
    exp(2*Pi*I/(p-1)/p^(e-1)*t*znlog(n,znprimroot(p^e)))
  )
}

charconductor(eps)=
{
  if(istrivial(eps),
    1
  ,
    prod(i=1,matsize(eps[2])[1],localconductor(eps[2][i,]))
  )
}

localconductor(v)=
{
  local(p,e,t,tp);
  p=v[1];e=v[2];t=v[3];
  if(p==2, \\ deal with this case first
    if(e==1||t==0,
      1
    ,
      if(e==2||t==-1,
        4
      ,
        \\ p=2 and e>=3 and im(5) isn't trivial
        tp=if(t>=0,t,-1-t); \\ 1<=tp<2^(e-2)
        2^(e-valuation(tp,2))
      )
    )
  ,
  \\ p odd
    if(t==0,
      1
    ,
      p^(e-valuation(t,p))
    )
  )
}

z2log(t)=
{
  \\ This returns the 2-adic log of t.
  if(padicprec(t,2)<3,
    1
  ,
    truncate(log(lift(t)+O(2^padicprec(t,2)))/log(5+O(2^padicprec(t,2))))
  )
}

idxG0(N)=
{
  local(facN);
  facN=factor(N);
  prod(i=1,matsize(facN)[1],facN[i,1]^facN[i,2]+facN[i,1]^(facN[i,2]-1))
}

istrivial(eps)=(eps[2][,3]==vectorv(matsize(eps[2])[1],i,0))

TrivialCharacter(N)=
{
  local(facN);
  facN=factor(N);

  [N,matrix(matsize(facN)[1],3,i,j,if(j<3,facN[i,j],0))]
}

DirichletCharacter(N,v)=
{
  local(facN);
  facN=factor(N);
  [N,matrix(matsize(facN)[1],3,i,j,if(j<3,facN[i,j],v[i]))]
}

\\
\\ Syntax for DirichletCharacter: takes (N,v) where N>0 is an integer
\\ and v is a vector of length equal to the number of factors of N;
\\ p_i will get assigned v[i]. Note that we *CAN* have factors
\\ in a non-standard order. Equality of chars will be difficult but
\\ who cares!

charmult(chi1,chi2)=
{
  if(gcd(chi1[1],chi2[1])>1,
    error("Too lazy to multiply Dirichlet characters together properly")
  );

  [chi1[1]*chi2[1],concat(chi1[2]~,chi2[2]~)~]
}


