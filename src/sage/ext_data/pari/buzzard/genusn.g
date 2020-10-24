\\ genusn.g, V1.3, by Bruce Kaskel, tinkered with by Kevin Buzzard.
\\
\\ Please report any problems to Kevin Buzzard (buzzard@ic.ac.uk)
\\
\\ These are Bruce Kaskel's programs for computing various things like
\\ dimensions of S_k(Gamma_0(N)) and so on.
\\
\\ Version history:
\\
\\ V1.3, 10/8/2001. Changed the syntax of local variables to
\\ correspond with current pari-2 philosophy.
\\
\\ V1.2: a minor correction by KMB, unfortunately I forgot what it was. I
\\ think it was to do with a "\" sign being interpreted at slightly the
\\ wrong time.
\\
\\ V1.1: original, by Bruce Kaskel.

mu0(n)=local(m,npd);m=factor(n);npd=matsize(m)[1];prod(x=1,npd,1+1/m[x,1],n);
mu20(n)=
{
  local(m,npd);
  if(n%4==0,0,m=factor(n);npd=matsize(m)[1];
    prod(x=1,npd,1+kronecker(-4,m[x,1]),1));
}
mu30(n)=
{
  local(m,npd);
  if(n%2==0 || n%9==0,0,
   m=factor(n);npd=matsize(m)[1];prod(x=1,npd,1+kronecker(-3,m[x,1]),1));
}
c0(n)=sumdiv(n,x,eulerphi(gcd(x,n/x)));
g0(n)=1+mu0(n)/12-mu20(n)/4-mu30(n)/3-c0(n)/2;
\\
mu1(n)=if(n<=2,mu0(n),eulerphi(n)*mu0(n)/2);
mu21(n)=if(n<4,mu20(n),0);
mu31(n)=if(n<4,mu30(n),0);
c1(n)=if(n<=3,c0(n),if(n==4,3,sumdiv(n,x,eulerphi(x)*eulerphi(n/x))/2));
g1(n)=1+mu1(n)/12-mu21(n)/4-mu31(n)/3-c1(n)/2;
\\
ss0(n,p)=
{
  if(isprime(p)==0 || n%p==0,
   print("Function format: ss0(N,p) where p is prime and doesn't divide N.")
  ,g0(n*p)-2*g0(n)+1);
}
\\
muXNp(n,p)=mu1(n)*mu0(p);
mu2XNp(n,p)=0;
mu3XNp(n,p)=0;
cXNp(n,p)=2*c1(n);
gXNp(n,p)=
if(n<4,g0(N*p),1+muXNp(n,p)/12-mu2XNp(n,p)/4-mu3XNp(n,p)/3-cXNp(n,p)/2);
ss1(n,p)=
{
  if(isprime(p)==0 || n%p==0,
   print("Function format: ss1(N,p) where p is prime and doesn't divide N.")
  ,gXNp(n,p)-2*g1(n)+1);
}
\\
eisen(p)=
{
  if(isprime(p)==0,print("Function format: eisen(N) where N is prime.")
  ,numer((p-1)/12));
}
\\
S0(n,k)=
{
  if(n<1,print("S0(N,k), N must be a positive integer")
  ,if(k<=0 || k%2,0,
  if(k==2,g0(n),(k-1)*(g0(n)-1)+(k/2-1)*c0(n)+mu20(n)*(k\4)+mu30(n)*(k\3))));
}
\\
M0(n,k)=
{
  if(n<1,print("M0(N,k), N must be a positive integer")
  ,if(k<=0 || k%2,0,S0(n,k)+c0(n)-(k==2)));
}
S1(n,k)=
{
  if(n<1,print("S1(N,k), N must be a positive integer")
  ,if(k<=0 || (n<3 && k%2),0,
  if(k==1,print("S1(N,k), k=1 not programmed for N>2.")
  ,if(k==2,g1(n)
  ,if(n<3,S0(n,k)
  ,(k-1)*(g1(n)-1)+(k/2-1)*c1(n)+(n==4 && k%2)/2+(n==3)*k\3)))));
}
M1(n,k)=
{
  if(n<1,print("M1(N,k), N must be a positive integer")
  ,if(k<=1,print("M1(N,k), k must be >=2"),S1(n,k)+c1(k)-(k==2)))
}
\\  coefficients of zeta(s)^-2
mumu(n,x,k,p)=
{
  if(type(n)!="t_INT" || n<1,print("mumu(N): N must be a positive integer.")
  ,x=factor(n);k=matsize(x)[1];p=1;
  for(y=1,k,if(x[y,2]>2,p=0,if(x[y,2]==1,p=-2*p,)));p);
}
\\
nf0(n,k)=
{
  if(n<1 || type(n)!="t_INT",print("nf0(N,k): N must be a positive integer.")
  ,if(type(k)!="t_INT",print("nf0(N,k): k must be an integer.")
  ,sumdiv(n,x,S0(x,k)*mumu(n/x))));
}
\\
nf1(n,k)=
{
  if(n<1 || type(n)!="t_INT",print("nf1(N,k): N must be a positive integer.")
  ,if(type(k)!="t_INT",print("nf1(N,k): k must be an integer.")
  ,sumdiv(n,x,S1(x,k)*mumu(n/x))));
}
\\
print("c0(N) -- calculates the number of cusps on X_0(N).");
print("g0(N) -- calculates the genus of X_0(N).");
print("ss0(N,p) -- calculates the number of s.s. points at p of X_0(N).");
print("S0(N,k) -- calculates the dim of the wt. k cusp forms on Gamma_0(N).");
print("nf0(N,k) -- calculates the dim of the wt. k newforms on Gamma_0(N).");
print("eisen(N) -- (N PRIME!) calcs. the order of (0)-(\infty) on J_0(N).");
print("c1(N) -- calculates the number of cusps on X_1(N).");
print("g1(N) -- calculates the genus of X_1(N).");
print("ss1(N,p) -- calculates the number of s.s. points at p of X_1(N).");
print("S1(N,k) -- calculates the dim of the wt. k cusp forms on Gamma_1(N).");
print("nf1(N,k) -- calculates the dim of the wt. k newforms on Gamma_1(N).");
print("Please notify me (buzzard@ic.ac.uk) if any of these functions seem");
print("to be giving erroneous answers.");
print("Original version by Bruce Kaskel.");

