##############################################################
##
##           Helps compute the decomposition of
##           the ramification module and the equivariant
##           degree of a multiple of a simple orbit for
##           the Hurwitz curves X with automorphism group
##           G = PSL(2,q), q a prime.
##
##           Using these computations and Borne's formula, the
##           G-module structure of the RR spaces of equivariant
##           divisors can be determined explicitly.
##
## Reference: David Joyner, Amy Ksir, Roger Vogeler,
##            "Group representations on Riemann-Roch spaces
##            of some Hurwitz curves," preprint, 2006.
##
## 6-9-2006, wdj
##############################################################

## To read this in, use for example:
## Read("/home/wdj/gapfiles/curves/hurwitz_crv_rr_sp.gap");
## To log your results, use for example:
## LogTo("/home/wdj/gapfiles/curves/hurwitz_crv_rr_sp1.log");

hurwitz_primes:=function(n)
# Returns the primes < n s.t. PSL(2,p) is a Hurwitz gp
# Assumes n<1000.
local p,L;
L:=[];
for p in Primes do
 if Int((p+1)/7)=((p+1)/7) and p<n then
   L:=Concatenation(L,[p]);
 fi;
 if Int((p-1)/7)=((p-1)/7) and p<n then
   L:=Concatenation(L,[p]);
 fi;
od;
return L;
end;

hurwitz_prime_powers:=function(n)
## Returns the prime powers q=p^3 < n s.t. PSL(2,q)
## is a Hurwitz gp
## Assumes n<1000.
local p,L;
L:=[];
for p in Primes do
 if Int((p+2)/7)=((p+2)/7) and p^3<n then
   L:=Concatenation(L,[p^3]);
 fi;
 if Int((p-2)/7)=((p-2)/7) and p^3<n then
   L:=Concatenation(L,[p^3]);
 fi;
 if Int((p+3)/7)=((p+3)/7) and p^3<n then
   L:=Concatenation(L,[p^3]);
 fi;
 if Int((p-3)/7)=((p-3)/7) and p^3<n then
   L:=Concatenation(L,[p^3]);
 fi;
od;
return L;
end;


ram_module_hurwitz:=function(p)
##
## input: p is a Hurwitz prime; output: [m1,...,mn],
## where n = # conj classes of G=PSL(2,p)
## and mi = mult of pi_i in ram mod of Hurwitz crv with
## aut gp G.
## Here Irr(G) = [pi_1,...,pi_n] (in that order).
##
local A,B,C,D,pi,G,i,j,n,n0,H,G1,H_chars,CG,G_chars,w,m,theta,pi_theta_to_the;
  if not(p in hurwitz_primes(500)) then
      Print("Input is not a (small) Hurwitz prime.\n");
      return [];
  fi;
  G:=PSL(2,p);
  H:=[];
  H_chars:=[];
  n0:=[ [ Z(p)^0, Z(p) ], [ Z(p)*0, Z(p)^0 ] ];
  CG:=ConjugacyClassesSubgroups(G);
  H[1]:=Representative(CG[2]); # size 2
  H[2]:=Representative(CG[3]); # size 3
  n :=Size(CG);
  for i in [1..n] do
    if Size(Representative(CG[i]))=7 then
      H[3]:=Representative(CG[i]);
    fi;
  od;
## H[3]:=Group(n0);       # size 7
  for i in [1..Size(H)] do
   H_chars[i]:=Irr(H[i]);
  od;
  G_chars:=Irr(G);
  m:=[];
  m[1]:=[];m[2]:=[];m[3]:=[];
  theta:=List([1..3],i->H_chars[i][2]);
  pi_theta_to_the:=List([1..3],i->Sum([1..(Size(H[i])-1)],
      j->j*InducedClassFunction(theta[i]^(j),G)));
#Print("\n\n pi_theta_to_the = ",Sum(pi_theta_to_the),"\n\n");
  for i in [1..3] do
    m[i]:=List(G_chars, pi->ScalarProduct(pi_theta_to_the[i],pi))/Size(H[i]);
  od;
#Print("\n\n m = ",m,"\n\n");
#Print("Multiplicities (by definition):\n ",Sum(m),"\n");
return Sum(m);
####
n:=Length(G_chars);
A:=[];
B:=[];
C:=[];
D:=[];
for i in [1..n] do
  pi:=G_chars[i];
  A[i]:=ScalarProduct(H_chars[1][1],RestrictedClassFunction(pi,H[1]));
  B[i]:=ScalarProduct(H_chars[2][1],RestrictedClassFunction(pi,H[2]));
  C[i]:=ScalarProduct(H_chars[3][1],RestrictedClassFunction(pi,H[3]));
  D[i]:=DegreeOfCharacter(pi);
od;
#Print("Multiplicities (predicted by JK - assumes all R_\ell=1):\n ",(1/2)*(3*D-A-B-C),"\n");
end;

#########example
#gap> ram_module_hurwitz(13); time;
#
#m = [ [0,2,2,3,3,3,3,4,3],[0,2,2,4,4,4,4,5,5],[0,3,3,5,5,5,6,6,6] ]
#
#Multiplicities (by definition):
# [ 0, 7, 7, 12, 12, 12, 13, 15, 14 ]
#Multiplicities (predicted by JK - assumes all R_ell=1):
# [ 0, 7, 7, 12, 12, 12, 13, 15, 14 ]
#1640
#

ram_module_X:=function(p)
## p is a Hurwitz prime
## output is (m1,...,mn)
## where n = # conj classes of G=PSL(2,p)
## and mi = mult of pi_i in ram mod of G
##
local G,i,j,n,n0,H,G1,H_chars,CG,G_chars,w,m,theta,pi_theta_to_the;
  if not(p in hurwitz_primes(500)) then
      Print("Input is not a (small) Hurwitz prime.\n");
      return [];
  fi;
  G:=PSL(2,p);
  H:=[];
  H_chars:=[];
  n0:=[ [ Z(p)^0, Z(p) ], [ Z(p)*0, Z(p)^0 ] ];
  CG:=ConjugacyClassesSubgroups(G);
  H[1]:=Representative(CG[2]); # size 2
  H[2]:=Representative(CG[3]); # size 3
  n :=Size(CG);
  for i in [1..n] do
    if Size(Representative(CG[i]))=7 then
      H[3]:=Representative(CG[i]);
    fi;
  od;
## H[3]:=Group(n0);       # size 7
  for i in [1..Size(H)] do
   H_chars[i]:=Irr(H[i]);
  od;
  G_chars:=Irr(G);
  m:=[];
  m[1]:=[];m[2]:=[];m[3]:=[];
  theta:=List([1..3],i->H_chars[i][2]);
  pi_theta_to_the:=List([1..3],i->Sum([1..(Size(H[i])-1)],
      j->j*InducedClassFunction(theta[i]^(j),G)));
#Print("\n\n pi_theta_to_the = ",Sum(pi_theta_to_the),"\n\n");
  for i in [1..3] do
    m[i]:=List(G_chars, pi->ScalarProduct(pi_theta_to_the[i],pi))/Size(H[i]);
  od;
#Print("\n\n m = ",m,"\n\n");
return Sum(m);
end;


ram_module_X_JK:=function(p)
## p is a Hurwitz prime
## output is (m1,...,mn)
## where n = # conj classes of G=PSL(2,p)
## and mi = "mult of pi_i in ram mod of G" using JK formula
##
local G,i,j,n,n0,H,G1,H_chars,CG,G_chars,A,B,C,D,pi;
G:=PSL(2,p);
H:=[];
H_chars:=[];
n0:=[ [ Z(p)^0, Z(p) ], [ Z(p)*0, Z(p)^0 ] ];
CG:=ConjugacyClassesSubgroups(G);
H[1]:=Representative(CG[2]); # size 2
H[2]:=Representative(CG[3]); # size 3
n:=Size(CG);
for i in [1..n] do
  if Size(Representative(CG[i]))=7 then
    H[3]:=Representative(CG[i]);
  fi;
od;
## H[3]:=Group(n0);       # size 7
for i in [1..Size(H)] do
 H_chars[i]:=Irr(H[i]);
od;
G_chars:=Irr(G);
n:=Length(G_chars);
A:=[];
B:=[];
C:=[];
D:=[];
for i in [1..n] do
  pi:=G_chars[i];
  A[i]:=ScalarProduct(H_chars[1][1],RestrictedClassFunction(pi,H[1]));
  B[i]:=ScalarProduct(H_chars[2][1],RestrictedClassFunction(pi,H[2]));
  C[i]:=ScalarProduct(H_chars[3][1],RestrictedClassFunction(pi,H[3]));
  D[i]:=DegreeOfCharacter(pi);
od;
return (1/2)*(3*D-A-B-C); #assumes all R_\ell=1
end;


equiv_deg_module_hurwitz:=function(p,ii,r)
## p is a Hurwitz prime
## ii  = 1 for H[1], size 2
## ii  = 2 for H[2], size 3
## ii  = 3 for H[3], size p
## output is (m1,...,mn)
## where n = # conj classes of G=PSL(2,p)
## and mi = mult of pi_i in deg_equiv module of G
##
local G,i,j,n,n0,H,G1,H_chars,CG,G_chars,w,m,theta,pi_theta,pi_theta_to_the;
G:=PSL(2,p);
H:=[]; # 3 decomp gps
H_chars:=[];
n0:=[ [ Z(p)^0, Z(p) ], [ Z(p)*0, Z(p)^0 ] ];
CG:=ConjugacyClassesSubgroups(G);
H[1]:=Representative(CG[2]); # size 2
H[2]:=Representative(CG[3]); # size 3
n:=Size(CG);
for i in [1..n] do
  if Size(Representative(CG[i]))=7 then
    H[3]:=Representative(CG[i]);
  fi;
od;
## H[3]:=Group(n0);       # size 7
for i in [1..3] do
 H_chars[i]:=Irr(H[i]);
od;
G_chars:=Irr(G);
m:=[];
m[1]:=[];m[2]:=[];m[3]:=[];
theta:=List([1..3],i->H_chars[i][2]);
pi_theta_to_the:=List([1..3],i->List([1..r],j->InducedClassFunction(theta[i]^(-j),G)));
for i in [1..3] do
 for j in [1..r] do
  m[i][j]:=List(G_chars, pi->ScalarProduct(pi_theta_to_the[i][j],pi));
 od;
od;
return Sum(m[ii]);
end;


# equiv_deg_module_X(7,2,1);
#[ 0, 1, 1, 2, 2, 3 ]
# equiv_deg_module_X(7,1,1);
#[ 0, 2, 2, 2, 4, 4 ]
# equiv_deg_module_X(7,3,1);
#[ 0, 1, 0, 1, 1, 1 ]
# equiv_deg_module_X(7,3,7);
#[ 1, 3, 3, 6, 7, 8 ]
# equiv_deg_module_X(7,1,2);
#[ 1, 3, 3, 6, 7, 8 ]


cuberoots:=function(q)
local x,L;
 L:=[];
 for x in GF(q^2) do
   if x<>Zero(GF(q^2)) and x<>One(GF(q^2)) and x^3=One(GF(q^2)) then
     L:=Concatenation(L,[x]);
   fi;
 od;
 return L;
end;

#gap> cuberoots(41);
#[ Z(41^2)^560, Z(41^2)^1120 ]
#gap> cuberoots(43);
#[ Z(43)^14, Z(43)^28 ]
#gap> fourthroots(43);
#[ Z(43^2)^462, Z(43^2)^1386 ]
#gap> fourthroots(41);
#[ Z(41)^10, Z(41)^30 ]
#gap> fourteenthroots(41);
#[ Z(41^2)^120, Z(41^2)^600, Z(41^2)^1560, Z(41^2)^1320, Z(41^2)^1080, Z(41^2)^360 ]
#gap> fourteenthroots(43);
#[ Z(43)^3, Z(43)^9, Z(43)^15, Z(43)^27, Z(43)^33, Z(43)^39 ]


fourteenthroots:=function(q)
local x,L;
 L:=[];
 for x in GF(q^2) do
   if x<>Zero(GF(q^2)) and x<>One(GF(q^2)) and x^2<>One(GF(q^2))
                     and x^7<>One(GF(q^2)) and x^(14)=One(GF(q^2)) then
     L:=Concatenation(L,[x]);
   fi;
 od;
 return L;
end;

fourthroots:=function(q)
local x,L;
 L:=[];
 for x in GF(q^2) do
   if x<>Zero(GF(q^2)) and x<>One(GF(q^2)) and
      x^2<>One(GF(q^2)) and x^(4)=One(GF(q^2)) then
     L:=Concatenation(L,[x]);
   fi;
 od;
 return L;
end;

Nalpha:=function(alpha,q)
 # alpha a character of GF(q)^x
 local root3,root4,root14,count;
 count := 0;
 root3:=cuberoots(q)[1];
 root4:=fourthroots(q)[1];
 root14:=fourteenthroots(q)[1];
 if root3 in GF(q) and root3^alpha<>1 then count:=count+1; fi;
 if root4 in GF(q) and root4^alpha<>1 then count:=count+1; fi;
 if root14 in GF(q) and root14^alpha<>1 then count:=count+1; fi;
 return count;
end;

Nbeta:=function(beta,q)
 # beta a character of GF(q^2)^x
 local root3,root4,root14,count;
 count := 0;
 root3:=cuberoots(q)[1];
 root4:=fourthroots(q)[1];
 root14:=fourteenthroots(q)[1];
 if root3 in GF(q) and root3^beta<>1 then count:=count+1; fi;
 if root4 in GF(q) and root4^beta<>1 then count:=count+1; fi;
 if root14 in GF(q) and root14^beta<>1 then count:=count+1; fi;
 return count;
end;
