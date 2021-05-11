##############################################################
##
##           Helps compute the decomposition of
##           the ramification module and the equivariant
##           degree of a multiple of a simple orbit for
##           the modular curves X(p) with automorphism group
##           G = PSL(2,p), p a prime.
##
##           Using these computations and Borne's formula, the
##           G-module structure of the RR spaces of equivariant
##           divisors can be determined explicitly.
##
## Reference: D. Joyner and A. Ksir, "Modular representations
##            on some Riemann-Roch spaces of modular curves
##            $X(N)$, Computational Aspects of Algebraic Curves,
##            (Editor: T. Shaska) Lecture Notes in Computing,
##            WorldScientific, 2005.)
##
##
## 12-30-2004, wdj
###########################################################

## To read this in, use for example:
### Read("/home/wdj/gapfiles/curves/mod_crv_aut_gp3.gap");
## To log your results, use for example:
### LogTo("/home/wdj/gapfiles/mod_crv_aut_gp1.log");

ram_module_X:=function(p)
## p is a prime
## output is [m1,...,mn]
## where n = # conj classes of G=PSL(2,p)
## and mi = mult of pi_i in ram mod of modular curve
## X with AutGp(X) = G.
## Here Irr(G) = [pi_1,...,pi_n] (in that order).
##
local G,i,j,n,n0,H,G1,H_chars,CG,G_chars,w,m,theta,pi_theta_to_the;
  G:=PSL(2,p);
  H:=[];
  H_chars:=[];
  n0:=[ [ Z(p)^0, Z(p) ], [ Z(p)*0, Z(p)^0 ] ];
  CG:=ConjugacyClassesSubgroups(G);
  H[1]:=Representative(CG[2]); # size 2
  H[2]:=Representative(CG[3]); # size 3
  n :=Size(CG);
  for i in [1..n] do
    if Size(Representative(CG[i]))=p then
      H[3]:=Representative(CG[i]);
    fi;
  od;
## H[3]:=Group(n0);       # size p
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

#ram_module_X(5);
#[ 0, 3, 3, 4, 5 ]
#ram_module_X(7);
#[ 0, 4, 3, 6, 7, 8 ]
#ram_module_X(11);
#[ 0, 5, 6, 11, 10, 12, 12, 12 ]
#ram_module_X(13);
#[ 0, 7, 7, 13, 13, 13, 13, 14, 15 ]
#ram_module_X(17);
#[ 0, 9, 9, 18, 17, 17, 17, 18, 18, 19, 19 ]
#ram_module_X(19);
#[ 0, 9, 10, 20, 20, 19, 19, 20, 20, 21, 21, 21 ]
#ram_module_X(23);
#[ 0, 14, 11, 24, 24, 24, 23, 23, 25, 25, 25, 25, 25, 25 ]
#ram_module_X(29);
#[ 0, 16, 16, 30, 31, 31, 30, 30, 30, 30, 31, 32, 32, 32, 31, 31, 31 ]
#ram_module_X(31);
#[ 0, 18, 15, 33, 33, 33, 32, 32, 32, 32, 33, 34, 33, 33, 34, 34, 34, 34 ]
#ram_module_X(37);
#[ 0, 20, 20, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 40, 39, 41, 41, 41, 40, 40, 40 ]

ram_module_X_JK:=function(p)
## p is a prime
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
    if Size(Representative(CG[i]))=p then
      H[3]:=Representative(CG[i]);
    fi;
  od;
## H[3]:=Group(n0);       # size p
  for i in [1..Size(H)] do
   H_chars[i]:=Irr(H[i]);
  od;
  G_chars:=Irr(G);
  n:=Length(G_chars);
  A:=[];
  for i in [1..n] do
    pi:=G_chars[i];
    A[i]:=ScalarProduct(H_chars[1][1],RestrictedClassFunction(pi,H[1]));
  od;
  B:=[];
  for i in [1..n] do
    pi:=G_chars[i];
    B[i]:=ScalarProduct(H_chars[2][1],RestrictedClassFunction(pi,H[2]));
  od;
  C:=[];
  for i in [1..n] do
    pi:=G_chars[i];
    C[i]:=ScalarProduct(H_chars[3][1],RestrictedClassFunction(pi,H[3]));
  od;
  D:=[];
  for i in [1..n] do
    pi:=G_chars[i];
    D[i]:=DegreeOfCharacter(pi);
  od;
 return (1/2)*(3*D-A-B-C);
end;

#ram_module_X_JK(5);
#[ 0, 3, 3, 4, 5 ]
#ram_module_X_JK(7);
#[ 0, 7/2, 7/2, 6, 7, 8 ]
#ram_module_X_JK(11);
#[ 0, 11/2, 11/2, 11, 10, 12, 12, 12 ]
#ram_module_X_JK(13);
#[ 0, 7, 7, 13, 13, 13, 13, 15, 14 ]
#ram_module_X_JK(17);
#[ 0, 9, 9, 18, 17, 17, 17, 18, 18, 19, 19 ]
#ram_module_X_JK(19);
#[ 0, 19/2, 19/2, 20, 20, 19, 19, 20, 20, 21, 21, 21 ]
#ram_module_X_JK(23);
#[ 0, 25/2, 25/2, 24, 24, 24, 23, 23, 25, 25, 25, 25, 25, 25 ]
#ram_module_X_JK(29);
#[ 0, 16, 16, 30, 31, 31, 30, 30, 30, 30, 31, 32, 32, 32, 31, 31, 31 ]

pieces_ram_module_X:=function(p)
## p is a prime
## output is (m1,...,mn)
## where n = # conj classes of G=PSL(2,p)
## and mi = mult of pi_i in ram mod of G
##
local G,i,j,n,n0,H,G1,H_chars,CG,G_chars,w,m,theta,pi_theta_to_the;
  G:=PSL(2,p);
  H:=[];
  H_chars:=[];
  n0:=[ [ Z(p)^0, Z(p) ], [ Z(p)*0, Z(p)^0 ] ];
  CG:=ConjugacyClassesSubgroups(G);
  H[1]:=Representative(CG[2]); # size 2
  H[2]:=Representative(CG[3]); # size 3
  n :=Size(CG);
  for i in [1..n] do
    if Size(Representative(CG[i]))=p then
      H[3]:=Representative(CG[i]);
    fi;
  od;
## H[3]:=Group(n0);       # size p
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
    m[i]:=List(G_chars, pi->ScalarProduct(pi_theta_to_the[i],pi));
  od;
#Print("\n\n m = ",m,"\n\n");
return m;
end;

equiv_deg_module_X:=function(p,ii,r)
## p is a prime
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
  if Size(Representative(CG[i]))=p then
    H[3]:=Representative(CG[i]);
  fi;
od;
## H[3]:=Group(n0);       # size p
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
