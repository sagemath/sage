\\ Tpprog.g, V2.2---a program to compute conjectural slopes.
\\ By Kevin Buzzard (buzzard@ic.ac.uk)
\\ Please report any problems to buzzard@ic.ac.uk
\\
\\ Version history:
\\
\\ V2.2, 10/8/01, minor tinkering.
\\
\\ Version 2.1, 9/8/01. Complete rewrite motivated by desire to publish
\\ a paper.
\\
\\ Version history:
\\
\\ Version 1.2, 7/8/00 minor edit for p=2 because I guess
\\ that kmin should be 6 not 4, judging by N=5 case; also I
\\ redefined "ordinary" for p=2 by saying "as many 0s as possible
\\ in wt 4 and then the rest 1s".

\\ Version 1.1, 26/7/00 (post France, just found it in the mess.)
\\                      I don't yet know what it does but we live in hope!
\\
\\ Uses: function c0() from genusn.g
\\       function DimensionCuspForms() from DimensionSk.g


print();
print("tpslopes(p,N,kmax) returns a vector of length kmax, whose k'th entry");
print("(1<=k<=kmax) is the conjectural sequence of valuations of eigenvalues");
print("of T_p on forms of level N, weight k, and trivial character.");
print("The conjecture is only valid if p doesn't divide N and if p is");
print("Gamma_0(N)-regular.");
print()
print("Example of use:");
print("(12:48) gp > c=tpslopes(2,1,100);");
print("(12:48) gp > c[100]");
print("%2 = [3, 9, 9, 17, 17, 26, 26, 32]");
print()
print("Hence I conjecture that the 2-adic valuations of the eigenvalues of");
print("T_2 on cusp forms of level 1 and weight 100 are [3,9,9,...,32],");
print("which indeed they are, as one can verify by an explicit computation.");


tpslopes(p,N,kmax)=
{
  local(t,s,kmin,a,b,c,alg,k1,k2,k3,v1,v2,v3,B,e,S,s2,
        e2,epsilon,V1,V,V1tmp,w,w0,kmaxtemp);

  kmaxtemp=if(p==2,5,p+2);
  if(kmax>kmaxtemp,kmaxtemp=kmax);

  t=vector(kmaxtemp,i,[]);

  \\ initialise t vector in regular case.

  t[2]=vector(S0(N*p,2)-S0(N,2),i,0);
  if(p==2,
    t[4]=concat(t[2],vector(S0(N,4)-length(t[2]),i,1));
    kmin=6
  ,
    forstep(k=4,p+1,2,
      t[k]=vector(S0(N,k),i,0)
    );
    kmin=p+3
  );

  s=t;
  s[2]=vector(S0(N,2),i,0);

  \\ Are we done yet?
  if(kmax<kmin,
    vector(kmax,i,s[i])
  ,
    \\ no we're not!

    m=c0(N); \\ for non-trivial (even) character at *N*, I seem to have once
             \\ thought that the analogue should be dim(mod forms wt 4)
             \\ minus dim(cusp forms wt 4).

    forstep(k=kmin,kmax,2,

      a=1;while(p^a<k-1,a++);a=a-1;

      b=1;while(p^a*b<k-1,b++);b=b-1;

      c=1+floor((k-2-p^a*b)/(p^(a-1)));

      if(b+c<=p-1,alg=1,
        if(b<p-1,alg=2,
          alg=3
        )
      );

      if(alg==1,
        k1=k-b*(p-1)*p^(a-1);
        k2=k-(b-1)*(p-1)*p^(a-1)-2*(b+c-1)*p^(a-1);
        v1=t[k1];
        v2=t[k2];
        B=p^a*b+p^(a-1)*(c-1)+1;
        e=k-B;
        epsilon=DirichletCharacter(p,[(B-1)%(p-1)]);
        S=1+DimensionCuspForms(charmult(epsilon,TrivialCharacter(N)),1+e);

        if(length(v1)>=S-1,
          V1=vector(S-1,i,v1[i])
        ,
          V1=concat(v1,vector(S-1-length(v1),i,e-v2[S-length(v1)-i]))
        );

        V=concat(V1,vector(m,i,e));
      ,
        if(alg==2,
          k1=k-(b+1)*p^(a-1)*(p-1);
          k2=k-p^(a-1)*(p-1);
          v1=t[k1];
          v2=t[k2];
          B=(b+1)*p^(a-1)*(p-1)+1;
          e=k-B;
          S=1+S0(N*p,1+e);
          s2=(S-1)\2;
          e2=e\2;

          if(length(v1)>=S-1,
            V1=vector(S-1,i,v1[i])
          ,
            if(S-1<=2*length(v1),
              V1=concat(v1,vector(S-1-length(v1),i,e-v1[S-length(v1)-i]))
            ,
              w=vector(s2-length(v1),i,v2[length(v1)+i]);
              V1=concat(v1,w);
              if(S%2==0,V1=concat(V1,[e2]));
              V1=concat(V1,vector(length(w),i,e-1-w[length(w)+1-i]));
              V1=concat(V1,vector(length(v1),i,e-v1[length(v1)+1-i]));

            )
          );
          V=concat(V1,vector(m-(e==1),i,e));
        ,
          \\ alg=3.
          k1=k-p^a*(p-1);
          k2=k-p^(a-1)*(p-1);
          v1=t[k1];
          v2=t[k2];
          B=p^a*(p-1)+1;
          e=k-B;
          S=1+S0(N*p,1+e);
          s2=(S-1)\2;
          e2=e\2;
          if(length(v1)>=S-1,
            V1=vector(S-1,i,v1[i])
          ,
            if(S-1<=2*length(v1),
              V1=concat(v1,vector(S-1-length(v1),i,e-v1[S-length(v1)-i]))
            ,
              \\ case I don't understand :-*
              w0=vector(s2-length(v1),i,v2[length(v1)+i]);
              w=vector(s2-length(v1),i,if(w0[i]<e2,w0[i]+1,e2));
              V1=concat(v1,w);
              if(S%2==0,V1=concat(V1,[e2]));
              V1=concat(V1,vector(length(w),i,e-1-w[length(w)+1-i]));
              V1=concat(V1,vector(length(v1),i,e-v1[length(v1)+1-i]));
            )
          );
        V=concat(V1,vector(m-(e==1),i,e));
        )
      );
      if(length(V)>=S0(N,k),
        t[k]=vector(S0(N,k),i,V[i])
      ,
        k3=2*B-k;
        v3=t[k3];
        t[k]=concat(V,vector(S0(N,k)-length(V),i,e+v3[i]))
      );

      s[k]=t[k];

    );
  s
  )
}
