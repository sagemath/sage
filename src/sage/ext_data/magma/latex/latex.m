// Latex printing for MAGMA objects.

/***************************************************************

       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
                     2006 Jennifer Balakrishnan <jenb@mit.edu>

  Distributed under the terms of the GNU General Public License (GPL)

 ***************************************************************/

/*
This converts MAGMA output to LaTeX. It's a work-in-progress that
    currently handles a few basic types, matrices, polynomials,
    power series, binary quadratic forms, elements of number fields,
    finite fields, certain p-adic rings/fields, points, and elliptic curves.
*/

intrinsic Latex(x::RngIntElt) -> MonStgElt
{}
    return Sprint(x);
end intrinsic;

intrinsic Latex(x::FldReElt) -> MonStgElt
{}
    return Sprint(x);
end intrinsic;

intrinsic Latex(x::FldRatElt) -> MonStgElt
{}
    if Denominator(x) eq 1 then
        return Latex(Numerator(x));
    end if;
    return Sprintf("\\frac{%o}{%o}", Numerator(x), Denominator(x));
end intrinsic;

Letters:={@
"$.1",
"alpha",
"beta",
"gamma",
"delta",
"epsilon",
"varepsilon",
"zeta",
"eta",
"theta",
"theta",
"vartheta",
"iota",
"kappa",
"lambda",
"mu",
"nu",
"xi",
"pi",
"varpi",
"rho",
"varrho",
"sigma",
"varsigma",
"tau",
"upsilon",
"phi",
"varphi",
"chi",
"psi",
"omega",
"Gamma",
"Delta",
"Theta",
"Lambda",
"Xi",
"Pi",
"Sigma",
"Upsilon",
"Phi",
"Psi" @};


intrinsic Latex(x::RngPadElt) -> MonStgElt
{}
   z := Integers()!x;
   p := Prime(Parent(x));
   l := p^(Degree(Parent(x)));
   i := 0;
   j :=AbsolutePrecision(x);
   s := "";
   if IsPrime(l) then

  while z ne 0 and i eq 0 do
   c := z mod p;
        if c ne 0 then
           s *:= Sprintf("%o+", c);
        end if;
   z := z div p;
   i := i + 1;
   end while;

   while z ne 0 and i eq 1 do
    c := z mod p;
         if c ne 0 then
           if c ne 1 then
         s *:= Sprintf("%o\\cdot{}", c);
   end if;
   s *:= Sprintf("%o^{%o} + ", p, i);
   end if;
   z := z div p;
   i := i + 1;
   end while;

   while z ne 0 and i lt Precision(Parent(x)) do
      c := z mod p;
        if c ne 0 then
           if c ne 1 then
        s *:= Sprintf("%o\\cdot{}", c);
   end if;
   s *:= Sprintf("%o^{%o} + ", p, i);
   end if;
   z := z div p;
   i := i + 1;
   end while;

   else return Sprintf("\\mbox{\\rm %o}", x);
   end if;
   s *:= Sprintf("O(%o^{%o})", p, j);
   return s;
end intrinsic;

intrinsic Latex(x::FldPadElt) -> MonStgElt
{}
   v := Valuation(x);
   z := RationalField()!x;
   p := Prime(Parent(x));
   l := p^(Degree(Parent(x)));
   i := 0;
   j :=AbsolutePrecision(x);
   s := "";
   if IsPrime(l) then

        z:=Numerator(z);
         while z ne 0 and i eq 0 do
            c := z mod p;
            if c ne 0 then
               if c ne 1 then
                  s *:= Sprintf("%o", c);
               end if;
               if i+v ne 0 then
                  s*:= Sprintf("\\cdot{}%o^{%o} + ",p, i+v);
               else if c ne 0 and c ne 1 then
                  s*:= Sprintf("+ ");
               else if c eq 1 then
                  s*:= Sprintf("1 + ");
               end if;
               end if;
               end if;
            end if;
            z := z div p;
            i := i + 1;
         end while;

         while z ne 0 and i eq 1 do
            c := z mod p;
            if c ne 0 then
               if c ne 1 then
                   s *:= Sprintf("%o", c);
               end if;
               if i+v ne 0 then
                  s *:= Sprintf("\\cdot{}%o^{%o} + ", p, i+v);
               else if c ne 0 and c ne 1 then
                 s*:=Sprintf("+ ");
               else if c eq 1 then
                 s*:=Sprintf("1 + ");
               end if;
               end if;
               end if;
            end if;
            z := z div p;
            i := i + 1;
         end while;

         while z ne 0 and i lt Precision(Parent(x)) do
            c := z mod p;
            if c ne 0 then
               if c ne 1 then
                  s *:= Sprintf("%o", c);
               end if;
               if i+v ne 0 then
                  s *:= Sprintf("\\cdot{}%o^{%o} + ", p, i+v);
               else if c ne 0 and c ne 1 then
                  s*:= Sprintf("+ ");
               else if c eq 1 then
                  s*:= Sprintf("1 + ");
               end if;
               end if;
               end if;
          end if;
            z := z div p;
            i := i + 1;
         end while;
        else return Sprintf("\\text{%o}",x);

   end if;
   s *:= Sprintf("O(%o^{%o})", p, j);
   return s;
end intrinsic;


function S(x)
   if x gt 0 then
       return "+";
   end if;
   if x lt 0 then
       return "-";
   end if;
   if x eq 0 then
       return "";
   end if;
end function;

function Abs(x)
   return Latex(AbsoluteValue(x));
end function;

intrinsic Latex(f::RngUPolElt) -> MonStgElt
{}
    if Sprintf("%o",Parent(f).1) in Letters then
        x:=Sprintf("\\%o",Parent(f).1);
    else  x:=Sprintf("%o",Parent(f).1);
    end if;
    v := Eltseq(f);
    if v[1] ne 0 then
      s :=Abs(v[1]);
    else s:= "";
    end if;

    if s eq "" then

        if AbsoluteValue(v[2]) eq 1 then
            s:= Sprintf("%o",x) * s;
        else if v[2] eq 0 then
           s := S(v[1])*s;
        else if v[2] ne 0 then
            s := Abs(v[2]) * Sprintf("%o",x) *S(v[1])* s;
        end if;
        end if;
        end if;


    else   if AbsoluteValue(v[2]) eq 1 then
               s:=  Sprintf("%o",x) * S(v[1])* s;
           else if v[2] eq 0 then
               s := S(v[1])*s;
           else if v[2] ne 0 then
                s := Abs(v[2]) * Sprintf("%o",x) *S(v[1])* s;
           end if;
           end if;
           end if;

    end if;

    for i in [3..#v-1] do
      if s eq "" then

        if AbsoluteValue(v[i]) eq 1 then
            s:= Sprintf("%o",x)* Sprintf("^{%o}", i-1) * S(v[i-1]) * s;
        else

        if v[i] eq 0 then
           s := S(v[i-1])*s;

        else

        if v[i] ne 0 then
            s :=  Abs(v[i]) * Sprintf("%o",x) * Sprintf("^{%o}", i-1) * S(v[i-1]) * s;
        end if;
        end if;
        end if;

   else
        if AbsoluteValue(v[i]) eq 1 then
            s:=  Sprintf("%o",x) * Sprintf("^{%o}", i-1) * S(v[i-1]) * s;
        else

        if v[i] eq 0 then
           s := S(v[i-1])*s;

        else

        if v[i] ne 0 then
            s :=  Abs(v[i]) * Sprintf("%o",x)*Sprintf("^{%o}", i-1) *S(v[i-1])* s;
        end if;
        end if;
        end if;
end if;

    end for;

   if #v eq 2 then
       if S(v[2]) eq "-" then
          s := S(v[2])*s;
       end if;
   end if;

   if #v gt 2 then
        if AbsoluteValue(v[#v]) eq 1 then
             if S(v[#v]) eq "-" then
                 s:= Sprintf("-%o",x)* Sprintf("^{%o}", #v-1)*S(v[#v-1])*s;
             else s:=Sprintf("%o",x)* Sprintf("^{%o}", #v-1)*S(v[#v-1])*s;
             end if;
        else

        if v[#v] ne 0 then
             if S(v[#v]) eq "-" then
                 s := S(v[#v])*Abs(v[#v])*Sprintf("x")
                        *Sprintf("^{%o}", #v-1) *S(v[#v-1])* s;
             else s := Abs(v[#v])*Sprintf("%o",x)*Sprintf("^{%o}", #v-1) *S(v[#v-1])* s;
             end if;
        end if;
        end if;
  end if;

  return s;
end intrinsic;


intrinsic Latex(f::RngSerElt) -> MonStgElt
{}
  s:="";
  n:=AbsolutePrecision(f);
  d:=Degree(LeadingTerm(f));
  v:=ElementToSequence(f);
  m:=#v;
  if Sprintf("%o",Parent(f).1) in Letters then
        zn:=Sprintf("\\%o",Parent(f).1);
  else  zn:=Sprintf("%o",Parent(f).1);
  end if;

  if d eq 0 then
     if v[1] ne 0 then
        s:=s*Latex(v[1]);
     end if;

     if v[2] ne 0 then
         if s eq "" then
            if AbsoluteValue(v[2]) eq 1 then
                if S(v[2]) eq "-" then
                   s:= s*S(v[2])*Sprintf("%o",zn);
                else s:= s*Sprintf("%o",zn);
                end if;
            else if S(v[2]) eq "-" then
                   s:= s*S(v[2])*Abs(v[2])*Sprintf("%o",zn);
                else s:= s*Abs(v[2])*Sprintf("%o",zn);
                end if;
             end if;
         else  if AbsoluteValue(v[2]) eq 1 then
                   s:= s*S(v[2])*Sprintf("%o",zn);
            else   s:= s*S(v[2])*Abs(v[2])*Sprintf("%o",zn);
            end if;
        end if;
     end if;

     for i in [3..m] do
      if v[i] ne 0 then
         if s eq "" then
            if AbsoluteValue(v[i]) eq 1 then
                if S(v[i]) eq "-" then
                   s:= s*S(v[i])*Sprintf("%o^{%o}",zn,i-1);
                else s:= s*Sprintf("%o^{%o}",zn,i-1);
                end if;
            else if S(v[i]) eq "-" then
                   s:= s*S(v[i])*Abs(v[i])*Sprintf("%o^{%o}",zn,i-1);
                else s:= s*Abs(v[i])*Sprintf("%o^{%o-1}",zn,i-1);
                end if;
             end if;
         else  if AbsoluteValue(v[i]) eq 1 then
                   s:= s*S(v[i])*Sprintf("%o^{%o}",zn,i-1);
            else   s:= s*S(v[i])*Abs(v[i])*Sprintf("%o^{%o}",zn,i-1);
            end if;
        end if;
      end if;
   end for;

else if d eq 1 then

    if v[1] ne 0 then
         if s eq "" then
            if AbsoluteValue(v[1]) eq 1 then
                if S(v[1]) eq "-" then
                   s:= s*S(v[1])*Sprintf("%o",zn);
                else s:= s*Sprintf("%o",zn);
                end if;
            else if S(v[1]) eq "-" then
                   s:= s*S(v[1])*Abs(v[1])*Sprintf("%o",zn);
                else s:= s*Abs(v[1])*Sprintf("%o",zn);
                end if;
             end if;
         else  if AbsoluteValue(v[1]) eq 1 then
                   s:= s*S(v[1])*Sprintf("%o",zn);
            else   s:= s*S(v[1])*Abs(v[1])*Sprintf("%o",zn);
            end if;
        end if;
     end if;

     for i in [2..m] do
      if v[i] ne 0 then
         if s eq "" then
            if AbsoluteValue(v[i]) eq 1 then
                if S(v[i]) eq "-" then
                   s:= s*S(v[i])*Sprintf("%o^{%o}",zn,i);
                else s:= s*Sprintf("%o^{%o}",zn,i);
                end if;
            else if S(v[i]) eq "-" then
                   s:= s*S(v[i])*Abs(v[i])*Sprintf("%o^{%o}",zn,i);
                else s:= s*Abs(v[i])*Sprintf("%o^{%o}",zn,i);
                end if;
          end if;
      else  if AbsoluteValue(v[i]) eq 1 then
                s:= s*S(v[i])*Sprintf("%o^{%o}",zn,i);
         else   s:= s*S(v[i])*Abs(v[i])*Sprintf("%o^{%o}",zn,i);
         end if;
     end if;
   end if;
  end for;
else for i in [1..m] do
   if v[i] ne 0 then
      if s eq "" then
         if AbsoluteValue(v[i]) eq 1 then
             if S(v[i]) eq "-" then
                s:= s*S(v[i])*Sprintf("%o^{%o}",zn,d+i-1);
             else s:= s*Sprintf("%o^{%o}",zn,d+i-1);
             end if;
         else if S(v[i]) eq "-" then
                s:= s*S(v[i])*Abs(v[i])*Sprintf("%o^{%o}",zn,d+i-1);
             else s:= s*Abs(v[i])*Sprintf("%o^{%o}",zn,d+i-1);
             end if;
          end if;
      else  if AbsoluteValue(v[i]) eq 1 then
                s:= s*S(v[i])*Sprintf("%o^{%o}",zn,d+i-1);
         else   s:= s*S(v[i])*Abs(v[i])*Sprintf("%o^{%o}",zn,d+i-1);
         end if;
     end if;
   end if;
end for;
end if;
end if;

s:=s*Sprintf("+O(%o^{%o})",zn,n);
  return s;
end intrinsic;

intrinsic Latex(E::CrvEll) -> MonStgElt
{}
  v:=aInvariants(E);

  s:="y^2";

  if v[1] ne 0 then
     if AbsoluteValue(v[1]) eq 1 then
         s:=s*S(v[1])*Sprintf("xy");
     else
     s:= s*S(v[1])*Abs(v[1])*Sprintf("xy");
     end if;
  end if;

  if v[3] ne 0 then
     if AbsoluteValue(v[3]) eq 1 then
        s:=s*S(v[3])*Sprintf("y");
     else
        s:=s*S(v[3])*Abs(v[3])*Sprintf("y");
     end if;
  end if;

  s:=s*Sprintf("=x^3");

  if v[2] ne 0 then
    if AbsoluteValue(v[2])  eq 1 then
       s:= s*S(v[2])*Sprintf("x^2");
    else
       s:=s*S(v[2])*Abs(v[2])*Sprintf("x^2");
    end if;
  end if;

  if v[4] ne 0 then
    if AbsoluteValue(v[4]) eq 1 then
       s:=s*S(v[4])*Sprintf("x");
    else
       s:=s*S(v[4])*Abs(v[4])*Sprintf("x");
    end if;
  end if;

  if v[5] ne 0 then
     s:=s*S(v[5])*Abs(v[5]);
  end if;

return s;
end intrinsic;

intrinsic Latex(f::FldFunRatElt) -> MonStgElt
{}
    if Denominator(f) eq 1 then
        return Latex(Numerator(f));
    end if;
    return Sprintf("\\frac{%o}{%o}", Latex(Numerator(f)), Latex(Denominator(f)));
end intrinsic;

intrinsic Latex(f::QuadBinElt) -> MonStgElt
{}
  s:="";
  if AbsoluteValue(f[1]) eq 1 then
         if S(f[1]) eq "-" then
              s:= s*Sprintf("-x^2");
         else s:= s*Sprintf("x^2");
         end if;
  else
         if S(f[1]) eq "-" then
              s:= s*S(f[1])*Abs(f[1])*Sprintf("x^2");
         else s:= s*Abs(f[1])*Sprintf("x^2");
         end if;
  end if;

  if f[2] ne 0 then
     if AbsoluteValue(f[2]) eq 1 then
        s:=s*S(f[2])*Sprintf("xy");
     else
        s:=s*S(f[2])*Abs(f[2])*Sprintf("xy");
     end if;
  end if;

  if AbsoluteValue(f[3]) eq 1 then
         s:=s*S(f[3])*Sprintf("y^2");
     else
     s:= s*S(f[3])*Abs(f[3])*Sprintf("y^2");
  end if;
  return s;
end intrinsic;


intrinsic Latex(M::Mtrx) -> MonStgElt
{}
    m:=NumberOfRows(M);
    n:=NumberOfColumns(M);
    s:=Sprintf("\\left(\\begin{array}{");
    for i in [1..n] do
       s := s * Sprintf("c");
    end for;
    s:= s * Sprintf("}");
    for i in [1..m-1] do
        for j in [1..n-1] do
           s := s * Latex(M[i,j]) * Sprintf("&");
        end for;
        s := s * Latex(M[i,n]) * Sprintf("\\\\");
    end for;

    for j in [1..n-1] do
           s := s * Latex(M[m,j]) * Sprintf("&");
        end for;
    s := s * Latex(M[m,n]) *  "\\end{array}\\right)";
    return s;
end intrinsic;

intrinsic Latex(P::PtEll) -> MonStgElt
{}
   return Sprintf("(%o,%o)",Latex(P[1]),Latex(P[2]));

end intrinsic;

intrinsic Latex(P::Pt) -> MonStgElt
{}
  n:=#Coordinates(P);
  s:="(";
  for i in [1..n-1] do
      s:=s*Latex(P[i])*Sprintf(",");
  end for;
  s:=s*Latex(P[n])*Sprintf(")");
  return s;
end intrinsic;

intrinsic Latex(a::FldQuadElt) -> MonStgElt
{}
  s:="";
  v:=ElementToSequence(a);
  D:=Discriminant(Parent(a))/4;

  if v[1] ne 0 then
      s:=s*Latex(v[1]);
  end if;

  if v[2] ne 0 then
      if s eq "" then
         if AbsoluteValue(v[2]) eq 1 then
             if S(v[2]) eq "-" then
                s:= s*S(v[2])*Sprintf("\\sqrt{%o}",D);
             else s:= s*Sprintf("\\sqrt{%o}",D);
             end if;
         else if S(v[2]) eq "-" then
                s:= s*S(v[2])*Abs(v[2])*Sprintf("\\sqrt{%o}",D);
             else s:= s*Abs(v[2])*Sprintf("\\sqrt{%o}",D);
             end if;
          end if;
      else  if AbsoluteValue(v[2]) eq 1 then
                s:= s*S(v[2])*Sprintf("\\sqrt{%o}",D);
         else   s:= s*S(v[2])*Abs(v[2])*Sprintf("\\sqrt{%o}",D);
         end if;
     end if;
  end if;
  return s;
end intrinsic;

intrinsic Latex(a::FldCycElt) -> MonStgElt
{}
  s:="";
  v:=ElementToSequence(a);
  n:=CyclotomicOrder(Parent(a));
  zn:=Sprintf("\\zeta_{%o}",n);
  m:=Degree(Parent(a));

  if v[1] ne 0 then
      s:=s*Latex(v[1]);
  end if;

  if v[2] ne 0 then
      if s eq "" then
         if AbsoluteValue(v[2]) eq 1 then
             if S(v[2]) eq "-" then
                s:= s*S(v[2])*Sprintf("%o",zn);
             else s:= s*Sprintf("%o",zn);
             end if;
         else if S(v[2]) eq "-" then
                s:= s*S(v[2])*Abs(v[2])*Sprintf("%o",zn);
             else s:= s*Abs(v[2])*Sprintf("%o",zn);
             end if;
          end if;
      else  if AbsoluteValue(v[2]) eq 1 then
                s:= s*S(v[2])*Sprintf("%o",zn);
         else   s:= s*S(v[2])*Abs(v[2])*Sprintf("%o",zn);
         end if;
     end if;
  end if;

  for i in [3..m-1] do
   if v[i] ne 0 then
      if s eq "" then
         if AbsoluteValue(v[i]) eq 1 then
             if S(v[i]) eq "-" then
                s:= s*S(v[i])*Sprintf("%o^{%o}",zn,i-1);
             else s:= s*Sprintf("%o^{%o}",zn,i-1);
             end if;
         else if S(v[i]) eq "-" then
                s:= s*S(v[i])*Abs(v[i])*Sprintf("%o^{%o}",zn,i-1);
             else s:= s*Abs(v[i])*Sprintf("%o^{%o-1}",zn,i-1);
             end if;
          end if;
      else  if AbsoluteValue(v[i]) eq 1 then
                s:= s*S(v[i])*Sprintf("%o^{%o}",zn,i-1);
         else   s:= s*S(v[i])*Abs(v[i])*Sprintf("%o^{%o}",zn,i-1);
         end if;
     end if;
   end if;
  end for;
  return s;
end intrinsic;


intrinsic Latex(a::FldNumElt) -> MonStgElt
{}
  s:="";
  v:=ElementToSequence(a);
  n:=Degree(Parent(a));

  if Sprintf("%o",Parent(a).1) in Letters then
     zn:=Sprintf("\\%o",Parent(a).1);
  else  zn:=Sprintf("%o",Parent(a).1);
  end if;

  if v[1] ne 0 then
      s:=s*Latex(v[1]);
  end if;

  if v[2] ne 0 then
      if s eq "" then
         if AbsoluteValue(v[2]) eq 1 then
             if S(v[2]) eq "-" then
                s:= s*S(v[2])*Sprintf("%o",zn);
             else s:= s*Sprintf("%o",zn);
             end if;
         else if S(v[2]) eq "-" then
                s:= s*S(v[2])*Abs(v[2])*Sprintf("%o",zn);
             else s:= s*Abs(v[2])*Sprintf("%o",zn);
             end if;
          end if;
      else  if AbsoluteValue(v[2]) eq 1 then
                s:= s*S(v[2])*Sprintf("%o",zn);
         else   s:= s*S(v[2])*Abs(v[2])*Sprintf("%o",zn);
         end if;
     end if;
  end if;

  for i in [3..n] do
   if v[i] ne 0 then
      if s eq "" then
         if AbsoluteValue(v[i]) eq 1 then
             if S(v[i]) eq "-" then
                s:= s*S(v[i])*Sprintf("%o^{%o}",zn,i-1);
             else s:= s*Sprintf("%o^{%o}",zn,i-1);
             end if;
         else if S(v[i]) eq "-" then
                s:= s*S(v[i])*Abs(v[i])*Sprintf("%o^{%o}",zn,i-1);
             else s:= s*Abs(v[i])*Sprintf("%o^{%o-1}",zn,i-1);
             end if;
          end if;
      else  if AbsoluteValue(v[i]) eq 1 then
                s:= s*S(v[i])*Sprintf("%o^{%o}",zn,i-1);
         else   s:= s*S(v[i])*Abs(v[i])*Sprintf("%o^{%o}",zn,i-1);
         end if;
     end if;
   end if;
  end for;
  return s;
end intrinsic;

intrinsic Latex(K::FldFin) -> MonStgElt
{}
   p := Characteristic(K);
   n := Degree(K);
   s := Sprintf("\\mathbf{F}_{{%o}",p);
   if n gt 1 then
      s *:= Sprintf("^{%o}}",n);
   else
      s *:= "}";
   end if;
   return s;
end intrinsic;

intrinsic Latex(a::FldFinElt) -> MonStgElt
{}
  s:="";
  v:=ElementToSequence(a);
  m:=#v;
  n:=#(Parent(a));
  if IsPrime(n) then
    return Sprint(a);
  else
      if Sprintf("%o",Parent(a).1) in Letters then
          zn:=Sprintf("\\%o",Parent(a).1);
      else  zn:=Sprintf("%o",Parent(a).1);
      end if;
  end if;

  if v[1] ne (Parent(a))!0 then
      s:=s*Sprintf("%o", v[1]);
  end if;

  if v[2] ne (Parent(a))!0 then
      if s eq "" then
         if v[2] eq Parent(a)!1 then
               s:= s*Sprintf("%o",zn);
          else
               s:= s*Sprintf("%o%o",v[2],zn);
         end if;
      else
      if v[2] eq Parent(a)!1 then
         s:= s*Sprintf("+%o",zn);
      else
         s:= s*Sprintf("+%o%o",v[2],zn);
      end if;
      end if;
  end if;

  for i in [3..m] do
   if v[i] ne Parent(a)!0 then
      if s eq "" then
         if v[i] eq Parent(a)!1 then
             s:= s*Sprintf("%o^{%o}",zn,i-1);
         else
             s:= s*Sprintf("%o%o^{%o-1}",v[i],zn,i-1);
         end if;
      else
      if v[i] eq Parent(a)!1 then
                s:= s*Sprintf("+%o^{%o}",zn,i-1);
         else   s:= s*Sprintf("+%o%o^{%o}",v[i],zn,i-1);
      end if;
   end if;
   end if;
   end for;
   return s;
end intrinsic;

intrinsic Latex(x::.) -> MonStgElt
{}
    return Sprintf("\\mbox{\\rm %o}", x);
end intrinsic;


/*

These are the MAGMA types:
     (*) indicates done

AlgAssElt
AlgAssVElt
AlgBasElt
AlgChtrElt
AlgClffElt
AlgExtElt
AlgFPElt
AlgFPEltOld
AlgFPGElt
AlgFPLieElt
AlgFinDElt
AlgFrElt
AlgGenElt
AlgGrpElt
AlgHckElt
AlgIUEElt
AlgInfDElt
AlgLieElt
AlgMatElt
AlgMatLieElt
AlgMatVElt
AlgPBWElt
AlgQUEElt
AlgQuatElt
AlgQuatOrdElt
AlgSymElt
AlgUEElt
AutCrvEll
CopElt
(*)CrvEll
DiffCrvElt
DiffFunElt
DivCrvElt
DivFunElt
DivNumElt
ExtReElt
FldACElt
FldAlgElt
FldComElt
(*)FldCycElt
FldElt
(*)FldFin
(*)FldFinElt
FldFracElt
FldFunElt
FldFunFracSchElt
FldFunFracSchEltOld
FldFunGElt
FldFunOrdElt
(*)FldFunRatElt
FldFunRatMElt
FldFunRatUElt
(*)FldNumElt
FldNumGElt
FldOrdElt
FldPadElt
FldPrElt
(*)FldQuadElt
(*)FldRatElt
(*)FldReElt
FldResLst
FldResLstElt
FldTimeElt
GenMPolBElt
GenMPolBGElt
GenMPolElt
GenMPolGElt
GenMPolResElt
GrpAbElt
GrpAbGenElt
GrpAtcElt
GrpAutoElt
GrpBBElt
GrpBrdElt
GrpCaygElt
GrpDrchElt
GrpDrchEltNew
GrpElt
GrpFPCosElt
GrpFPCoxElt
GrpFPDcosElt
GrpFPElt
GrpFrmlElt
GrpGPCElt
GrpGenElt
GrpLieAutoElt
GrpLieElt
GrpMatElt
GrpMatProjElt
GrpPCElt
GrpPSL2
GrpPSL2Elt
GrpPermDcosElt
GrpPermElt
GrpPermLcosElt
GrpPermRcosElt
GrpRWSElt
GrpSLPElt
HilbSpcElt
Infty
IsoCrvEll
LatElt
List
MPolElt
MapCrvEll
MapIsoCrvEll
ModAbVarElt
ModAlgBasElt
ModAlgElt
ModAltElt
ModBrdtElt
ModCycElt
ModDedElt
ModEDElt
ModExtElt
ModFldElt
ModFrmElt
ModGrpElt
ModHgnElt
ModHrmElt
ModLatElt
ModMPolElt
ModMatFldElt
ModMatGrpElt
ModMatRngElt
ModRngElt
ModRngMPolRedElt
ModSSElt
ModSymElt
ModThetaElt
ModTupAlgElt
ModTupFldElt
ModTupRngElt
MonAbElt
MonFPElt
MonOrdElt
MonPlcElt
MonRWSElt
MonStgElt
MonStgGenElt
(*)Mtrx
MtrxSpcElt
OFldComElt
OFldReElt
PicCrvElt
PicHypSngElt
PlcCrvElt
PlcFunElt
PlcNumElt
(*)Pt
(*)PtEll
PtGrp
PtHyp
(*)QuadBinElt
Rec
RecField
RecFrmt
Ref
RegExp
RegExpAlg
Rel
RelElt
RngCycElt
RngDiffElt
RngDiffOpElt
RngElt
RngFracElt
RngFrmElt
RngFunFracElt
RngFunFracSchElt
RngFunFracSchEltOld
RngFunFracSchOld
RngFunGElt
RngFunOrdElt
RngFunOrdGElt
RngFunOrdIdl
RngGalElt
RngHckElt
RngHckIdl
(*)RngIntElt
RngIntResElt
RngMPolElt
RngMPolResElt
RngMSerElt
RngOrdElt
RngOrdFracIdl
RngOrdIdl
RngOrdResElt
RngPadElt
RngPadResElt
RngPadResExtElt
RngPowLazElt
RngQuadElt
RngQuadFracIdl
RngQuadIdl
RngReSubElt
RngRelKElt
RngSLPolElt
(*)RngSerElt
RngSerLaurElt
RngSerPowElt
RngSerPuisElt
(*)RngUPolElt
RngUPolResElt
RngValElt
RngWittElt
SchGrpEll
SeqEnum
Set
SetCspElt
SetPtEll
SgpFPElt
SpcFldElt
SpcHypElt
SpcRngElt
SubFldLatElt
SubGrpLatElt
SubModLatElt
SymCrvEll
UnusedMapCrvEll

*/
