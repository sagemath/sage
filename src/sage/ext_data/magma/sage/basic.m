// -*- magma -*-

function PreparseElts(R)
    if Type(R) eq RngInt then
       return false;
    end if;
    return true;
end function;

intrinsic Sage(X::.) -> MonStgElt, BoolElt
{Default way to convert a Magma object to Sage if we have not
written anything better.}
    return Sprintf("%o", X), true;
end intrinsic;

intrinsic Sage(X::SeqEnum) -> MonStgElt, BoolElt
{Convert an enumerated sequence to Sage.}
    Y := [Sage(z) : z in X];
    return Sprintf("%o", Y), true;
end intrinsic;

intrinsic Sage(X::Tup) -> MonStgElt, BoolElt
{Return a Magma Tuple as a Sage Tuple}
    if #X eq 0 then
        return "()", true;
    elif #X eq 1 then
        return Sprintf("(%o,)",Sage(X[1])),true;
    end if;
    r := Sprintf("%o",[Sage(x) : x in X]);
    return "(" cat Substring(r,2,#r-2) cat ")",true;
end intrinsic;

intrinsic Sage(X::SetEnum) -> MonStgElt, BoolElt
{Convert an enumerated set to Sage.}
    Y := [Sage(z) : z in X];
    return Sprintf("Set(%o)", Y), true;
end intrinsic;

intrinsic Sage(X::SetIndx) -> MonStgElt, BoolElt
{Convert an indexed set to Sage.
 WARNING: Sage does not have an analogue of indexed sets (yet!),
 so we just return a Python list.}
    Y := [Sage(z) : z in X];
    return Sprintf("%o", Y), true;
end intrinsic;

intrinsic Sage(X::SetMulti) -> MonStgElt, BoolElt
{Convert a multiset to Sage.
 WARNING: Sage does not have an analogue of multisets yet, so we return a Python list.}
    Y := [Sage(z) : z in X];
    return Sprintf("%o", Y), true;
end intrinsic;

intrinsic Sage(X::RngInt) -> MonStgElt, BoolElt
{Conver the ring of integers to Sage.}
    return "ZZ", false;
end intrinsic;

intrinsic Sage(X::FldRat) -> MonStgElt, BoolElt
{}
    return "QQ", false;
end intrinsic;

intrinsic Sage(X::RngIntElt) -> MonStgElt, BoolElt
{}
    return Sprintf("Integer('%h')", X), false;
end intrinsic;

/* Matrices */

function convert_matrix(X, preparse_entries)
    if preparse_entries then
        return Sprintf("matrix(%o, %o, %o, %o)", Sage(BaseRing(X)),
                    Nrows(X), Ncols(X), [Sage(y) : y in Eltseq(X)]);
    else
        return Sprintf("matrix(%o, %o, %o, %o)", Sage(BaseRing(X)),
                    Nrows(X), Ncols(X), Eltseq(X));
    end if;
end function;

intrinsic Sage(X::AlgMatElt) -> MonStgElt, BoolElt
{}
    pp := PreparseElts(BaseRing(X));
    return convert_matrix(X, pp), pp;
end intrinsic;

intrinsic Sage(X::ModMatRngElt) -> MonStgElt, BoolElt
{}
    pp := PreparseElts(BaseRing(X));
    return convert_matrix(X, pp), pp;
end intrinsic;


intrinsic SageCreateWithNames(X::., names::.) -> .
{Assign the given names to the object X, then return X.}
    AssignNames(~X, names);
    return X;
end intrinsic;

/* Finite fields */

intrinsic Sage(X::FldFin) -> MonStgElt, BoolElt
{}
  if IsPrimeField(X) then
    return Sprintf("GF(%o)", Characteristic(X)), false;
  else
    return Sprintf("GF(%o, '%o'.replace('$.', 'x').replace('.', ''), modulus=%o)", #X, X.1, Sage(DefiningPolynomial(X))), false;
  end if;
end intrinsic;

intrinsic Sage(X::FldFinElt) -> MonStgElt, BoolElt
{}
    P := Parent(X);
    if IsPrimeField(P) then
        return Sprintf("%o(%o)", Sage(P), Integers()!X), false;
    else
        return Sprintf("%o(%o)", Sage(Parent(X)), Sage(Polynomial(Eltseq(X)))), false;
    end if;
end intrinsic;

/* Finite quotients of ZZ */

intrinsic Sage(X::RngIntRes) -> MonStgElt, BoolElt
{}
  return Sprintf("Zmod(%o)", Characteristic(X)), false;
end intrinsic;

/* Approximate real and complex fields */

intrinsic Sage(X::FldRe) -> MonStgElt, BoolElt
{}
    return Sprintf("RealField(%o)", Precision(X : Bits := true)), false;
end intrinsic;

intrinsic Sage(X::FldCom) -> MonStgElt, BoolElt
{}
    return Sprintf("ComplexField(%o)", Precision(X : Bits := true)), false;
end intrinsic;

intrinsic Sage(X::FldReElt) -> MonStgElt, BoolElt
{}
    return Sprintf("%o(%o)", Sage(Parent(X)), X), true;
end intrinsic;

intrinsic Sage(X::FldComElt) -> MonStgElt, BoolElt
{}
  return Sprintf("%o([%o, %o])", Sage(Parent(X)), Sage(Real(X)), Sage(Imaginary(X))), true;
end intrinsic;

/* p-adic rings and fields */

intrinsic Sage(X::RngPad) -> MonStgElt, BoolElt
{p-adic rings, either free precision model or exact model}
    prec := Precision(X);
    if Type(prec) eq Infty then
      return Sprintf("Zp(%o, %o, 'relaxed')", Sage(Prime(X)), Sage(prec), false;
    else
      return Sprintf("Zp(%o, %o, 'capped-rel')", Sage(Prime(X)), Sage(prec)), false;
    end if;
end intrinsic;

intrinsic Sage(X::FldPad) -> MonStgElt, BoolElt
{p-adic fields, either free precision model or exact model}
    prec := Precision(X);
    if Type(prec) eq Infty then
      return Sprintf("Qp(%o, %o, 'relaxed')", Sage(Prime(X)), Sage(prec)), false;
    else
      return Sprintf("Qp(%o, %o, 'capped-rel')", Sage(Prime(X)), Sage(prec)), false;
    end if;
end intrinsic;

intrinsic Sage(X::RngPadRes) -> MonStgElt, BoolElt
{fixed precision model}
    return Sprintf("Zp(%o, %o, 'fixed-mod')", Sage(Prime(X)), Sage(Precision(X))), false;
end intrinsic;


/* Polynomials */

intrinsic SageNamesHelper(X::.) -> MonStgElt
{}
  /* XXX */
  i := NumberOfNames(X);
  if "$" in Sprint(X.i) then
    /* unnamed variables */
    return "(" * (&* [ Sprintf("'x%o', ", j) : j in [ 1..i ] ]) * ")";
  else
    /* named variables */
    return "(" * (&* [ Sprintf("'%o'.replace('.', ''), ", X.j) : j in [ 1..i ] ]) * ")";

end if;
end intrinsic;

intrinsic Sage(X::RngUPol) -> MonStgElt, BoolElt
{}
  txt := "PolynomialRing(%o, %o)";
  return Sprintf(txt, Sage(BaseRing(X)), SageNamesHelper(X)), false;
end intrinsic;

intrinsic Sage(X::RngUPolElt) -> MonStgElt, BoolElt
{}
  pp := PreparseElts(BaseRing(X));
  return Sprintf("%o(%o)", Sage(Parent(X)), Sage(Coefficients(X))), pp;
end intrinsic;

intrinsic Sage(X::RngMPol) -> MonStgElt, BoolElt
{}
  txt := "PolynomialRing(%o, %o)";
  return Sprintf(txt, Sage(BaseRing(X)), SageNamesHelper(X)), false;
end intrinsic;

intrinsic Sage(X::RngMPolElt) -> MonStgElt, BoolElt
{}
  Y := Sage([ < <e : e in Exponents(t)>, Coefficients(t)[1]> : t in Terms(X)]);
  R := Sage(Parent(X));
  return Sprintf("%o(dict(%o))",R,Y),true;
end intrinsic;

/* Number fields  */

intrinsic Sage(K::FldNum) -> MonStgElt, BoolElt
{}
  gens := GeneratorsSequence(K);
  if "$" in Sprint(gens[1]) then
    /* unnamed variables */
    names := "(" * (&* [ Sprintf("'a%o', ", j) : j in [ 1..#gens ] ]) * ")";
  else
    /* named variables */
    names := "(" * (&* [ Sprintf("'%o'.replace('.', ''), ", a) : a in gens]) * ")";
  end if;
  polynomials := DefiningPolynomial(K);
  return Sprintf("NumberField(%o, %o)", Sage(polynomials), names), false;
end intrinsic;

intrinsic Sage(A::FldNumElt) -> MonStgElt, BoolElt
{Converts a number field element to Sage.
 Only number fields generated by a single element over
 the base field are supported.
 It seems that FldNum is always represented by
 a power basis. Just in case, this function
 checks whether this is true.}
    K := Parent(A);
    gens := GeneratorsSequence(K);
    if #gens ne 1 then return Sprint(A), true; end if;
    gen := gens[1];
    deg := Degree(K);
    bas := Basis(K);
    for a in [0..deg-1] do
        if gen^a ne bas[a+1] then
            return Sprint(A), true;
        end if; end for;
    seq := Eltseq(A);
    return Sprintf("%o(%o)", Sage(K), Sage(seq)), true;
end intrinsic;

intrinsic Sage(O::RngOrd) -> MonStgElt, BoolElt
{Converts an order of a number field to sage.}
K:=NumberField(O);
if IsMaximal(O) then
    return Sprintf("%o.maximal_order()",Sage(K)), true;
end if;
B:=Basis(O);
seq := [K!B[i] : i in [1..#B]];
return Sprintf("%o.order(%o)", Sage(K),Sage(seq)), true;
end intrinsic;

intrinsic Sage(I::RngOrdIdl) -> MonStgElt, BoolElt
{Converts an ideal of a number field to sage.}
O:=Order(I);
K:=NumberField(O);
gens:=Generators(I);
seq := [K!gens[i] : i in [1..#gens]];
return Sprintf("%o.ideal(%o)", Sage(K),Sage(seq)), true;
end intrinsic;

/* Symmetric functions */

intrinsic Sage(X::AlgSym) -> MonStgElt, BoolElt
{}
if HasSchurBasis(X) then
  return Sprintf("SymmetricFunctions(%o).s()", Sage(BaseRing(X))), false;
elif HasHomogeneousBasis(X) then
  return Sprintf("SymmetricFunctions(%o).h()", Sage(BaseRing(X))), false;
elif HasElementaryBasis(X) then
  return Sprintf("SymmetricFunctions(%o).e()", Sage(BaseRing(X))), false;
elif HasPowerSumBasis(X) then
  return Sprintf("SymmetricFunctions(%o).p()", Sage(BaseRing(X))), false;
elif HasMonomialBasis(X) then
  return Sprintf("SymmetricFunctions(%o).m()", Sage(BaseRing(X))), false;
end if;
end intrinsic;

intrinsic Sage(X::AlgSymElt) -> MonStgElt, BoolElt
{}
PA := Parent(X);
SF := Sage(PA);
BR := Sage(BaseRing(PA));
parts, coeffs := Support(X);
dict := (&* [ Sprintf("Partition(%o):%o(%o),", Sage(parts[i]), BR, Sage(coeffs[i])) : i in [1..#parts] ]);
return Sprintf("%o._from_dict({%o})", SF, dict), false;
end intrinsic;

/* Elliptic curves */

intrinsic Sage(X::CrvEll) -> MonStgElt, BoolElt
{}
  as := aInvariants(X);
  return Sprintf("EllipticCurve(%o)", Sage(as)), true;
end intrinsic;

/* Hyperelliptic curves */

intrinsic Sage(X::CrvHyp) -> MonStgElt, BoolElt
{}
  f, g := HyperellipticPolynomials(X);
  return Sprintf("HyperellipticCurve(%o, %o)", Sage(f), Sage(g)), true;
end intrinsic;

/* Modules and vector spaces */

intrinsic Sage(X::ModTupRng) -> MonStgElt, BoolElt
{}
  if IsIdentity(InnerProductMatrix(X)) then
      return Sprintf("FreeModule(%o, %o)", Sage(BaseRing(X)), Sage(Rank(X))), true;
  else
      return Sprintf("FreeModule(%o, %o, inner_product_matrix=%o)", Sage(BaseRing(X)), Sage(Rank(X)), Sage(InnerProductMatrix(X))), true;
  end if;
end intrinsic;

intrinsic Sage(X::ModTupRngElt) -> MonStgElt, BoolElt
{}
    return Sprintf("%o(%o)", Sage(Parent(X)), Sage(ElementToSequence(X))), true;
end intrinsic;

/* Power series rings */

intrinsic Sage(X::RngSerPow) -> MonStgElt, BoolElt
{}
  txt := "PowerSeriesRing(%o, %o)";
  var := Sprintf("['%o']", X.1);
  return Sprintf(txt, Sage(BaseRing(X)), var), false;
end intrinsic;

intrinsic Sage(X::RngSerLaur) -> MonStgElt, BoolElt
{}
  txt := "LaurentSeriesRing(%o, %o)";
  var := Sprintf("['%o']", X.1);
  return Sprintf(txt, Sage(BaseRing(X)), var), false;
end intrinsic;

intrinsic Sage(X::RngSerPuis) -> MonStgElt, BoolElt
{}
  txt := "PuiseuxSeriesRing(%o, %o)";
  var := Sprintf("['%o']", X.1);
  return Sprintf(txt, Sage(BaseRing(X)), var), false;
end intrinsic;
