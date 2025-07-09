(* ::Package:: *)

(* Created with the Wolfram Language : www.wolfram.com *)
(*-----------------------Source Coefficients---------------------------------------------------*)
(*Here, quntities like "Subsuperscript[A, 0,lm, (1,0)]" are represented as A0lm, 
and quntities like "Subsuperscript[A, 0,lm, (1,1)]" are represented as \[Delta]A0lm*)
(*--------------------Master variables------------------------------------------------------*)
(*Here, master variables are defined as follows:
 R[r]=Subsuperscript[\[CapitalPsi], R, (1,0)] ;
 \[Delta]R[r]=Subsuperscript[\[CapitalPsi], R, (1,1)];
 Kh[r]=Subsuperscript[\[CapitalPsi], Z, (1,0)] ;
 \[Delta]Kh[r]=Subsuperscript[\[CapitalPsi], Z, (1,1)];
  XF[r]=Subsuperscript[\[CapitalPsi], F, (1,1)] ;
n1=n;
r01v1=Subscript[r, (0,1)]*)
(*-----------------Description of the matrix---------------------------------------------------------*)
(*The matrix below is as follows
 {{"Source term at order \[Epsilon] for axial gravity case: Subsuperscript[S, R, (1,0)]",
 "Source term at order \[Epsilon]\[Xi] for axial gravity case: Subsuperscript[S, R, (1,1)]",
 "Source term at order \[Epsilon] for polar gravity case: Subsuperscript[S, Z, (1,0)]",
 "Source term at order \[Epsilon]\[Xi] for polar fluid case: Subsuperscript[S, F, (1,1)]"
 ,"Subsuperscript[C, 1,lm, (1,0)]","Subsuperscript[C, 2,lm, (1,0)]","Subsuperscript[C, 1,lm, (1,1)]","Subsuperscript[C, 2,lm, (1,1)]",
 "Source term at order \[Epsilon]\[Xi] for polar gravity case: Subsuperscript[S, Z, (1,1)]"}}*)
{((8*I)*Sqrt[2]*Pi*(-2 + r)*Dlm[r])/(Sqrt[n1*(1 + n1)]*r^2) + 
  ((8*I)*Pi*(1 - 2/r)^2*Qlm[r])/Sqrt[1 + n1] + 
  ((4*I)*Sqrt[2]*Pi*(-2 + r)^2*Derivative[1][Dlm][r])/(Sqrt[n1*(1 + n1)]*r), 
 ((8*I)*Sqrt[2]*Pi*(-2 + r)*\[Delta]Dlm[r])/(Sqrt[n1*(1 + n1)]*r^2) + 
  Qlm[r]*(((8*I)*Pi*(-2 + r)^2*r01v1)/(Sqrt[1 + n1]*r^3) + 
    ((16*I)*Pi*(-2 + r)*\[Delta]m[r])/(Sqrt[1 + n1]*r^2)) + 
  ((8*I)*Pi*(1 - 2/r)^2*\[Delta]Qlm[r])/Sqrt[1 + n1] + 
  ((8*I)*Pi*(1 - 2/r)^2*r01v1*Derivative[1][Qlm][r])/Sqrt[1 + n1] + 
  Derivative[1][Dlm][r]*(((8*I)*Sqrt[2]*Pi*(-2 + r)*\[Delta]m[r])/
     (Sqrt[n1*(1 + n1)]*r) + ((4*I)*Sqrt[2]*Pi*(2 - r)*
      (2*(1 - r)*r01v1 + (2 - r)*r*Derivative[1][r01][r]))/
     (Sqrt[n1*(1 + n1)]*r^2)) + ((4*I)*Sqrt[2]*Pi*(-2 + r)^2*
    Derivative[1][\[Delta]Dlm][r])/(Sqrt[n1*(1 + n1)]*r) + 
  Derivative[1][R][r]*((4*(-2 + r)*\[Delta]m[r])/(r^2*(2*r - r^2)) + 
    (4*(2 - 3*r)*\[Delta]m[r] - 3*(2 - r)*r*
       ((2 - r)*Derivative[1][\[Delta]f][r] + 4*Derivative[1][\[Delta]m][r]))/
     (2*r^3)) + Dlm[r]*(((-16*I)*Sqrt[2]*Pi*(-2 + r)*\[Delta]m[r])/
     (Sqrt[n1*(1 + n1)]*(2 - r)*r^2) + 
    ((4*I)*Sqrt[2]*Pi*(-2*(8 - 6*r + r^2)*r01v1 + 
       r*(2*r*\[Delta]m[r] + (2 - r)*((2 - r)*(2*Derivative[1][r01][r] + 
             r*Derivative[1][\[Delta]f][r]) + 2*r*Derivative[1][\[Delta]m][
             r]))))/(Sqrt[n1*(1 + n1)]*r^3)) + 
  ((4*I)*Sqrt[2]*Pi*(-2 + r)^2*r01v1*Derivative[2][Dlm][r])/
   (Sqrt[n1*(1 + n1)]*r) + R[r]*(\[Omega]^2*\[Delta]f[r] + 
    (2*(2*(-2 + r)*(-3 + r + n1*r) - r^4*\[Omega]^2)*\[Delta]m[r])/
     ((-2 + r)*r^4) + (2*(12 - 10*r + r^2 + r^4*\[Omega]^2)*\[Delta]m[r])/
     ((2 - r)*r^4) + ((-4 + r^2)*Derivative[1][\[Delta]f][r] - 
      4*(1 + r)*Derivative[1][\[Delta]m][r] + (-2 + r)^2*r*
       Derivative[2][\[Delta]f][r])/(2*r^3)), 
 (8*Sqrt[2]*Pi*(2 - r)*(3 - (3 + n1)*r)*A1lm[r])/(r*(3 + n1*r)^2*\[Omega]) + 
  (8*Pi*(-2 + r)^2*Alm[r])/(3 + n1*r) + 
  (8*Pi*(2 - r)*(12 + 3*(-2 + n1)*r + n1^2*r^2)*B0lm[r])/
   (Sqrt[1 + n1]*r*(3 + n1*r)^2*\[Omega]) + (8*Pi*(-2 + r)^2*Blm[r])/
   (Sqrt[1 + n1]*(3 + n1*r)) + (8*Sqrt[2]*Pi*(2 - r)*Flm[r])/
   Sqrt[n1*(1 + n1)] + (4*Sqrt[2]*Pi*(-2 + r)^2*Derivative[1][A1lm][r])/
   ((3 + n1*r)*\[Omega]) + (8*Pi*(-2 + r)^2*Derivative[1][B0lm][r])/
   (Sqrt[1 + n1]*(3 + n1*r)*\[Omega]), 
 -1/8*((1 - 2/r)^((3 + Csr^(-2))/4)*r^(-3 - Cst^2/Csr^2)*
     (-4*(1 + n1) + r*(2 + 2*n1 + r^2*\[Omega]^2))*
     ((-8*Sqrt[2]*Pi*A1lm[r])/(2*n1*\[Omega] + (6*\[Omega])/r) + 
      (16*Pi*(1 - 2/r)*r^2*Alm[r])/(2*n1 + 6/r) + 
      (16*Pi*(1 - 2/r)*r^2*Blm[r])/(Sqrt[1 + n1]*(2*n1 + 6/r)))*
     Derivative[1][\[Delta]m][r])/(Csr^2*Pi) + 
  ((1 - 2/r)^((-1 + Csr^(-2))/4)*(-2 + r)^2*r^(-5 - Cst^2/Csr^2)*
    Derivative[1][Kh][r]*
    (-((3*(4 + n1) - 3*n1*r - 5*n1^2*r + 6*r^3*\[Omega]^2 + 
        2*n1*r^2*(1 + n1 + r^2*\[Omega]^2))*Derivative[1][\[Delta]m][r]) + 
     r*(3*(2 + n1) - n1^2*r + 3*r^3*\[Omega]^2 + 
       n1*r^2*(1 + n1 + r^2*\[Omega]^2))*Derivative[2][\[Delta]m][r]))/
   (8*Csr^2*Pi*(3 + n1*r)) + ((1 - 2/r)^((-1 + Csr^(-2))/4)*
    r^(-6 - Cst^2/Csr^2)*Kh[r]*((36*(4 + n1) - 2*n1^2*(1 + n1)*r^5 + 
       12*r*(-12 + n1^2 + 6*r^2*\[Omega]^2) + 
       n1*r^4*(2*n1^3 + 6*r^2*\[Omega]^2 + n1*(17 - 3*r^2*\[Omega]^2) + 
         n1^2*(19 - 2*r^2*\[Omega]^2)) + 3*r^2*(12 - 6*n1^2 + 4*n1^3 - 
         33*r^2*\[Omega]^2 + n1*(-15 + 2*r^2*\[Omega]^2)) - 
       2*r^3*(18*n1^3 + 2*n1^4 - 9*r^2*\[Omega]^2 + 
         2*n1^2*(5 + 3*r^2*\[Omega]^2) + 3*n1*(-3 + 7*r^2*\[Omega]^2)))*
      Derivative[1][\[Delta]m][r] - (-2 + r)^2*r*(9*(2 + n1) + 9*n1*r + 
       6*n1^2*r - n1^2*(1 + n1)*r^3 + 9*r^3*\[Omega]^2 + 
       3*n1*r^2*(n1 + n1^2 + r^2*\[Omega]^2))*Derivative[2][\[Delta]m][r]))/
   (8*Csr^2*Pi*(3 + n1*r)^2) - ((1 - 2/r)^((-1 + Csr^(-2))/4)*(-2 + r)*
    r^(-3 - Cst^2/Csr^2)*((4*Sqrt[2]*Pi*A1lm[r])/\[Omega] + 
     ((-8*Sqrt[2]*Pi*A1lm[r])/(2*n1*\[Omega] + (6*\[Omega])/r) + 
       (16*Pi*(1 - 2/r)*r^2*Alm[r])/(2*n1 + 6/r) + 
       (16*Pi*(1 - 2/r)*r^2*Blm[r])/(Sqrt[1 + n1]*(2*n1 + 6/r)))/r - 
     (8*Sqrt[2]*Pi*r*Flm[r])/Sqrt[n1*(1 + n1)])*
    ((3*n1 + n1*r^3*(-2 - 2*n1 + r^2*\[Omega]^2) - 
       r*(6 + 16*n1 + 9*n1^2 + 9*r^2*\[Omega]^2) + 
       r^2*(7*n1 + 7*n1^2 + 6*r^2*\[Omega]^2 - n1*r^2*\[Omega]^2))*
      Derivative[1][\[Delta]m][r] - (2 - r)*r*(3 + 6*n1 - 4*n1*r + 
       6*r^3*\[Omega]^2 + n1*r^2*(1 + 2*r^2*\[Omega]^2))*
      Derivative[2][\[Delta]m][r]))/(8*Csr^2*Pi*(3 + n1*r)^2) - 
  ((I/8)*(1 - 2/r)^((-1 + Csr^(-2))/4)*(-2 + r)^2*r^(-5 - Cst^2/Csr^2)*
    (((8*I)*Pi*r*B0lm[r])/(Sqrt[1 + n1]*(1 - 2/r)*\[Omega]) - 
     (I*((-8*Sqrt[2]*Pi*A1lm[r])/(2*n1*\[Omega] + (6*\[Omega])/r) + 
        (16*Pi*(1 - 2/r)*r^2*Alm[r])/(2*n1 + 6/r) + 
        (16*Pi*(1 - 2/r)*r^2*Blm[r])/(Sqrt[1 + n1]*(2*n1 + 6/r))))/
      (1 - 2/r) + ((8*I)*Sqrt[2]*Pi*r^2*Flm[r])/(Sqrt[n1*(1 + n1)]*
       (1 - 2/r)))*((6*(-1 + 2*n1) + r*(3 + 11*n1 + 21*n1^2 - 
         12*r^2*\[Omega]^2) + n1*r^3*(2 + n1 - n1^2 - n1*r^2*\[Omega]^2) + 
       r^2*(-15*n1^2 + n1^3 - 3*r^2*\[Omega]^2 - n1*(14 + 9*r^2*\[Omega]^2)))*
      Derivative[1][\[Delta]m][r] + (2 - r)*r*
      (3 - r*(n1 + 2*n1^2 - 6*r^2*\[Omega]^2) + 
       n1*r^2*(1 + n1 + 2*r^2*\[Omega]^2))*Derivative[2][\[Delta]m][r]))/
   (Csr^2*Pi*(3 + n1*r)^2) + ((1 - 2/r)^((3 + Csr^(-2))/4)*
    r^(-3 - Cst^2/Csr^2)*(-8*Sqrt[2]*Pi*r^2*(3 + n1*r)*
      (1 + n1 - r^3*\[Omega]^2)*Flm[r]*Derivative[1][\[Delta]m][r] - 
     Sqrt[n1*(1 + n1)]*(2 - r)*(Derivative[1][\[Delta]m][r]*
        (r^2*(1 - n1 + n1*r - r^3*\[Omega]^2)*
          (-(((-8*Sqrt[2]*Pi*A1lm[r])/(2*n1*\[Omega] + (6*\[Omega])/r) + 
              (16*Pi*(1 - 2/r)*r^2*Alm[r])/(2*n1 + 6/r) + (16*Pi*(1 - 2/r)*
                r^2*Blm[r])/(Sqrt[1 + n1]*(2*n1 + 6/r)))/r^2) - 
           (8*Sqrt[2]*Pi*Flm[r])/Sqrt[n1*(1 + n1)] + 
           (4*Sqrt[2]*Pi*Derivative[1][A1lm][r])/\[Omega] + 
           ((-48*Sqrt[2]*Pi*\[Omega]*A1lm[r])/(r^2*(2*n1*\[Omega] + 
                 (6*\[Omega])/r)^2) + (96*Pi*(1 - 2/r)*Alm[r])/
              (2*n1 + 6/r)^2 + (32*Pi*Alm[r])/(2*n1 + 6/r) + 
             (32*Pi*(1 - 2/r)*r*Alm[r])/(2*n1 + 6/r) + 
             (96*Pi*(1 - 2/r)*Blm[r])/(Sqrt[1 + n1]*(2*n1 + 6/r)^2) + 
             (32*Pi*Blm[r])/(Sqrt[1 + n1]*(2*n1 + 6/r)) + 
             (32*Pi*(1 - 2/r)*r*Blm[r])/(Sqrt[1 + n1]*(2*n1 + 6/r)) - 
             (8*Sqrt[2]*Pi*Derivative[1][A1lm][r])/(2*n1*\[Omega] + 
               (6*\[Omega])/r) + (16*Pi*(1 - 2/r)*r^2*Derivative[1][Alm][r])/
              (2*n1 + 6/r) + (16*Pi*(1 - 2/r)*r^2*Derivative[1][Blm][r])/
              (Sqrt[1 + n1]*(2*n1 + 6/r)))/r - 
           (8*Sqrt[2]*Pi*r*Derivative[1][Flm][r])/Sqrt[n1*(1 + n1)]) + 
         I*(2 - r)*((1 + n1 - r^3*\[Omega]^2)*(((8*I)*Pi*B0lm[r])/
              (Sqrt[1 + n1]*(1 - 2/r)*\[Omega]) - ((16*I)*Pi*B0lm[r])/
              (Sqrt[1 + n1]*(1 - 2/r)^2*r*\[Omega]) + 
             ((2*I)*((-8*Sqrt[2]*Pi*A1lm[r])/(2*n1*\[Omega] + (6*\[Omega])/
                   r) + (16*Pi*(1 - 2/r)*r^2*Alm[r])/(2*n1 + 6/r) + 
                (16*Pi*(1 - 2/r)*r^2*Blm[r])/(Sqrt[1 + n1]*(2*n1 + 6/r))))/
              ((1 - 2/r)^2*r^2) - ((16*I)*Sqrt[2]*Pi*Flm[r])/
              (Sqrt[n1*(1 + n1)]*(1 - 2/r)^2) + ((16*I)*Sqrt[2]*Pi*r*Flm[r])/
              (Sqrt[n1*(1 + n1)]*(1 - 2/r)) + ((8*I)*Pi*r*Derivative[1][B0lm][
                r])/(Sqrt[1 + n1]*(1 - 2/r)*\[Omega]) - 
             (I*((-48*Sqrt[2]*Pi*\[Omega]*A1lm[r])/(r^2*(2*n1*\[Omega] + 
                    (6*\[Omega])/r)^2) + (96*Pi*(1 - 2/r)*Alm[r])/
                 (2*n1 + 6/r)^2 + (32*Pi*Alm[r])/(2*n1 + 6/r) + 
                (32*Pi*(1 - 2/r)*r*Alm[r])/(2*n1 + 6/r) + (96*Pi*(1 - 2/r)*
                  Blm[r])/(Sqrt[1 + n1]*(2*n1 + 6/r)^2) + (32*Pi*Blm[r])/
                 (Sqrt[1 + n1]*(2*n1 + 6/r)) + (32*Pi*(1 - 2/r)*r*Blm[r])/
                 (Sqrt[1 + n1]*(2*n1 + 6/r)) - (8*Sqrt[2]*Pi*Derivative[1][
                    A1lm][r])/(2*n1*\[Omega] + (6*\[Omega])/r) + 
                (16*Pi*(1 - 2/r)*r^2*Derivative[1][Alm][r])/(2*n1 + 6/r) + 
                (16*Pi*(1 - 2/r)*r^2*Derivative[1][Blm][r])/(Sqrt[1 + n1]*
                  (2*n1 + 6/r))))/(1 - 2/r) + ((8*I)*Sqrt[2]*Pi*r^2*
               Derivative[1][Flm][r])/(Sqrt[n1*(1 + n1)]*(1 - 2/r))) + 
           I*r*(3 + n1*r)*((-576*Sqrt[2]*Pi*\[Omega]^2*A1lm[r])/
              (r^4*(2*n1*\[Omega] + (6*\[Omega])/r)^3) + 
             (96*Sqrt[2]*Pi*\[Omega]*A1lm[r])/(r^3*(2*n1*\[Omega] + 
                 (6*\[Omega])/r)^2) + (32*Pi*(1 - 2/r)*Alm[r])/(2*n1 + 6/r) + 
             (1152*Pi*(1 - 2/r)*Alm[r])/((2*n1 + 6/r)^3*r^2) + 
             (384*Pi*Alm[r])/((2*n1 + 6/r)^2*r^2) + (192*Pi*(1 - 2/r)*Alm[r])/
              ((2*n1 + 6/r)^2*r) + (64*Pi*Alm[r])/((2*n1 + 6/r)*r) + 
             (32*Pi*(1 - 2/r)*Blm[r])/(Sqrt[1 + n1]*(2*n1 + 6/r)) + 
             (1152*Pi*(1 - 2/r)*Blm[r])/(Sqrt[1 + n1]*(2*n1 + 6/r)^3*r^2) + 
             (384*Pi*Blm[r])/(Sqrt[1 + n1]*(2*n1 + 6/r)^2*r^2) + 
             (192*Pi*(1 - 2/r)*Blm[r])/(Sqrt[1 + n1]*(2*n1 + 6/r)^2*r) + 
             (64*Pi*Blm[r])/(Sqrt[1 + n1]*(2*n1 + 6/r)*r) - 
             (96*Sqrt[2]*Pi*\[Omega]*Derivative[1][A1lm][r])/
              (r^2*(2*n1*\[Omega] + (6*\[Omega])/r)^2) + 
             (192*Pi*(1 - 2/r)*Derivative[1][Alm][r])/(2*n1 + 6/r)^2 + 
             (64*Pi*Derivative[1][Alm][r])/(2*n1 + 6/r) + 
             (64*Pi*(1 - 2/r)*r*Derivative[1][Alm][r])/(2*n1 + 6/r) + 
             (192*Pi*(1 - 2/r)*Derivative[1][Blm][r])/(Sqrt[1 + n1]*
               (2*n1 + 6/r)^2) + (64*Pi*Derivative[1][Blm][r])/
              (Sqrt[1 + n1]*(2*n1 + 6/r)) + (64*Pi*(1 - 2/r)*r*Derivative[1][
                 Blm][r])/(Sqrt[1 + n1]*(2*n1 + 6/r)) - 
             (8*Sqrt[2]*Pi*Derivative[2][A1lm][r])/(2*n1*\[Omega] + 
               (6*\[Omega])/r) + (16*Pi*(1 - 2/r)*r^2*Derivative[2][Alm][r])/
              (2*n1 + 6/r) + (16*Pi*(1 - 2/r)*r^2*Derivative[2][Blm][r])/
              (Sqrt[1 + n1]*(2*n1 + 6/r))))) + 
       (3 + n1*r)*((-48*Sqrt[2]*Pi*\[Omega]*A1lm[r])/
          (r^2*(2*n1*\[Omega] + (6*\[Omega])/r)^2) + (96*Pi*(1 - 2/r)*Alm[r])/
          (2*n1 + 6/r)^2 + (32*Pi*Alm[r])/(2*n1 + 6/r) + 
         (32*Pi*(1 - 2/r)*r*Alm[r])/(2*n1 + 6/r) + (96*Pi*(1 - 2/r)*Blm[r])/
          (Sqrt[1 + n1]*(2*n1 + 6/r)^2) + (32*Pi*Blm[r])/
          (Sqrt[1 + n1]*(2*n1 + 6/r)) + (32*Pi*(1 - 2/r)*r*Blm[r])/
          (Sqrt[1 + n1]*(2*n1 + 6/r)) - (8*Sqrt[2]*Pi*Derivative[1][A1lm][r])/
          (2*n1*\[Omega] + (6*\[Omega])/r) + (16*Pi*(1 - 2/r)*r^2*
           Derivative[1][Alm][r])/(2*n1 + 6/r) + 
         (16*Pi*(1 - 2/r)*r^2*Derivative[1][Blm][r])/(Sqrt[1 + n1]*
           (2*n1 + 6/r)))*(4*Derivative[1][\[Delta]m][r] + 
         (-2 + r)*r*Derivative[2][\[Delta]m][r]))))/
   (8*Csr^2*Sqrt[n1*(1 + n1)]*Pi*(3 + n1*r)), 
 (4*Sqrt[2]*Pi*(2 + n1*r)*A1lm[r])/((3 + n1*r)*\[Omega]) - 
  (8*Pi*(2 - r)*r*Alm[r])/(3 + n1*r) - (8*Pi*(2 - r)*r*Blm[r])/
   (Sqrt[1 + n1]*(3 + n1*r)) - (8*Sqrt[2]*Pi*r*Flm[r])/Sqrt[n1*(1 + n1)], 
 ((-4*I)*Sqrt[2]*Pi*r^2*A1lm[r])/((2 - r)*(3 + n1*r)*\[Omega]) - 
  ((8*I)*Pi*r^3*Alm[r])/(3 + n1*r) + ((8*I)*Pi*r*B0lm[r])/
   (Sqrt[1 + n1]*(\[Omega] - (2*\[Omega])/r)) - 
  ((8*I)*Pi*r^3*Blm[r])/(Sqrt[1 + n1]*(3 + n1*r)) - 
  ((8*I)*Sqrt[2]*Pi*r^3*Flm[r])/(Sqrt[n1*(1 + n1)]*(2 - r)), 
 (4*Pi*(1 - 2/r)^((-3 - Csr^(-2))/4)*r^(-4 + Cst^2/Csr^2)*
    (2 - 2*Csr^2 + 8*Cst^2 + (4*Cst^2*(-1 + n1) + n1 - Csr^2*n1)*r - 
     2*Cst^2*n1*r^2 + 2*Csr^2*r^4*\[Omega]^2)*XF[r])/
   ((3 + n1*r)*\[Omega]^2) + (4*Sqrt[2]*Pi*(2 + n1*r)*\[Delta]A1lm[r])/
   ((3 + n1*r)*\[Omega]) - (8*Pi*(2 - r)*r*\[Delta]Alm[r])/(3 + n1*r) - 
  (8*Pi*(2 - r)*r*\[Delta]Blm[r])/(Sqrt[1 + n1]*(3 + n1*r)) - 
  (8*Sqrt[2]*Pi*r*\[Delta]Flm[r])/Sqrt[n1*(1 + n1)] - 
  (8*Sqrt[2]*Pi*r*r01v1*Derivative[1][Flm][r])/Sqrt[n1*(1 + n1)] + 
  (8*Csr^2*Pi*(1 - 2/r)^(1/4 - 1/(4*Csr^2))*r^(-2 + Cst^2/Csr^2)*(2 + n1*r)*
    Derivative[1][XF][r])/((3 + n1*r)*\[Omega]^2) - 
  (8*Pi*(2 - r)*Derivative[1][Alm][r]*(r^2*r01v1*\[Omega]^2 + 
     (2 - r)*Derivative[1][\[Delta]m][r]))/(r*(3 + n1*r)*\[Omega]^2) - 
  (8*Pi*(2 - r)*Derivative[1][Blm][r]*(r^2*r01v1*\[Omega]^2 + 
     (2 - r)*Derivative[1][\[Delta]m][r]))/(Sqrt[1 + n1]*r*(3 + n1*r)*
    \[Omega]^2) + (4*Sqrt[2]*Pi*Derivative[1][A1lm][r]*
    (r^2*(2 + n1*r)*r01v1*\[Omega]^2 + (-2 + r)*Derivative[1][\[Delta]m][r]))/
   (r^2*(3 + n1*r)*\[Omega]^3) - 
  (2*Sqrt[2]*Pi*A1lm[r]*(2*n1*r^4*\[Omega]^2*\[Delta]m[r] + 
     n1*(-2 + r)*r^5*\[Omega]^2*Derivative[1][\[Delta]f][r] + 
     2*(6 + (-5 + n1)*r - n1*r^2 - n1*r^5*\[Omega]^2)*
      Derivative[1][\[Delta]m][r]))/(r^3*(3 + n1*r)^2*\[Omega]^3) + 
  (C1lm1[r]*((-2*r^3*(6*(-1 + n1) - 2*r*(7*n1 + n1^2 - 3*r^2*\[Omega]^2) + 
        n1*r^2*(3 + 2*r^2*\[Omega]^2))*\[Delta]m[r])/(2 - r) + 
     (2 - r)*r^3*(6 + 8*n1*r - 3*n1*r^2)*Derivative[1][\[Delta]f][r] + 
     (2*(3*(1 + 5*n1) + n1*r^2*(2 - 5*n1 + 19*r^2*\[Omega]^2) + 
        r*(-9*n1 + 7*n1^2 + 21*r^2*\[Omega]^2) + 
        n1*r^3*(n1 - 3*r^2*\[Omega]^2 + 2*n1*r^2*\[Omega]^2))*
       Derivative[1][\[Delta]m][r])/\[Omega]^2))/(2*r^2*(3 + n1*r)^3) - 
  (8*Sqrt[2]*Pi*Flm[r]*(2*(3 + n1*r)*r01v1 + 
     r^2*((2 - r)*Derivative[1][\[Delta]f][r] + 
       2*Derivative[1][\[Delta]m][r])))/(Sqrt[n1*(1 + n1)]*(3 + n1*r)) - 
  (4*Pi*Alm[r]*(2*(3 + 2*n1)*r^4*\[Omega]^2*\[Delta]m[r] + 
     (2 - r)*(3*(2 - r)*r^4*\[Omega]^2*Derivative[1][\[Delta]f][r] + 
       2*(9 + 2*(-4 + n1)*r - 2*n1*r^2 + 3*r^4*\[Omega]^2)*
        Derivative[1][\[Delta]m][r])))/(r^2*(3 + n1*r)^2*\[Omega]^2) - 
  (4*Pi*Blm[r]*(2*(2 - r)*r^2*(3 + n1*r)*r01v1*\[Omega]^2 + 
     2*(3 + 2*n1)*r^4*\[Omega]^2*\[Delta]m[r] + 
     (2 - r)*(3*(2 - r)*r^4*\[Omega]^2*Derivative[1][\[Delta]f][r] + 
       2*(12 + (-8 + 3*n1)*r - 2*n1*r^2 + 3*r^4*\[Omega]^2)*
        Derivative[1][\[Delta]m][r])))/(Sqrt[1 + n1]*r^2*(3 + n1*r)^2*
    \[Omega]^2) + (Derivative[1][Kh][r]*
    (-2*r^3*\[Omega]^2*(6*(-1 + n1) - 2*r*(7*n1 + n1^2 - 3*r^2*\[Omega]^2) + 
       n1*r^2*(3 + 2*r^2*\[Omega]^2))*\[Delta]m[r] + 
     (2 - r)*((2 - r)*r^3*(6 + 8*n1*r - 3*n1*r^2)*\[Omega]^2*
        Derivative[1][\[Delta]f][r] + 
       2*(6*(2 + n1) + r*(6*n1 + n1^2 + 12*r^2*\[Omega]^2) + 
         n1*r^2*(2 + 2*n1 - n1^2 + 13*r^2*\[Omega]^2) + 
         n1*r^3*(n1 + n1^2 - 3*r^2*\[Omega]^2 + n1*r^2*\[Omega]^2))*
        Derivative[1][\[Delta]m][r])))/(2*(-2 + r)*r^4*(3 + n1*r)^2*
    \[Omega]^2) + 
  ((I/2)*C2lm1[r]*(-2*r^3*\[Omega]^2*(6*(-1 + n1) - 
       2*r*(7*n1 + n1^2 - 3*r^2*\[Omega]^2) + n1*r^2*(3 + 2*r^2*\[Omega]^2))*
      \[Delta]m[r] + (2 - r)*((2 - r)*r^3*(6 + 8*n1*r - 3*n1*r^2)*\[Omega]^2*
        Derivative[1][\[Delta]f][r] + 2*(-3*(-1 + n1) - 5*n1^2*r + 
         21*r^3*\[Omega]^2 + n1*r^2*(2 + n1 - 2*n1^2 + 19*r^2*\[Omega]^2) + 
         n1*r^3*(n1 + n1^2 - 3*r^2*\[Omega]^2 + 2*n1*r^2*\[Omega]^2))*
        Derivative[1][\[Delta]m][r])))/(r^4*(3 + n1*r)^3*\[Omega]^2) + 
  (Kh[r]*(-2*(2 - r)*r^6*(3 + n1*r)^3*\[Omega]^4*\[Delta]f[r] - 
     2*r^3*\[Omega]^2*(36*(-1 + n1) + 6*r*(-14*n1 + 4*n1^2 + 
         15*r^2*\[Omega]^2) + 3*r^2*(n1 - 13*n1^2 + 4*n1^3 - 
         9*r^2*\[Omega]^2 + 22*n1*r^2*\[Omega]^2) - 
       n1*r^3*(20*n1^2 + 12*r^2*\[Omega]^2 + n1*(11 - 18*r^2*\[Omega]^2)) - 
       n1^2*r^4*(-3 + n1^2 + r^2*\[Omega]^2 - 2*n1*(1 + r^2*\[Omega]^2)))*
      \[Delta]m[r] + (2 - r)*((2 - r)*r^3*\[Omega]^2*(36 + 66*n1*r + 
         n1^2*r^4*(-3 - 2*n1 + n1^2 + 3*r^2*\[Omega]^2) + 
         3*r^2*(-n1 + 13*n1^2 + 9*r^2*\[Omega]^2) + 
         n1*r^3*(5*n1 + 14*n1^2 + 18*r^2*\[Omega]^2))*
        Derivative[1][\[Delta]f][r] + 
       2*(36*(2 + n1) + 6*r*(-6 + 9*n1 + 7*n1^2 + 12*r^2*\[Omega]^2) + 
         n1^2*r^5*(n1 + n1^2 - 3*r^2*\[Omega]^2 - 2*n1*r^2*\[Omega]^2 + 
           n1^2*r^2*\[Omega]^2 + 3*r^4*\[Omega]^4) + 
         3*r^2*(3*n1^2 + 8*n1^3 - 6*r^2*\[Omega]^2 + 
           4*n1*(-3 + 8*r^2*\[Omega]^2)) + n1*r^4*(-5*n1^3 + 
           18*r^4*\[Omega]^4 + 2*n1*(1 + r^2*\[Omega]^2) + 
           n1^2*(-3 + 14*r^2*\[Omega]^2)) + r^3*(-10*n1^3 + 6*n1^4 - 
           18*n1*r^2*\[Omega]^2 + 27*r^4*\[Omega]^4 + 
           n1^2*(-19 + 45*r^2*\[Omega]^2)))*Derivative[1][\[Delta]m][r])))/
   (2*(-2 + r)^2*r^5*(3 + n1*r)^3*\[Omega]^2), 
 ((-4*I)*Pi*(1 - 2/r)^((-3 - Csr^(-2))/4)*r^(-2 + Cst^2/Csr^2)*
    (-1 + Csr^2 + 20*Cst^2 + 2*Cst^2*(-5 + 4*n1)*r - 4*Cst^2*n1*r^2 + 
     2*Csr^2*r^4*\[Omega]^2)*XF[r])/((-2 + r)*(3 + n1*r)*\[Omega]^2) - 
  ((4*I)*Sqrt[2]*Pi*r^2*\[Delta]A1lm[r])/((2 - r)*(3 + n1*r)*\[Omega]) - 
  ((8*I)*Pi*r^3*\[Delta]Alm[r])/(3 + n1*r) + ((8*I)*Pi*r*\[Delta]B0lm[r])/
   (Sqrt[1 + n1]*(\[Omega] - (2*\[Omega])/r)) - 
  ((8*I)*Pi*r^3*\[Delta]Blm[r])/(Sqrt[1 + n1]*(3 + n1*r)) - 
  ((8*I)*Sqrt[2]*Pi*r^3*\[Delta]Flm[r])/(Sqrt[n1*(1 + n1)]*(2 - r)) + 
  ((8*I)*Pi*r*B0lm[r]*((-2 + r)*r01v1 + 2*r*\[Delta]m[r]))/
   (Sqrt[1 + n1]*(-2 + r)^2*\[Omega]) - 
  ((4*I)*Sqrt[2]*Pi*r^2*r01v1*Derivative[1][A1lm][r])/
   ((2 - r)*(3 + n1*r)*\[Omega]) - ((8*I)*Pi*r^3*r01v1*Derivative[1][Alm][r])/
   (3 + n1*r) - ((8*I)*Pi*r^3*r01v1*Derivative[1][Blm][r])/
   (Sqrt[1 + n1]*(3 + n1*r)) - 
  ((8*I)*Pi*r^2*r01v1*(Sqrt[n1*(1 + n1)]*Derivative[1][B0lm][r] + 
     Sqrt[2]*Sqrt[1 + n1]*r*\[Omega]*Derivative[1][Flm][r]))/
   (Sqrt[n1]*(1 + n1)*(2 - r)*\[Omega]) + 
  ((8*I)*Csr^2*Pi*(1 - 2/r)^((-3 - Csr^(-2))/4)*r^(-1 + Cst^2/Csr^2)*
    Derivative[1][XF][r])/((3 + n1*r)*\[Omega]^2) + 
  ((4*I)*Pi*Blm[r]*(-2*r^2*(3 + n1*r)*r01v1*\[Omega]^2 + 
     6*r^3*\[Omega]^2*\[Delta]m[r] - 6*r^4*\[Omega]^2*
      Derivative[1][\[Delta]f][r] + 3*r^5*\[Omega]^2*Derivative[1][\[Delta]f][
       r] + 24*Derivative[1][\[Delta]m][r] - 
     14*r*Derivative[1][\[Delta]m][r] + 6*n1*r*Derivative[1][\[Delta]m][r] - 
     4*n1*r^2*Derivative[1][\[Delta]m][r] - 6*r^4*\[Omega]^2*
      Derivative[1][\[Delta]m][r]))/(Sqrt[1 + n1]*(3 + n1*r)^2*\[Omega]^2) + 
  ((4*I)*Pi*Alm[r]*(6*r^3*\[Omega]^2*\[Delta]m[r] + 
     3*(-2 + r)*r^4*\[Omega]^2*Derivative[1][\[Delta]f][r] + 
     2*(15 + (-7 + 4*n1)*r - 2*n1*r^2 - 3*r^4*\[Omega]^2)*
      Derivative[1][\[Delta]m][r]))/((3 + n1*r)^2*\[Omega]^2) + 
  (C2lm1[r]*((2*(54 + 3*(r + 18*n1*r) + n1*r^3*(3 - n1 + 2*r^2*\[Omega]^2) + 
        2*r^2*(-7*n1 + 4*n1^2 + 3*r^2*\[Omega]^2))*\[Delta]m[r])/(2 - r) + 
     (2 - r)*r*(12 + 4*n1*r + n1*(3 + 2*n1)*r^2)*Derivative[1][\[Delta]f][
       r] + (2*(36 + 3*(-4 + 17*n1)*r + n1*r^3*(7 - 9*n1 + 
          11*r^2*\[Omega]^2) + r^2*(-36*n1 + 11*n1^2 + 24*r^2*\[Omega]^2) + 
        n1*r^4*(3*r^2*\[Omega]^2 + n1*(2 + 3*r^2*\[Omega]^2)))*
       Derivative[1][\[Delta]m][r])/(r^3*\[Omega]^2)))/(2*(3 + n1*r)^3) + 
  ((8*I)*Sqrt[2]*Pi*Flm[r]*(-2*(2 - r)*r^2*(3 + n1*r)*r01v1*\[Omega]^2 + 
     2*r^3*(3 + n1*r)*\[Omega]^2*\[Delta]m[r] + 
     (2 - r)*((-2 + r)*r^4*\[Omega]^2*Derivative[1][\[Delta]f][r] + 
       (3 + n1*r - 2*r^4*\[Omega]^2)*Derivative[1][\[Delta]m][r])))/
   (Sqrt[n1*(1 + n1)]*(-2 + r)^2*(3 + n1*r)*\[Omega]^2) + 
  ((2*I)*Sqrt[2]*Pi*A1lm[r]*(2*r^3*(6 + n1*r^2)*\[Omega]^2*\[Delta]m[r] + 
     (2 - r)*(n1*(2 - r)*r^5*\[Omega]^2*Derivative[1][\[Delta]f][r] + 
       2*(12 + (-4 + 3*n1)*r - n1*r^2 + n1*r^5*\[Omega]^2)*
        Derivative[1][\[Delta]m][r])))/((-2 + r)^2*r*(3 + n1*r)^2*
    \[Omega]^3) - 
  ((I/2)*C1lm1[r]*(2*r^3*\[Omega]^2*(54 + 3*(r + 18*n1*r) + 
       n1*r^3*(3 - n1 + 2*r^2*\[Omega]^2) + 
       2*r^2*(-7*n1 + 4*n1^2 + 3*r^2*\[Omega]^2))*\[Delta]m[r] + 
     (2 - r)*((2 - r)*r^4*(12 + 4*n1*r + n1*(3 + 2*n1)*r^2)*\[Omega]^2*
        Derivative[1][\[Delta]f][r] + 2*(36 + 3*(-4 + 17*n1)*r + 
         n1*r^3*(7 - 9*n1 + 11*r^2*\[Omega]^2) + 
         r^2*(-36*n1 + 11*n1^2 + 24*r^2*\[Omega]^2) + 
         n1*r^4*(2*n1 + 3*r^2*\[Omega]^2 + 3*n1*r^2*\[Omega]^2))*
        Derivative[1][\[Delta]m][r])))/((-2 + r)^2*r*(3 + n1*r)^3*
    \[Omega]^2) + ((I/2)*Derivative[1][Kh][r]*
    (2*r^3*\[Omega]^2*(54 + 3*(r + 18*n1*r) + 
       n1*r^3*(3 - n1 + 2*r^2*\[Omega]^2) + 
       2*r^2*(-7*n1 + 4*n1^2 + 3*r^2*\[Omega]^2))*\[Delta]m[r] + 
     (2 - r)*((2 - r)*r^4*(12 + 4*n1*r + n1*(3 + 2*n1)*r^2)*\[Omega]^2*
        Derivative[1][\[Delta]f][r] + 2*(36 + 3*(-4 + 17*n1)*r + 
         n1*r^3*(7 - 9*n1 + 11*r^2*\[Omega]^2) + 
         r^2*(-36*n1 + 11*n1^2 + 24*r^2*\[Omega]^2) + 
         n1*r^4*(2*n1 + 3*r^2*\[Omega]^2 + 3*n1*r^2*\[Omega]^2))*
        Derivative[1][\[Delta]m][r])))/((-2 + r)^2*r^3*(3 + n1*r)^2*
    \[Omega]^2) - ((I/2)*Kh[r]*(2*(2 - r)*r^7*(3 + n1*r)^3*\[Omega]^4*
      \[Delta]f[r] + 2*r^3*\[Omega]^2*(324 + 18*(-11 + 24*n1)*r + 
       3*r^2*(-91*n1 + 76*n1^2 + 48*r^2*\[Omega]^2) + 
       n1^2*r^5*(3 + n1 - 4*n1^2 - r^2*\[Omega]^2 + 4*n1*r^2*\[Omega]^2) + 
       n1*r^4*(-39*n1^2 + 8*n1^3 - 12*r^2*\[Omega]^2 + 
         n1*(-11 + 36*r^2*\[Omega]^2)) + 3*r^3*(-38*n1^2 + 24*n1^3 - 
         9*r^2*\[Omega]^2 + n1*(4 + 40*r^2*\[Omega]^2)))*\[Delta]m[r] + 
     (2 - r)*(-((2 - r)*r^4*\[Omega]^2*(9 + 48*n1*r + 
          n1^2*r^4*(-3 - 4*n1 + 3*r^2*\[Omega]^2) + 
          3*r^2*(-7*n1 + 9*n1^2 + 9*r^2*\[Omega]^2) + 
          n1*r^3*(-7*n1 + 8*n1^2 + 18*r^2*\[Omega]^2))*
         Derivative[1][\[Delta]f][r]) + 2*(216 + 18*(-10 + 17*n1)*r + 
         3*r^2*(12 - 93*n1 + 52*n1^2 + 30*r^2*\[Omega]^2) + 
         3*r^3*(-40*n1^2 + 18*n1^3 - 21*r^2*\[Omega]^2 + 
           7*n1*(3 + r^2*\[Omega]^2)) + n1^2*r^6*(2*n1^2 + 3*r^2*\[Omega]^2 - 
           3*r^4*\[Omega]^4 + n1*(2 + 3*r^2*\[Omega]^2)) - 
         n1*r^5*(8*n1^3 + 18*r^4*\[Omega]^4 + n1*(-7 + 5*r^2*\[Omega]^2) + 
           n1^2*(-5 + 7*r^2*\[Omega]^2)) - r^4*(45*n1^3 - 8*n1^4 + 
           27*n1*r^2*\[Omega]^2 + 27*r^4*\[Omega]^4 + 
           n1^2*(-7 + 12*r^2*\[Omega]^2)))*Derivative[1][\[Delta]m][r])))/
   ((2 - r)^3*r^4*(3 + n1*r)^3*\[Omega]^2), 
 (3*(2 - r)*(1 - (2 + n1)*r)*\[Delta]C1lm1[r])/(r*(3 + n1*r)^2) + 
  (I*(-2 + r)^2*(-3 - n1 + n1*(1 + n1)*r)*\[Delta]C2lm1[r])/
   (r^2*(3 + n1*r)^2) + ((-2 + r)^2*Derivative[1][\[Delta]C1lm1][r])/
   (3 + n1*r) + (I*(2 - r)^3*Derivative[1][\[Delta]C2lm1][r])/
   (r^2*(3 + n1*r))}
