# Flow quantities for hemodynamics in a vessel with elliptic cross section
# Copyright (c) Robert M. Corless 2023  MIT Release license
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# -------------------------------------------------------------------------------


# The active case has constant area, same as a circle
# with radius recorded in the global "radius".

# Input: 
# f_e = fraction the minor axis gets squished to.  Must be in 0 < f_e < 1
#      (the coordinate change is singular at f_e = 1 and at f_e = 0 the ellipse would 
#      be infinitely long and have zero width so that's singular too).
# omeg = applied frequency in Hz.  Typically between 0.1 and 10, 
#        but can be in (0,infinity)

Active := proc( f_e, omeg )
  local alpha, asp_ratio, beta, E, eccentri;
  global confracf, Refine, radius, rho, mu;
  alpha := radius/f_e;
  beta  := f_e*radius;
  asp_ratio := beta/alpha;
  eccentri := sqrt( 1 - asp_ratio^2 );
  E := EllipticE( eccentri );
  # Now get the flow quantities.  Final boolean ensures "active=true" in answer
  FlowQuantities( f_e, omeg, alpha, beta, asp_ratio, eccentri, E, true );
end proc:

# The passive case has constant circumference, same as a circle
# with radius recorded in the global "radius".
# Input: 
# f_e = fraction the minor axis gets squished to.  Must be in 0 < f_e < 1
#      (the coordinate change is singular at f_e = 1 and at f_e = 0 the ellipse would 
#      be infinitely long and have zero width so that's singular too).
# omeg = applied frequency in Hz.  Typically between 0.1 and 10, 
#        but can be in (0,infinity)

Passive := proc( f_e, omeg )
  local alpha, asp_ratio, beta, E, eccentri, t, g;
  global confracf, Refine, radius, rho, mu;
  eccentri := fsolve( Pi*sqrt(1-t^2)/(2*EllipticE(t))- f_e, t );
  E := EllipticE(eccentri);
  g := Pi/(2*E);
  beta  := f_e*radius;
  alpha := g*radius;
  asp_ratio := beta/alpha;
  # Now get the flow quantities. Final boolean ensures "active=false" in answer
  FlowQuantities( f_e, omeg, alpha, beta, asp_ratio, eccentri, E, false );
end proc:


# The following code is not particularly elegant.
# In particular, it uses several global variables:
# confracf = a function to evaluate Blanch's continued fraction
# Refine   = a function to use Newton's method on Blanch's continued fraction
# radius = the radius of the blood vessel
# rho = density of blood
# mu  = viscosity of blood
# N = number of terms to take for the solution

# The code also uses Digits.  For the eigenvalues internally, it temporarily triples
# Digits, but then returns it to the user setting.

# The input variables (as preprocessed by either Active or Passive) are as below.

# The output is a Record with the following fields:
#
# Ce -- Array of blendstrings for Ce[0] through Ce[N] 
# Eigenvalues -- Vector of eigenvalues of the Mathieu equation
# ExpansionCoeffs -- Array(0..N) of expansion coefficients for the solution
# I2 -- Array(0..N) of integrals of the squares of ce[k] over [0,2*Pi]
# OscillatoryFlowRate -- procedure for maximum oscillatory flow rate
# active -- boolean; true if Active case, false if not
# aspect_ratio -- beta/alpha
# ce -- Array of blendstrings for ce[0] through ce[N]
# eccentricity -- the eccentricity of the ellipse
# focus -- ellipse foci are at -focus and +focus 
# fraction -- beta = fraction*radius, amount of compression
# lambdae -- dimensionless frequency parameter
# max_flow_rate -- maximum oscillatory flow rate
# max_velocity -- maximum velocity
# max_wall_shear_stress -- minimum shear stress at wall
# min_wall_shear_stress -- maximum shear stress at wall
# omega -- imposed frequency
# q -- the Mathieu parameter for this solution
# shearstress -- procedure to compute wall shear stress (depends on eta)
# stresses -- sampled stresses at 21 equally-spaced points
# w -- procedure to compute velocity
# xi0 -- parameter value at edge of the ellipse

# Notes: the procedured returned use various local variables and 
# the parameters have to be passed in from the Record, appropriately.


FlowQuantities := proc( f_e, omeg, alpha, beta, asp_ratio, eccentri, E, activeboolean )
  local A, AE, AS, ASE, As, Ces, I2loca, Lambda_e, OsclltryFlwRt, XpansionCoeffs, 
        answer, ces, cesnosquareIntegrated, cesq, cesqIntegrated, foc, j, k, 
        lambda_e, ndegtayl, qlocal, rtol, shrstrss, sigma, strsss, w, wa, wma, xi;

  global N, Refine, confracf, mu, radius, rho;

  Digits := max( 15, Digits ); # Default 10 Digits is too low.

  # Container for all the flow quantities computed.
  #
  answer := Record( 'Ce', 'Eigenvalues', 'ExpansionCoeffs', 'I2', 
                    'OscillatoryFlowRate', 'active', 'aspect_ratio', 
                    'ce', 'eccentricity', 'focus', 'fraction', 'lambdae',
                    'max_flow_rate', 'max_velocity', 'max_wall_shear_stress', 
                    'min_wall_shear_stress', 'omega', 'q', 'shearstress', 
                    'stresses', 'w', 'xi0'
                  );

  answer:-fraction := f_e;
  answer:-omega := omeg;
  answer:-active := activeboolean;
  answer:-aspect_ratio := asp_ratio;
  answer:-eccentricity := eccentri;
  foc := alpha*eccentri;
  answer:-focus := foc;
  sigma := sqrt(2/(1/alpha^2 + 1/beta^2));
  lambda_e := rho*omeg*sigma^2/mu;
  answer:-lambdae := lambda_e;
  Lambda_e := rho*omeg*foc^2/mu;
  qlocal := -Lambda_e*I/4;
  xi[0] := arccosh(alpha/foc);
  answer:-xi0 := xi[0];
  answer:-q := qlocal;

  # Get eigenvalues; temporarily increase Digits to do so.
  Digits := 3*Digits;
  A := GetMat(4*N, qlocal);
  AE := LinearAlgebra:-Eigenvalues(evalf(A));
  ASE := convert(AE, list);
  ASE := sort(AE, (a, b) -> abs(a) <= abs(b));
  # Don't refine the eigenvalues if q is too small
  for k to N + 1 while abs(qlocal) > 0.2 do
    for j to 20 while Float(1, 4 - Digits) < abs(confracf(ASE[k],qlocal)) do
      ASE[k] := Refine(ASE[k],qlocal, 0);
    end do;
  end do;
  if abs(qlocal) > 0.2 and 
     0.1e-8 < max(seq(abs(confracf(ASE[k],qlocal)),k=1..N+1)) then
		error "Eigenvalue computation didn't wholly succeed", asp_ratio, omeg;
  end if;
  Digits := Digits/3;
  As := Array(0 .. N, [seq(ASE[j], j = 1 .. N + 1)]);
  answer:-Eigenvalues := AS;

  # Now compute the solution.
  ces := Array(0 .. N);
  Ces := Array(0 .. N);
  rtol := 0.5*Float(1,-Digits);
  ndegtayl := Digits;
  for k from 0 to N do
    wa := BothIandII(As[k], qlocal, evalf(2*Pi), 
                     'residualTolerance' = rtol, 
                     'nDeg' = ndegtayl 
                     );
    ces[k] := eval(wa[1]);
    wma := BothIandII(As[k], qlocal, xi[0]*I, 
                      'residualTolerance' = rtol,
                      'nDeg' = ndegtayl 
                      );
    Ces[k] := eval(wma[1]);
  end do;
  answer:-Ce := Ces;
  answer:-ce := ces;
        
  # Find the norms
  cesq := Array(0 .. N);
  cesqIntegrated := Array(0 .. N);
  cesnosquareIntegrated := Array(0 .. N);
  for k from 0 to N do
    cesq[k] := zipBlendstrings((a, b) -> a*b, ces[k], ces[k]);
    cesqIntegrated[k] := intBlend(cesq[k]);
    cesnosquareIntegrated[k] := intBlend(ces[k]);
  end do;
  XpansionCoeffs := Array(0 .. N);
  I2loca := Array(0 .. N);


  # Expansion coefficients in the simple eigenvalue case
        
  for k from 0 to N do
    I2loca[k] := deval(cesqIntegrated[k], 'x' = evalf(2*Pi));
    XpansionCoeffs[k] := deval(cesnosquareIntegrated[k],'x' = evalf(2*Pi))
                               /(I2loca[k]*deval(Ces[k], 'x' = xi[0]*I));
  end do;
  answer:-I2 := I2loca;
  answer:-ExpansionCoeffs := XpansionCoeffs;
  w := proc(xi, eta, Cesp, cesp, Xpansionp, lambdae )
    local k, ans, c, C;
    option remember;
	ans := -1.0;
    for k from 0 to N do
      C := deval(Cesp[k], 'x' = xi*I);
      c := deval(cesp[k], 'x' = eta);
      ans := Xpansionp[k]*C*c + ans;
    end do;
    return 4/(lambdae*I)*ans;
  end proc;
 
  answer:-w := eval(w);
  answer:-max_velocity := abs(w(0, evalf(Pi/2), Ces, ces, XpansionCoeffs, lambda_e ));
  # Should pass parameters in if we're going to return the procedure out of this scope
  OsclltryFlwRt := proc(  Cesp, cesp, Xpansionp, lambdae, I2p, xi0 )
    local  k, ans, C;
    option remember;
    ans := 0.;
    for k from 0 to N do
      C := deval(Cesp[k], 'x' = xi0*I, ':-nder' = 1);
      ans := Xpansionp[k]^2*I2loca[k]*C[0]*C[1]*I + ans; # CHAIN RULE
    end do;
    return 8*(1 - tanh(2*xi0)*ans/evalf(Pi*lambdae*I))/(lambdae*I);
  end proc;
  
  answer:-OscillatoryFlowRate :=  eval(OsclltryFlwRt);
  answer:-max_flow_rate := abs(OsclltryFlwRt( Ces, ces, XpansionCoeffs, 
                                              lambda_e, I2loca, xi[0] ));
  shrstrss := proc(eta,  Cesp, cesp, Xpansionp, lambdae, xi0, focu )
    local delta, k, ans, c, C;
    option remember;
    ans := 0.;
    delta := focu*sqrt(cosh(xi0)^2 - cos(eta)^2);
    for k from 0 to N do
      C := deval(Cesp[k], 'x' = xi0*I, ':-nder' = 1);
      c := deval(cesp[k], 'x' = eta);
      ans := Xpansionp[k]*C[1]*c*I + ans;  # CHAIN RULE
    end do;
    return 2/(lambdae*I)*ans/delta;
  end proc;
  
  answer:-shearstress := eval(shrstrss);
  # Probably not necessary.  Max always occurs at pi/2
  strsss := seq(abs(shrstrss(evalf(k*Pi/20), Ces, ces, XpansionCoeffs, lambda_e, xi[0], foc)), k = 0 .. 20);
  answer:-max_wall_shear_stress := abs( shrstrss(evalf(Pi/2), Ces, ces, XpansionCoeffs, lambda_e, xi[0], foc) );
  answer:-min_wall_shear_stress := abs( shrstrss(0.0, Ces, ces, XpansionCoeffs, lambda_e, xi[0], foc) );
  answer:-stresses := strsss;
  
  return answer;

end proc:
