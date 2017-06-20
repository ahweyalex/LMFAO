% Paul Winniford
% taken from the engaging classic, "Characteristics of Some Asymmetrical
% Coupled Transmission Lines" by SS Bedair


function [Casym]=asymmetric_line(W_L,W_R,S,epsr,h)
% notation matches that of the paper. left and right line widths, spacing,
% permittivity, and substrate height are the inputs. A structure with the
% capacitances associated with two lines of asymmetric width is the output.
% Casym.LR is mutual capacitance, Casym.LL is self capacitance (with the
% ground plane) and the "a" postscript indicates values associated with an
% air dielectric structure


[C_pLa,C_pRa,C_fsLa,C_fsRa,C_primefsLa,C_primefsRa,C_gd,C_fLa,C_fRa,C_cps1]...
    =asymmetrical_sub(W_L,W_R,S,1,h);
%air gap structure
Casym.LLa=C_fLa+C_pLa+C_primefsLa+2*C_gd;
Casym.RRa=C_fRa+C_pRa+C_primefsRa+2*C_gd;
Casym.LRa=2*C_gd;

del_C_f=(C_fLa-C_fsLa)*(C_fRa-C_fsRa)/(C_fLa+C_fRa-C_fsLa-C_fsRa); %eqn20
[C_pL,C_pR,~,~,C_primefsL,C_primefsR,C_gd,C_fL,C_fR,C_cps2]...
    =asymmetrical_sub(W_L,W_R,S,epsr,h);
% dielectric structure
C_cps=C_cps1+C_cps2;
C_ga=0.5*C_cps-del_C_f; %eqn21
Casym.LL=C_fL+C_pL+C_primefsL+C_gd+C_ga; %eqn 15
Casym.RR=C_fR+C_pR+C_primefsR+C_gd+C_ga;
Casym.LR=C_gd+C_ga;
end

function [C_pL,C_pR,C_fsL,C_fsR,C_primefsL,C_primefsR,C_gd,C_fL,C_fR,...
    C_cps]=asymmetrical_sub(W_L,W_R,S,epsr,h)
% called twice by asymmetrical master, with dielectric and with air gap
eps0=8.85e-12;
c_v=3e8;
% "Some definitions of the complete elliptic integrals use the modulus
%     k instead of the parameter M.  They are related by M = k^2."
% eqn 12, I use the notation of pi and c rather than odd and even

% k is probably the problem? [Alex]
% when w/h is higher than 6.5 tanh(pi/4*w/h) becomes 1
k_piL=tanh(pi/4*W_L/h)./tanh(pi/4*(S+W_L)/h);   %pi mode (odd), left line
k_cL=tanh(pi/4*W_L/h)*tanh(pi/4*(S+W_L)/h);     %c mode (even), left line
k_piR=tanh(pi/4*W_R/h)./tanh(pi/4*(S+W_R)/h);   %pi mode (odd),right line
k_cR=tanh(pi/4*W_R/h)*tanh(pi/4*(S+W_R)/h);     %c mode (even),right line
k_sL=tanh(pi/4*W_L/h);                          %single line,left line, eqn 13
k_sR=tanh(pi/4*W_R/h);                          %single line,right line

keyboard

% The following six lines have a nasty tendency towards Infinity. Not
% everyone should shoot for the stars. Specifically, all of the numerator
% elliptic functions have this tendency. But that's because their inputs
% are 1! The problem is really in the previous lines defining k
C_OsL=1/(c_v*60*pi)*ellipke(k_piL)./ellipke(sqrt(1-k_piL.^2));
C_OsR=1/(c_v*60*pi)*ellipke(k_piR)./ellipke(sqrt(1-k_piR.^2));
C_EsL=1/(c_v*60*pi)*ellipke(k_cL)./ellipke(sqrt(1-k_cL.^2));
C_EsR=1/(c_v*60*pi)*ellipke(k_cR)./ellipke(sqrt(1-k_cR.^2));
C_sL=1/(c_v*60*pi)*ellipke(k_sL)./ellipke(sqrt(1-k_sL.^2));%based on stripline
C_sR=1/(c_v*60*pi)*ellipke(k_sR)./ellipke(sqrt(1-k_sR.^2));%
%keyboard
if(isinf(C_OsL)|isinf(C_OsR))
    sprintf('Error: W/h ratio out of bounds!')
    return
end
% C_sL=single_line(W_L,epsr,h);%which ones!
% C_sR=single_line(W_R,epsr,h);%Pozar microstrip formulations
C_pL=eps0*epsr*W_L/h; %eqn 8
C_pR=eps0*epsr*W_R/h;
C_fsL=0.5*(C_sL-C_pL); %eqn 9
C_fsR=0.5*(C_sR-C_pR);
C_primefsL=C_EsL-C_pL-C_fsL;
C_primefsR=C_EsR-C_pR-C_fsR;
C_gdL=0.5*(C_OsL-C_EsL);
C_gdR=0.5*(C_OsR-C_EsR);
C_gd=sqrt(C_gdL.*C_gdR);

C_mL=epsre(W_L,h,epsr)/(c_v*Z_m(W_L,h)); %eqn 17
C_mR=epsre(W_R,h,epsr)/(c_v*Z_m(W_R,h));
C_fL=0.5*(C_mL-C_pL); %eqn 16
C_fR=0.5*(C_mR-C_pR);

k_squared=(1+W_L./S+W_R./S)./((1+W_L./S).*(1+W_R./S));
C_cps=eps0*epsr*ellipke(1-k_squared)./ellipke(k_squared);
end

function [Z]=Z_m(W_Y,h)
% equation 18, consider replacing with Pozar
if (W_Y/h<=1)
    Z=60*log(8*h/W_Y+W_Y/(4*h));
else
    Z=120*pi/(W_Y/h+1.393+0.677*log(W_Y/h+1.444));
end
end

function [eps]=epsre(W_Y,h,epsr)
% equation 19, consider replacing with Pozar
if(W_Y/h<=1)
    F=(1+12*h/W_Y)^(-1/2)+0.04*(1-W_Y/h)^2;
else
    F=(1+12*h/W_Y)^(-1/2);
end
eps=0.5*((epsr+1)+(epsr-1)*F);
end