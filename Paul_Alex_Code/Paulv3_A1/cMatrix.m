% Paul Winniford
% function file formerly known as 'meowzillaHarbingerOfDeath.m'

function [C,L]=cMatrix(w,h,s,epsr,f,nonadjacent,fdependent)
% this function combines capacitance calculations from "Accurate Wide Range
% Design Equations for..." (Kirschining and Jansen) with "Characteristics
% of some asymmetrical coupled Transmission Lines" (Bedair) to calculate 
% the capacitance matrix for an arbitrary n (parallel) line structure
% A vector of line widths (dimension n), substrate height, vector of line 
% spacings (n-1), and substrate permittivity are inputs
% Outputs are the PUL C matrix and L matrix, where L is calculated from the
% C matrix when the dielectric is replaced with air
% The stated bounds on accurate capacitance calculation is from 0.5-20 of
% the width/spacing normalized to the substrate height.
% "nonadjacent" is a logical that specifies taking nonadjacent coupling
% into account, which is what I recommend
% frequency dependent calculations can result in NaN output from kirschning

mu0=pi*4e-7;
eps0=8.85e-12;
n=numel(w);
g=s/h; %normalized spacing
u=w'/h; %normalized width
Cs=zeros(n); %self capacitance
Cas=zeros(n); %air dielectric
Cm=zeros(n); %mutual capacitance
Cam=zeros(n); 
if nonadjacent
    for idx=1:n-1 %over reference line
        spacing=cumsum(s(idx:(n-1))+[0 w(idx+1:(n-1))]);
        for jdx=(idx+1):n %cross coupling with other lines
            if w(idx)==w(jdx)
                % Alex 06-19-2017
                % g=s/h; u=w/h; %
                % Valid from 0.1<=u<=10, 0.1<=g<=10, 1<=er<=18
                if(u(idx)>=10) % out of bounds
                    u(idx) = 10;
                end
                g = spacing(jdx-idx)/h;
                if(g>=10)
                    g=10;
                elseif(g<=0.1)
                    g=0.1;
                end 
                [epsre,Z0m]=single_ustrip(epsr,w(idx),h);
                [C0]=kirschning(epsr,epsre,u(idx),g,Z0m,h,f,fdependent); 
                kir = 'kirschning' %alex
                Cs((jdx-1),idx)=C0.s;
                Cm(jdx,idx)=C0.m;
                Cas((jdx-1),idx)=C0.sa;
                Cam(jdx,idx)=C0.ma;
                if idx==(n-1)
                    Cs(n,n)=C0.s; %fills in just the last element
                    Cas(n,n)=C0.sa;
                end
                %keyboard %alex
            else
                %keyboard
                [Casym]=asymmetric_line(w(idx),w(jdx),spacing(jdx-idx),epsr,h);
                a_L = 'asymmetric_line'
                Cs((jdx-1),idx)=Casym.LL; 
                Cm(jdx,idx)=Casym.LR;
                Cs(jdx,idx+1)=Casym.RR;
                Cas((jdx-1),idx)=Casym.LLa;
                Cam(jdx,idx)=Casym.LRa;
                Cas(jdx,idx+1)=Casym.RRa;
% this is often overwritten but is necessary to fill the last index.
                %keyboard %alex
            end
        end
    end
else %using only adjacent line capacitance
    for idx=1:n-1 %over reference line
    spacing=cumsum(s(idx:(n-1))+[0 w(idx+1:(n-1))]);
    jdx=(idx+1); %cross coupling with other lines
        if w(idx)==w(jdx)
            % Alex 06-19-2017
            % g=s/h; u=w/h; %
            % Valid from 0.1<=u<=10, 0.1<=g<=10, 1<=er<=18
            if(u(idx)>=10) % out of bounds
                u(idx) = 10;
            end
            g = spacing(jdx-idx)/h;
            if(g>=10)
                g=10;
            elseif(g<=0.1)
                g=0.1;
                end
            [epsre,Z0m]=single_ustrip(epsr,w(idx),h);
            [Z0,C0]=kirschning(epsr,epsre,u(idx),spacing(jdx-idx)/h,Z0m,h,f);
            kir = 'kirschning' %alex
            Cs((jdx-1),idx)=C0.s;
            Cm(jdx,idx)=C0.m;
            Cas((jdx-1),idx)=C0.sa;
            Cam(jdx,idx)=C0.ma;
            if idx==(n-1)
                Cs(n,n)=C0.s; %fills in just the last element
                Cas(n,n)=C0.sa;
            end
            %keyboard%alex
        else
            [Casym]=asymmetric_line(w(idx),w(jdx),spacing(jdx-idx),epsr,h);
            a_l = 'asymmetric_line'
            Cs((jdx-1),idx)=Casym.LL;
            Cm(jdx,idx)=Casym.LR;
            Cs(jdx,idx+1)=Casym.RR;
            Cas((jdx-1),idx)=Casym.LLa;
            Cam(jdx,idx)=Casym.LRa;
            Cas(jdx,idx+1)=Casym.RRa;
%             this is often overwritten but is necessary to fill the last
%             index.
        end
    end
end
%keyboard
Cm=Cm'.*(1-eye(n))+Cm;
Cam=Cam'.*(1-eye(n))+Cam;
C=Cs.*eye(n)+ones(n,1)*sum(Cm,1).*eye(n)-Cm;
Ca=Cas.*eye(n)+ones(n,1)*sum(Cam,1).*eye(n)-Cam;
% C=Cs.*eye(n)-Cm;
% Ca=Cas.*eye(n)-Cam;
%keyboard
L=mu0*eps0*Ca^-1; 
%L1=mu0*eps0*pinv(Ca);  %alex
%L2=Ca.\(mu0*eps0);     %alex
% gmres     % Alex
%keyboard
end

function [epsre,Z0]=single_ustrip(epsr,w,h)
% calculates effective permittivity and characteristic impedance for a
% single microstrip line, Pozar section 3.8 of the 4th edition
epsre=(epsr+1)/2+(epsr-1)/2*1/sqrt(1+12*h/w);
if w/h<=1
    Z0=60/sqrt(epsre)*log(8*d/w+w/(4*d));
else
    Z0=120*pi/(sqrt(epsre)*(w/h+1.393+0.667*log(w/h+1.444)));
end
end
    
function [C0]= kirschning(epsr,epsre,u,g,Z0m,h,varargin)
% Kirschning and Jansen's "Accurate Wide Range Design Equations for...".
% Cell outputs are self and mutual capacitance and impedance, with the
% capacitance given both in terms of air and dielectric values
% g=s/h; u=w/h;
% Valid from 0.1<=u<=10, 0.1<=g<=10, 1<=er<=18
%frequency dependent calculations can result in NaN output, so I strongly
%disrecommend using freq dependent calc
c=3e8; 
eps0=8.85e-12;
mu0=pi*4e-7;
v=u*(20+g.^2)./(10+g.^2)+g.*exp(-g);
ae=1 + log((v.^4+(v/52).^2)./(v.^4+0.432))/49 + log(1+(v/18.1).^3)/18.7;
be=0.564*((epsr-0.9)/(epsr+3))^0.053;
epsre_e=(epsr+1)/2+(epsr-1)/2*(1+10./v).^(-ae*be); 
% effective static even mode permittivity
a0=0.7287*(epsre-0.5*(epsr+1))*(1-exp(-0.179*u));
b0=0.747*epsr/(0.15+epsr);
c0=b0-(b0-0.207)*exp(-0.414*u);
d0=0.593+0.694*exp(-0.562*u);
epsre_o= epsre +(0.5*(epsr+1)-epsre+a0)*exp(-c0*(g.^d0));
% effective static odd mode permittivity
% %{
if varargin{2} %optional, buggy, frequency dependent epsilon calculation
    f=varargin{1};
    fn=f/1e9*h/1e-3; %normalized freq, f and h in GHz and mm
    P1=0.27488+(0.6315+0.525/(1+0.0157*fn)^20)*u-0.065683*exp(-8.7513*u);
    P2=0.33622*(1-exp(-0.03442*epsr));
    P3=0.0363*exp(-4.6*u)*(1-exp(-(fn/38.7)^4.97));
    P4=1+2.751*(1-exp(-(epsr/15.916)^8));
    P5=0.334*exp(-3.3*(epsr/15)^3)+0.746;
    P6=P5*exp(-(fn/18)^0.368);
    P7=1+4.069*P6*g^0.479*exp(-1.347*g^0.595-0.17*g^2.5);
    Fe=P1*P2*((P3*P4+0.1844*P7)*fn)^1.5763;
    eps_effe=epsr-(epsr-epsre_e)/(1+Fe);
    P8=0.7168*(1+1.076/(1+0.0576*(epsr-1)));
    P9=P8-0.7913*(1-exp(-(fn/20)^1.424))*atan(2.481*(epsr/8)^0.946);
    P10=0.242*(epsr-1)^0.55;
    P11=0.6366*(exp(-0.3401*fn)-1)*atan(1.263*(u/3)^1.629);
    P12=P9+(1-P9)/(1+1.183*u^1.376);
    P13=1.695*P10/(0.414+1.605*P10);
    P14=0.8928+0.1072*(1-exp(-0.42*(fn/20)^3.215));
    P15=abs(1-0.8928*(1+P11)*P12*exp(-P13*g^1.092)/P14);
    Fo=P1*P2*((P3*P4+0.1844)*fn*P15)^1.5763;
    eps_effo=epsr-(epsr-epsre_o)/(1+Fo);
end
% %}
Q1=0.8695*u^0.194;
Q2=1+0.7519*g+0.189*g.^2.31;
Q3=0.1975+(16.6+(8.4./g).^6).^-0.387+1/241*log(g.^10./(1+(g/3.4).^10));
Q4=2*Q1./Q2.*1./(u.^Q3.*exp(-g)+(2-exp(-g)).*u.^-Q3);
ZLe=Z0m*sqrt(epsre./epsre_e)./(1-Q4.*sqrt(epsre)*Z0m/377); %Z0m same as Z0
%characteristic line impedance, even mode, static
Q5=1.794+1.14*log(1+0.638./(g+0.517*g.^2.43));
Q6=0.2305+1/281.3*log(g.^10./(1+(g/5.8).^10))+1/5.1*log(1+0.598*g.^1.154);
Q7=(10+190*g.^2)./(1+82.3*g.^3);
Q8=exp(-(6.5+0.95*log(g)+(g/0.15).^5)); %this expression must be wrong, always=0
Q9=log(Q7).*(Q8+1/16.5);
Q10=Q4-Q5./Q2.*exp(Q6.*log(u)./u.^Q9);
ZLo=Z0m*sqrt(epsre./epsre_o)./(1-Q10*sqrt(epsre)*Z0m/377);
%Z0 is static characteristic impedance for ustrip line
Ce=1/c*sqrt(epsre_e)./ZLe;
Co=1/c*sqrt(epsre_o)./ZLo;
Ce_a=1./(c*ZLe.*sqrt(epsre_e));
Co_a=1./(c*ZLo.*sqrt(epsre_o));

C0.s=Ce;
C0.m=(Co-Ce)/2;
C0.sa=Ce_a; 
C0.ma=(Co_a-Ce_a)/2;
Z0.e=ZLe;
Z0.o=ZLo;

% %{
if varargin{2} %optional frequency-dependent calc
    Q11=0.893*(1-0.3/(1+0.7*(epsr-1)));
    Q12=2.121*((fn/20)^4.91/(1+Q11*(fn/20)^4.91))*exp(-2.87*g)*g^0.902;
    Q13=0.038*(epsr/8)^5.1; %this can yield low values ~10^-5, which can cause Q15 to be NaN
    Q14=1+1.203*(epsr/15)^4/(1+(epsr/15)^4);
    Q15=1.887*exp(-1.5*g^0.84)*g^Q14*(1+0.41*(fn/15)^3*u^(2/Q13)/(0.15+u^(1.626/Q13)))^(-1); %NaN, maybe because of Q13
    Q16=(1+9/(1+0.403*(epsr-1)^2))*Q15; %will be NaN if Q15 is NaN
    Q17=0.394*(1-exp(-1.47*(u/7)^0.672))*(1-exp(-4.25*(fn/20)^1.87));
    Q18=0.61*(1-exp(-2.13*(u/8)^1.593))/(1+6.554*g^4.17);
    Q19=0.21*g^4*((1+0.18*g^4.9)*(1+0.1*u^2)*(1+(fn/24)^3))^(-1);
    Q20=(0.09+1/(1+0.1*(epsr-1)^2.7))*Q19;
    Q21=abs(1-42.52*g^0.133*exp(-0.812*g)*u^2.5/(1+0.033*u^2.5));
    re=(fn/28.843)^12; %can be very small, ~10^-23
    qe=0.016+(0.0514*epsr*Q21)^4.524;
    pe=4.766*exp(-3.228*u^0.641);
    Ce=1+1.275*(1-exp(-0.004625*pe*epsr^1.674*(fn/18.365)^2.745))-Q12+Q16-Q17+Q18+Q20; %can be NaN due to Q16
    de=5.086*qe*(re/(0.3838+0.386*qe))*(exp(-22.2*u^1.92)/(1+1.2992*re))*((epsr-1)^6/(1+10*(epsr-1)^6)); %can be very small because of re, ~10^-100
% ZLf is the frequency dependent characteristic impedance of a single
% microstrip line
    [Q0,ZLf]=freqDependentuStripZ0(epsr,epsre,u,Z0m,h,f);
    Zfe=ZLe*(0.9408*epsre^Ce-0.9603)^Q0*1/((0.9408-de)*epsre^Ce-0.9603)^Q0; %Ce and de can make this NaN
    Q29=15.16/(1+0.196*(epsr-1)^2);
    Q26=30-22.2*(((epsr-1)/13)^12/(1+3*((epsr-1)/13)^12))-Q29;
    Q27=0.4*g^0.84*(1+2.5*(epsr-1)^1.5/(5+(epsr-1)^1.5));
    Q28=0.149*(epsr-1)^3/(94.5+0.038*(epsr-1)^3);
    Q22=0.925*(fn/Q26)^1.536/(1+0.3*(fn/30)^1.536);
    Q23=1+0.005*fn*Q27*((1+0.812*(fn/15)^1.9)*(1+0.025*u^2))^(-1);
    Q24=2.506*Q28*u^0.894*((1+1.3*u)*fn/99.25)^4.29*(3.575+u^0.894)^(-1);
    Q25=(0.3*fn^2/(10+fn^2))*(1+2.333*(epsr-1)^2/(5+(epsr-1)^2));
% Zl(fn) needs to be included, expressed as Zlf
    Zfo=ZLf+(ZLo*(eps_effo/epsre_o)^Q22-ZLf*Q23)*(1+Q24+(0.46*g)^2.2*Q25)^(-1);

    Ce=1/c*sqrt(eps_effe)./Zfe; %Zfe can make this NaN
    Co=1/c*sqrt(eps_effo)./Zfo;
    Ce_a=1./(c*Zfe.*sqrt(eps_effe)); %Zfe here as well
    Co_a=1./(c*Zfo.*sqrt(eps_effo));

    C0.s=Ce; %NaN errors will likely cause all C0 to be NaN
    C0.m=(Co-Ce)/2;
    C0.sa=Ce_a; 
    C0.ma=(Co_a-Ce_a)/2;
    Z0.e=ZLe; %Z0 is likely unaffected by NaN errors
    Z0.o=ZLo;
end
% %}
end

function [R17,ZLf]=freqDependentuStripZ0(epsr,epsre,u,Zl,h,f)
% Taken from "Arguments and an Accurate Model for the Power-Current
% Formulation of Microstrip Characteristic Impedance"
[epsref]=freqDependentuStripEpsilon(epsr,epsre,h,f,u);
fn=f/1e9*h/1e-3; %normalized freq, f and h in GHz and mm
R1=0.03891*epsr^1.4;
R2=0.267*u^7;
R3=4.766*exp(-3.228*u^0.641);
R4=0.016+(0.0514*epsr)^4.524;
R5=(fn/28.843)^12;
R6=22.2*u^1.92;
R7=1.206-0.3144*exp(-R1)*(1-exp(-R2));
R8=1+1.275*(1-exp(-0.004625*R3*epsr^1.674*(fn/18.365)^2.745));
R9=5.086*R4*R5/(0.3838+0.386*R4)*exp(-R6)/(1+1.2992*R5)*(epsr-1)^6/(1+10*(epsr-1)^6);
R10=0.00044*epsr^2.136+0.0184;
R11=(fn/19.47)^6/(1+0.0962*(fn/19.47)^6);
R12=1/(1+0.00245*u^2);
R13=0.9408*epsref^R8-0.9603;
R14=(0.9408-R9)*epsre^R8-0.9603;
R15=0.707*R10*(fn/12.3)^1.097;
R16=1+0.0503*epsr^2*R11*(1-exp(-0.026*fn^1.15656-R15));
R17=R7*(1-1.1241*R12/R16*exp(-0.026*fn^1.15656-R15));
ZLf=Zl*(R13/R14)^R17;
end

function [epsref]=freqDependentuStripEpsilon(epsr,epsre,h,f,u)
% taken from "Accurate Model for Effective Dielectric Constant of
% Microstrip with Validity up to Millimetre-Wave Frequencies"
% The output is the frequency dependent effective permittivity of a single
% microstrip line
%h*f is normalized frequency, given in GHz*cm
h=h*100; %convert to cm
P1=0.27488+(0.6315+0.525/(1+0.157*f*h)^20)*u-0.065683*exp(-8.7513*u);
P2=0.33622*(1-exp(-0.0344*2*epsr));
P3=0.0363*exp(-4.6*u)*(1-exp(-(f*h/3.87)^4.97));
P4=1+2.751*(1-exp(-(epsr/15.916)^8));
P=P1*P2*((0.1844+P3*P4)*10*f*h)^1.5763;
epsref=epsr-(epsr-epsre)/(1+P);
end

% function [cnotlight]=aStupidFunctionforAndy(w,h)
% 
% c=0:0.01:20;
% stupid1=floor(pi/2*w/h*1e2); %varies from  0.1571 to 15.7080
% stupid2=floor((c-asinh(c))*1e2);
% stupidIndex=stupid1==stupid2;
% creturn=c(stupidIndex); 
% if(isempty(creturn))
%     sprintf('Error: Out of bounds w/h ratio!')
%     return
% else
% cnotlight=creturn(1); %in case there are multiple values
% end
% end