% Alexander Moreno
% 2-13-2017
% 

function [VSWR] = findOpitmalSerrations_Paul(f,h,er,Zc,Lpp,maxL,maxW,dimserr)
%% initializing varaibles
Dserr = dimserr(1); % [m] (middle) serration depth
Wserr = dimserr(2); % [m] (middle) serratopm width
Lgpw = maxW;        % [m] parallel plate
M  = 50;            % [scalar] Number of segments
Zg = 50;            % [Ohms] Gen Impedance
eta0 = 377;         % [Ohms] Impedance of free space
c  = 3e8;           % [m/s] speed of light
mu0 = 1.256637e-6;  % [scalar] mu0 
e0  = 8.8542e-12;   % [scalar] e0
lambda0 = c./f;     % [1D vector] wavelength(s)
%% Compute edge serrations parameters
wdim = maxW/Wserr;
N1   = floor(wdim); % [integer] number of middle serrations
%
if(wdim==N1)
    if(N1==1)
    else
        N1=N1-1;
    end
end
% check if there are edge serrations
wm1=Wserr; % middle serration width
if(mod(Lgpw,Wser)==0)
    % there is no edge serration(s)
    Me  = 0;    % [integer] Number of segments for edge serration
    wm2 = 0;    % [scalar] (edge) serration width
    De  = 0;    % [scalar] (edge) serration depth 
else
    if(wdim==N1)
        if(N1==1)
        else
            N1=N1-1;
        end
    end
    wm2 = (Lgpw - N1*Wserr)/2;  % [scalar] (edge) serration width 
    De = ((2*Dserr)/wm1)*wm2;   % [scalar] (edge) serration depth
    N = N1+2;                   % [integer] total number of serrations
    ratio = De/Dserr;           % [scalar] ratio between middle and edge serration
    Me = ceil(ratio*M);         % [scalar] number of semegents for edge serrations 
end
Dim = [Dserr,Wserr,Lpp,maxW,maxL,h,er,De,wm2];
Var = [M,Me,lambda0];
%% Compute Input Impedance and Segement Parameters
[Zin, Zserr] = Impedance

%%
end