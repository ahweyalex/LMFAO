% Alexander Moreno
% 01-04-2016
%
%{
clear all; close all; clc;
%% Define values
Wm = 30e-3;    % middle serration width [m]
We = 3e-3;    % edge serration width [m]
Ds = 70e-3;     % middle serration depth [m]
M = 50;         % number of segments
w  = [We Wm We];
h =  0.001524;  % substrate height [m] 5mil
% h=0.762e-3; %substrate height
s0 = 21.5e-3;
s  = [s0 s0 s0];   % seperation between serrations [m], n-1 
er = 2.2;   %relative permittivity
f0 = 433e6; % frequency [Hz]
l  = Ds/M;  % line length

%% Obtain Z-Parameters for each segment
[Zpaul] = ZParameters(w,h,s,er,f0,l);
%}

%%
clear all; close all; clc;
% example(s)
Ds = 70e-3;     % middle serration depth [m]
M = 50;
er=2.2; %relative permittivity
f0= 433e6; %frequency
%f0 = 18e9;

%w=[2e-3 5e-3 2e-3]; % error 1
%w=[20e-3 20e-3 20e-3]; % error 2
%w=[4e-3 5e-3 5e-3 4e-3];

w=[41e-3 50e-3 41e-3]; % A
s=[4.6e-3 4.6e-3]; %A

%w=[4e-3 5e-3 4e-3]; %B, error 1
%s=[4e-3 4e-3]; %B, error 1

%w=[4e-3 5e-3 4.1e-3]; %C, error 1
%s=[4e-3 4e-3]; %C, error 1

%s=[50.1e-3 40e-3 50.1e-3]; %E, error 1
%w=[4e-3 5e-3 5e-3 4e-3];	%E, error 1

%w=[2e-3 5e-3 2e-3]; %E 
%s=[4.6e-3 4.6e-3]; %E

%w=[10e-3 10e-3 10e-3].*2;  % x, error 1
%s=[5e-3 5e-3 5e-3 5e-3];    % x, error 1
%w = [3e-3 3e-3 3e-3];
%s = [30e-3 30e-3];

l= Ds/M;
%h = 0.00254; % 100mil
%h= 0.001524; %substrate height % 60 mil
%h=0.001016; % 40mil
%h=0.000889; % 35mil
%h=0.000762; % 30mil
%h=0.000381; %15mil
h = 0.000254; %10mil
%h=0.000127; % 5mil
%h=0.0002

Z0 = 50;
ZL = 377;
[Zpaul] = ZParameters(w,h,s,er,f0,l);
%% convert to s-parameter (ver 1.0)
% this assumes that you loading the last ports from the S-Parameter
% Correction* it assumes you are loading the from the bottom and right
% example
s = z2s(Zpaul);
zn = numel(Zpaul(1,:));
N = 1; % number of ports
n = abs(N-zn);
Gamma = eye(N)*(ZL-Z0)/(ZL+Z0);
P11 = s(1:n,1:n);
P12 = s(1:n,(n+1):end);
P21 = s((n+1):end,1:n);
P22 = s((n+1):end,(n+1):end);
sR  = P11 + (P12*Gamma)*((eye(N)*P22*Gamma)^-1)*P21;
%% convert to T-Parameters (ver 1.0)
% assume: balanced T-Parameters i.e. same number of in and out ports
% and S-Parameter is even by even large
% example
tn = zn/2;

SI_I = s(1:tn,1:tn);
SI_II = s(1:tn,(tn+1):end);
SII_I = s((tn+1):end,1:tn);
SII_II = s((tn+1):end,(tn+1):end);

TI_I = SI_II - (SI_I*(SII_I^-1)*SII_II);
TI_II = SI_I*(SII_I^-1);
TII_I = -(SII_I^-1)*SII_II;
TII_II = SII_I^-1;
T = [TI_I, TI_II;TII_I, TII_II];

