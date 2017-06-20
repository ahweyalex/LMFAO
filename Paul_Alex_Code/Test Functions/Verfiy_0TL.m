clear all; close all; clc;

S = magic(6);
S(1,:)=0; S(:,1)=0; S(end,:)=0; S(:,end)=0; 
%S = S/30 + 1j*S/20;
S(1,end)=1; S(end,1)=1;
S0 = S;
S = S(2:end-1,2:end-1);

P = [1,3,4,5,6];
ZL = 377*ones(1,numel(P));
S_0TL = loadSPorts(S0,P,ZL);
P = [2,3,4];
ZL = 377*ones(1,numel(P));
S_TL = loadSPorts(S,P,ZL);

comp = S_0TL == S_TL