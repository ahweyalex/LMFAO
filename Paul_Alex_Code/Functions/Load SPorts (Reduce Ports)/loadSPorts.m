% Alexander Moreno
% 01-19-2017
% ver 1.0
% Reduces an N-Port (S-Parameter matrix) to an M-Port (S-Parameter matrix),
% where N>M, by loading x-amount of ports.
% 
% Inputs:
%   [nS]:  (NxN S-Parameter matrix)
%   [P]:  1-D vector indicating which ports to load, entries should be
%   increasing order example: P = [1 4 6] refering to ports 1,4,and 6.
%   [ZL]: [DxD unity matrix] Load used upon the speficied ports in [P],
%   these vectors should be one to one. refering that first entry [P]
%   corresponds to first entry of [ZL]
%
% Outputs:
%   [mS]: (MxM S-Parameter matrix) 

%% convert to s-parameter (ver 1.0)
% this assumes that you loading the last ports from the S-Parameter
% Correction* it assumes you are loading the from the bottom and right
% example
% function [mS] = loadSPorts(S,P,ZL)
function [mS] = loadSPorts(S,P,ZL)
%{
[ZLrow, ZLcol] = size(ZL);
if(ZLrow~=ZLcol)
    msg = 'ZL matrix is not a square matrix';
    error(msg);
end
%}
pSize = numel(P); % number of ports to be loaded
N = numel(ZL);
if(N~=pSize)
    msg = 'Number of Ports to be loaded does not correspond with ZL dimensions';
    error(msg)
end
nSize = numel(S(1,:));
mSize = abs(nSize - pSize); % total number of ports after loading
ZL = diag(ZL);
Z0 = eye(pSize)*50; % assumes Z0 is 50 ohms
Gamma = (ZL-Z0)/(ZL+Z0);

% fill values for submatrices
nSize = numel(S(:,1));
pSize = numel(P);
mSize = abs(pSize - nSize);
M = 1:nSize;
M(P) = []; % remove entrices 

P11 = S(M,M); % mxm (row x col)
P12 = S(M,P); % mxp
P21 = S(P,M); % pxm
P22 = S(P,P); % pxp

% compute reduce S-Parameter matrix
mS = P11 + (P12*Gamma)*((eye(N)*P22*Gamma)^-1)*P21;
end