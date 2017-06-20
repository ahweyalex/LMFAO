% Alexander Moreno
% 01-14-2017
% ver 1.0 
%
% Converts T-Parameter matrix into a S-Parameter matrix.
% Reference: Multiport S-Parameter and T-Parameter Conversion With Symmetry 
%            Extension by James Frei, Xiao-Ding Cai, and Stephen Muller
% Input Parameters:
%   s: (T-Parameter) [nxm matrix]
%
% Output Parameters:
%   [S] is [NxN matrix] 

function [S]=TtoS(T)

[rowSize,colSize] = size(T); % dimension of T-Parameter
rw=rowSize/2; 
cl=colSize/2;

TI_I   = T(1:rw,1:cl);
TI_II  = T(1:rw,(cl+1):end);
TII_I  = T((rw+1):end,1:cl);
TII_II = T((rw+1):end,(cl+1):end);

SI_I   = TI_II*pinv(TII_II);  
SI_II  = TI_I-(TI_II)*pinv(TII_II)*TII_I;
SII_I  = pinv(TII_II);
SII_II = -pinv(TII_II)*TII_I;

S = [SI_I, SI_II; SII_I, SII_II];
end