% Alexander Moreno
% 01-13-2017
% 
% Converts S-Parameter matrix into a T-Parameter matrix.
% Reference: Multiport S-Parameter and T-Parameter Conversion With Symmetry 
%            Extension by James Frei, Xiao-Ding Cai, and Stephen Muller
% Input Parameters:
%   s: (S-Parameter) [NxN matrix]
%   m: number of entrance ports [scalar]   
%   n: number of exit ports [scalar]   
%
% Output Parameters:
%   [T] is [nxm matrix] 
%
function [T]=StoT(s,m,n)
%
% initial error handling
[sm,sn] = size(s);
if(sm~=sn)
    msg = 'S-Parameter variable is not NxN';
    error(msg)
elseif (sm<m+n)
    msg = 'Number of T-entrance and T-exit exceed number of ports.';
    error(msg)
elseif(sm<m)
    msg = 'Number of T-entrance exceed number of S-ports.';
    error(msg)
elseif(sn<n)
    msg = 'Number of T-exit exceed number of S-ports.';
    error(msg)
end

if(sm==sn && m==n && mod(sm,2)==0) % balanced, ie S-Parameter NxN,N is even
    tn=sm/2;
    SI_I = s(1:tn,1:tn);
    SI_II = s(1:tn,(tn+1):end);
    SII_I = s((tn+1):end,1:tn);
    SII_II = s((tn+1):end,(tn+1):end);

    TI_I = SI_II - (SI_I*(SII_I^-1)*SII_II);
    TI_II = SI_I*(SII_I^-1);
    TII_I = -(SII_I^-1)*SII_II;
    TII_II = SII_I^-1;
    T = [TI_I, TI_II;TII_I, TII_II];
    
    keyboard
elseif(m>n) 
    % entrance ports(m)>(n) exit ports
    % [SI,I] and [SII,II] are square
    % [SI,I] largest
    % [SII,II] smallest
    
    % working on, this isnt correct atm (5/10/17)
    SI_I   = s(1:m,1:m);  % square, largest
    SI_II  = s(1:m,m+1:end);
    SII_I  = s(m+1:end,1:m);
    SII_II = s(m+1:end,m+1:end);   % square, smallest
    
    keyboard
    
    % 05-10-2017
    TI_I   = SI_II - (SI_I*pinv(SII_I)*SII_II);
    TI_II  = SI_I*pinv(SII_I);
    TII_I  = -pinv(SII_I)*SII_II;
    TII_II = pinv(SII_I);
    T = [TI_I, TI_II;TII_I, TII_II];    
    
    % 01-13-2017
    %{
    tn=sm;
    SI_I   = s(1:(tn-1),1:(tn-1));  % square, largest
    SI_II  = s(1:(tn-1),end);
    SII_I  = s(tn,1:(tn-1));
    SII_II =  s(tn,tn);   % square, smallest
    
    TI_I   = SI_II - (SI_I*(SII_I^-1)*SII_II);
    TI_II  = SI_I*(SII_I^-1);
    TII_I  = -(SII_I^-1)*SII_II;
    TII_II = SII_I^-1;
    T = [TI_I, TI_II;TII_I, TII_II];
    %}
elseif(m<n) 
    % exit ports > entrance ports,S-Parameter is NxN, N is odd 
    % this may need to be double checked if done correctly
    % refering to m and n should maybe play into how submatrices are
    % created
    tn=sm;
    SI_I   = s(1:m,1:m)  % square, smallest
    SI_II  = s(1:m,m+1:end)
    SII_I  = s(m+1:end,1:m)
    SII_II = s(end-m:end,end-m:end)   % square, largest
    
    % 05-10-2017
    TI_I   = SI_II - (SI_I*pinv(SII_I)*SII_II);
    TI_II  = SI_I*pinv(SII_I);
    TII_I  = -pinv(SII_I)*SII_II;
    TII_II = pinv(SII_I);
    T = [TI_I, TI_II;TII_I, TII_II];

    % 01-13-2017
    %{
    %TI_I   = SI_II - (SI_I*(SII_I^-1)*SII_II);
    %TI_II  = SI_I*(SII_I^-1);
    %TII_I  = -(SII_I^-1)*SII_II;
    %TII_II = SII_I^-1;
    %T = [TI_I, TI_II;TII_I, TII_II];
    %}
end


end