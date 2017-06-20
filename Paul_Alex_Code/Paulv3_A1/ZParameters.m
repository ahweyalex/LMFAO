% Paul Winniford
% Composed on the eighteenth and nineteenth of the eleventh month,
% November, in the year of our Lord, two-thousand and fifteen

function [Z_paul]=ZParameters(w,h,s,epsr,f0,l)
% calculates Z parameters for a 2n port network consisting of n-coupled
% lines based on my own work. 

% epsr=2.2; %relative permittivity
% f0=0.5e9:0.001e9:5e9; %frequency
% w=0.381e-3*[1 1 1 1 1 1 1 1 1 1]; %width, defning n lines
% s=0.381e-3*[1 1 1 1 1 1 1 1 1]; %spacing, defining n-1 spaces
% l=50e-3; %line length
% h=0.762e-3; %substrate height
n=numel(w);
Z_paul=zeros(2*n,2*n,numel(f0));
for idx=1:length(f0)
    [Z_paul(:,:,idx)]=singleFrequency(w,h,s,epsr,f0(idx),l);
end

end

function [Z]=singleFrequency(w,h,s,epsr,f0,l)
% -eigval is gamma squared, eigvec is the modal voltages (eigenpotentials).
% If requested, the left eigenvectors would give the modal currents. Since
% both the eigenpotentials and eigencurrents are normalized quantities, it
% is better to only use one, and then use either the characteristic
% admittance or characteristic impedance to convert to the other, in order
% to keep a meaningful magnitude relationship between the two
% -Characteristic admittance. This step took effort to derive. Please treat
% it with respect. The characteristic impedance is the elementwise inverse
% of the characteristic impedance.
% -The "phase" terms below are the exponential terms of the modal forward
% (negative gamma) and backward (positive gamma) travelling waves for the 
% multiline system, evaluated at x=0 and x=l
n=numel(w);
nonadjacent=1; % 1=calculate non-adjacent line coupling, 0=adjacent only
fdependent=0; %1=use frequency dependent calculations for even and odd-mode
% permittivity. This corresponds to the Kirschning and Jansen code. I would
% strongly recommend not using because it is prone to fatal NaN errors.
[C,L]=cMatrix(w,h,s,epsr,f0,nonadjacent,fdependent); %calculates all capacitances
% [C,L]=adjacentCMatrix(w,h,s,epsr,f0); %only adjacent line coupling
z=1i*2*pi*f0*L; %series impedance (pul)
y=1i*2*pi*f0*C; %shunt admittance (pul)
[V_hat,eigval,I_hat]=eig(z*y); 
% z*y is known as the system matrix in eigenvalue problems, V_hat and I_hat
% are modal voltages and currents
gam=sqrt(eigval);
Yc=z^(-1)*V_hat*gam./V_hat; 
Zc=y^(-1)*I_hat*gam./I_hat; 
% Characteristic admittance and impedance, Zc=Yc.^-1=V0+/I0+ on TL
phase1=exp(kron(gam.*0, [-1 1])).*kron(eye(n),[1 1]);
phase2=exp(kron(gam.*l, [-1 1])).*kron(eye(n),[1 1]);
% Network voltage matrix, similar to the form given by Vijai Tripathi
V=[V_hat*phase1;flipud(V_hat)*phase2];
% The network currents need an alternating negative sign, so I redefine phase to
% add that
phase1=exp(kron(gam.*0, [-1 1])).*kron(eye(n),[1 -1]); 
phase2=exp(kron(gam.*l, [-1 1])).*kron(eye(n),[1 -1]);
I=[V_hat.*Yc*phase1;-1*flipud(V_hat.*Yc)*phase2];
Z=V*I^(-1);

%keyboard %alex

% Check out the following if you're confused by "kron". It is the tensor
% product of two matrices, and is used here to very concisely double each
% column
% n=2;
% a=[1 0; 0 2];%eye(n);
% b=[-1,1];
% phase=kron(a,b)
% V_hat=[1 2;3 4];
% V_hat*c
end