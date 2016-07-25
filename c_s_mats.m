%% Matrices for Li Diffusion in Solid Phase, c_s(r,t)
%   Created May 21, 2012 by Scott Moura
%   Modified Jan 10, 2014 by Hector Perez 

function [An,Bn,Ap,Bp,varargout] = c_s_mats(p)

% Load Pade Approximations
N = p.PadeOrder;
% load(['ws/pade_cs' num2str(N) '.mat']);

% Replace symbolic expressions with numerical values
b_vecn = symsubsnum(p.R_s_n, p.D_s_n); %New Jan 10, 2014
a_vecn = symsubsden(p.R_s_n, p.D_s_n); %New Jan 10, 2014

b_vecp = symsubsnum(p.R_s_p, p.D_s_p); %New Jan 10, 2014
a_vecp = symsubsden(p.R_s_p, p.D_s_p); %New Jan 10, 2014

bp_vecn = b_vecn / a_vecn(end);
ap_vecn = a_vecn / a_vecn(end);

bp_vecp = b_vecp / a_vecp(end);
ap_vecp = a_vecp / a_vecp(end);

% Create state-space representation
An = diag(ones(N-1,1),1);
Ap = diag(ones(N-1,1),1);
for idx = 1:N
    An(N,idx) = -ap_vecn(idx);
    Ap(N,idx) = -ap_vecp(idx);
end

Bn = zeros(N,1);
Bp = zeros(N,1);
Bn(N) = 1;
Bp(N) = 1;

C1n = bp_vecn';
C2n = b_vecn(1) * ap_vecn(2:end)';
Cn = [C1n; C2n];

C1p = bp_vecp';
C2p = b_vecp(1) * ap_vecp(2:end)';
Cp = [C1p; C2p];

%% Convert to Jordan-Form   NEW Apr 24, 2014
% Similarity Matrix, computed symbolically
Vn = [-(p.R_s_n^4*(9*sqrt(2429)-457))/(381150*p.D_s_n^2), (p.R_s_n^4*(9*sqrt(2429)+457))/(381150*p.D_s_n^2), 1;...
       p.R_s_n^2*(sqrt(2429)-63)/(2310*p.D_s_n), -p.R_s_n^2*(sqrt(2429)+63)/(2310*p.D_s_n), 0;...
       1, 1, 0];
Vp = [-(p.R_s_p^4*(9*sqrt(2429)-457))/(381150*p.D_s_p^2), (p.R_s_p^4*(9*sqrt(2429)+457))/(381150*p.D_s_p^2), 1;...
       p.R_s_p^2*(sqrt(2429)-63)/(2310*p.D_s_p), -p.R_s_p^2*(sqrt(2429)+63)/(2310*p.D_s_p), 0;...
       1, 1, 0];

An1 = Vn\An*Vn;
Ap1 = Vp\Ap*Vp;

Bn1 = Vn\Bn;
Bp1 = Vp\Bp;

Cn1 = Cn*Vn;
Cp1 = Cp*Vp;

% Perform additional transformation that scales third state, such that it's
% exactly \bar{c}_s^\pm
Vn2 = diag([1 1 1/Cn1(2,3)]);
Vp2 = diag([1 1 1/Cp1(2,3)]);

An2 = Vn2\An1*Vn2;
Ap2 = Vp2\Ap1*Vp2;

Bn2 = Vn2\Bn1;
Bp2 = Vp2\Bp1;

Cn2 = Cn1*Vn2;
Cp2 = Cp1*Vp2;

An = An2;
An(3,:)=[0 0 0]; %reduce numerical error
An(1,2)=0; %reduce numerical error
An(2,1)=0; %reduce numerical error
An_normalized_D=An/p.D_s_n; %added by Federico for sensitivity analysis (equivalent to having D_s_n=1)
An_normalized_R=An*p.R_s_n^2; %added by Scott for sensitivity analysis (equivalent to having R_s_n=1)
Ap = Ap2;
Ap(3,:)=[0 0 0]; %reduce numerical error
Ap(1,2)=0; %reduce numerical error
Ap(2,1)=0; %reduce numerical error
Ap_normalized_D=Ap/p.D_s_p; %added by Federico for sensitivity analysis (equivalent to having D_s_p=1)
Ap_normalized_R=Ap*p.R_s_p^2; %added by Scott for sensitivity analysis (equivalent to having R_s_n=1)
Bn = Bn2;
Bp = Bp2;
Cn = Cn2;
Cn(2,1)=0; %reduce numerical error
Cn(2,2)=0; %reduce numerical error
Cp = Cp2;
Cp(2,1)=0; %reduce numerical error
Cp(2,2)=0; %reduce numerical error
   
%% Set Varargout
varargout{1} = Cn;
varargout{2} = Cp;
varargout{3} = An_normalized_D;
varargout{4} = Ap_normalized_D;
varargout{5} = An_normalized_R;
varargout{6} = Ap_normalized_R;

% Crank-Nicolson dicretization
% F1n = eye(N) - p.delta_t/2 * An;
% F2n = -p.delta_t/2 * Bn;
% 
% F1p = eye(N) - p.delta_t/2 * Ap;
% F2p = -p.delta_t/2 * Bp;
% 
% A1n = eye(N) + p.delta_t/2 * An;
% A2n = p.delta_t/2 * Bn;
% 
% A1p = eye(N) + p.delta_t/2 * Ap;
% A2p = p.delta_t/2 * Bp;
