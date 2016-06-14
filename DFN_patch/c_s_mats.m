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

% b_vecn = subs(pade_cs.num(:),{'R';'D'},[p.R_s_n; p.D_s_n]); %Changed Jan 10, 2014, was pade_cs.num
% a_vecn = subs(pade_cs.den(:),{'R';'D'},[p.R_s_n; p.D_s_n]); %Changed Jan 10, 2014, was pade_cs.den

% b_vecp = subs(pade_cs.num(:),{'R';'D'},[p.R_s_p; p.D_s_p]); %Changed Jan 10, 2014, was pade_cs.num
% a_vecp = subs(pade_cs.den(:),{'R';'D'},[p.R_s_p; p.D_s_p]); %Changed Jan 10, 2014, was pade_cs.den

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

varargout{1} = Cn;
varargout{2} = Cp;

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
