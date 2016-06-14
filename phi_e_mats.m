%% Matrices for Electric Potential in Electrolyte Phase, phi_e(x,t)
%   Created May 8, 2012 by Scott Moura
%   Modified Apr 26, 2016 by Saehong Park
%   Modified May 7, 2016 by Scott Moura - FDM scheme completely updated
%   Modified May 30, 2016 by Scott Moura - Precalculate matrices for speed

%   Input: electrolyte concentration c_e (entire vector, w/ BCs)

function [M1,M2,M3,M4,C] = phi_e_mats(p)

% The electrolyte concentration vector is arranged as follows
% [c_0, c_1, ...,    c_ns, c_Nxn+2, ...,  c_sp, c_Nxn+Nxp      , c_pN]
% [------- c_e_n ------|------- c_e_s ------|------- c_e_p --------]
% length(c_e_n) = Nxn+1 11
% length(c_e_s) = Nxs+1 9
% length(c_e_p) = Nsp+1 11
% c_e_n(end) = c_e_s(1)
% c_e_s(end) = c_e_p(1)
% length(c_ex) = Nx+1 31

% Conductivity and FD Length Coefficients
alpha_n = 1 / (2 * p.L_n * p.delta_x_n);
alpha_ns = 1 / (p.L_n * p.delta_x_n + p.L_s * p.delta_x_s);
alpha_s = 1 / (2 * p.L_s * p.delta_x_s);
alpha_sp = 1 / (p.L_s * p.delta_x_s + p.L_p * p.delta_x_p);
alpha_p = 1 / (2 * p.L_p * p.delta_x_p);

% % Electrolyte Conductivity
% kappa = electrolyteCond(c_ex);
% 
% kappa0 = kappa(1);                              % BC1
% kappa_n = kappa(2:p.Nxn);
% kappa_ns = kappa(p.Nxn+1);
% kappa_s = kappa(p.Nxn+2 : p.Nxn+2+p.Nxs-2);
% kappa_sp = kappa(p.Nxn+2+p.Nxs-1);
% kappa_p = kappa(p.Nxn+2+p.Nxs : end-1);
% kappaN = kappa(end);                            % BC2
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Effective conductivity - multiply by p.epsilon_e_x ^ (p.brug) % Apr.22 by Saehong Park
% 
% kappa_eff0 = kappa0 .* p.epsilon_e_n.^(p.brug);
% kappa_eff_n = kappa_n .* p.epsilon_e_n.^(p.brug);
% 
% kappa_eff_ns = kappa_ns .* ((p.epsilon_e_n + p.epsilon_e_s)/2).^(p.brug);
% 
% kappa_eff_s = kappa_s .* p.epsilon_e_s.^(p.brug);
% 
% kappa_eff_sp = kappa_sp .* ((p.epsilon_e_s + p.epsilon_e_p)/2).^(p.brug);
% 
% kappa_eff_p = kappa_p .* p.epsilon_e_p.^(p.brug);
% kappa_effN = kappaN .* p.epsilon_e_p.^(p.brug);
% 
% Kap_eff_n = sparse(diag(kappa_eff_n));
% Kap_eff_s = sparse(diag(kappa_eff_s));
% Kap_eff_p = sparse(diag(kappa_eff_p));
% 
% % Diffusional Conductivity
% bet = (2*p.R*p.T_amp)/(p.Faraday) * (p.t_plus - 1) * (1 + p.dactivity);


%% Block Matrices
% M1 : phi_e x = (1/(2*Delta_x)) * kappa_eff_j(c_e) * [ones up diag, -ones low diag]

% eqns for anode, n
M1n_c = alpha_n * p.M1_pen_skel;
M1n_r = zeros(p.Nxn-1,1);
M1n_r(end) = alpha_n;

M1_n = [M1n_c, M1n_r, zeros(p.Nxn-1, (p.Nxs-1 + 1 + p.Nxp-1))];


% eqns for ns interface
M1ns_c = [-alpha_ns, 0, alpha_ns];
M1_ns = [zeros(1,p.Nxn-2), M1ns_c, zeros(1,p.Nxs-2 + 1 + p.Nxp-1)];


% eqns for separator, s
M1s_c = alpha_s * p.M1_pes_skel;
M1s_l = zeros(p.Nxs-1,1);
M1s_l(1) = -alpha_s;
M1s_r = zeros(p.Nxs-1,1);
M1s_r(end) = alpha_s;

M1_s = [zeros(p.Nxs-1,p.Nxn-1), M1s_l, M1s_c, M1s_r, zeros(p.Nxs-1,p.Nxp-1)];


% eqns for sp interface
M1sp_c = [-alpha_sp, 0, alpha_sp];
M1_sp = [zeros(1,p.Nxn-1 + 1 + p.Nxs-2), M1sp_c, zeros(1,p.Nxp-2)];


% eqns for cathode, p
M1p_c = alpha_p * p.M1_pep_skel;
M1p_l = zeros(p.Nxp-1,1);
M1p_l(1) = -alpha_p;

M1_p = [zeros(p.Nxp-1, (p.Nxn-1 + 1 + p.Nxs-1)), M1p_l, M1p_c];

% assemble submatrices
M1 = [M1_n; M1_ns; M1_s; M1_sp; M1_p];
M1 = sparse(M1);


% M2 : phi_e z
M2 = zeros(size(M1,1), 2);
M2(1,1) = -alpha_n;
M2(end,end) = alpha_p;
M2 = sparse(M2);


% M3 : i_e
M3 = speye(size(M1));


% M4 : ln(c_ex)
M4 = [M2(:,1), M1, M2(:,2)];
M4 = sparse(M4);

%% Boundary Conditions
N1 = zeros(2,size(M1,1));
N1(1,1) = 0;
N1(1,2) = 0;
N1(2,end-1) = 1;
N1(2,end) = -4;
N1 = sparse(N1);

N2inv = sparse(diag([1,1/3]));

%% Form F matrices
% F1 = M1 - M2*(N2inv*N1);
% F2 = M3;
% F3 = M4;

C = -(N2inv*N1);

