%% Matrices for Electrolyte Current, i_e(x,t)
%   Created May 8, 2012 by Scott Moura

function [F1n,F1p,F2n,F2p,F3n,F3p] = i_e_mats(p)

% Conductivity and FD Length Coefficients
alpha_n = 1 / (2 * p.L_n * p.delta_x_n * p.a_s_n);
alpha_p = 1 / (2 * p.L_p * p.delta_x_p * p.a_s_p);

beta_n = p.Faraday;
beta_p = p.Faraday;

%% Block Matrices
% M1 : i_e x
M1n = alpha_n * (diag(ones(p.Nxn-2,1),1) + diag(-ones(p.Nxn-2,1),-1));
M1p = alpha_p * (diag(ones(p.Nxp-2,1),1) + diag(-ones(p.Nxp-2,1),-1));

% M2 : i_e z
M2n = zeros(p.Nxn-1,2);
M2n(1,1) = -alpha_n;
M2n(end,end) = alpha_n;

M2p = zeros(p.Nxp-1,2);
M2p(1,1) = -alpha_p;
M2p(end,end) = alpha_p;

% M3 : J
M3n = -beta_n*diag(ones(p.Nxn-1,1));
M3p = -beta_p*diag(ones(p.Nxp-1,1));

% N1 : i_e x
N1n = zeros(2,p.Nxn-1);

N1p = zeros(2,p.Nxp-1);

% N2 : i_e z
N2n = eye(2);
N2p = eye(2);

% N3 : jn
N3n = zeros(2,p.Nxn-1);

N3p = zeros(2,p.Nxp-1);

% N4 : I
N4n = [0; 1];
N4p = [1; 0];

%% Form F matrices
F1n = M1n - M2n*(N2n\N1n);
F2n = M3n - M2n*(N2n\N3n);
F3n = M2n*(N2n\N4n);

F1p = M1p - M2p*(N2p\N1p);
F2p = M3p - M2p*(N2p\N3p);
F3p = M2p*(N2p\N4p);

