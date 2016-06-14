%% Matrices for Electric Potential in Solid Phase, phi_s(x,t)
%   Created May 6, 2016 by Scott Moura
%   An update of phi_s_mats with 2nd order accurate boundary conditions

function [F1n,F1p,F2n,F2p,Gn,Gp,varargout] = phi_s_mats(p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update BRUGGEMAN RELATION % Apr.26 2016 by Saehong Park

sigma_n_eff = p.sig_n * (p.epsilon_s_n + p.epsilon_f_n) ^ p.brug; % Match DUALFOIL (Sig1)
sigma_p_eff = p.sig_p * (p.epsilon_s_p + p.epsilon_f_p) ^ p.brug; % Match DUALFOIL (Sig3)

alpha_n = 1 / (2 * p.L_n * p.delta_x_n);
alpha_p = 1 / (2 * p.L_p * p.delta_x_p);

%% Block Matrices
% M1 : phi_s x
M1n = alpha_n * (diag(ones(p.Nxn-2,1),1) + diag(-1*ones(p.Nxn-2,1),-1));
M1p = alpha_p * (diag(ones(p.Nxp-2,1),1) + diag(-1*ones(p.Nxp-2,1),-1));

% M2 : phi_s z
M2n = zeros(p.Nxn-1,2);
M2n(1,1) = -alpha_n;
M2n(end,end) = alpha_n;

M2p = zeros(p.Nxp-1,2);
M2p(1,1) = -alpha_p;
M2p(end,end) = alpha_p;

% M3 : i_e
M3n = [zeros(p.Nxn-1,1), diag(-1*ones(p.Nxn-1,1)), zeros(p.Nxn-1,1)] ./ sigma_n_eff;
M3p = [zeros(p.Nxp-1,1), diag(-1*ones(p.Nxp-1,1)), zeros(p.Nxp-1,1)] ./ sigma_p_eff;

% M4 : I
M4n = 1/sigma_n_eff * ones(p.Nxn-1,1);
M4p = 1/sigma_p_eff * ones(p.Nxp-1,1);

% N1 : phi_s x
N1n = zeros(2,p.Nxn-1);
N1n(1,1) = 4*alpha_n;
N1n(1,2) = -1*alpha_n;
N1n(2,end) = -4*alpha_n;
N1n(2,end-1) = 1*alpha_n;

N1p = zeros(2,p.Nxp-1);
N1p(1,1) = 4*alpha_p;
N1p(1,2) = -1*alpha_p;
N1p(2,end) = -4*alpha_p;
N1p(2,end-1) = 1*alpha_p;

% N2 : phi_s z
N2n = [-3*alpha_n, 0; 0, 3*alpha_n];
N2p = [-3*alpha_p, 0; 0, 3*alpha_p];

% N3 : I
N3n = [1/sigma_n_eff; 0];
N3p = [0; 1/sigma_p_eff];

%% Form F, G matrices
F1n = M1n - M2n*(N2n\N1n);
F2n = M3n;
Gn = M4n - M2n*(N2n\N3n);

F1p = M1p - M2p*(N2p\N1p);
F2p = M3p;
Gp = M4p - M2p*(N2p\N3p);

%% Compute C,D matrices for boundary values
Cn = -(N2n\N1n);
Dn = -(N2n\N3n);

Cp = -(N2p\N1p);
Dp = -(N2p\N3p);

varargout{1} = Cn;
varargout{2} = Cp;
varargout{3} = Dn;
varargout{4} = Dp;
