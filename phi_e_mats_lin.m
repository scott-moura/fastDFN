%% Matrices for Electric Potential in Electrolyte Phase, phi_e(x,t)
%   Created May 8, 2012 by Scott Moura

function [F1,F2,F3,varargout] = phi_e_mats_lin(p,c_ex)

% The electrolyte concentration vector is arranged as follows
% [c_0, c_1, ...,    c_sn, c_Nxn+2, ...,  c_sp, c_Nxn+Nxp      , c_pN]
% [------- c_e_n ------|------- c_e_s ------|------- c_e_p --------]
% length(c_e_n) = Nxn+1 11
% length(c_e_s) = Nxs+1 9
% length(c_e_p) = Nsp+1 11
% c_e_n(end) = c_e_s(1)
% c_e_s(end) = c_e_p(1)
% length(c_ex) = Nx+1 31

% Conductivity and FD Length Coefficients
alpha_n = 1 / (2 * p.L_n * p.delta_x_n);
alpha_s = 1 / (2 * p.L_s * p.delta_x_s);
alpha_p = 1 / (2 * p.L_p * p.delta_x_p);

% Electrolyte Conductivity
kappa = electrolyteCond(1e3*ones(size(c_ex)));
kappa0 = kappa(1);
kappa_n = kappa(2:p.Nxn);
kappa_ns = kappa(p.Nxn+1);
kappa_s = kappa(p.Nxn+2 : p.Nxn+2+p.Nxs-2);
kappa_sp = kappa(p.Nxn+2+p.Nxs-1);
kappa_p = kappa(p.Nxn+2+p.Nxs : end-1);
kappaN = kappa(end);

kappa_eff0 = kappa0; % .* p.epsilon_e_n.^(p.brug);
kappa_eff_n = kappa_n; % .* p.epsilon_e_n.^(p.brug);
kappa_eff_ns = kappa_ns; % .* p.epsilon_e_n.^(p.brug);
kappa_eff_s = kappa_s; % .* p.epsilon_e_s.^(p.brug);
kappa_eff_sp = kappa_sp; % .* p.epsilon_e_p.^(p.brug);
kappa_eff_p = kappa_p; % .* p.epsilon_e_p.^(p.brug);
kappa_effN = kappaN; % .* p.epsilon_e_p.^(p.brug);

% Diffusional Conductivity
beta_n = (p.R*p.T_amp)/(p.alph*p.Faraday) * (1 - p.t_plus) * (1 + 0) / (2 * p.L_n * p.delta_x_n);
beta_s = (p.R*p.T_amp)/(p.alph*p.Faraday) * (1 - p.t_plus) * (1 + 0) / (2 * p.L_s * p.delta_x_s);
beta_p = (p.R*p.T_amp)/(p.alph*p.Faraday) * (1 - p.t_plus) * (1 + 0) / (2 * p.L_p * p.delta_x_p);

gam = (p.R*p.T_amp)/(p.alph*p.Faraday) * (1 - p.t_plus) * (1 + 0);

%% Block Matrices
% M1 : phi_e x
M1n = alpha_n * (diag(ones(p.Nxn-2,1),1) + diag(-ones(p.Nxn-2,1),-1));
M1s = alpha_s * (diag(ones(p.Nxs-2,1),1) + diag(-ones(p.Nxs-2,1),-1));
M1p = alpha_p * (diag(ones(p.Nxp-2,1),1) + diag(-ones(p.Nxp-2,1),-1));

rsM1 = [p.Nxn-1; p.Nxs-1; p.Nxp-1];
csM1 = rsM1';
M1 = blkdiagFast(rsM1,csM1, M1n, M1s, M1p);

% M2 : phi_e z
M2n = zeros(p.Nxn-1,2);
M2n(1,1) = -alpha_n;
M2n(end,end) = alpha_n;

M2s = zeros(p.Nxs-1,2);
M2s(1,1) = -alpha_s;
M2s(end,end) = alpha_s;

M2p = zeros(p.Nxp-1,2);
M2p(1,1) = -alpha_p;
M2p(end,end) = alpha_p;

M2 = [M2n, zeros(p.Nxn-1,2); ...
      zeros(p.Nxs-1,1), M2s, zeros(p.Nxs-1,1);...
      zeros(p.Nxp-1,2), M2p];

% M3 : i_e
M3n = diag(1./kappa_eff_n);
M3s = diag(1./kappa_eff_s);
M3p = diag(1./kappa_eff_p);
rs = [p.Nxn-1; p.Nxs-1; p.Nxp-1];
cs = [p.Nxn, p.Nxs+1, p.Nxp];
M3 = blkdiagFast(rs,cs,[zeros(p.Nxn-1,1) M3n],...
                    [zeros(p.Nxs-1,1) M3s zeros(p.Nxs-1,1)],...
                    [M3p, zeros(p.Nxp-1,1)]);
                
% M4 : ln(c_ex)
% M4 = zeros(p.Nx-3,p.Nx+1);
% 
% for idx = 1:(p.Nx-3)
%     
%     if(idx <= (p.Nxn-1))
%         M4(idx,idx) = -beta_n;
%         M4(idx,idx+2) = beta_n;
%     elseif(idx <= (p.Nxn+p.Nxs-2))
%         M4(idx,idx+1) = -beta_s;
%         M4(idx,idx+3) = beta_s;
%     else
%         M4(idx,idx+2) = -beta_p;
%         M4(idx,idx+4) = beta_p;
%     end
%     
% end

M4n_1 = beta_n * (diag(-ones(p.Nxn-2,1),-1) + diag(ones(p.Nxn-2,1),1));
M4s_1 = beta_s * (diag(-ones(p.Nxs-2,1),-1) + diag(ones(p.Nxs-2,1),1));
M4p_1 = beta_p * (diag(-ones(p.Nxp-2,1),-1) + diag(ones(p.Nxp-2,1),1));

M4n_2 = [zeros(p.Nxn-1,1), M4n_1, zeros(p.Nxn-1,1)];
M4n_2(1,1) = -beta_n;
M4n_2(end,end) = beta_n;

M4p_2 = [zeros(p.Nxp-1,1), M4p_1, zeros(p.Nxp-1,1)];
M4p_2(1,1) = -beta_p;
M4p_2(end,end) = beta_p;

rsM4 = [p.Nxn-1; p.Nxs-1; p.Nxp-1];
csM4 = [p.Nxn+1, p.Nxs-1, p.Nxp+1];
M4 = blkdiagFast(rsM4, csM4, M4n_2, M4s_1, M4p_2);
M4(p.Nxn,p.Nxn+1) = -beta_s;
M4(p.Nxn-1+p.Nxs-1,p.Nxn+1+p.Nxs) = beta_s;

%% Boundary Conditions
N1 = zeros(4,p.Nx-3);   % phi_e x
N2 = zeros(4);          % phi_e z
N3 = zeros(4,p.Nx+1);   % i_e
N4 = zeros(4,p.Nx+1);   % ln(c_ex)

% BC1
N1(1,1) = 2*alpha_n;    % phi_e_{n,1}
N2(1,1) = -2*alpha_n;   % phi_e_{n,0}
N3(1,1) = 1/kappa_eff0; % i_e_{n,0}
N4(1,1) = 2*beta_n(1);  % ln(c_e)_{n,0}
N4(1,2) = -2*beta_n(1); % ln(c_e)_{n,1}

% BC2
N1(2,p.Nxn-1) = -alpha_n; % phi_e_{n,N-1}
N1(2,p.Nxn) = alpha_s;    % phi_e_{s,1}

N2(2,2) = alpha_n - alpha_s; % phi_e_{ns}

N3(2,p.Nxn+1) = 1/kappa_eff_ns; % i_e_{ns}

N4(2,p.Nxn) = gam*alpha_n;                  % ln(c_e)_{n,N-1}
N4(2,p.Nxn+1) = -gam*(alpha_n - alpha_s);   % ln(c_e)_{ns}
N4(2,p.Nxn+2) = -gam*alpha_s;               % ln(c_e)_{s,1}

% BC3
N1(3,p.Nxn-1+p.Nxs-1) = -alpha_s;   % phi_e_{s,N-1}
N1(3,p.Nxn-1+p.Nxs) = alpha_p;      % phi_e_{p,1}

N2(3,3) = alpha_s - alpha_p;        % phi_e_{sp}

N3(3,p.Nxn+2+p.Nxs-1) = 1/kappa_eff_sp; % i_e_{sp}

N4(3,p.Nxn+2+p.Nxs-2) = gam*alpha_s; % ln(c_e)_{s,N-1}
N4(3,p.Nxn+2+p.Nxs-1) = -gam*(alpha_s - alpha_p); % ln(c_e)_{sp}
N4(3,p.Nxn+2+p.Nxs) = -gam*alpha_p; % ln(c_e)_{p,1}

% BC4
N2(4,4) = alpha_p; % phi_e_{p,N}

% Invert N2
N2inv = diag(1./diag(N2));

%% Form F matrices
F1f = M1 - M2*N2inv*N1;
F2f = M3 - M2*N2inv*N3;
F3f = M4 - M2*N2inv*N4;

F1 = sparse(F1f);
F2 = sparse(F2f);
F3 = sparse(F3f);

C1 = sparse(-N2inv*N1);
C2 = sparse(-N2inv*N3);
C3 = sparse(-N2inv*N4);
varargout{1} = C1;
varargout{2} = C2;
varargout{3} = C3;

