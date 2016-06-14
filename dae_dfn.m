%% DAEs for Doyle-Fuller-Newman Model
%   Created May 22, 2012 by Scott Moura

% Modified June 19, 2015 by Scott Moura
% Electrolyte concentration dynamics re-written with
%   - 2nd order boundary conditions
%   - rigorous implementation of nonlinear dynamics with Finite Difference

% Modified Apr 26, 2016 by Saehong Park
% Bruggeman relationship
%   - 'Effective' Electrolyte diffusion coefficients in here
%   - 'Effective' Electrolyte conductivity (Kappa) in phi_e_mats
%   - 'Effective' Solid-phase conductivity (Sigma) in phi_s_mats
%   - Modified electrolyte governing equation with volume fraction
%   - update phi_e_mats.m / phi_s_mats.m
%   - Added p.epsilon_f_n, p.epsilon_f_p in params_dualfoil


function [f, g, varargout] = dae_dfn(x,z,Cur,p)
%% Parse out states
Ncsn = p.PadeOrder * (p.Nxn-1);
Ncsp = p.PadeOrder * (p.Nxp-1);
Nce = p.Nx - 3;
Nc = Ncsn+Ncsp+Nce;
Nn = p.Nxn - 1;
Np = p.Nxp - 1;
Nnp = Nn+Np;
Nx = p.Nx - 3;

% Solid Concentration
c_s_n = x(1:Ncsn);
c_s_p = x(Ncsn+1:Ncsn+Ncsp);

% Reformat into matrices
c_s_n_mat = reshape(c_s_n,p.PadeOrder,p.Nxn-1);
c_s_p_mat = reshape(c_s_p,p.PadeOrder,p.Nxp-1);

% Electrolyte concentration
c_e = x((Ncsn+Ncsp + 1):(Nc));

% Temperature
T = x(end);

% Solid Potential
phi_s_n = z(1:Nn);
phi_s_p = z(Nn+1:Nnp);

% Electrolyte Current
i_en = z(Nnp+1 : Nnp+Nn);
i_ep = z(Nnp+Nn+1 : 2*Nnp);

% Electrolyte Potential
phi_e = z(2*Nnp+1:2*Nnp+Nx+2);

% Molar ionic flux
jn = z(2*Nnp+Nx+3 : 2*Nnp+Nx+2+Nn);
jp = z(2*Nnp+Nx+2+Nn+1 : end);

%% Li Diffusion in Solid Phase: c_s(x,r,t)

% Anode LTI System
c_s_n_dot_mat = p.A_csn * c_s_n_mat + p.B_csn * jn';
y_csn = p.C_csn * c_s_n_mat;

% Anode Parse LTI outputs
c_ss_n = y_csn(1,:)';
c_avg_n = y_csn(2,:)';
c_s_n_dot = reshape(c_s_n_dot_mat,numel(c_s_n_dot_mat),1);

% Cathode LTI System
c_s_p_dot_mat = p.A_csp * c_s_p_mat + p.B_csp * jp';
y_csp = p.C_csp * c_s_p_mat;

% Cathode Parse LTI outputs
c_ss_p = y_csp(1,:)';
c_avg_p = y_csp(2,:)';
c_s_p_dot = reshape(c_s_p_dot_mat,numel(c_s_p_dot_mat),1);

%% Li Diffusion in Electrolyte Phase: c_e(x,t)
% Compute Boundary Conditions
c_e_bcs = p.ce.C*c_e;

% Separate and aggregate
c_en = c_e(1:Nn);
c_es = c_e(Nn+1:Nn+p.Nxs-1);
c_ep = c_e(Nn+p.Nxs : end);
c_ex = [c_e_bcs(1); c_en; c_e_bcs(2); c_es; c_e_bcs(3); c_ep; c_e_bcs(4)];

% Compute Electrolyte Diffusion Coefficient and Derivative
[D_en0,dD_en0] = electrolyteDe(c_en);
[D_es0,dD_es0] = electrolyteDe(c_es);
[D_ep0,dD_ep0] = electrolyteDe(c_ep);

% Adjustment for Arrhenius
Arrh_De = exp(p.E.De/p.R*(1/p.T_ref - 1/T));
D_en = D_en0 * Arrh_De;
D_es = D_es0 * Arrh_De;
D_ep = D_ep0 * Arrh_De;
dD_en = dD_en0 * Arrh_De;
dD_es = dD_es0 * Arrh_De;
dD_ep = dD_ep0 * Arrh_De;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD BRUGGEMAN RELATION % Apr.22 2016 by Saehong Park
D_en_eff = D_en .* p.epsilon_e_n.^(p.brug-1);
dD_en_eff = dD_en .* p.epsilon_e_n.^(p.brug-1);

% DO same s,p 
D_es_eff = D_es .* p.epsilon_e_s.^(p.brug-1);
dD_es_eff = dD_es .* p.epsilon_e_s.^(p.brug-1);

D_ep_eff = D_ep .* p.epsilon_e_p.^(p.brug-1);
dD_ep_eff = dD_ep .* p.epsilon_e_p.^(p.brug-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO change D_en, D_es D_ep as D_eff... % Apr.22 2016 by Saehong Park

% System Matrices have all been precomputed & stored in param struct "p"

% Compute derivative
c_en_dot = dD_en_eff.*(p.ce.M1n*c_en + p.ce.M2n*c_e_bcs(1:2)).^2 ...
    + D_en_eff.*(p.ce.M3n*c_en + p.ce.M4n*c_e_bcs(1:2)) + p.ce.M5n*jn;

c_es_dot = dD_es_eff.*(p.ce.M1s*c_es + p.ce.M2s*c_e_bcs(2:3)).^2 ...
    + D_es_eff.*(p.ce.M3s*c_es + p.ce.M4s*c_e_bcs(2:3));

c_ep_dot = dD_ep_eff.*(p.ce.M1p*c_ep + p.ce.M2p*c_e_bcs(3:4)).^2 ...
    + D_ep_eff.*(p.ce.M3p*c_ep + p.ce.M4p*c_e_bcs(3:4)) + p.ce.M5p*jp;

% Assemble c_e_dot
c_e_dot = [c_en_dot; c_es_dot; c_ep_dot];

%% Potential in Solid Phase: phi_s(x,t)
% Algebraic eqns (semi-explicit form)
i_enn = [0; i_en; Cur];
i_epp = [Cur; i_ep; 0];

phi_sn_dot = p.F1_psn*phi_s_n + p.F2_psn*i_enn + p.G_psn*Cur;
phi_sp_dot = p.F1_psp*phi_s_p + p.F2_psp*i_epp + p.G_psp*Cur;

% Terminal Voltage
phi_s_n_bcs = p.C_psn * phi_s_n + p.D_psn * Cur;
phi_s_p_bcs = p.C_psp * phi_s_p + p.D_psp * Cur;

Volt = phi_s_p_bcs(2) - phi_s_n_bcs(1) - p.R_c*Cur;% IR drop added by Satadru, November 2015


%% Electrolyte Current: i_e(x,t)
i_en_dot = p.F1_ien*i_en + p.F2_ien*jn + p.F3_ien*Cur;
i_ep_dot = p.F1_iep*i_ep + p.F2_iep*jp + p.F3_iep*Cur;

% Electrolyte current across all three regions
i_e_in = [i_en; Cur*ones(p.Nxs+1,1); i_ep];

%% Potential in Electrolyte Phase: phi_e(x,t)

% Electrolyte Conductivity
kappa_ref = electrolyteCond(c_ex);

% Adjustment for Arrhenius
kappa = kappa_ref * exp(p.E.kappa_e/p.R*(1/p.T_ref - 1/T));

kappa0 = kappa(1);                              % BC1
kappa_n = kappa(2:p.Nxn);
kappa_ns = kappa(p.Nxn+1);
kappa_s = kappa(p.Nxn+2 : p.Nxn+2+p.Nxs-2);
kappa_sp = kappa(p.Nxn+2+p.Nxs-1);
kappa_p = kappa(p.Nxn+2+p.Nxs : end-1);
kappaN = kappa(end);                            % BC2

% Effective conductivity - multiply by p.epsilon_e_x ^ (p.brug) % Apr.22 by Saehong Park
kappa_eff0 = kappa0 .* p.epsilon_e_n.^(p.brug);
kappa_eff_n = kappa_n .* p.epsilon_e_n.^(p.brug);

kappa_eff_ns = kappa_ns .* ((p.epsilon_e_n + p.epsilon_e_s)/2).^(p.brug);

kappa_eff_s = kappa_s .* p.epsilon_e_s.^(p.brug);

kappa_eff_sp = kappa_sp .* ((p.epsilon_e_s + p.epsilon_e_p)/2).^(p.brug);

kappa_eff_p = kappa_p .* p.epsilon_e_p.^(p.brug);
kappa_effN = kappaN .* p.epsilon_e_p.^(p.brug);

% Form into vector
kappa_eff = [kappa_eff_n; kappa_eff_ns; kappa_eff_s; kappa_eff_sp; kappa_eff_p];
Kap_eff = sparse(diag(kappa_eff));

% Diffusional Conductivity
bet = (2*p.R*T)/(p.Faraday) * (p.t_plus - 1) * (1 + p.dactivity);

% Modified effective conductivity
Kap_eff_D = bet*Kap_eff;

% Form Matrices
M2_pe = p.M2_pe;
M2_pe(1,1) = M2_pe(1,1) * kappa_eff0;
M2_pe(end,end) = M2_pe(end,end) * kappa_effN;

F1_pe = Kap_eff*p.M1_pe + M2_pe*p.C_pe;
F2_pe = p.M3_pe;
F3_pe = Kap_eff_D*p.M4_pe;

% Algebraic eqns (semi-explicit form)
phi_e_dot = F1_pe*phi_e + F2_pe*i_e_in + F3_pe*log(c_ex);

%% Butler-Volmer Equation
aFRT = (p.alph*p.Faraday)/(p.R*T);

% Exchange Current Density, i_0^{\pm}
[i_0n,i_0p] = exch_cur_dens(p,c_ss_n,c_ss_p,c_e);

% Equilibrium Potential, U^{\pm}(c_ss)
theta_n = c_ss_n / p.c_s_n_max;
theta_p = c_ss_p / p.c_s_p_max;
Unref = refPotentialAnode(p, theta_n);
Upref = refPotentialCathode(p, theta_p);

% Overpotential, \eta
eta_n = phi_s_n - phi_e(1:Nn) - Unref - p.Faraday*p.R_f_n*jn;
eta_p = phi_s_p - phi_e(end-Np+1:end) - Upref - p.Faraday*p.R_f_p*jp;

% Algebraic eqns (semi-explicit form)
jn_dot = 2/p.Faraday * i_0n .* sinh(aFRT * eta_n) - jn;
jp_dot = 2/p.Faraday * i_0p .* sinh(aFRT * eta_p) - jp;

%% Temperature
% Equilibrium Potential and Gradient wrt bulk concentration
theta_avg_n = c_avg_n / p.c_s_n_max;
theta_avg_p = c_avg_p / p.c_s_p_max;
[Unb,~,dUnbdT] = refPotentialAnode(p, theta_avg_n);
[Upb,~,dUpbdT] = refPotentialCathode(p, theta_avg_p);

% Heat generated from intercalation (w/o boundaries for NOW)
Q_nx = p.a_s_n*p.Faraday * jn .* (Unb - T*dUnbdT);
Q_n = sum(Q_nx) * p.delta_x_n * p.L_n;

Q_px = p.a_s_p*p.Faraday * jp .* (Upb - T*dUpbdT);
Q_p = sum(Q_px) * p.delta_x_p * p.L_p;

Q_inter = Q_n + Q_p;

% Temperature ODE
T_dot = (p.h*(p.T_amb - T) - Cur*Volt - Q_inter) / (p.rho_avg*p.C_p);

%% Concatenate Time Derivatives
f = [c_s_n_dot; c_s_p_dot; c_e_dot; T_dot];
g = [phi_sn_dot; phi_sp_dot; i_en_dot; i_ep_dot; ...
     phi_e_dot; jn_dot; jp_dot];

%% Optional Output Vars

% Conservation of Li-ion matters
n_Li_s = sum(c_avg_n) * p.L_n*p.delta_x_n * p.epsilon_s_n * p.Area ...
     + sum(c_avg_p) * p.L_p*p.delta_x_p * p.epsilon_s_p * p.Area;
 
n_Li_e = sum(c_e(1:Nn)) * p.L_n*p.delta_x_n * p.epsilon_e_n * p.Area ...
     + sum(c_e(Nn+1:end-Np)) * p.L_s*p.delta_x_s * p.epsilon_e_s * p.Area ...
     + sum(c_e(end-Np+1:end)) * p.L_p*p.delta_x_p * p.epsilon_e_p * p.Area;

% nLidot = sum(jn) * p.delta_x_n * Nn + sum(jp) * p.delta_x_p * Np;

c_e0n = c_ex(1);
c_e0p = c_ex(end);

eta_s_Ln = phi_s_n_bcs(2) - phi_e(Nn+1);

% Aggregate Outputs
varargout{1} = [c_ss_n; c_ss_p; c_avg_n; c_avg_p; c_ex; eta_n; eta_p;...
    c_e0n; c_e0p; eta_s_Ln;...
    Volt; n_Li_s; n_Li_e];


