%% Jacobian w.r.t. Parameters for Doyle-Fuller-Newman Model
%   Created Feb 18, 2014 by Scott Moura
%   Rebooted Jul 23, 2016 by Scott Moura for Samsung GRO project
%
%   INPUTS
%   x   : States      c_s_n, c_s_p, c_e, T
%   z   : Alg. vars   phi_s_n, phi_s_p, i_en, i_ep, phi_e, jn, jp
%   Cur : Applied Current
%   p   : Model parameter structure
%
%   OUTPUTS
%   JF = \partial F / \partial \theta    Jacobian of x ODEs w.r.t. theta
%   JFb = \partial Fb / \partial \theta    Jacobian of x BCs w.r.t. theta
%   JG = \partial G / \partial \theta    Jacobian of z eqns w.r.t. theta
%   JGb = \partial Gb / \partial \theta    Jacobian of z BCs w.r.t. theta
%
%   UNCERTAIN PARAMETERS, theta
%   1  : D_s_n
%   2  : D_s_p
%   3  : R_s_n
%   4  : R_s_p
%   5  : epsilon_s_n
%   6  : epsilon_s_p
%   7  : 1/sig_n
%   8  : 1/sig_p
%   9  : D_e
%   10 : epsilon_e_n
%   11 : epsilon_e_s
%   12 : epsilon_e_p
%   13 : kappa
%   14 : (1+t_plus)(1 + d ln f_ca / d ln c_e)
%   15 : k_n
%   16 : k_p
%   17 : R_f_n
%   18 : R_f_p
%   19 : L_n
%   20 : L_s
%   21 : L_p

function [JF, JFb, JG, JGb] = jac_p_dfn(x,z,Cur,p)


%% Parse out states

% Vector Lengths
Ncsn = p.PadeOrder * (p.Nxn-1);
Ncsp = p.PadeOrder * (p.Nxp-1);
Nce = p.Nx - 3;
Nc = Ncsn+Ncsp+Nce;
Nn = p.Nxn - 1;
Ns = p.Nxs - 1;
Np = p.Nxp - 1;
Nnp = Nn+Np;
Nx = p.Nx - 3;
Nz = 3*Nnp + Nx + 2;
Nt = 21;    % Number of uncertain params

% State Indices
ind_csn = 1:Ncsn;
ind_csp = Ncsn+1:Ncsn+Ncsp;

ind_cs = 1:Ncsn+Ncsp;
ind_ce = Ncsn+Ncsp+1:Nc;

ind_phi_s_n = 1:Nn;
ind_phi_s_p = Nn+1:Nnp;

ind_ien = Nnp+1:Nnp+Nn;
ind_iep = Nnp+Nn+1:2*Nnp;

ind_phi_e = 2*Nnp+1 : 2*Nnp+Nx+2;

ind_jn = 2*Nnp+Nx+3 : 2*Nnp+Nx+2+Nn;
ind_jp = 2*Nnp+Nx+2+Nn+1 : Nz;

% Parameter indices
ind_Dsn = 1;        %   1  : D_s_n
ind_Dsp = 2;        %   2  : D_s_p
ind_Rsn = 3;        %   3  : R_s_n
ind_Rsp = 4;        %   4  : R_s_p
ind_epsilonsn = 5;  %   5  : epsilon_s_n
ind_epsilonsp = 6;  %   6  : epsilon_s_p
ind_sign = 7;       %   7  : 1/sig_n
ind_sigp = 8;       %   8  : 1/sig_p
ind_De = 9;         %   9  : D_e
ind_epsilonen = 10; %   10 : epsilon_e_n
ind_epsilones = 11; %   11 : epsilon_e_s
ind_epsilonep = 12; %   12 : epsilon_e_p
ind_kappa = 13;     %   13 : kappa
ind_tplusfca = 14;  %   14 : (1+t_plus)(1 + d ln f_ca / d ln c_e)
ind_kn = 15;        %   15 : k_n
ind_kp = 16;        %   16 : k_p
ind_Rfn = 17;       %   17 : R_f_n
ind_Rfp = 18;       %   18 : R_f_p
ind_Ln = 19;        %   19 : L_n
ind_Ls = 20;        %   20 : L_s
ind_Lp = 21;        %   21 : L_p


% PARSE OUT THE STATES
% Solid Concentration
c_s_n = x(1:Ncsn);
c_s_p = x(Ncsn+1:Ncsn+Ncsp);

% Reformat into matrices
c_s_n_mat = reshape(c_s_n,p.PadeOrder,p.Nxn-1);
c_s_p_mat = reshape(c_s_p,p.PadeOrder,p.Nxp-1);

y_csn = p.C_csn * c_s_n_mat;
c_ss_n = y_csn(1,:)';
c_avg_n = y_csn(2,:)';

y_csp = p.C_csp * c_s_p_mat;
c_ss_p = y_csp(1,:)';
c_avg_p = y_csp(2,:)';

% Electrolyte concentration
c_e = x(ind_ce);
c_en = c_e(1:Nn);
c_es = c_e(Nn+1:Nn+p.Nxs-1);
c_ep = c_e(Nn+p.Nxs : end);
c_e_bcs = p.ce.C*c_e;
c_ex = [c_e_bcs(1); c_en; c_e_bcs(2); c_es; c_e_bcs(3); c_ep; c_e_bcs(4)];

% Temperature
T = x(end);

% Solid Potential
phi_s_n = z(ind_phi_s_n);
phi_s_p = z(ind_phi_s_p);

% Terminal Voltage
phi_s_n_bcs = p.C_psn * phi_s_n + p.D_psn * Cur;
phi_s_p_bcs = p.C_psp * phi_s_p + p.D_psp * Cur;
Volt = phi_s_p_bcs(2) - phi_s_n_bcs(1) - p.R_c*Cur;

% Electrolyte current
i_en = z(ind_ien);
i_ep = z(ind_iep);

% Electrolyte current across all three regions
i_e_in = [i_en; Cur*ones(p.Nxs+1,1); i_ep];

% Electrolyte potential
phi_e = z(ind_phi_e);

% Molar ionic flux
jn = z(ind_jn);
jp = z(ind_jp);

%% Preallocate Jacobian
JF = zeros(Nc+1,Nt);
JFb = zeros(10,Nt);
JG = zeros(Nz,Nt);
JGb = zeros(8,Nt);

%% [IN PROGRESS] Terms of JF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [DONE] Li Diffusion in Solid Phase: c_s(x,r,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% [DONE] Jacobian of c_s w.r.t. D_s_n, D_s_p
%%% [IN-PROGRESS] Jacobian of c_s w.r.t. R_s_n, R_s_p ???

%%%%%%%%%%%%%%%%%% Code written by Federico %%%%%%%%%%%%%%%%
% Loop through each "comb tooth" in anode
JF_csn_D = zeros(3,Nn);
JF_csn_R = zeros(3,Nn);
for idx = 1:Nn
    JF_csn_D(:,idx) = p.A_csn_normalized_D*c_s_n_mat(:,idx);
    JF_csn_R(:,idx) = p.A_csn_normalized_R*c_s_n_mat(:,idx);
end
JF(ind_csn,ind_Dsn) = reshape(JF_csn_D,[numel(JF_csn_D),1]);
JF(ind_csn,ind_Rsn) = reshape(JF_csn_R,[numel(JF_csn_R),1]);

% Loop through each "comb tooth" in cathode
JF_csp_D = zeros(3,Np);
JF_csp_R = zeros(3,Np);
for idx = 1:Np
    JF_csp_D(:,idx) = p.A_csp_normalized_D*c_s_p_mat(:,idx);
    JF_csp_R(:,idx) = p.A_csp_normalized_R*c_s_p_mat(:,idx);
end
JF(ind_csp,ind_Dsp) = reshape(JF_csp_D,[numel(JF_csp_D),1]);
JF(ind_csp,ind_Rsp) = reshape(JF_csp_R,[numel(JF_csp_R),1]);






%% [IN PROGRESS] Li Diffusion in Electrolyte Phase: c_e(x,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% BRUGGEMAN RELATION
D_en_eff = D_en .* p.epsilon_e_n.^(p.brug-1);
dD_en_eff = dD_en .* p.epsilon_e_n.^(p.brug-1);

D_es_eff = D_es .* p.epsilon_e_s.^(p.brug-1);
dD_es_eff = dD_es .* p.epsilon_e_s.^(p.brug-1);

D_ep_eff = D_ep .* p.epsilon_e_p.^(p.brug-1);
dD_ep_eff = dD_ep .* p.epsilon_e_p.^(p.brug-1);

%%% [DONE] Jacobian of c_e w.r.t. D_e

JF(ind_ce(1:Nn),ind_De) = ...
    dD_en_eff.*(p.ce.M1n*c_en + p.ce.M2n*c_e_bcs(1:2)).^2 ...
    + D_en_eff.*(p.ce.M3n*c_en + p.ce.M4n*c_e_bcs(1:2));

JF(ind_ce(Nn+1:Nn+p.Nxs-1),ind_De) = ...
    dD_es_eff.*(p.ce.M1s*c_es + p.ce.M2s*c_e_bcs(2:3)).^2 ...
    + D_es_eff.*(p.ce.M3s*c_es + p.ce.M4s*c_e_bcs(2:3));

JF(ind_ce(Nn+p.Nxs : end),ind_De) = ...
    dD_ep_eff.*(p.ce.M1p*c_ep + p.ce.M2p*c_e_bcs(3:4)).^2 ...
    + D_ep_eff.*(p.ce.M3p*c_ep + p.ce.M4p*c_e_bcs(3:4));



%%% [DONE] Jacobian of c_e w.r.t. epsilon_e_n, epsilon_e_s, epsilon_e_p

% BRUGGEMAN RELATION
D_en_eff_depse = D_en .* p.epsilon_e_n.^(p.brug-2) * (p.brug-1);
dD_en_eff_depse = dD_en .* p.epsilon_e_n.^(p.brug-2) * (p.brug-1);

D_es_eff_depse = D_es .* p.epsilon_e_s.^(p.brug-2) * (p.brug-1);
dD_es_eff_depse = dD_es .* p.epsilon_e_s.^(p.brug-2) * (p.brug-1);

D_ep_eff_depse = D_ep .* p.epsilon_e_p.^(p.brug-2) * (p.brug-1);
dD_ep_eff_depse = dD_ep .* p.epsilon_e_p.^(p.brug-2) * (p.brug-1);


JF(ind_ce(1:Nn),ind_epsilonen) = ...
    dD_en_eff_depse.*(p.ce.M1n*c_en + p.ce.M2n*c_e_bcs(1:2)).^2 ...
    + D_en_eff_depse.*(p.ce.M3n*c_en + p.ce.M4n*c_e_bcs(1:2)) + p.ce.M5n*jn*(-1/p.epsilon_e_n);

JF(ind_ce(Nn+1:Nn+p.Nxs-1),ind_epsilones) = ...
    dD_es_eff_depse.*(p.ce.M1s*c_es + p.ce.M2s*c_e_bcs(2:3)).^2 ...
    + D_es_eff_depse.*(p.ce.M3s*c_es + p.ce.M4s*c_e_bcs(2:3));

JF(ind_ce(Nn+p.Nxs : end),ind_epsilonep) = ... 
    dD_ep_eff_depse.*(p.ce.M1p*c_ep + p.ce.M2p*c_e_bcs(3:4)).^2 ...
    + D_ep_eff_depse.*(p.ce.M3p*c_ep + p.ce.M4p*c_e_bcs(3:4)) + p.ce.M5p*jp*(-1/p.epsilon_e_p);


%%% [DONE] Jacobian of c_e w.r.t. (1-t_plus)*(1+0)
JF(ind_ce(1:Nn),ind_tplusfca) = p.ce.M5n*jn / (1-p.t_plus);
JF(ind_ce(Nn+p.Nxs : end),ind_tplusfca) = p.ce.M5p*jp / (1-p.t_plus);


%%% [DONE] Jacobian of c_e w.r.t. epsilon_s_n, epsilon_s_p (via a_s)
JF(ind_ce(1:Nn),ind_epsilonsn) = p.ce.M5n*jn / (p.epsilon_s_n);
JF(ind_ce(Nn+p.Nxs : end),ind_epsilonsp) = p.ce.M5p*jp / (p.epsilon_s_p);


%%% [DONE] Jacobian of c_e w.r.t. R_s_n, R_s_p (via a_s)
JF(ind_ce(1:Nn),ind_Rsn) = p.ce.M5n*jn * (-1/p.R_s_n);
JF(ind_ce(Nn+p.Nxs : end),ind_Rsp) = p.ce.M5p*jp * (-1/p.R_s_p);


%%% [IN PROGRESS] Jacobian of c_e w.r.t. L_n, L_s, L_p
% ignore for now


%% [IN PROGRESS] Temperature: T(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ignore for now

%% Terms of JFb
% For now, disregard parameter variations in the BCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [IN PROGRESS] Terms of JG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [IN PROGRESS] Potential in Solid Phase: phi_s(x,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_enn = [0; i_en; Cur];
i_epp = [Cur; i_ep; 0];

%%% [DONE] Jacobian of phi_s w.r.t. 1/sig_n, 1/sig_p
JG(ind_phi_s_n, ind_sign) = (p.F2_psn*i_enn + p.G_psn*Cur) .* p.sig_n;
JG(ind_phi_s_p, ind_sigp) = (p.F2_psp*i_epp + p.G_psp*Cur) .* p.sig_p;


%%% [DONE] Jacobian of phi_s w.r.t. epsilon_s_n, epsilon_s_p
JG(ind_phi_s_n, ind_epsilonsn) = (p.F2_psn*i_enn + p.G_psn*Cur) .* (-p.brug/(p.epsilon_s_n + p.epsilon_f_n));
JG(ind_phi_s_p, ind_epsilonsp) = (p.F2_psp*i_epp + p.G_psp*Cur) .* (-p.brug/(p.epsilon_s_p + p.epsilon_f_p));


%%% [IN PROGRESS] Jacobian of phi_s w.r.t. L_n, L_s, L_p



%% [IN PROGRESS] Electrolyte Current: i_e(x,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% [DONE] Jacobian of i_e w.r.t. epsilon_s_n, epsilon_s_p (via a_s)
JG(ind_ien, ind_epsilonsn) = p.F1_ien*i_en * (-1/p.epsilon_s_n);
JG(ind_iep, ind_epsilonsp) = p.F1_iep*i_ep * (-1/p.epsilon_s_p);


%%% [DONE] Jacobian of i_e w.r.t. R_s_n, R_s_p (via a_s)
JG(ind_ien, ind_Rsn) = p.F1_ien*i_en * (1/p.R_s_n);
JG(ind_iep, ind_Rsp) = p.F1_iep*i_ep * (1/p.R_s_p);


%%% [IN PROGRESS] Jacobian of i_e w.r.t. L_n, L_s, L_p



%% [IN PROGRESS] Potential in Electrolyte Phase: phi_e(x,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%% [DONE] Jacobian of phi_e w.r.t. kappa (actually, unity coeff)
JG(ind_phi_e, ind_kappa) = (Kap_eff*p.M1_pe)*phi_e + (Kap_eff_D*p.M4_pe)*log(c_ex);


%%% [DONE] Jacobian of phi_e w.r.t. epsilon_e_n, epsilon_e_s, epsilon_e_p

% Effective conductivity - multiply by p.epsilon_e_x ^ (p.brug) % Apr.22 by Saehong Park
kappa_eff0_epsen = kappa0 .* p.brug * p.epsilon_e_n.^(p.brug-1);
kappa_eff_n_epsen = kappa_n .* p.brug * p.epsilon_e_n.^(p.brug-1);

kappa_eff_ns_epsen = kappa_ns .* p.brug/2 * ((p.epsilon_e_n + p.epsilon_e_s)/2).^(p.brug-1);
kappa_eff_ns_epses = kappa_ns .* p.brug/2 * ((p.epsilon_e_n + p.epsilon_e_s)/2).^(p.brug-1);

kappa_eff_s_epses = kappa_s .* p.brug * p.epsilon_e_s.^(p.brug-1);

kappa_eff_sp_epses = kappa_sp .* p.brug/2 * ((p.epsilon_e_s + p.epsilon_e_p)/2).^(p.brug-1);
kappa_eff_sp_epsep = kappa_sp .* p.brug/2 * ((p.epsilon_e_s + p.epsilon_e_p)/2).^(p.brug-1);

kappa_eff_p_epsep = kappa_p .* p.brug * p.epsilon_e_p.^(p.brug-1);
kappa_effN_epsep = kappaN .* p.brug * p.epsilon_e_p.^(p.brug-1);

% Form into vector
kappa_eff_epsen = [kappa_eff_n_epsen; kappa_eff_ns_epsen; 0*kappa_eff_s; 0*kappa_eff_sp; 0*kappa_eff_p];
kappa_eff_epses = [0*kappa_eff_n; kappa_eff_ns_epses; kappa_eff_s_epses; kappa_eff_sp_epses; 0*kappa_eff_p];
kappa_eff_epsep = [0*kappa_eff_n; 0*kappa_eff_ns; 0*kappa_eff_s; kappa_eff_sp_epsep; kappa_eff_p_epsep];

Kap_eff_epsen = sparse(diag(kappa_eff_epsen));
Kap_eff_epses = sparse(diag(kappa_eff_epses));
Kap_eff_epsep = sparse(diag(kappa_eff_epsep));

% Modified effective conductivity
Kap_eff_D_epsen = bet*Kap_eff_epsen;
Kap_eff_D_epses = bet*Kap_eff_epses;
Kap_eff_D_epsep = bet*Kap_eff_epsep;


F1_pe_epsen = Kap_eff_epsen*p.M1_pe;
F3_pe_epsen = Kap_eff_D_epsen*p.M4_pe;

F1_pe_epses = Kap_eff_epses*p.M1_pe;
F3_pe_epses = Kap_eff_D_epses*p.M4_pe;

F1_pe_epsep = Kap_eff_epsep*p.M1_pe;
F3_pe_epsep = Kap_eff_D_epsep*p.M4_pe;


JG(ind_phi_e, ind_epsilonen) = F1_pe_epsen*phi_e + F3_pe_epsen*log(c_ex);
JG(ind_phi_e, ind_epsilones) = F1_pe_epses*phi_e + F3_pe_epses*log(c_ex);
JG(ind_phi_e, ind_epsilonep) = F1_pe_epsep*phi_e + F3_pe_epsep*log(c_ex);


%%% [DONE] Jacobian of phi_e w.r.t. (1-p.t_plus)*(1+p.activity)
JG(ind_phi_e, ind_tplusfca) = -F3_pe*log(c_ex) ./ ((p.t_plus - 1) * (1 + p.dactivity));


%%% [IN PROGRESS] Jacobian of phi_e w.r.t. L_n, L_s, L_p



%% [DONE] Butler-Volmer Equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% [DONE] Jacobian of jn w.r.t. k_n
JG(ind_jn,ind_kn) = 2/p.Faraday .* sinh(aFRT * eta_n) .* ((p.c_s_n_max - c_ss_n) .* c_ss_n .* c_en).^p.alph;

% [DONE] Jacobian of jp w.r.t. k_p
JG(ind_jp,ind_kp) = 2/p.Faraday .* sinh(aFRT * eta_p) .* ((p.c_s_p_max - c_ss_p) .* c_ss_p .* c_ep).^p.alph;

% [DONE] Jacobian of jn w.r.t. R_f_n
JG(ind_jn,ind_Rfn) = 2/p.Faraday * i_0n .* cosh(aFRT * eta_n) .* aFRT .* (-p.Faraday*jn);

% [DONE] Jacobian of jn w.r.t. R_f_n
JG(ind_jp,ind_Rfp) = 2/p.Faraday * i_0p .* cosh(aFRT * eta_p) .* aFRT .* (-p.Faraday*jp);



%% [IN PROGRESS] Terms of JGb
%  All terms are zero.

%% Sparsify Jacobian
JF = sparse(JF); %HEP
%JF = sparse(JF * diag(p.theta0)); %not required
JFb = sparse(JFb); 
JG = sparse(JG); %HEP
%JG = sparse(JG * diag(p.theta0)); %not required
JGb = sparse(JGb);

