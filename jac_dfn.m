%% Jacobian for Doyle-Fuller-Newman Model
%   Created May 22, 2012 by Scott Moura
%   State-dependent elements

function [f_x_full, f_z_full, g_x_full, g_z_full, varargout] = jac_dfn(x,z,Cur,f_x,f_z,g_x,g_z,p)


%% Parse out states
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

% Solid Concentration
c_s_n = x(1:Ncsn);
c_s_p = x(Ncsn+1:Ncsn+Ncsp);

% Reformat into matrices
c_s_n_mat = reshape(c_s_n,p.PadeOrder,p.Nxn-1);
c_s_p_mat = reshape(c_s_p,p.PadeOrder,p.Nxp-1);

% Electrolyte concentration
c_e = x(ind_ce);
c_en = c_e(1:Nn);
c_es = c_e(Nn+1:Nn+p.Nxs-1);
c_ep = c_e(end-Np+1:end);

% c_e across entire sandwich
c_e_bcs = p.ce.C*c_e;
c_ex = [c_e_bcs(1); c_en; c_e_bcs(2); c_es; c_e_bcs(3); c_ep; c_e_bcs(4)];

% Solid Potential
phi_s_n = z(ind_phi_s_n);
phi_s_p = z(ind_phi_s_p);

% Electrolyte current
i_en = z(ind_ien);
i_ep = z(ind_iep);

% Electrolyte potential
phi_e = z(ind_phi_e);

% Molar ionic flux
jn = z(ind_jn);
jp = z(ind_jp);

% Temperature
T = x(end);

% Preallocate B matrices
B1 = zeros(Nc+1,1);
B2 = zeros(Nz,1);

%% Li Diffusion in Electrolyte Phase: c_e(x,t)
% I'm going to ignore dependence of D_e on c_e in these Jacobians...

% Compute Electrolyte Diffusion Coefficient and Derivative
[D_en0] = electrolyteDe(c_en);
[D_es0] = electrolyteDe(c_es);
[D_ep0] = electrolyteDe(c_ep);

% Adjustment for Arrhenius
Arrh_De = exp(p.E.De/p.R*(1/p.T_ref - 1/T));
D_en = D_en0 * Arrh_De;
D_es = D_es0 * Arrh_De;
D_ep = D_ep0 * Arrh_De;

% ADD BRUGGEMAN RELATION % Apr.22 2016 by Saehong Park
D_en_eff = D_en .* p.epsilon_e_n.^(p.brug-1);
D_es_eff = D_es .* p.epsilon_e_s.^(p.brug-1);
D_ep_eff = D_ep .* p.epsilon_e_p.^(p.brug-1);

% rM1 = [Nn; Ns; Np];
% cM1 = rM1';
% M1 = sparse(blkdiagFast(rM1, cM1, p.ce.M3n, p.ce.M3s, p.ce.M3p));
M2 = [p.ce.M4n, zeros(p.Nxn-1,2); ...
      zeros(p.Nxs-1,1), p.ce.M4s, zeros(p.Nxs-1,1);...
      zeros(p.Nxp-1,2), p.ce.M4p];
G = p.ce.M3 + M2*p.ce.C;

D_e_eff = sparse(diag([D_en_eff; D_es_eff; D_ep_eff]));

dcedce = D_e_eff * G;

% Save Jacobians
f_x(ind_ce,ind_ce) = dcedce;

f_z(ind_ce(1:Nn),ind_jn) = p.ce.M5n;
f_z(ind_ce(Nn+Ns+1:end),ind_jp) = p.ce.M5p;

% if(nargout > 4)
%     diexdI = [0; zeros(Nn,1); ones(p.Nxs+1,1); zeros(Np,1); 0];
%     B1(ind_ce) = B_ce * diexdI;
% end

%% Temperature: T(t)
% Loop through each "comb tooth" in anode
y_csn = p.C_csn * c_s_n_mat;
c_avg_n = y_csn(2,:)';
% for idx = 1:Nn
%     y_csn = p.C_csn * c_s_n_mat(:,idx);
%     c_avg_n(idx) = y_csn(2);
% end

% Loop through each "comb tooth" in cathode
y_csp = p.C_csp * c_s_p_mat;
c_avg_p = y_csp(2,:)';
% for idx = 1:Np
%     y_csp = p.C_csp * c_s_p_mat(:,idx);
%     c_avg_p(idx) = y_csp(2);
% end

% Equilibrium Potential and Gradient wrt bulk concentration
[Unb,dUnb] = refPotentialAnode(p, c_avg_n / p.c_s_n_max);
[Upb,dUpb] = refPotentialCathode(p, c_avg_p / p.c_s_p_max);

% Derivatives wrt c_s
dfTdcsn_int = -(p.a_s_n*p.Faraday * jn .* dUnb) / (p.rho_avg*p.C_p);
dfTdcsn = dfTdcsn_int * p.C_csn(2,:) * p.delta_x_n * p.L_n;

dfTdcsp_int = -(p.a_s_p*p.Faraday * jp .* dUpb) / (p.rho_avg*p.C_p);
dfTdcsp = dfTdcsp_int * p.C_csp(2,:) * p.delta_x_p * p.L_p;

% Derivatives wrt T
dfTdT = -p.h/(p.rho_avg*p.C_p);

% Derivatives wrt jn
dfTdjn_int = -(p.a_s_n*p.Faraday * Unb)/(p.rho_avg*p.C_p);
dfTdjn = dfTdjn_int * p.delta_x_n * p.L_n;

dfTdjp_int = -(p.a_s_p*p.Faraday * Upb)/(p.rho_avg*p.C_p);
dfTdjp = dfTdjp_int * p.delta_x_p * p.L_p;

f_x(end,ind_csn) = reshape(dfTdcsn',1,numel(dfTdcsn));
f_x(end,ind_csp) = reshape(dfTdcsp',1,numel(dfTdcsp));
f_x(end,end) = dfTdT;
f_z(end,ind_jn) = dfTdjn';
f_z(end,ind_jp) = dfTdjp';

if(nargout > 4)
    
    % Terminal Voltage
    phi_s_n_bcs = p.C_psn * phi_s_n + p.D_psn * Cur;
    phi_s_p_bcs = p.C_psp * phi_s_p + p.D_psp * Cur;
    Volt = phi_s_p_bcs(2) - phi_s_n_bcs(1);
    
    B1(end) = (Volt + Cur*(p.D_psp(end)-p.D_psn(1))) / (p.rho_avg*p.C_p);
end

%% Potential in Solid Phase: phi_s(x,t)
if(nargout > 4)
    B2(ind_phi_s_n) = p.F2_psn(:,end) + p.G_psn;
    B2(ind_phi_s_p) = p.F2_psp(:,1) + p.G_psp;
end

%% Electrolyte Current: i_e(x,t)
if(nargout > 4)
    B2(ind_ien) = p.F3_ien;
    B2(ind_iep) = p.F3_iep;
end

%% Potential in Electrolyte Phase: phi_e(x,t)... 
% 2016.10.25
% Disregarding the kappa(c_e) and dactivity(c_e) jacobians for now.

% % System matrices
% [F1_pe,F2_pe,F3_pe] = phi_e_mats_new(p,c_ex);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diffusional Conductivity - electrolyteAct % Oct.25 by Saehong Park

%%%bet_org = (2*p.R*T)/(p.Faraday) * (p.t_plus - 1) * (1 + p.dactivity); %
%When dactivity is constant

[dactivity, ddactivity] = electrolyteAct(c_ex,T,p);
dActivity0 = dactivity(1);                              % BC1
dActivity_n = dactivity(2:p.Nxn);
dActivity_ns = dactivity(p.Nxn+1);
dActivity_s = dactivity(p.Nxn+2 : p.Nxn+2+p.Nxs-2);
dActivity_sp = dactivity(p.Nxn+2+p.Nxs-1);
dActivity_p = dactivity(p.Nxn+2+p.Nxs : end-1);
dActivityN = dactivity(end);                            % BC2

dActivity = [dActivity_n; dActivity_ns; dActivity_s; dActivity_sp; dActivity_p];

bet = (2*p.R*T)/(p.Faraday) * (p.t_plus - 1) * (1 + dActivity);
bet_mat = sparse(diag(bet));

% Modified effective conductivity
%Kap_eff_D_org = bet_org*Kap_eff;% When dactivity is constant
Kap_eff_D = bet_mat*Kap_eff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Form Matrices
M2_pe = p.M2_pe;
M2_pe(1,1) = M2_pe(1,1) * kappa_eff0;
M2_pe(end,end) = M2_pe(end,end) * kappa_effN;

F1_pe = Kap_eff*p.M1_pe + M2_pe*p.C_pe;
F2_pe = p.M3_pe;
F3_pe = Kap_eff_D*p.M4_pe;

% Save Jacobians
g_z(ind_phi_e,ind_phi_e) = F1_pe;

dpedie = F2_pe(:,[1:Nn, (end-Np+1):end]);
g_z(ind_phi_e, [ind_ien, ind_iep]) = dpedie;

% Derivatives w.r.t. c_e here
pe_ce_on_n = [diag(ones(Nn,1)), zeros(Nn, p.Nxs-1+Np)];
pe_ce_on_s = [zeros(p.Nxs-1, Nn), diag(ones(p.Nxs-1,1)), zeros(p.Nxs-1, Np)];
pe_ce_on_p = [zeros(Np, p.Nxs-1+Nn), diag(ones(Np,1))];

F3_pe_bcs = [p.ce.C(1,:); pe_ce_on_n; p.ce.C(2,:); pe_ce_on_s; p.ce.C(3,:); pe_ce_on_p; p.ce.C(4,:)];
c_ex_inv = sparse(diag(1./c_ex));
dpedce = F3_pe * c_ex_inv * F3_pe_bcs;
g_x(ind_phi_e,ind_ce) = dpedce;

% F3_pe_sub = [F3_pe(:, 2:(Nn+1)), F3_pe(:, (Nn+3:Nn+3+(p.Nxs-2))), F3_pe(:, (end-Np):(end-1))];
% 
% g_x(ind_phi_e,ind_ce) = F3_pe_sub * diag(1./c_e);

% g_x(ind_phi_e,ind_ce) = F3_pe(:,[2:Nn+1,Nn+3:Nn+p.Nxs+1,Nn+p.Nxs+3:end-1]) * diag(1./c_e);

if(nargout > 4)
    B2(ind_phi_e) = F2_pe * diexdI;
end

%% Butler-Volmer Equation
% Surface concentration
c_ss_n = (p.C_csn(1,:) * c_s_n_mat)';
c_ss_p = (p.C_csp(1,:) * c_s_p_mat)';

% Param
aFRT = (p.alph*p.Faraday)/(p.R*T);

% Exchange Current Density
[i_0n,i_0p] = exch_cur_dens(p,c_ss_n,c_ss_p,c_e);

di0dcssn = p.k_n*c_en.*(p.c_s_n_max - 2*c_ss_n) ./ ...
            (2 * sqrt(c_en.*c_ss_n.*(p.c_s_n_max - c_ss_n)));
di0dcssp = p.k_p*c_ep.*(p.c_s_p_max - 2*c_ss_p) ./ ...
            (2 * sqrt(c_ep.*c_ss_p.*(p.c_s_p_max - c_ss_p)));

di0dcen = p.k_n*sqrt(c_en.*c_ss_n.*(p.c_s_n_max - c_ss_n)) ./ (2*c_en);
di0dcep = p.k_p*sqrt(c_ep.*c_ss_p.*(p.c_s_p_max - c_ss_p)) ./ (2*c_ep);

% Equilibrium Potential
[Unref,dUnref] = refPotentialAnode(p, c_ss_n / p.c_s_n_max);
[Upref,dUpref] = refPotentialCathode(p, c_ss_p / p.c_s_p_max);

% Overpotential
eta_n = phi_s_n - phi_e(1:Nn) - Unref - p.Faraday*p.R_f_n*jn;
eta_p = phi_s_p - phi_e(end-Np+1:end) - Upref - p.Faraday*p.R_f_p*jp;

% Components of Jacobian
df6dcssn = 2/p.Faraday * di0dcssn .* sinh(aFRT*eta_n) ...
          - 2/p.Faraday * i_0n .* cosh(aFRT*eta_n) * aFRT .* dUnref;
df6dcssp = 2/p.Faraday * di0dcssp .* sinh(aFRT*eta_p) ...
          - 2/p.Faraday * i_0p .* cosh(aFRT*eta_p) * aFRT .* dUpref;
      
df6dcen = 2/p.Faraday * di0dcen .* sinh(aFRT*eta_n);
df6dcep = 2/p.Faraday * di0dcep .* sinh(aFRT*eta_p);

df6psn = 2/p.Faraday*i_0n.*cosh(aFRT*eta_n) * aFRT;
df6psp = 2/p.Faraday*i_0p.*cosh(aFRT*eta_p) * aFRT;

df6pen = -2/p.Faraday*i_0n.*cosh(aFRT*eta_n) * aFRT;
df6pep = -2/p.Faraday*i_0p.*cosh(aFRT*eta_p) * aFRT;

df6Jn = -2/p.Faraday*i_0n.*cosh(aFRT*eta_n) * aFRT*p.R_f_n*p.Faraday - 1;
df6Jp = -2/p.Faraday*i_0p.*cosh(aFRT*eta_p) * aFRT*p.R_f_p*p.Faraday - 1;

df6Tn = 2/p.Faraday * i_0n .* cosh(aFRT * eta_n) * (p.alph*p.Faraday)/(p.R) .* eta_n * (-1/T^2);
df6Tp = 2/p.Faraday * i_0p .* cosh(aFRT * eta_p) * (p.alph*p.Faraday)/(p.R) .* eta_p * (-1/T^2);

% Input into Jacobian

% Loop through each "comb tooth" in anode
%%%  SCOTT, YOU CAN RE-WRITE THIS TO AVOID LOOPS AND BLKDIAGFAST!!!
% Cell_df6n = cell(Nn,1);
% Cell_df6p = cell(Np,1);
% for idx = 1:Nn
%     Cell_df6n{idx} = df6dcssn(idx) * p.C_csn(1,:);
% end
% for idx = 1:Np
%     Cell_df6p{idx} = df6dcssp(idx) * p.C_csp(1,:);
% end
% 
% rsn = ones(Nn,1);
% rsp = ones(Np,1);
% csn = p.PadeOrder * ones(1,Nn);
% csp = p.PadeOrder * ones(1,Np);
% g_x(ind_jn,ind_csn) = blkdiagFast(rsn, csn, Cell_df6n{:});
% g_x(ind_jp,ind_csp) = blkdiagFast(rsp, csp, Cell_df6p{:});

df6dcssn_diag = sparse(diag(df6dcssn));
df6dcssp_diag = sparse(diag(df6dcssp));

g_x(ind_jn,ind_csn) = df6dcssn_diag * g_x(ind_jn,ind_csn);
g_x(ind_jp,ind_csp) = df6dcssp_diag * g_x(ind_jp,ind_csp);



g_x(ind_jn,ind_ce(1:Nn)) = diag(df6dcen);
g_x(ind_jp,ind_ce(end-Np+1:end)) = diag(df6dcep);

g_x(ind_jn,end) = df6Tn;
g_x(ind_jp,end) = df6Tp;



g_z(ind_jn,ind_phi_s_n) = diag(df6psn);
g_z(ind_jp,ind_phi_s_p) = diag(df6psp);

g_z(ind_jn,ind_phi_e(1:Nn)) = diag(df6pen);
g_z(ind_jp,ind_phi_e(end-Np+1:end)) = diag(df6pep);

g_z(ind_jn,ind_jn) = diag(df6Jn);
g_z(ind_jp,ind_jp) = diag(df6Jp);

%%
f_x_full = f_x;
f_z_full = f_z;
g_x_full = g_x;
g_z_full = g_z;

% spy(Jac)
% pause;

if(nargout > 4)
    varargout{1} = B1;
    varargout{2} = B2;
end
