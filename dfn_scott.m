%% Doyle-Fuller-Newman Model
%   Published June 14, 2016 by Professor Scott Moura
%   Energy, Controls, and Applications Lab (eCAL)
%   University of California, Berkeley
%   http://ecal.berkeley.edu/
clc;
clear;
tic;

disp('Fast DFN')
disp('%%%%%%%%')

%% Electrochemical Model Parameters
% Load Lithium Cobolt Oxide Params, adopted from DUALFOIL
run param/params_LCO


%% Input charge/discharge Current Data %%
% % Current | Positive <=> Discharge, Negative <=> Charge

% Calculate C-rate in terms of [A/m^2] using low/high voltage cutoffs
[cn_low,cp_low] = init_cs(p,p.volt_min);
[cn_high,cp_high] = init_cs(p,p.volt_max);
Delta_cn = cn_high-cn_low;
Delta_cp = cp_low-cp_high;
p.OneC = min(p.epsilon_s_n*p.L_n*Delta_cn*p.Faraday/3600, p.epsilon_s_p*p.L_p*Delta_cp*p.Faraday/3600);

%%%%%%%%%%%%%%% CONSTANT C-RATE %%%%%%%%%%%%%%%%%%%%%%%%%
p.delta_t = 1;
t = 0:p.delta_t:(120);
I = 0*p.OneC*ones(size(t));

%%%%%%%%%%%%%%% DYNAMIC CHARGE/DISCHARGE CYCLES FROM EXPERIMENTS %%%%%%%%%%%%%%%
% load('data/Int_Obs/UDDS_data_Oct_26_2015_Sample_05sec');
% 
% I = -current_exp'/p.Area;
% t = time_exp';
% p.delta_t = t(2)-t(1);

NT = length(t);


%% Initial Conditions & Preallocation
% Solid concentration
V0 = 4.1;
[csn0,csp0] = init_cs(p,V0);

% Electrolyte concentration
ce0 = p.c_e;

% Temperature
T0 = p.T_amb;



% Vector lengths
Ncsn = p.PadeOrder * (p.Nxn-1);
Ncsp = p.PadeOrder * (p.Nxp-1);
Nce = p.Nx - 3;
Nc = Ncsn+Ncsp+Nce;
Nn = p.Nxn - 1;
Ns = p.Nxs - 1;
Np = p.Nxp - 1;
Nnp = Nn+Np;
Nx = p.Nx - 3;
Nz = 3*Nnp + Nx;

% Output Discretization params
disp('Discretization Params');
fprintf(1,'No. of FDM nodes in Anode | Separator | Cathode : %1.0f | %1.0f | %1.0f\n',p.Nxn,p.Nxs,p.Nxp);
fprintf(1,'Order of Pade Approx for Solid Concentration : %1.0f\n',p.PadeOrder);
fprintf(1,'Time Step : %2.2f sec\n',p.delta_t);
disp(' ');

c_s_n0 = zeros(p.PadeOrder,1);
c_s_p0 = zeros(p.PadeOrder,1);

%%%%% Initial condition based on Jordan form
c_s_n0(3) = csn0;
c_s_p0(3) = csp0;

c_s_n = zeros(Ncsn,NT);
c_s_p = zeros(Ncsp,NT);

c_s_n(:,1) = repmat(c_s_n0, [Nn 1]);
c_s_p(:,1) = repmat(c_s_p0, [Nn 1]);

% Electrolyte concentration
c_e = zeros(Nx,NT);
c_e(:,1) = ce0 * ones(Nx,1);

c_ex = zeros(Nx+4,NT);
c_ex(:,1) = c_e(1,1) * ones(Nx+4,1);

% Temperature
T = zeros(NT,1);
T(1) = T0;

% Solid Potential
Uref_n0 = refPotentialAnode(p, csn0(1)*ones(Nn,1) / p.c_s_n_max);
Uref_p0 = refPotentialCathode(p, csp0(1)*ones(Np,1) / p.c_s_p_max);

phi_s_n = zeros(Nn,NT);
phi_s_p = zeros(Np,NT);
phi_s_n(:,1) = Uref_n0;
phi_s_p(:,1) = Uref_p0;

% Electrolyte Current
i_en = zeros(Nn,NT);
i_ep = zeros(Np,NT);

% Electrolyte Potential
phi_e = zeros(Nx+2,NT);

% Molar Ionic Flux
jn = zeros(Nn,NT);
jp = zeros(Np,NT);

% Surface concentration
c_ss_n = zeros(Nn,NT);
c_ss_p = zeros(Np,NT);
c_ss_n(:,1) = repmat(csn0, [Nn 1]);
c_ss_p(:,1) = repmat(csp0, [Np 1]);

% Volume average concentration
c_avg_n = zeros(Nn,NT);
c_avg_p = zeros(Np,NT);
c_avg_n(:,1) = repmat(csn0, [Nn 1]);
c_avg_p(:,1) = repmat(csp0, [Np 1]);

% SOC (Bulk Anode SOC)
SOC = zeros(NT,1);
SOC(1) = (mean(c_avg_n(:,1)) - cn_low) / (cn_high - cn_low);

% Overpotential
eta_n = zeros(Nn,NT);
eta_p = zeros(Np,NT);

% Constraint Outputs
c_e_0n = zeros(NT,1);
c_e_0n(1) = c_ex(1,1);

c_e_0p = zeros(NT,1);
c_e_0p(1) = c_ex(end,1);

eta_s_Ln = zeros(NT,1);
eta_s_Ln(1) = phi_s_p(1,1) - phi_e(1,1);

% Voltage
Volt = zeros(NT,1);
Volt(1) = phi_s_p(end,1) - phi_s_n(1,1) - p.R_c*I(1);

% Conservation of Li-ion matter
n_Li_s = zeros(NT,1);
n_Li_e = zeros(NT,1);

n_Li_e(1) = sum(c_e(1:Nn,1)) * p.L_n*p.delta_x_n * p.epsilon_e_n * p.Area ...
     + sum(c_e(Nn+1:end-Np,1)) * p.L_s*p.delta_x_s * p.epsilon_e_s * p.Area ...
     + sum(c_e(end-Np+1:end,1)) * p.L_p*p.delta_x_p * p.epsilon_e_p * p.Area;

% Stats
newtonStats.iters = zeros(NT,1);
newtonStats.relres = cell(NT,1);
newtonStats.condJac = zeros(NT,1);

% Initial Conditions
x0 = [c_s_n(:,1); c_s_p(:,1); c_e(:,1); T(1)];

z0 = [phi_s_n(:,1); phi_s_p(:,1); i_en(:,1); i_ep(:,1);...
      phi_e(:,1); jn(:,1); jp(:,1)];

%% Preallocate states
x = zeros(length(x0), NT);
z = zeros(length(z0), NT);

x(:,1) = x0;
z(:,1) = z0;

%% Precompute data
% % Solid concentration matrices
% [A_csn,B_csn,A_csp,B_csp,C_csn,C_csp,A_csn_normalized, A_csp_normalized] = c_s_mats(p);
% p.A_csn = A_csn;
% p.A_csn_normalized= A_csn_normalized;
% p.B_csn = B_csn;
% p.A_csp = A_csp;
% p.A_csp_normalized=A_csp_normalized;
% p.B_csp = B_csp;
% p.C_csn = C_csn;
% p.C_csp = C_csp;
% 
% clear A_csn B_csn A_csp B_csp C_csn C_csp A_csn_normalized A_csp_normalized;

% Adjust Temperature Dependent Parameters, based on present temperaure
% Solid concentration matrices
p.D_s_n = p.D_s_n0 * exp(p.E.Dsn/p.R*(1/p.T_ref - 1/T(1)));
p.D_s_p = p.D_s_n0 * exp(p.E.Dsp/p.R*(1/p.T_ref - 1/T(1)));

[A_csn,B_csn,A_csp,B_csp,C_csn,C_csp] = c_s_mats(p);
p.A_csn = A_csn;
p.B_csn = B_csn;
p.A_csp = A_csp;
p.B_csp = B_csp;
p.C_csn = C_csn;
p.C_csp = C_csp;

clear A_csn B_csn A_csp B_csp C_csn C_csp;

% Electrolyte concentration matrices
[M1n,M2n,M3n,M4n,M5n, M1s,M2s,M3s,M4s, M1p,M2p,M3p,M4p,M5p, C_ce] = c_e_mats(p);

p.ce.M1n = M1n;
p.ce.M2n = M2n;
p.ce.M3n = M3n;
p.ce.M4n = M4n;
p.ce.M5n = M5n;

p.ce.M1s = M1s;
p.ce.M2s = M2s;
p.ce.M3s = M3s;
p.ce.M4s = M4s;

p.ce.M1p = M1p;
p.ce.M2p = M2p;
p.ce.M3p = M3p;
p.ce.M4p = M4p;
p.ce.M5p = M5p;

p.ce.C = C_ce;

rM3 = [Nn; Ns; Np];
cM3 = rM3';
p.ce.M3 = sparse(blkdiagFast(rM3, cM3, p.ce.M3n, p.ce.M3s, p.ce.M3p));

clear M1n M2n M3n M4n M5n M1s M2s M3s M4s M1p M2p M3p M4p M5p C_ce rM1 cM1;

% Solid Potential
[F1_psn,F1_psp,F2_psn,F2_psp,G_psn,G_psp,...
    C_psn,C_psp,D_psn,D_psp] = phi_s_mats(p);
p.F1_psn = F1_psn;
p.F1_psp = F1_psp;
p.F2_psn = F2_psn;
p.F2_psp = F2_psp;
p.G_psn = G_psn;
p.G_psp = G_psp;
p.C_psn = C_psn;
p.C_psp = C_psp;
p.D_psn = D_psn;
p.D_psp = D_psp;

clear F1_psn F1_psp F2_psn F2_psp G_psn G_psp C_psn C_psp D_psn D_psp;

% Electrolyte Current
[F1_ien,F1_iep,F2_ien,F2_iep,F3_ien,F3_iep] = i_e_mats(p);
p.F1_ien = F1_ien;
p.F1_iep = F1_iep;
p.F2_ien = F2_ien;
p.F2_iep = F2_iep;
p.F3_ien = F3_ien;
p.F3_iep = F3_iep;

clear F1_ien F1_iep F2_ien F2_iep F3_ien F3_iep;

% Electrolyte Potential
p.M1_pen_skel = sparse(diag(ones(p.Nxn-2,1),1) + diag(-ones(p.Nxn-2,1),-1));
p.M1_pes_skel = sparse(diag(ones(p.Nxs-2,1),1) + diag(-ones(p.Nxs-2,1),-1));
p.M1_pep_skel = sparse(diag(ones(p.Nxp-2,1),1) + diag(-ones(p.Nxp-2,1),-1));

[M1_pe,M2_pe,M3_pe,M4_pe,C_pe] = phi_e_mats(p);
p.M1_pe = M1_pe;
p.M2_pe = M2_pe;
p.M3_pe = M3_pe;
p.M4_pe = M4_pe;
p.C_pe = C_pe;

clear M1_pe M2_pe M3_pe M4_pe C_pe

% Jacobian
[f_x, f_z, g_x, g_z] = jac_dfn_pre(p);
p.f_x = f_x;
p.f_z = f_z;
p.g_x = g_x;
p.g_z = g_z;
clear f_x f_z g_x g_z

%% Integrate!
disp('Simulating DFN Model...');

for k = 1:(NT-1)
    
    % Current
    if(k == 1)
        Cur_vec = [I(k), I(k), I(k+1)];
    else
        Cur_vec = [I(k-1), I(k), I(k+1)];
    end
    
    % Adjust Temperature Dependent Parameters, based on present temperaure
    p.k_n = p.k_n0 * exp(p.E.kn/p.R*(1/p.T_ref - 1/T(k)));
    p.k_p = p.k_p0 * exp(p.E.kp/p.R*(1/p.T_ref - 1/T(k)));
    
    % Solid concentration matrices
    p.D_s_n = p.D_s_n0 * exp(p.E.Dsn/p.R*(1/p.T_ref - 1/T(k)));
    p.D_s_p = p.D_s_n0 * exp(p.E.Dsp/p.R*(1/p.T_ref - 1/T(k)));
    
    [A_csn,B_csn,A_csp,B_csp,C_csn,C_csp] = c_s_mats(p);
    p.A_csn = A_csn;
    p.B_csn = B_csn;
    p.A_csp = A_csp;
    p.B_csp = B_csp;
    p.C_csn = C_csn;
    p.C_csp = C_csp;
    
    clear A_csn B_csn A_csp B_csp C_csn C_csp;
    
    % Step-forward in time
    [x(:,k+1), z(:,k+1), stats] = cn_dfn(x(:,k),z(:,k),Cur_vec,p);

    % Parse out States
    c_s_n(:,k+1) = x(1:Ncsn, k+1);
    c_s_p(:,k+1) = x(Ncsn+1:Ncsn+Ncsp, k+1);
    c_e(:,k+1) = x(Ncsn+Ncsp+1:Nc, k+1);
    T(k+1) = x(end, k+1);
    
    phi_s_n(:,k+1) = z(1:Nn, k+1);
    phi_s_p(:,k+1) = z(Nn+1:Nnp, k+1);
    i_en(:,k+1) = z(Nnp+1:Nnp+Nn, k+1);
    i_ep(:,k+1) = z(Nnp+Nn+1:2*Nnp, k+1);
    phi_e(:,k+1) = z(2*Nnp+1:2*Nnp+Nx+2, k+1);
    jn(:,k+1) = z(2*Nnp+Nx+3:2*Nnp+Nx+Nn+2, k+1);
    jp(:,k+1) = z(2*Nnp+Nx+Nn+3:end, k+1);
    
    newtonStats.iters(k+1) = stats.iters;
    newtonStats.relres{k+1} = stats.relres;
%     newtonStats.condJac(k+1) = stats.condJac;
    
    % Output data
    [~, ~, y] = dae_dfn(x(:,k+1),z(:,k+1),I(k+1),p);
    
    c_ss_n(:,k+1) = y(1:Nn);
    c_ss_p(:,k+1) = y(Nn+1:Nnp);
    
    c_avg_n(:,k+1) = y(Nnp+1:Nnp+Nn);
    c_avg_p(:,k+1) = y(Nnp+Nn+1 : 2*Nnp);
    SOC(k+1) = (mean(c_avg_n(:,k+1)) - cn_low) / (cn_high - cn_low);
    
    c_ex(:,k+1) = y(2*Nnp+1:2*Nnp+Nx+4);
    
    eta_n(:,k+1) = y(2*Nnp+Nx+4+1 : 2*Nnp+Nx+4+Nn);
    eta_p(:,k+1) = y(2*Nnp+Nx+4+Nn+1 : 2*Nnp+Nx+4+Nn+Np);
    
    c_e_0n(k+1) = y(end-5);
    c_e_0p(k+1) = y(end-4);
    eta_s_Ln(k+1) = y(end-3);
    
    Volt(k+1) = y(end-2);
    n_Li_s(k+1) = y(end-1);
    n_Li_e(k+1) = y(end);
    
    eta_s_n = phi_s_n - phi_e(1:Nn,:);
    eta_s_p = phi_s_p - phi_e(end-Np+1:end, :);
    
    fprintf(1,'Time : %3.2f sec | C-rate : %2.2f | Temp : %2.1fdegC | SOC : %1.3f | Voltage : %2.3fV | Newton Iters : %2.0f\n',...
        t(k),I(k+1)/p.OneC,T(k+1)-273.15,SOC(k+1),Volt(k+1),stats.iters);
    
    if(Volt(k+1) < p.volt_min)
        fprintf(1,'Min Voltage of %1.1fV exceeded\n',p.volt_min);
        beep;
        %break;
    elseif(Volt(k+1) > p.volt_max)
        fprintf(1,'Max Voltage of %1.1fV exceeded\n',p.volt_max);
        beep;
        break;
    elseif(any(c_ex(:,k) < 1))
        fprintf(1,'c_e depleted below 1 mol/m^3\n');
        beep;
        %break;
    end

end


%% Outputs
disp('Simulating Output Vars...');
simTime = toc;
fprintf(1,'Simulation Time : %3.2f min\n',simTime/60);
disp('To plots results, run...');
disp(' plot_dfn')
disp(' animate_dfn')


%% Save Output Data for Plotting
% out.date=date;
% out.time=t;
% out.cur=I;
% out.volt=Volt;
% out.soc=SOC;
% out.c_ss_n=c_ss_n;
% out.c_ss_p=c_ss_p;
% out.temp = T;
% % out.eta_s_Ln=eta_s_Ln;
% out.ce0n=c_e_0n;
% out.ce0p=c_e_0p;
% out.simtime=simTime;

% save('data/Int_Obs/dfn_UDDS_NCM20Q.mat', '-struct', 'out');
% save('data/new/dfn_ce_new.mat', '-struct', 'out'); %10C Discharge LiCoO2

%% Save Output Data for Sensitivity Analysis
out.date=date;
out.p = p;
out.time=t;
out.cur=I;
out.volt=Volt;
out.soc=SOC;
out.temp=T;
out.x = x;
out.z = z;
out.simtime=simTime;

save('data/sensitivity/zero_dfn.mat', 'out');



