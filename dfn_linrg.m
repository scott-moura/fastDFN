%% Doyle-Fuller-Newman Model with Linear Reference Governor
%   Created September 11, 2012 by Scott Moura
clc;
clear;
tic;

%% Model Construction
% Electrochemical Model Parameters
run ../params_dualfoil

% Vector lengths
Ncsn = p.PadeOrder * (p.Nxn-1);
Ncsp = p.PadeOrder * (p.Nxp-1);
Nce = p.Nx - 3;
Nc = Ncsn+Ncsp+Nce;
Nn = p.Nxn - 1;
Np = p.Nxp - 1;
Nnp = Nn+Np;
Nx = p.Nx - 3;
Nz = 3*Nnp + Nx;

%% Input Signal

% Manual Data
% t = -2:p.delta_t:(60*60);
% Iamp = zeros(length(t),1);
% Iamp(t >= 0) = 35;
% Iamp(t >= (20*60)) = 0;
% Iamp(t >= 20) = 0;
% Iamp(t >= 30) = 10;
% Iamp(t >= 40) = 5;
% I = Iamp;

% Pulse Data
t = -10:p.delta_t:(120);
% Ir(t >= 0) = -35;
Ir(mod(t,20) < 10) = 350; %350
% Ir(end) = 0;
Iamp = Ir;

% Experimental Data
% load('data/UDDSx2_batt_ObsData.mat');
% tdata = t;
% Tfinal = tdata(end);
% t = -2:p.delta_t:Tfinal;
% Iamp = interp1(tdata,I,t,'spline',0);
% Ah_amp = trapz(tdata,I)/3600;
% I = Iamp * (0.4*35)/Ah_amp;

NT = length(t);

%% Initial Conditions & Preallocation
% Reference Govorner Current
I = zeros(NT,1);
beta = zeros(NT,1);

% Solid concentration
csn0 = 0.6 * p.c_s_n_max; % [mols/cm^3] 0.9
csp0 = 0.74 * p.c_s_p_max; % [mols/cm^3] 0.55

c_s_n0 = zeros(p.PadeOrder,1);
c_s_p0 = zeros(p.PadeOrder,1);

c_s_n0(1) = csn0 * (-p.R_s_n/3) * (p.R_s_n^4 / (3465 * p.D_s_n^2));
c_s_p0(1) = csp0 * (-p.R_s_p/3) * (p.R_s_p^4 / (3465 * p.D_s_p^2));

c_s_n = zeros(Ncsn,NT);
c_s_p = zeros(Ncsp,NT);

c_s_n(:,1) = repmat(c_s_n0, [Nn 1]);
c_s_p(:,1) = repmat(c_s_p0, [Nn 1]);

% Electrolyte concentration
c_e = zeros(Nx,NT);
c_e(:,1) = 1e3 * ones(Nx,1);

c_ex = zeros(Nx+4,NT);
c_ex(:,1) = c_e(1,1) * ones(Nx+4,1);

% Temperature
T = zeros(NT,1);
T(1) = p.T_amp;

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
phi_e = zeros(Nx,NT);

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

SOC = zeros(NT,1);
SOC(1) = mean(c_avg_n(:,1)) / p.c_s_n_max;

% Overpotential
eta_n = zeros(Nn,NT);
eta_p = zeros(Np,NT);

% Constraint Outputs
c_e_0p = zeros(NT,1);
eta_s_Ln = zeros(NT,1);

% Voltage
Volt = zeros(NT,1);
Volt(1) = phi_s_p(end,1) - phi_s_n(1,1);

% Conservation of Li-ion matter
nLi = zeros(NT,1);
nLidot = zeros(NT,1);

% Stats
newtonStats.iters = zeros(NT,1);
newtonStats.relres = cell(NT,1);
newtonStats.condJac = zeros(NT,1);

% Initial Conditions
x0 = [c_s_n(:,1); c_s_p(:,1); c_e(:,1); T(1)];

z0 = [phi_s_n(:,1); phi_s_p(:,1); i_en(:,1); i_ep(:,1);...
      phi_e(:,1); jn(:,1); jp(:,1)];

%% Preallocate
x = zeros(length(x0), NT);
z = zeros(length(z0), NT);

I(1) = Ir(1);
x(:,1) = x0;
z(:,1) = z0;

%% Precompute data
% Solid concentration matrices
[A_csn,B_csn,A_csp,B_csp,C_csn,C_csp] = c_s_mats(p);
p.A_csn = A_csn;
p.B_csn = B_csn;
p.A_csp = A_csp;
p.B_csp = B_csp;
p.C_csn = C_csn;
p.C_csp = C_csp;

% Electrolyte concentration matrices
[~,~,C_ce] = c_e_mats(p,c_ex);
p.C_ce = C_ce;

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

% Electrolyte Current
[F1_ien,F1_iep,F2_ien,F2_iep,F3_ien,F3_iep] = i_e_mats(p);
p.F1_ien = F1_ien;
p.F1_iep = F1_iep;
p.F2_ien = F2_ien;
p.F2_iep = F2_iep;
p.F3_ien = F3_ien;
p.F3_iep = F3_iep;

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
    
    % Modified Reference Governor
    if(k > 1 && Ir(k) ~= 0)
        beta(k) = linrg(p, x(:,k-1),z(:,k-1),I(k-1), x(:,k),z(:,k),Ir(k));
    else
        beta(k) = 1;
    end


    I(k+1) = beta(k) * Ir(k);
    
    % Current
    if(k == 1)
        Cur_vec = [I(k), I(k), I(k+1)];
    else
        Cur_vec = [I(k-1), I(k), I(k+1)];
    end
    
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
    phi_e(:,k+1) = z(2*Nnp+1:2*Nnp+Nx, k+1);
    jn(:,k+1) = z(2*Nnp+Nx+1:2*Nnp+Nx+Nn, k+1);
    jp(:,k+1) = z(2*Nnp+Nx+Nn+1:end, k+1);
    
%     i_en(:,k+1)
%     i_ep(:,k+1)
%     
%     phi_s_n(:,k+1)
%     phi_s_p(:,k+1)
%     
%     phi_e(:,k+1)
%     
%     jn(:,k+1)
%     jp(:,k+1)
    
    newtonStats.iters(k+1) = stats.iters;
    newtonStats.relres{k+1} = stats.relres;
%     newtonStats.condJac(k+1) = stats.condJac;
    
    % Output data
    [~, ~, y] = dae_dfn(x(:,k+1),z(:,k+1),I(k+1),p);
    
    c_ss_n(:,k+1) = y(1:Nn);
    c_ss_p(:,k+1) = y(Nn+1:Nnp);
    
    c_avg_n(:,k+1) = y(Nnp+1:Nnp+Nn);
    c_avg_p(:,k+1) = y(Nnp+Nn+1 : 2*Nnp);
    SOC(k+1) = mean(c_avg_n(:,k+1)) / p.c_s_n_max;
    
    c_ex(:,k+1) = y(2*Nnp+1:2*Nnp+Nx+4);
    
    eta_n(:,k+1) = y(2*Nnp+Nx+4+1 : 2*Nnp+Nx+4+Nn);
    eta_p(:,k+1) = y(2*Nnp+Nx+4+Nn+1 : 2*Nnp+Nx+4+Nn+Np);
    
    c_e_0p(k) = y(end-4);
    eta_s_Ln(k) = y(end-3);
    
    Volt(k+1) = y(end-2);
    nLi(k+1) = y(end-1);
    nLidot(k+1) = y(end);
    
    eta_s_n = phi_s_n - phi_e(1:Nn,:);
    eta_s_p = phi_s_p - phi_e(end-Np+1:end, :);
    
    fprintf(1,'Time : %3.2f sec | Current : %2.4f A/m^2 | SOC : %1.3f | Voltage : %2.4fV\n',...
        t(k),I(k+1),SOC(k+1),Volt(k+1));
    
    if(Volt(k+1) < p.volt_min)
        fprintf(1,'Min Voltage of %1.1fV exceeded\n',p.volt_min);
        beep;
        break;
    elseif(Volt(k+1) > p.volt_max)
        fprintf(1,'Max Voltage of %1.1fV exceeded\n',p.volt_max);
        beep;
        break;
    elseif(any(c_ex(:,k) < 1))
        fprintf(1,'c_e depleted below 1 mol/m^3\n');
        beep;
        break;
    end

end


%% Outputs
disp('Simulating Output Vars...');
simTime = toc;
fprintf(1,'Simulation Time : %3.2f min\n',simTime/60);

%% Plot Results

