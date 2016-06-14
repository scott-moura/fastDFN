%% Sensitivity Eqn Analysis for DFN Model
%   Created Feb 20, 2014 by Scott Moura
%   Modified Apr 15, 2014 by Federico Bribiesca
%
%   This code analyzes the sensitivity of all states in the DFN model to
%   perturbations in the uncertain parameters. These parameters include
%
%   UNCERTAIN PARAMETERS, theta
%   1  : D_s_n
%   2  : D_s_p
%   3  : D_e_n
%   4  : D_e_s
%   5  : D_e_p
%   6  : (1-t_plus)
%   7  : 1/sig_n
%   8  : 1/sig_p
%   9  : 1/kappa
%   10 : (1 + d ln f_ca / d ln c_e)
%   11 : k_n
%   12 : k_p
%   13 : R_f_n
%   14 : R_f_p
%   15 : epsilon_e_n
%   16 : epsilon_e_s
%   17 : epsilon_e_p
%   18 : c_s_n_max
%   19 : c_s_p_max
%   20 : h
%   21 : 1/(rho_avg * C_p)
%
%   OUTPUTS OF INTEREST, y
%   1  : Volt
%   2  : SOC
%   3  : T

clc;
clear;
tic;

%% Load DFN Data

% DFN Data filename
%fn = 'data/sensitivity/zero_dfn.mat';
%fn = 'data/sensitivity/Federico_test_newBC.mat';
%fn='data/sensitivity/Federico_test_newBC_UDDS_500s.mat';
fn='data/sensitivity/Federico_test_newcs_UDDS_500s.mat';
%fn='data/sensitivity/Federico_test_newcs_constant_discharge.mat';
load(fn);
disp(['Loaded DFN data file:  ' fn]);

% Parse output data
t = out.time;
p = out.p;
Cur = out.cur;
x = out.x;
z = out.z;
SOC = out.soc;
Volt = out.volt;

% Vector Lengths
NT = length(t); % Length of time
Nt = 21;        % Number of params
Ny = 3;         % Number of outputs

Ncsn = p.PadeOrder * (p.Nxn-1);
Ncsp = p.PadeOrder * (p.Nxp-1);
Nce = p.Nx - 3;
Nc = Ncsn+Ncsp+Nce;
Nn = p.Nxn - 1;
Np = p.Nxp - 1;
Nnp = Nn+Np;

% Indices
ind_csn = 1:Ncsn;
ind_phi_s_n = 1:Nn;
ind_phi_s_p = Nn+1:Nnp;
ind_T = Nc+1;

clear out;

% Plot DFN Response
fs = 15;
figure(1); clf;
set(gcf,'Position',[311, 16, 613, 684]);

subplot(311)
plot(t,Cur,'LineWidth',2)
ylabel('Current [A]','FontSize',fs);
xlim([0 t(end)]);
set(gca,'Fontsize',fs);

subplot(312)
plot(t,SOC,'LineWidth',2)
ylabel('SOC','FontSize',fs);
xlim([0 t(end)]);
set(gca,'Fontsize',fs);

subplot(313)
plot(t,Volt,'LineWidth',2)
ylabel('Volt','FontSize',fs);
xlabel('Time [sec]','FontSize',fs)
xlim([0 t(end)]);
set(gca,'Fontsize',fs);

%% Nominal Parameters
theta0 = [p.D_s_n;...
    p.D_s_p;
    1; 1; 1;
    (1-p.t_plus);
    1/p.sig_n;
    1/p.sig_p;
    1/1;
    1 + 0;
    p.k_n;
    p.k_p;
    p.R_f_n;
    p.R_f_p;
    p.epsilon_e_n;
    p.epsilon_e_s;
    p.epsilon_e_p;
    p.c_s_n_max;
    p.c_s_p_max;
    p.h;
    (1/p.rho_avg*p.C_p)];
p.theta0 = theta0;

%% Precompute Data
% Solid concentration matrices
[A_csn,B_csn,A_csp,B_csp,C_csn,C_csp, A_csn_normalized, A_csp_normalized] = c_s_mats(p);
p.A_csn = A_csn;
p.A_csn_normalized=A_csn_normalized;
p.B_csn = B_csn;
p.A_csp = A_csp;
p.A_csp_normalized=A_csp_normalized;
p.B_csp = B_csp;
p.C_csn = C_csn;
p.C_csp = C_csp;

% Electrolyte concentration matrices  
c_ex = x(Ncsn+Ncsp+1,1) * ones(p.Nx+4,1);
[trash_var,trash_var,C_ce] = c_e_mats_federico(p,c_ex);
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

% Jacobian w.r.t. states
[A11, A12, A21, A22] = jac_dfn_pre(p);
p.A11 = A11;
p.A12 = A12;
p.A21 = A21;
p.A22 = A22;
clear A11 A12 A21 A22

%% Matrices for output equations
CCx = zeros(Ny, size(x,1));
CCz = zeros(Ny, size(z,1));
DD = zeros(Ny, Nt);

% Voltage
CCz(1,ind_phi_s_p) = p.C_psp(2,:);
CCz(1,ind_phi_s_n) = -p.C_psn(1,:);


% SOC
for i = 1:Nn
    CCx(2,(1:3)+3*(i-1)) = 1/(Nn*p.c_s_n_max) * p.C_csn(2,:);
end

% Temperautre
CCx(3,ind_T) = 1;

%% Preallocate Sensitivity Vars
S1 = zeros(size(x,1), Nt, NT);  % Sensitivity of x vars
S2 = zeros(size(z,1), Nt, NT);  % Sensitivity of z vars
S3 = zeros(Ny, Nt, NT);         % Sensitivity of outputs

%% Solve Sensitivity Eqns
disp('Solving Sensitivity Eqns...');

for k = 1:(NT-1)
        
    %%% Jacobian w.r.t. states, @ current and nxt time step
    [A11, A12, A21, A22] = jac_dfn_federico(x(:,k),z(:,k),Cur(k),p.A11,p.A12,p.A21,p.A22,p);
    [A_corr_p, A_corr_n, A_corr_s] = jac_dfn_sensitivity_federico(x(:,k),z(:,k),Cur(k),p.A11,p.A12,p.A21,p.A22,p);
    
    [A11_nxt, A12_nxt, A21_nxt, A22_nxt] = ...
        jac_dfn_federico(x(:,k+1),z(:,k+1),Cur(k+1),p.A11,p.A12,p.A21,p.A22,p);
    [A_corr_p_nxt, A_corr_n_nxt, A_corr_s_nxt] = jac_dfn_sensitivity_federico(x(:,k+1),z(:,k+1),Cur(k+1),p.A11,p.A12,p.A21,p.A22,p);
    
    %%% Jacobian w.r.t. params, @ current and nxt time step
    [B1, trash_var, B2, trash_var] = jac_p_dfn_federico(x(:,k),z(:,k),Cur(k),p);
    [B1_nxt, trash_var, B2_nxt, trash_var] = jac_p_dfn_federico(x(:,k+1),z(:,k+1),Cur(k+1),p);
    
    %%% Assemble Matrices %% Sensitivity BCs enter as a correction to the
    %%% A*S term
    A = A11 - A12*(A22\A21);
    B = B1 - A12*(A22\B2); 
    
    A_nxt = A11_nxt - A12_nxt*(A22_nxt\A21_nxt);
    B_nxt = B1_nxt - A12_nxt*(A22_nxt\B2_nxt);
    
    C_nxt = -(A22_nxt\A21_nxt);
    D_nxt = -(A22_nxt\B2_nxt);
    
    MM = speye(size(x,1)) - p.delta_t/2 * A_nxt;
    AA = speye(size(x,1)) + p.delta_t/2 * A;
    BB = p.delta_t/2 * (B + B_nxt);
    AA_corr_p=p.delta_t/2 * (A_corr_p + A_corr_p_nxt);
    AA_corr_n=p.delta_t/2 * (A_corr_n + A_corr_n_nxt);
    AA_corr_s=p.delta_t/2 * (A_corr_s + A_corr_s_nxt);
    
    %%% Output sensitivity matrices
    % Jacobian of SOC (y=2) w.r.t. c_s_n_max (theta=18)
    c_s_n =x(1:Ncsn,k+1);
    c_s_n_mat = reshape(c_s_n,p.PadeOrder,p.Nxn-1);
    DD(2,18) = -1/(p.c_s_n_max^2) * mean(p.C_csn(2,:) * c_s_n_mat);
    
    % Jacobian of SOC (y=2) w.r.t. D_s_n (theta=1)
    dGdD =[0, 0 ,0]; %[3/p.R_s_n*6930*p.D_s_n/p.R_s_n^4, 3/p.R_s_n*189/p.R_s_n^2, 0];
    DD(2,1) = 1/(p.c_s_n_max) * mean(dGdD * c_s_n_mat);
    
    %%% Time-stepping

%%%%%%%%%%%%%%%%%%%%%%%% Commented by Federico %%%%%%%%%%%%%%%%%%%%%%%%
%    % Not sure about stability of this method
%    % Forward Euler (Faster, Less Accurate)
%   S1(:,:,k+1) = AA*S1(:,:,k) + p.delta_t/2 * B;
%   S2(:,:,k+1) = C_nxt*S1(:,:,k+1) + D_nxt;
%   S3(:,:,k+1) = CCx*S1(:,:,k+1) + CCz*S2(:,:,k+1) + DD*diag(p.theta0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Central Difference in time (More Accurate, Slower)

    S1(:,:,k+1) = (MM\AA)*S1(:,:,k) + (MM\BB);
    
    %%% Correct for the effect of BCs in sensitivity
     S1(:,3,k+1) = S1(:,3,k+1) + (MM\AA_corr_n);
     S1(:,15,k+1) = S1(:,15,k+1) + (MM\AA_corr_n);
     S1(:,4,k+1) = S1(:,4,k+1) + (MM\AA_corr_s);
     S1(:,16,k+1) = S1(:,16,k+1) + (MM\AA_corr_s);
     S1(:,5,k+1) = S1(:,5,k+1) + (MM\AA_corr_p);
     S1(:,17,k+1) = S1(:,17,k+1) + (MM\AA_corr_p);
    
     S2(:,:,k+1) = C_nxt*S1(:,:,k+1) + D_nxt;
     %%%%%%%%%%%%%%%%%% Addition by Federico %%%%%%%%%%%%%%%%%%%%%%%
     S3(:,:,k+1) = CCx*S1(:,:,k+1) + CCz*S2(:,:,k+1) + DD*diag(p.theta0);
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Error Checking
    if(isnan(any(S1(:,:,k+1))))
        disp('Nan Error');
        disp(S1(:,:,k+1));
    end
    
    %%% Output Progress to Command Prompt
    fprintf(1,'Time : %3.2f sec | Current : %2.4f A/m^2 | SOC : %1.3f | Voltage : %2.4fV\n',...
    t(k),Cur(k+1),SOC(k+1),Volt(k+1));
    
end

%%
simTime = toc;
fprintf(1,'Simulation Time : %3.2f min\n',simTime/60);






