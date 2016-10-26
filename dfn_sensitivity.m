%% Sensitivity Eqn Analysis for DFN Model
%   Originally Created Feb 20, 2014 by Scott Moura
%   Modified Apr 15, 2014 by Federico Bribiesca
%   Rebooted Jul 23, 2016 by Scott Moura for Samsung GRO project
%   Updated Sep 22, 2016 by Scott Moura for Bosch RTC project
%
%   This code analyzes the sensitivity of all states in the DFN model to
%   perturbations in the uncertain parameters. These parameters include
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
%   14 : t_plus
%   15 : d ln f_ca / d ln c_e
%   16 : k_n
%   17 : k_p
%   18 : R_f_n
%   19 : R_f_p
%   20 : n_Li_s
%   21 : c_e0

%
%   OUTPUTS OF INTEREST, y
%   1  : Volt (ONLY!!! for Samsung GRO project)

clc;
clear;
tic;

%% Load DFN Data

% DFN Data filename
fn = 'data/sensitivity/0C_dfn.mat';
load(fn);
disp(['Loaded DFN data file:  ' fn]);

% Parse output data
t = out.time;
p = out.p;
Cur = out.cur;
x = out.x;
z = out.z;
SOC = out.soc;
T = out.temp;
Volt = out.volt;

clear out;



% Vector Lengths
NT = length(t); % Length of time
Nt = 21;        % Number of params
Ny = 1;         % Number of outputs

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

% Indices
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
    p.R_s_n;
    p.R_s_p;
    p.epsilon_s_n;
    p.epsilon_s_p;
    1/p.sig_n;
    1/p.sig_p;
    1;
    p.epsilon_e_n;
    p.epsilon_e_s;
    p.epsilon_e_p;
    1;
    p.t_plus;
    1;
    p.k_n;
    p.k_p;
    p.R_f_n;
    p.R_f_p;
    p.n_Li_s;
    p.c_e];
p.theta0 = theta0;

% Normalization Factor
N=diag(theta0)^2;

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
ind_tplus = 14;     %   14 : t_plus
ind_fca = 15;       %   15 : d ln f_ca / d ln c_e
ind_kn = 16;        %   16 : k_n
ind_kp = 17;        %   17 : k_p
ind_Rfn = 18;       %   18 : R_f_n
ind_Rfp = 19;       %   19 : R_f_p
ind_nLis = 20;      %   20 : n_Li_s
ind_ce0 = 21;       %   21 : c_e0

%% Precompute Data
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

[A_csn,B_csn,A_csp,B_csp,C_csn,C_csp, ...
    A_csn_normalized_D,A_csp_normalized_D,A_csn_normalized_R,A_csp_normalized_R] = c_s_mats(p);
p.A_csn = A_csn;
p.B_csn = B_csn;
p.A_csp = A_csp;
p.B_csp = B_csp;
p.C_csn = C_csn;
p.C_csp = C_csp;

p.A_csn_normalized_D = A_csn_normalized_D;
p.A_csp_normalized_D = A_csp_normalized_D;
p.A_csn_normalized_R = A_csn_normalized_R;
p.A_csp_normalized_R = A_csp_normalized_R;

clear A_csn B_csn A_csp B_csp C_csn C_csp A_csn_normalized_D A_csp_normalized_D A_csn_normalized_R A_csp_normalized_R

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


%% Matrices for output equations
CCx = zeros(Ny, size(x,1));
CCz = zeros(Ny, size(z,1));
DD = zeros(Ny, Nt);

% Voltage Output
CCz(1,ind_phi_s_p) = p.C_psp(2,:);
CCz(1,ind_phi_s_n) = -p.C_psn(1,:);



%% Preallocate Sensitivity Vars
S1 = zeros(size(x,1), Nt, NT);  % Sensitivity of x vars
S2 = zeros(size(z,1), Nt, NT);  % Sensitivity of z vars
S3 = zeros(Ny, Nt, NT);         % Sensitivity of outputs

%% Initialize Sensitivities S_{1,0}
% preallocate IC for sensitivity diff eqns
S10 = zeros(size(x,1), Nt);

% Jacobian of c_s w.r.t. n_Li_s
J_csn0_nLis = zeros(p.PadeOrder,1);
J_csp0_nLis = zeros(p.PadeOrder,1);

J_csn0_nLis(3) = 1/(p.epsilon_s_n*p.L_n*p.Area);
J_csp0_nLis(3) = 1/(p.epsilon_s_p*p.L_p*p.Area);

S10(ind_csn,ind_nLis) = repmat(J_csn0_nLis, [Nn 1]);
S10(ind_csp,ind_nLis) = repmat(J_csp0_nLis, [Np 1]);

% Jacobian of c_e(0) w.r.t. p.c_e;
S10(ind_ce,ind_ce0) = 1;

% save IC into S1 array
S1(:,:,1) = S10;

%% Solve Sensitivity Eqns
disp('Solving Sensitivity Eqns...');

for k = 1:(NT-1)
    
    %%% Jacobian w.r.t. states, @ current and nxt time step
    [A11, A12, A21, A22] = jac_dfn(x(:,k),z(:,k),Cur(k), p.f_x,p.f_z,p.g_x,p.g_z, p);
    [A11_nxt, A12_nxt, A21_nxt, A22_nxt] = jac_dfn(x(:,k+1),z(:,k+1),Cur(k+1), p.f_x,p.f_z,p.g_x,p.g_z, p);
    
    % do i need correction terms for the BCs???
    
    
    %%% Jacobian w.r.t. params, @ current and nxt time step
    [B1, ~, B2, ~] = jac_p_dfn(x(:,k),z(:,k),Cur(k),p);
    [B1_nxt, ~, B2_nxt, ~] = jac_p_dfn(x(:,k+1),z(:,k+1),Cur(k+1),p);
    
    
    
    
    %%% Assemble Matrices %% Sensitivity BCs enter as a correction to the
    %%% A*S term
    A = A11 - A12*(A22\A21);
    B = (B1 - A12*(A22\B2))*N; %Normalized, HEP;
    
    A_nxt = A11_nxt - A12_nxt*(A22_nxt\A21_nxt);
    B_nxt = (B1_nxt - A12_nxt*(A22_nxt\B2_nxt))*N; %Normalized, HEP;
    
    C_nxt = -(A22_nxt\A21_nxt);
    D_nxt = -(A22_nxt\B2_nxt)*N; %Normalized, HEP;
    
    MM = speye(size(x,1)) - p.delta_t/2 * A_nxt;
    AA = speye(size(x,1)) + p.delta_t/2 * A;
    BB = p.delta_t/2 * (B + B_nxt);
%     AA_corr_p=p.delta_t/2 * (A_corr_p + A_corr_p_nxt);
%     AA_corr_n=p.delta_t/2 * (A_corr_n + A_corr_n_nxt);
%     AA_corr_s=p.delta_t/2 * (A_corr_s + A_corr_s_nxt);
    
    
    
    % Central Difference in time (More Accurate, Slower)
    S1(:,:,k+1) = (MM\AA)*S1(:,:,k) + (MM\BB);
    
    %%% Correct for the effect of BCs in sensitivity
%     S1(:,3,k+1) = S1(:,3,k+1) + (MM\AA_corr_n);
%     S1(:,15,k+1) = S1(:,15,k+1) + (MM\AA_corr_n);
%     S1(:,4,k+1) = S1(:,4,k+1) + (MM\AA_corr_s);
%     S1(:,16,k+1) = S1(:,16,k+1) + (MM\AA_corr_s);
%     S1(:,5,k+1) = S1(:,5,k+1) + (MM\AA_corr_p);
%     S1(:,17,k+1) = S1(:,17,k+1) + (MM\AA_corr_p);
    
    S2(:,:,k+1) = C_nxt*S1(:,:,k+1) + D_nxt;
     %%%%%%%%%%%%%%%%%% Addition by Federico %%%%%%%%%%%%%%%%%%%%%%%
    S3(:,:,k+1) = CCx*S1(:,:,k+1) + CCz*S2(:,:,k+1) + DD*N; %Normalized, HEP;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%% Error Checking
    if(isnan(any(S1(:,:,k+1))))
        disp('Nan Error');
        disp(S1(:,:,k+1));
    end
    
    %%% Output Progress to Command Prompt
    fprintf(1,'Time : %3.2f sec | C-rate : %2.2f | Temp : %2.1fdegC | SOC : %1.3f | Voltage : %2.3fV\n',...
        t(k),Cur(k+1)/p.OneC,T(k+1)-273.15,SOC(k+1),Volt(k+1));
    
end


%% Computations on Sensitivity Vector
S3_sqeeze = squeeze(S3);
S3sq = S3_sqeeze*(S3_sqeeze');

sensitivity = diag(S3sq).^(0.5);

