%% Compare DFN (Padé) versus SPM (FDM)
%   Created May 30, 2016 by Scott Moura

clear;
clc;
close all;
fs = 16;

%% Input Data

% Electrochemical Model Parameters
run params_dualfoil_Saehong

% Calculate C-rate in terms of [A/m^2] using low/high voltage cutoffs
[cn_low,cp_low] = init_cs(p,p.volt_min);
[cn_high,cp_high] = init_cs(p,p.volt_max);
Delta_cn = cn_high-cn_low;
Delta_cp = cp_low-cp_high;
OneC = min(p.epsilon_s_n*p.L_n*Delta_cn*p.Faraday/3600, p.epsilon_s_p*p.L_p*Delta_cp*p.Faraday/3600);

% % Constant Current Data
% p.delta_t = 1;
% t = 0:p.delta_t:(360);
% I = 1*OneC*ones(size(t));
% % I(t>=180) = 0;

% % Pulsed Current Data
% p.delta_t = 1;
% t = 0:p.delta_t:(300);
% I = 0*OneC*ones(size(t));
% 
% I((t>=30) & (t<60)) = 1*OneC;
% I((t>=90) & (t<120)) = -1*OneC;
% I((t>=150) & (t<180)) = 2*OneC;
% I((t>=210) & (t<240)) = -2*OneC;

%%%%%%%%%%%%%%%%% INPUT and STATES FROM DFN MODEL %%%%%%%%%%%%%%%%
datafile = 'data/UDDSx2/UDDSx2_DFN_Co_1sec.mat';
% datafile = 'data/Discharge_NEW/Discharge_DFN_Co_1C.mat';
load(datafile);

t = out.time;

% Current | Positive <=> Discharge, Negative <=> Charge
I = out.cur;

NT = length(t);

%% Initial Conditions & Preallocation
% Solid concentration
V0 = 4.0;
[csn0,csp0] = init_cs(p,V0);



%% Generate SPM (FDM) Data

%%% Finite difference for spherical particle
p.Nr = 100; % 100 Make this very large so it closely approximates the true model
Nr = p.Nr;
p.delta_r_n = 1/p.Nr;
p.delta_r_p = 1/p.Nr;
r_vec = (0:p.delta_r_n:1)';
r_vecx = r_vec(2:end-1);

% Construct (A,B) matrices for Li diffusion
[A_n,A_p,B_n,B_p,C_n,C_p,D_n,D_p,MN_mats] = spm_plant_obs_mats(p);

sys_n = ss(A_n,B_n,C_n,D_n);
sys_p = ss(A_p,B_p,C_p,D_p);

disp('Simulating SPMe Plant');

% Initial conditions for SPM
csn0_spm = csn0 * ones(Nr-1,1);
csp0_spm = csp0 * ones(Nr-1,1);

% SIMULATE!
[c_ss_n,t,c_nx] = lsim(sys_n,I,t,csn0_spm);     % Anode Particle
[c_ss_p,t,c_px] = lsim(sys_p,I,t,csp0_spm);     % Cathode Particle

% Normalize
theta_n_spm = c_ss_n/p.c_s_n_max;
theta_p_spm = c_ss_p/p.c_s_p_max;

%% Generate DFN (Pade) Data

csn0_dfn = zeros(p.PadeOrder,1);
csp0_dfn = zeros(p.PadeOrder,1);

%%%%% Initial condition based on controllable canonical form
% c_s_n0(1) = csn0 * (-p.R_s_n/3) * (p.R_s_n^4 / (3465 * p.D_s_n^2));
% c_s_p0(1) = csp0 * (-p.R_s_p/3) * (p.R_s_p^4 / (3465 * p.D_s_p^2));

%%%%% Initial condition based on Jordan form
csn0_dfn(3) = csn0;
csp0_dfn(3) = csp0;

% Solid concentration matrices
[A_csn,B_csn,A_csp,B_csp,C_csn,C_csp] = c_s_mats(p);

sys_n = ss(A_csn,B_csn,C_csn,zeros(2,1));
sys_p = ss(A_csp,B_csp,C_csp,zeros(2,1));

% Input for DFN
jn = I/(p.Faraday*p.a_s_n*p.L_n);
jp = -I/(p.Faraday*p.a_s_p*p.L_p);

% SIMULATE!
[yn,t,c_nx] = lsim(sys_n,jn,t,csn0_dfn);     % Anode Particle
[yp,t,c_px] = lsim(sys_p,jp,t,csp0_dfn);     % Cathode Particle

% Normalize
theta_n_dfn = yn(:,1)/p.c_s_n_max;
theta_p_dfn = yp(:,1)/p.c_s_p_max;

%% Compare / Plot Results

% Compute Root Mean Sqaure Error
rmse_n = rms(theta_n_dfn - theta_n_spm);
rmse_p = rms(theta_p_dfn - theta_p_spm);

% Plot Results
figure(1); clf;

subplot(3,1,1);
plot(t,I/OneC,'LineWidth',2);
ylabel('Current [C-rate]','FontSize',fs)
set(gca,'FontSize',fs);
xlim([t(1) t(end)])

subplot(3,1,2);
plot(t,theta_n_spm,'b-','LineWidth',2); hold on;
plot(t,theta_n_dfn,'r-','LineWidth',2);
legend('SPM (100 FDM nodes)','DFN (3rd order Pade)',0)
ylabel('Anode Conc [-]','FontSize',fs)
title_rms = sprintf('\\bf Anode RMS Error : %0.4f',rmse_n);
title(title_rms,'FontSize',fs)
set(gca,'FontSize',fs);
xlim([t(1) t(end)])

subplot(3,1,3);
plot(t,theta_p_spm,'b-','LineWidth',2); hold on;
plot(t,theta_p_dfn,'r-','LineWidth',2);
legend('SPM (100 FDM nodes)','DFN (3rd order Pade)',0)
ylabel('Cathode Conc [-]','FontSize',fs)
xlabel('Time [sec]','FontSize', fs)
title_rms = sprintf('\\bf Cathode RMS Error : %0.4f',rmse_p);
title(title_rms,'FontSize',fs)
set(gca,'FontSize',fs);
xlim([t(1) t(end)])



% Plot Results
figure(2); clf;

subplot(3,1,1);
plot(t,I/OneC,'LineWidth',2);
ylabel('Current [C-rate]','FontSize',fs)
set(gca,'FontSize',fs);
xlim([t(1) t(end)])

subplot(3,1,2);
plot(t,theta_n_spm-theta_n_dfn,'b-','LineWidth',2);
% legend('SPM (100 FDM nodes)','DFN (3rd order Pade)',0)
ylabel('Anode Conc Err [-]','FontSize',fs)
title_rms = sprintf('\\bf Anode RMS Error : %0.4f',rmse_n);
title(title_rms,'FontSize',fs)
set(gca,'FontSize',fs);
xlim([t(1) t(end)])

subplot(3,1,3);
plot(t,theta_p_spm-theta_p_dfn,'b-','LineWidth',2);
% legend('SPM (100 FDM nodes)','DFN (3rd order Pade)',0)
ylabel('Cathode Conc Err [-]','FontSize',fs)
xlabel('Time [sec]','FontSize', fs)
title_rms = sprintf('\\bf Cathode RMS Error : %0.4f',rmse_p);
title(title_rms,'FontSize',fs)
set(gca,'FontSize',fs);
xlim([t(1) t(end)])