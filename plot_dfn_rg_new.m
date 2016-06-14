%% Plot Doyle-Fuller-Newman, Reference Governor
%   Created July 24, 2012 by Scott Moura
close all;

%% Motivation

% Load Data
load('data/new/dfn_ce_new.mat');
t = time;
NT = length(time);
%I_MRG = cur;
ce0p_DFN = ce0p;
t_DFN=t(1:69);
ce0p_DFN=ce0p_DFN(1:69);
clear out

load('data/new/rg_ce_new.mat');
I = refcur;
I_MRG = cur;
ce0p_MRG = ce0p;
ce0p_MRG = ce0p_MRG(1:end-1);
t_MRG  = t(1:end-1);
clear out

%Ir(mod(t,20) < 10) = 350;
Ir(mod(t,20) < 10) = 670; %670, 10C Discharge for LiCoO2 Cell 67Ah/m^2 for 1C

k = NT;
Ircrate = I/67;
Icrate  = I_MRG/67;

fs1=12;
as1=12;

figure(1)
set(gcf,'Position',[20 20 600 480],'PaperPositionMode','auto')

subplot(2,1,1)
cla
plot(t,Ircrate,t,Icrate,'--','LineWidth',2);
xlim([0 120])
set(gca,'XTick',[])
ylim([-0.5 8])
ylabel('Current [C-rate]','FontSize',fs1)
legend({'$$I^r(t)$$';'MRG, $$I(t)$$'},'interpreter','latex')
set(gca,'FontSize',as1)
set(gca,'Position',[0.13 0.57 0.85 0.4])


subplot(2,1,2)
cla

% plot(t,eta_s_n(end,:),'LineWidth',2); hold on;
% plot(t,zeros(size(t)),'k--','LineWidth',2)
% ce0p = c_ex(end,:);
plot(t_DFN,ce0p_DFN/1e3,'LineWidth',2); hold on;
plot(t_MRG,ce0p_MRG/1e3,'color',[0 .5 0],'linestyle','--','LineWidth',2);hold on;
plot(t,0.15*ones(size(t)),'k--','LineWidth',2)

xlim([0 120])

% ylabel('Side Rxn Overpotential [V]','FontSize',19)
ylabel('Electrolyte Concentration [kmol/m^3]','FontSize',fs1)

xlabel('Time [sec]','FontSize',fs1)

% legend({'$$\eta_s(L^-,t)$$'},'interpreter','latex')
%legend({'$$c_e(0^+,t)$$';'$$c_{e,\min}$$'},'interpreter','latex')
legend({'$$c_e(0^+,t)$$';'MRG, $c_e(0^+,t)$'},'interpreter','latex')

% ylim([0, 1.8])

set(gca,'FontSize',as1)
set(gca,'Position',[0.13 0.13 0.85 0.4])

%% RG Results
clear

% Load Data
load('data/new/dfn_etas_new.mat');
t = time;
NT = length(time);
%I_MRG = cur;
etas_DFN = eta_s_Ln;
t_DFN=t;
etas_DFN=etas_DFN;
clear out

load('data/new/rg_etas_new.mat');
I = refcur;
I_MRG = cur;
etas_MRG = eta_s_Ln;
etas_MRG = etas_MRG(1:end);
t_MRG  = t(1:end);
clear out

k = NT;
Ircrate = I/67;
Icrate  = I_MRG/67;

fs3=12;
as3=12;

figure(2)
set(gcf,'Position',[20 20 600 480],'PaperPositionMode','auto')

subplot(2,1,1)
cla
plot(t,Ircrate,t,Icrate,'--','LineWidth',2);
xlim([0 120])
set(gca,'XTick',[])
ylim([-3.5 0.5])
ylabel('Current [C-rate]','FontSize',fs3)
legend({'$$I^r(t)$$';'MRG, $$I(t)$$'},'interpreter','latex')
set(gca,'FontSize',as3)
set(gca,'Position',[0.13 0.57 0.85 0.4])

subplot(2,1,2)
cla
plot(t_DFN,etas_DFN,'LineWidth',2); hold on;
plot(t_MRG,etas_MRG,'color',[0 .5 0],'linestyle','--','LineWidth',2);hold on;
plot(t,zeros(size(t)),'k--','LineWidth',2)

legend({'$$\eta_s(L^-,t)$$';'MRG, $$\eta_s(L^-,t)$$'},'interpreter','latex')

xlim([0 120])
ylim([-0.2 0.25])

ylabel('Side Rxn Overpotential [V]','FontSize',fs3)
xlabel('Time [sec]','FontSize',fs3)

set(gca,'FontSize',as3)
set(gca,'Position',[0.13 0.13 0.85 0.4])

%% Charging RG Results

% Load Data
load('data/new/baseline_CCCV_new')
t_min_cccv = time/60;
Icrate_cccv = cur/67;
Volt_cccv = volt;
SOC_cccv = soc;
eta_s_cccv = eta_s_Ln;
clear out

load('data/new/rg_CCCV_new')
t_min_rg = time/60;
Icrate_rg = cur/67;
Volt_rg = volt;
SOC_rg = soc;
eta_s_rg = eta_s_Ln;
clear out

fs4=12;
as4=12;

figure(3)
clf
set(gcf,'Position',[20 20 600 640],'PaperPositionMode','auto')

subplot(4,1,1)
cla
plot(t_min_cccv,Icrate_cccv,t_min_rg,Icrate_rg,'--','LineWidth',2);
xlim([-0.2 45])
set(gca,'XTick',[])
%ylim([-0.2 1.2])
ylim([-1.2 0.2])
ylabel('Current [C-rate]','FontSize',fs4)
legend('CCCV','MRG')
set(gca,'FontSize',as4)
set(gca,'Position',[0.13 0.77 0.85 0.21])

subplot(4,1,2)
cla
plot(t_min_cccv,Volt_cccv,t_min_rg,Volt_rg,'--','LineWidth',2);
xlim([-0.2 45])
ylim([3.6 4.4])
set(gca,'XTick',[])
ylabel('Voltage [V]','FontSize',fs4)
% legend({'$$V(t)$$'},'interpreter','latex')
set(gca,'FontSize',as4)
set(gca,'Position',[0.13 0.54 0.85 0.21])

subplot(4,1,3)
cla
plot(t_min_cccv,SOC_cccv,t_min_rg,SOC_rg,'--','LineWidth',2);
xlim([-0.2 45])
ylim([0.55 1.05])
set(gca,'XTick',[])
ylabel('SOC','FontSize',fs4)
% legend({'$$RG$$'},'interpreter','latex')
set(gca,'FontSize',as4)
set(gca,'Position',[0.13 0.32 0.85 0.21])

subplot(4,1,4)
cla
plot(t_min_cccv,eta_s_cccv,t_min_rg,eta_s_rg,'--','LineWidth',2); hold on;
plot(t_min_cccv,zeros(size(t_min_cccv)),'k--','LineWidth',2)
xlim([-0.2 45])
ylim([-0.05, 0.2])
ylabel('Side Rxn Overpotential [V]','FontSize',fs4)
xlabel('Time [min]','FontSize',fs4)
% legend({'$$\eta_s(L^-,t)$$'},'interpreter','latex')
set(gca,'FontSize',as4)
set(gca,'Position',[0.13 0.1 0.85 0.21])
