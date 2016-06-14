%% Plot Doyle-Fuller-Newman, Linear Reference Governor
%   Created Sept 13, 2012 by Scott Moura
close all;
clear;

%% Electrolyte Concentration
load('data/new/rg_ce_new.mat');
t = time;
t_MRG=t(1:end-1);
NT = length(time);
I_MRG = cur;
ce0p_MRG = ce0p;
ce0p_MRG=ce0p_MRG(1:end-1);

load('data/new/linrg_ce_new.mat');
t_LMRG=t(1:end-1);
I_LMRG = cur;
ce0p_LMRG = ce0p;
ce0p_LMRG=ce0p_LMRG(1:end-1);

fs3=12;
as3=12;

%Ir(mod(t,20) < 10) = 350;
Ir(mod(t,20) < 10) = 7*67; %670, 10C Discharge for LiCoO2 Cell 67Ah/m^2 for 1C

figure(4)
set(gcf,'Position',[20 20 600 480],'PaperPositionMode','auto')

subplot(2,1,1)
cla
plot(t,I_MRG/67,t,I_LMRG/67,'-.',t,Ir/67,'k--','LineWidth',2);
xlim([0 t(NT)])
set(gca,'XTick',[])
ylim([-0.5 8])
% ylim([-10.5 0.5])
ylabel('Current [C-rate]','FontSize',fs3)
legend({'MRG, $$I(t)$$';'LMRG, $$I(t)$$';'$$I^r(t)$$'},'interpreter','latex')
set(gca,'FontSize',as3)
set(gca,'Position',[0.13 0.57 0.85 0.4])


subplot(2,1,2)
cla

plot(t_MRG,ce0p_MRG/1e3,t_LMRG,ce0p_LMRG/1e3,'-.','LineWidth',2); hold on;
plot(t,0.15*ones(size(t)),'k--','LineWidth',2)
xlim([0 t(NT)])
ylabel('Electrolyte Concentraiton [kmol/m^3]','FontSize',fs3)
xlabel('Time [sec]','FontSize',as3)
legend({'MRG $$c_e(0^+,t)$$';'LMRG $$c_e(0^+,t)$$'},'interpreter','latex')
% legend({'$$c_e(0^+,t)$$'},'interpreter','latex')

ylim([-0.1, 1])
set(gca,'FontSize',fs3)
set(gca,'Position',[0.13 0.13 0.85 0.4])

%% Electrolyte Concentration
load('data/new/rg_etas_new.mat');
t = time;
NT = length(time);
I_MRG = cur;
etas_MRG = eta_s_Ln;

load('data/new/linrg_etas_new.mat');
I_LMRG = cur;
etas_LMRG = eta_s_Ln;

%Ir(mod(t,20) < 10) = -105;
Ir(mod(t,20) < 10) = -3*67; %-201, 3C Charge for LiCoO2 Cell 67Ah/m^2 for 1C

fs4=12;
as4=12;

figure(5)
set(gcf,'Position',[20 20 600 480],'PaperPositionMode','auto')

subplot(2,1,1)
cla
plot(t,I_MRG/67,t,I_LMRG/67,'-.',t,Ir/67,'k--','LineWidth',2);
xlim([0 t(NT)])
set(gca,'XTick',[])
ylim([-3.5 0.5])
% ylim([-10.5 0.5])
ylabel('Current [C-rate]','FontSize',fs4)
legend({'MRG, $$I(t)$$';'LMRG, $$I(t)$$';'$$I^r(t)$$'},'interpreter','latex')
set(gca,'FontSize',as4)
set(gca,'Position',[0.13 0.57 0.85 0.4])


subplot(2,1,2)
cla

plot(t,etas_MRG,t,etas_LMRG,'-.','LineWidth',2); hold on;
plot(t,0*ones(size(t)),'k--','LineWidth',2)
xlim([0 t(NT)])
ylabel('Side Rxn Overpotential [V]','FontSize',fs4)
xlabel('Time [sec]','FontSize',as4)
legend({'MRG $$\eta_s(L^-,t)$$';'LMRG $$\eta_s(L^-,t)$$'},'interpreter','latex')
% legend({'$$c_e(0^+,t)$$'},'interpreter','latex')

ylim([-0.05, 0.24])
set(gca,'FontSize',as4)
set(gca,'Position',[0.13 0.13 0.85 0.4])
