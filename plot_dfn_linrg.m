%% Plot Doyle-Fuller-Newman, Linear Reference Governor
%   Created Sept 13, 2012 by Scott Moura
close all;
clear;

%% Electrolyte Concentration
load('data/rg_ce.mat');
t = out.time;
NT = length(out.time);
I_MRG = out.cur;
ce0p_MRG = out.ce0p;

load('data/linrg_ce.mat');
I_LMRG = out.cur;
ce0p_LMRG = out.ce0p;

Ir(mod(t,20) < 10) = 350;

figure(3)
set(gcf,'Position',[20 20 600 480],'PaperPositionMode','auto')

subplot(2,1,1)
cla
plot(t,-I_MRG/35,t,-I_LMRG/35,'-.',t,-Ir/35,'k--','LineWidth',2);
xlim([0 t(NT)])
set(gca,'XTick',[])
ylim([-10.5 0.5])
% ylim([-10.5 0.5])
ylabel('Current [C-rate]','FontSize',19)
legend({'MRG, $$I(t)$$';'LMRG, $$I(t)$$';'$$I^r(t)$$'},'interpreter','latex')
set(gca,'FontSize',17)
set(gca,'Position',[0.13 0.57 0.85 0.4])


subplot(2,1,2)
cla

plot(t,ce0p_MRG/1e3,t,ce0p_LMRG/1e3,'-.','LineWidth',2); hold on;
plot(t,0.15*ones(size(t)),'k--','LineWidth',2)
xlim([0 t(NT)])
ylabel('Electrolyte Concentraiton [kmol/m^3]','FontSize',19)
xlabel('Time [sec]','FontSize',20)
legend({'MRG $$c_e(0^+,t)$$';'LMRG $$c_e(0^+,t)$$';'$$c_{e,\min}$$'},'interpreter','latex')
% legend({'$$c_e(0^+,t)$$'},'interpreter','latex')

ylim([-0.1, 1])
set(gca,'FontSize',17)
set(gca,'Position',[0.13 0.13 0.85 0.4])

%% Electrolyte Concentration
load('data/rg_etas.mat');
t = out.time;
NT = length(out.time);
I_MRG = out.cur;
etas_MRG = out.eta_s_Ln;

load('data/linrg_etas.mat');
I_LMRG = out.cur;
etas_LMRG = out.eta_s_Ln;

Ir(mod(t,20) < 10) = -105;

figure(4)
set(gcf,'Position',[20 20 600 480],'PaperPositionMode','auto')

subplot(2,1,1)
cla
plot(t,-I_MRG/35,t,-I_LMRG/35,'-.',t,-Ir/35,'k--','LineWidth',2);
xlim([0 t(NT)])
set(gca,'XTick',[])
ylim([-0.5 3.5])
% ylim([-10.5 0.5])
ylabel('Current [C-rate]','FontSize',19)
legend({'MRG, $$I(t)$$';'LMRG, $$I(t)$$';'$$I^r(t)$$'},'interpreter','latex')
set(gca,'FontSize',17)
set(gca,'Position',[0.13 0.57 0.85 0.4])


subplot(2,1,2)
cla

plot(t,etas_MRG,t,etas_LMRG,'-.','LineWidth',2); hold on;
plot(t,0*ones(size(t)),'k--','LineWidth',2)
xlim([0 t(NT)])
ylabel('Side Rxn Overpotential [V]','FontSize',19)
xlabel('Time [sec]','FontSize',20)
legend({'MRG $$\eta_s(L^-,t)$$';'LMRG $$\eta_s(L^-,t)$$';'$$\eta_s = 0$$'},'interpreter','latex')
% legend({'$$c_e(0^+,t)$$'},'interpreter','latex')

ylim([-0.05, 0.15])
set(gca,'FontSize',17)
set(gca,'Position',[0.13 0.13 0.85 0.4])
