%% Plot Doyle-Fuller-Newman, Reference Governor
%   Created July 24, 2012 by Scott Moura
close all;

%% Motivation

k = NT;
Icrate = I/35;

figure(1)
set(gcf,'Position',[20 20 600 480],'PaperPositionMode','auto')

subplot(3,1,1)
cla
plot(t,-Icrate,'LineWidth',2);
xlim([0 NT])
set(gca,'XTick',[])
ylim([-10.5 0.5])
ylabel('Current [C-rate]','FontSize',19)
legend({'$$I(t)$$'},'interpreter','latex')
set(gca,'FontSize',17)
set(gca,'Position',[0.13 0.7 0.85 0.27])

subplot(3,1,2)
cla
plot(t,Volt,'LineWidth',2);
xlim([0 NT])
ylim([3.2 3.9])
set(gca,'XTick',[])
ylabel('Voltage [V]','FontSize',19)
legend({'$$V(t)$$'},'interpreter','latex')
set(gca,'FontSize',17)
set(gca,'Position',[0.13 0.4 0.85 0.27])


subplot(3,1,3)
cla

% plot(t,eta_s_n(end,:),'LineWidth',2); hold on;
% plot(t,zeros(size(t)),'k--','LineWidth',2)
% ce0p = c_ex(end,:);
plot(t,ce0p/1e3,'LineWidth',2); hold on;
plot(t,0.15*ones(size(t)),'k--','LineWidth',2)

xlim([0 NT])

% ylabel('Side Rxn Overpotential [V]','FontSize',19)
ylabel('Electrolyte Concentraiton [kmol/m^3]','FontSize',19)

xlabel('Time [sec]','FontSize',20)

% legend({'$$\eta_s(L^-,t)$$'},'interpreter','latex')
legend({'$$c_e(0^+,t)$$';'$$c_{e,\min}$$'},'interpreter','latex')

% ylim([0, 1.8])
set(gca,'FontSize',17)
set(gca,'Position',[0.13 0.1 0.85 0.27])

%% RG Results
k = NT;
Icrate = I/35;
% Ircrate = Ir/35;

figure(3)
set(gcf,'Position',[20 20 600 480],'PaperPositionMode','auto')

subplot(2,1,1)
cla
plot(t,-Icrate,'--','color',[0 0.5 0],'LineWidth',2);
xlim([0 NT])
set(gca,'XTick',[])
ylim([-0.5 3.5])
% ylim([-10.5 0.5])
ylabel('Current [C-rate]','FontSize',19)
legend({'$$I(t)$$';'$$I^r(t)$$'},'interpreter','latex')
set(gca,'FontSize',17)
set(gca,'Position',[0.13 0.57 0.85 0.4])


subplot(2,1,2)
cla

plot(t,eta_s_n(:),'--','color',[0 0.5 0],'LineWidth',2); hold on;
plot(t,zeros(size(t)),'k--','LineWidth',2)
% plot(t,c_e(end,:)/1e3,'LineWidth',2); hold on;
% plot(t,0.1*ones(size(t)),'k--','LineWidth',2)

xlim([0 NT])

ylabel('Side Rxn Overpotential [V]','FontSize',19)
% ylabel('Electrolyte Concentraiton [kmol/m^3]','FontSize',19)

xlabel('Time [sec]','FontSize',20)

legend({'$$\eta_s(L^-,t)$$'},'interpreter','latex')
% legend({'$$c_e(0^+,t)$$'},'interpreter','latex')

% ylim([0, 1.8])
set(gca,'FontSize',17)
set(gca,'Position',[0.13 0.13 0.85 0.4])

%% Charging RG Results

% Load Data
load('data/rg_chg_cccv')
t_min_cccv = out.time/60;
Icrate_cccv = -out.cur/35;
Volt_cccv = out.volt;
SOC_cccv = out.soc;
eta_s_cccv = out.eta_s_Ln;
clear out

load('data/rg_chg_rg')
t_min_rg = out.time/60;
Icrate_rg = -out.cur/35;
Volt_rg = out.volt;
SOC_rg = out.soc;
eta_s_rg = out.eta_s_Ln;
clear out

figure(4)
clf
set(gcf,'Position',[20 20 600 640],'PaperPositionMode','auto')

subplot(4,1,1)
cla
plot(t_min_cccv,Icrate_cccv,t_min_rg,Icrate_rg,'--','LineWidth',2);
xlim([-0.2 45])
set(gca,'XTick',[])
ylim([-0.2 1.2])
ylabel('Current [C-rate]','FontSize',19)
legend('CCCV','MRG')
set(gca,'FontSize',17)
set(gca,'Position',[0.13 0.77 0.85 0.21])

subplot(4,1,2)
cla
plot(t_min_cccv,Volt_cccv,t_min_rg,Volt_rg,'--','LineWidth',2);
xlim([-0.2 45])
ylim([3.6 4.4])
set(gca,'XTick',[])
ylabel('Voltage [V]','FontSize',19)
% legend({'$$V(t)$$'},'interpreter','latex')
set(gca,'FontSize',17)
set(gca,'Position',[0.13 0.54 0.85 0.21])

subplot(4,1,3)
cla
plot(t_min_cccv,SOC_cccv,t_min_rg,SOC_rg,'--','LineWidth',2);
xlim([-0.2 45])
ylim([0.55 1.05])
set(gca,'XTick',[])
ylabel('SOC','FontSize',19)
% legend({'$$RG$$'},'interpreter','latex')
set(gca,'FontSize',17)
set(gca,'Position',[0.13 0.32 0.85 0.21])

subplot(4,1,4)
cla
plot(t_min_cccv,eta_s_cccv,t_min_rg,eta_s_rg,'--','LineWidth',2); hold on;
plot(t_min,zeros(size(t)),'k--','LineWidth',2)
xlim([-0.2 45])
ylim([-0.05, 0.2])
ylabel('Side Rxn Overpotential [V]','FontSize',19)
xlabel('Time [min]','FontSize',20)
% legend({'$$\eta_s(L^-,t)$$'},'interpreter','latex')
set(gca,'FontSize',17)
set(gca,'Position',[0.13 0.1 0.85 0.21])
