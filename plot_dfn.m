%% Plot Doyle-Fuller-Newman Model Results
%   Created May 23, 2012 by Scott Moura
close all;
fs = 16;


%% Plot Current, SOC, Voltage
figure(1)
clf

subplot(4,1,1)
plot(t,I/p.OneC,'LineWidth',2)
ylabel('Current [C-rate, h^{-1}]','FontSize', fs)
set(gca,'FontSize', fs)
legend({'$$I(t)$$'},'interpreter','latex')
xlim([t(1), t(end)])

subplot(4,1,2)
plot(t,SOC,'LineWidth',2)
ylabel('SOC','FontSize',fs)
set(gca,'FontSize', fs)
legend({'$$SOC(t)$$'},'interpreter','latex')
xlim([t(1), t(end)])

subplot(4,1,3)
plot(t,Volt,'LineWidth',2)
ylabel('Voltage','FontSize', fs)
% xlabel('Time [sec]','FontSize', fs)
legend({'$$V(t)$$'},'interpreter','latex')
set(gca,'FontSize', fs)
xlim([t(1), t(end)])

subplot(4,1,4)
plot(t,T-273.15,'LineWidth',2)
ylabel('Temperature [deg C]','FontSize', fs)
xlabel('Time [sec]','FontSize', fs)
legend({'$$T(t)$$'},'interpreter','latex')
set(gca,'FontSize', fs)
xlim([t(1), t(end)])


%% Plot Current, c_ss(x,t) & c_e(0^\pm,t)
figure(2)
clf

subplot(3,1,1)
plot(t,I/p.OneC,'LineWidth',2)
ylabel('Current [C-rate, h^{-1}]','FontSize', fs)
set(gca,'FontSize', fs)
legend({'$$I(t)$$'},'interpreter','latex')
xlim([t(1), t(end)])

subplot(3,1,2)
plot(t,c_ss_n(1,:)/p.c_s_n_max,'b-','LineWidth',2); hold on;
plot(t,c_ss_n(end,:)/p.c_s_n_max,'b--','LineWidth',2); hold on;
plot(t,c_ss_p(1,:)/p.c_s_p_max,'r--','LineWidth',2); hold on;
plot(t,c_ss_p(end,:)/p.c_s_p_max,'r-','LineWidth',2); hold on;
ylabel('Solid Conc. [Stoic, -]','FontSize',fs)
set(gca,'FontSize', fs)
legend({'$$c_{ss}^-(0^-,t)$$','$$c_{ss}^-(L^-,t)$$',...
    '$$c_{ss}^+(L^+,t)$$','$$c_{ss}^+(0^+,t)$$'},...
    'interpreter','latex','Location','Best')
xlim([t(1), t(end)])

subplot(3,1,3)
plot(t,c_e_0n,'b-',t,c_e_0p,'r-','LineWidth',2)
ylabel('Elec Conc. [mol/m^3]','FontSize', fs)
xlabel('Time [sec]','FontSize', fs)
legend({'$$c_e(0^-,t)$$','$$c_e(0^+,t)$$'},'interpreter','latex')
set(gca,'FontSize', fs)
xlim([t(1), t(end)])