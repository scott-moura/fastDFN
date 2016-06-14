%% Plot Doyle-Fuller-Newman Model Results
%   Created July 12, 2012 by Scott Moura
close all;

%% Animation

k = NT;

figure(2)
set(gcf,'Position',[20 20 640 480],'PaperPositionMode','auto')

subplot(2,3,1)
cla
plot(1:Nn,c_avg_n(:,k)*1e-3,1:Nn,c_ss_n(:,k)*1e-3,'--','LineWidth',2);
xlim([1 Nn])
ylim([5 20])
set(gca,'XTick',[])
ylabel('Anode Li Concentration, [kmol/m^3]','FontSize',19)
legend({'$$c_{avg}^-(x,t)$$'; '$$c_{ss}^-(x,t)$$'},'interpreter','latex')
set(gca,'FontSize',17)
set(gca,'Position',[0.12 0.55 0.35 0.4])

subplot(2,3,3)
cla
plot(1:Np,c_avg_p(:,k)*1e-3,1:Np,c_ss_p(:,k)*1e-3,'--','LineWidth',2);
xlim([1 Np])
ylim([20 28])
set(gca,'XTick',[])
ylabel('Cathode Li Concentration, [kmol/m^3]','FontSize',19)
legend({'$$c_{avg}^+(x,t)$$'; '$$c_{ss}^+(x,t)$$'},'interpreter','latex')
set(gca,'FontSize',17)
set(gca,'Position',[0.62 0.55 0.35 0.4])


subplot(2,3,[4 5 6])
cla
plot(1:(p.Nx+1), c_ex(:,k)*1e-3,'LineWidth',2);
xlim([1 p.Nx+1])
set(gca,'XTick',[])
ylabel('Electrolyte Li Concentration, [kmol/m^3]','FontSize',19)
xlabel('Space, x','FontSize',20)
legend({'$$c_{e}(x,t)$$'},'interpreter','latex')
ylim([0, 1.8])
set(gca,'FontSize',17)
set(gca,'Position',[0.12 0.1 0.85 0.4])