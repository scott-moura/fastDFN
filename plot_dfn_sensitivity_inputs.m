%% Plot Sensitivity Analysis Results of DFN
%   Created February 23, 2014 by Scott Moura

clear;
close all;
clc;

fs = 15;

params = {'$D_s^{^{\_}}$','$D_s^+$','$R_s^-$','$R_s^+$','$\varepsilon_s^{^{\_}}$','$\varepsilon_s^+$',...
    '$1/\sigma^{^{\_}}$','$1/\sigma^+$','$D_e$','$\varepsilon_e^{^{\_}}$','$\varepsilon_e^{sep}$','$\varepsilon_e^+$',...
    '$\kappa$','$t^{+}_{0}$','\frac{d \ln f_{c/a}}{d \ln c_e}$','$k^{^{\_}}$','$k^+$',...
    '$R_f^{^{\_}}$','$R_f^+$','$n_{Li,s}$','c_{e,0}'};
Nt = length(params);

%% Load Data;

%%%%%%%%%%%%%%%%%% Start Commented by Federico %%%%%%%%%%%%%%%%%%%%%%
% Load Sensitivities
fn_cell = {'data/sensitivity/sensitivity_1C_shortdur.mat';...
           'data/sensitivity/sensitivity_Cby2_meddur.mat';...
           'data/sensitivity/sensitivity_5C_meddur.mat'};

S = zeros(3,9);
       
for idx = 1:length(fn_cell)

    load(fn_cell{idx});
    disp(['Loaded Sensitivity data file:  ' fn_cell{idx}]); 

    % Parse sensitivity data
    dfn_fn = out.fn;
    S3 = out.S3;
    clear out;

    % Load DFN Data
    load(dfn_fn);
    disp(['Loaded DFN data file:  ' dfn_fn]);

    % Compute Norms over time
    Snorm.volt.unsort = zeros(Nt,1);
    % Snorm.soc.unsort = zeros(Nt,1);
    % Snorm.temp.unsort = zeros(Nt,1);

    for idx2 = 1:Nt
        Snorm.volt.unsort(idx2) = norm(squeeze(S3(1,idx2,:))) / out.time(end);
    %     Snorm.soc.unsort(idx) = norm(squeeze(S3(2,idx,1:(NT-3))));
    %     Snorm.temp.unsort(idx) = norm(squeeze(S3(3,idx,1:(NT-3))));
    end

    % Save sensitivity data into matrix
    S(idx,:) = Snorm.volt.unsort([1,2,3,4,9,10,12,14,18]);

end

%% Plot DFN Response
figure(1); clf;
set(gcf,'Position',[311, 16, 613, 684]);

subplot(311)
plot(t,Cur,'LineWidth',2)
ylabel('Current [A/m^2]','FontSize',fs);
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

%% Plot Time Repsonse of Sensitivities
figure(2); clf;
set(gcf,'Position',[311, 16, 613, 684]);

subplot(311)
plot(t,out.cur,'LineWidth',2);
ylabel('Current [C-rate]','FontSize',fs);
xlim([0 t(NT)]);
set(gca,'Fontsize',fs);

subplot(312)
plot(t,Volt,'LineWidth',2)
ylabel('Voltage','FontSize',fs);
xlim([0 t(NT)]);
set(gca,'Fontsize',fs);

subplot(313)
plot(t,squeeze(S3(1,:,:)),'LineWidth',2)
ylabel('Voltage Sensitivity','FontSize',fs);
xlim([0 t(NT)]);
xlabel('Time [sec]')
set(gca,'Fontsize',fs);
legend(params,'Interpreter','latex')

% subplot(312)
% plot(t,squeeze(S3(2,:,:)),'LineWidth',2)
% ylabel('SOC Sensitivity','FontSize',fs);
% xlim([0 t(NT)]);
% set(gca,'Fontsize',fs);
% 
% subplot(313)
% plot(t,squeeze(S3(3,:,:)),'LineWidth',2)
% ylabel('Temp. Sensitivity','FontSize',fs);
% xlabel('Time [sec]','FontSize',fs)
% xlim([0 t(NT)]);
% set(gca,'Fontsize',fs);

%% Bar Chart of Sensitivity Ranking

% % SOC Sensitivity
% figure(3); clf;
% set(gcf,'Position',[234     3   564   695],'PaperPositionMode','auto');
% 
% %barh(log10(Snorm.soc.unsort(end:-1:1)));
% barh(log10(Snorm.soc.sort(end:-1:1)));
% xlim([-11 2])
% ylim([0 Nt+1])
% set(gca,'YTick',1:Nt);
% set(gca,'Position',[0.2 0.1 0.75 0.85])
% % set(gca, 'YTickLabel', params(end:-1:1));
% %params_inv = params(end:-1:1);
% params_sorted=params(Snorm.soc.ind(end:-1:1));
% 
% %[hx,hy] = format_ticks(gca,' ',params_inv,-11:2:1,[],0,0,0.02,'FontSize',fs,'FontWeight','Bold');
% [hx,hy] = format_ticks(gca,' ',params_sorted,-11:2:1,[],0,0,0.02,'FontSize',fs,'FontWeight','Bold');
% 
% set(gca,'FontSize',fs);
% 
% xlabel('Sensitivity Magnitude [log scale]','FontSize',fs)
% title('\bf Sensitivity of SOC','FontSize',fs+2);
% 
% 
% % Temperature Sensitivity
% figure(4); clf;
% set(gcf,'Position',[234     3   564   695],'PaperPositionMode','auto');
% 
% barh(log10(Snorm.temp.sort(end:-1:1)));
% xlim([-10.5 0.5])
% ylim([0 Nt+1])
% set(gca,'YTick',1:Nt);
% set(gca,'Position',[0.2 0.1 0.75 0.85])
% %set(gca, 'YTickLabel', params);
% [hx,hy] = format_ticks(gca,' ',params(Snorm.temp.ind(end:-1:1)),-10:1:0,[],0,0,0.02,'FontSize',fs,'FontWeight','Bold');
% set(gca,'FontSize',fs);
% xlabel('Sensitivity Magnitude [log scale]','FontSize',fs)
% title('\bf Sensitivity of Temp.','FontSize',fs+2);

% Voltage Sensitivity
figure(5); clf;
set(gcf,'Position',[234     3   564   695],'PaperPositionMode','auto');

barh(log10(Snorm.volt.unsort(end:-1:1)));
xlim([-10.5 0.5])
ylim([0 Nt+1])
set(gca,'YTick',1:Nt);
set(gca,'Position',[0.2 0.1 0.75 0.85])
%set(gca, 'YTickLabel', params);
[hx,hy] = format_ticks(gca,' ',params(Snorm.volt.ind(end:-1:1)),-10:1:0,[],0,0,0.02,'FontSize',fs,'FontWeight','Bold');
set(gca,'FontSize',fs);
xlabel('Sensitivity Magnitude [log scale]','FontSize',fs)
title('\bf Sensitivity of Voltage','FontSize',fs+2);


% Compare Sensitivities across inputs
% Voltage Sensitivity
figure(6); clf;
set(gcf,'Position',[234     3   564   695],'PaperPositionMode','auto');

barh(log10(S'));
xlim([-19.5 0.5])
ylim([0 10])
set(gca,'YTick',1:9);
set(gca,'Position',[0.2 0.1 0.75 0.85])
%set(gca, 'YTickLabel', params);
[hx,hy] = format_ticks(gca,' ',params([1,2,3,4,9,10,12,14,18]),-19:1:0,[],0,0,0.02,'FontSize',fs,'FontWeight','Bold');
set(gca,'FontSize',fs);
xlabel('Sensitivity Magnitude [log scale]','FontSize',fs)
title('\bf Sensitivity of Voltage','FontSize',fs+2);

legend('1C Short','C/2 Med','5C Med','Location','NorthWest')

