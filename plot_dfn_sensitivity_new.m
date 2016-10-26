%% Plot Sensitivity Analysis Results of DFN
%   Created February 23, 2014 by Scott Moura

%clear;
close all;
clc;

fs = 18;

params = {'$D_s^{^{\_}}$','$D_s^+$','$R_s^-$','$R_s^+$','$\varepsilon_s^{^{\_}}$','$\varepsilon_s^+$',...
    '$1/\sigma^{^{\_}}$','$1/\sigma^+$','$D_e$','$\varepsilon_e^{^{\_}}$','$\varepsilon_e^{sep}$','$\varepsilon_e^+$',...
    '$\kappa$','$t^{+}_{0}$','\frac{d \ln f_{c/a}}{d \ln c_e}$','$k^{^{\_}}$','$k^+$',...
    '$R_f^{^{\_}}$','$R_f^+$','$n_{Li,s}$','c_{e,0}'};
Nt = length(params);

%'$(1+t^{+}_{0})(1+\frac{d \ln f_{c/a}}{d \ln c_e})$'

%% Load Data;

%%%%%%%%%%%%%%%%%% Start Commented by Federico %%%%%%%%%%%%%%%%%%%%%%
% Load Sensitivities
fn= 'data/sensitivity/sensitivity_5C_meddur.mat';
load(fn);
disp(['Loaded Sensitivity data file:  ' fn]); 

% Parse sensitivity data
dfn_fn = out.fn;
S1 = out.S1;
S2 = out.S2;
S3 = out.S3;
clear out;

% Load DFN Data
load(dfn_fn);
disp(['Loaded DFN data file:  ' fn]);

% Parse output data
t = out.time;
p = out.p;
Cur = out.cur;
SOC = out.soc;
T = out.temp;
Volt = out.volt;
clear out;

Cur_Crate = Cur/p.OneC;

% Vector Sizes
Nt = 21;
NT = length(t);

%%%%%%%%%%%%%%%%%% End Commented by Federico %%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Start Added by Federico %%%%%%%%%%%%%%%%%%%%%%%
% Nt = 21;
% NT =length(t);


%%%%%%%%%%%%%%%%%%% End Added by Federico %%%%%%%%%%%%%%%%%%%%%%%

% Compute Norms over time
Snorm.volt.unsort = zeros(Nt,1);
% Snorm.soc.unsort = zeros(Nt,1);
% Snorm.temp.unsort = zeros(Nt,1);

for idx = 1:Nt
    Snorm.volt.unsort(idx) = norm(squeeze(S3(1,idx,:)));
%     Snorm.soc.unsort(idx) = norm(squeeze(S3(2,idx,1:(NT-3))));
%     Snorm.temp.unsort(idx) = norm(squeeze(S3(3,idx,1:(NT-3))));
end

% Sort Norms
[Snorm.volt.sort, Snorm.volt.ind] = sort(Snorm.volt.unsort,1,'descend');
% [Snorm.soc.sort, Snorm.soc.ind] = sort(Snorm.soc.unsort,1,'descend');
% [Snorm.temp.sort, Snorm.temp.ind] = sort(Snorm.temp.unsort,1,'descend');

S_rank_log = log10(Snorm.volt.sort(end:-1:1));
S_rank_log_unsort = log10(Snorm.volt.unsort);

%% Apply Linear Transformation S'S = D*C*D
S_volt = squeeze(S3(1,1:21,:))';

% Use QR decomposition by successive orthnomralization
[Q,R,E] = qr(S_volt,0);

% Extract Diagonal
D = diag(R);

S_rank_log_ortho = flipud(log10((abs(D))));

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
plot(t,Cur_Crate,'LineWidth',2);
ylabel('Current [C-rate]','FontSize',fs);
xlim([0 t(NT)]);
set(gca,'Fontsize',fs);

subplot(312)
plot(t,Volt,'LineWidth',2)
ylabel('Voltage','FontSize',fs);
xlim([0 t(NT)]);
set(gca,'Fontsize',fs);

subplot(313)
plot(t,squeeze(S3(1,1:18,:)),'LineWidth',2)
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

barh(S_rank_log);


xlim([-10.5 0.5])
ylim([0 Nt+1])
set(gca,'YTick',1:Nt);
set(gca,'Position',[0.2 0.1 0.75 0.85])
%set(gca, 'YTickLabel', params);
[hx,hy] = format_ticks(gca,' ',params(Snorm.volt.ind(end:-1:1)),-10:1:0,[],0,0,0.02,'FontSize',fs,'FontWeight','Bold');
set(gca,'FontSize',fs);
xlabel('Sensitivity Magnitude [log scale]','FontSize',fs)
title('\bf Sensitivity of Voltage','FontSize',fs+2);


% Orthogonalized Sensitivity
figure(6); clf;
set(gcf,'Position',[234  3   580   695],'PaperPositionMode','auto');

barh(S_rank_log_ortho(5:end));
hold on

E_original = fliplr(E);

barh(S_rank_log_unsort(E_original(5:end)),0.4,'FaceColor',[1 0.5 0])
hold off

xlim([-10.5 0.5])
ylim([0 18])
set(gca,'YTick',1:Nt-4);
set(gca,'Position',[0.29 0.1 0.68 0.86])
%set(gca, 'YTickLabel', params);
[hx,hy] = format_ticks(gca,' ',params(E_original(5:end)),-10:1:0,[],0,0,0.02,'FontSize',fs,'FontWeight','Bold');
set(gca,'FontSize',fs);
xlabel('Sensitivity Magnitude [log scale]','FontSize',fs)
title('\bf Sensitivity of Voltage on UDDS cycle','FontSize',fs+2);

legend('Orthogonalized','Original','location','NorthWest')

txt1 = '$(1+t^{+}_{0})(1+\frac{d \ln f_{c/a}}{d \ln c_e})$';
text(-10,16,txt1,'interpreter','latex','FontSize',fs)
