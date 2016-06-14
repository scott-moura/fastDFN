%% Plot Fast Charging Results with Reference Governor
%   Created Sep 19, 2013 by Scott Moura
%   For Motorola ATAP Proposal
clear;

%% Initialize Figure & Data
figure(1)
set(gcf,'Position',[360    99   589   599],'PaperPositionMode','auto');
clf; hold on;
fs = 16;

fn = {'data/0p7CCCV.mat';'data/1p8CCCV.mat';'data/2p5CCCV.mat'};
col = {'k-';'b--';'r-.'};

%% Loop Through Tests
for k = 1:length(fn)

    % Load Data
    load(fn{k});
    
    time = out.time/60;
    cur = -out.cur/35;
    volt = out.volt;
    soc = out.soc;
    temp = out.temp;

    clear out;

    % Plot Data
    subplot(311); hold on;
    plot(time,cur,col{k},'LineWidth',2)
    xlim([time(1), time(end)])
    ylim([-0.05, 2.7])

    subplot(312); hold on;
    plot(time,volt,col{k},'LineWidth',2)
    xlim([time(1), time(end)])
    ylim([3.3 4.4])

    subplot(313); hold on;
    plot(time,soc,col{k},'LineWidth',2)
    xlim([time(1), time(end)])
    
    % Output fast charge times
    ind10 = find(soc >= 0.1,1,'first');
    ind20 = find(soc >= 0.2,1,'first');
    ind50 = find(soc >= 0.5,1,'first');
    
    fprintf(1,[fn{k} '-----\n']);
    fprintf(1,'Time to 10per : %2.2f min\n',time(ind10));
    fprintf(1,'Time to 20per : %2.2f min\n',time(ind20));
    fprintf(1,'Time to 50per : %2.2f min\n',time(ind50));
    
    pause;

end

subplot(311);
legend('0.7C','1.8C','2.5C');
set(gca,'XTickLabel','');
set(gca,'FontSize',fs)
set(gca,'Position',[0.09 0.69 0.88 0.28]);
ylabel('C-rate','FontSize',fs)

subplot(312);
set(gca,'XTickLabel','');
set(gca,'FontSize',fs)
set(gca,'Position',[0.09 0.39 0.88 0.28]);
ylabel('Voltage','FontSize',fs);

subplot(313);
set(gca,'FontSize',fs);
set(gca,'Position',[0.09 0.09 0.88 0.28]);
ylabel('SOC','FontSize',fs)
xlabel('Time [min]','FontSize',fs)