%% Animate Doyle-Fuller-Newman Model Results for Seminars
%   Created May 12, 2016 by Scott Moura
close all;
fs = 16;

%% Animation for Seminars
Icrate = I/35;

[c_sr_n, c_sr_p] = sim_cs(p,t,jn,jp,csn0,csp0);
th_sr_n = c_sr_n / p.c_s_n_max;
th_sr_p = c_sr_p / p.c_s_p_max;

%%
figure(3)
clf;

% vidObj = VideoWriter('img/dfn_5Cdispulse.avi');
% vidObj.FrameRate = 10;
% vidObj.Quality = 100;
% open(vidObj);

for k = 3:NT
    
    figure(3)
    set(gcf,'Position',[160 -16 850 700],'PaperPositionMode','auto')
    
    % Anode concentration
    vn = linspace(0,0.9,101);
    subplot(3,3,1)
    contourf(fliplr(th_sr_n(:,:,k))',vn);
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    ylabel('Space, r','FontSize',18)
    title('$$c_s^-(x,r,t)$$','FontSize',18,'Interpreter','latex')
    set(gca,'FontSize',16)
    set(gca,'Position',[0.08 0.63 0.35 0.33])
    
    % Cathode concentration
    vp = linspace(0.5,1,101);
    subplot(3,3,3)
    contourf(fliplr(th_sr_p(:,:,k))',vp);
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    ylabel('Space, r','FontSize',18)
    title('$$c_s^+(x,r,t)$$','FontSize',18,'Interpreter','latex')
    set(gca,'FontSize',16)
    set(gca,'Position',[0.63 0.63 0.35 0.33])
    
    % Electrolyte concentration
    subplot(3,3,[4 5 6])
    plot(linspace(0,1,Nx+4),c_ex(:,k)/1e3,'LineWidth',2);
    xlim([0 1])
    ylim([0 1.5])
    xlabh2 = xlabel('Space, x/L','FontSize',18);
    set(xlabh2,'Position',[0.5 -0.15 1]);
    ylabel('$$c_e(x,t)$$','FontSize',18,'Interpreter','latex')
    set(gca,'FontSize',16)
    set(gca,'Position',[0.08 0.34 0.9 0.25])
    
    % Current
    subplot(3,3,[7 8 9])
    plot(t,Icrate,'LineWidth',2); hold on;
%     plot(t(k), Icrate(k), 'ro', 'MarkerSize',16, 'MarkerFaceColor','r'); hold off;
    plot(t,Volt,'Color',[0 0.5 0],'LineWidth',2);
    plot([t(k) t(k)], [0 6],'k--'); hold off;
    xlim([0 t(end)]);
    xlabh = xlabel('Time, t','FontSize',18);
    set(xlabh, 'Position',[120, -0.5, 1])
    ylabel('Current (B) | Voltage (G)','FontSize',18)
    set(gca,'FontSize',16)
    set(gca,'Position',[0.08 0.08 0.9 0.18])
    
%     F = getframe(gcf);
%     writeVideo(vidObj,F);
    
end

% close(vidObj);

%% FREEZE FOR 5C AT 30SEC FOR CDC15 TUTORIAL PAPER 
figure(4)
clf;

k = 32;
    
    figure(4)
    set(gcf,'Position',[162   107   787   474],'PaperPositionMode','auto')
    
    % Anode concentration
    vn = linspace(0,0.9,101);
    subplot(2,3,1)
    contourf(fliplr(th_sr_n(:,:,k))',vn);
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    ylabel('Space, r','FontSize',18)
    title('$$c_s^-(x,r,t) / c_{s,\max}^-$$','FontSize',24,'Interpreter','latex')
    set(gca,'FontSize',16)
    set(gca,'Position',[0.08 0.46 0.36 0.47])
    
    % Cathode concentration
    vp = linspace(0.5,1,101);
    subplot(2,3,3)
    contourf(fliplr(th_sr_p(:,:,k))',vp);
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    ylabel('Space, r','FontSize',18)
    title('$$c_s^+(x,r,t) / c_{s,\max}^+$$','FontSize',24,'Interpreter','latex')
    set(gca,'FontSize',16)
    set(gca,'Position',[0.62 0.46 0.36 0.47])
    
    % Electrolyte concentration
    subplot(2,3,[4 5 6])
    plot(linspace(0,1,Nx+4),c_ex(:,k)/1e3,'LineWidth',2);
    xlim([0 1])
    ylim([0 2])
    xlabh2 = xlabel('Space, x/L','FontSize',18);
    set(xlabh2,'Position',[0.5 -0.15 1]);
    ylabel('$$c_e(x,t)$$','FontSize',24,'Interpreter','latex')
    set(gca,'FontSize',16)
    set(gca,'Position',[0.08 0.1 0.9 0.3])
    
%     % Current
%     subplot(3,3,[7 8 9])
%     plot(t,Icrate,'LineWidth',2); hold on;
% %     plot(t(k), Icrate(k), 'ro', 'MarkerSize',16, 'MarkerFaceColor','r'); hold off;
%     plot(t,Volt,'Color',[0 0.5 0],'LineWidth',2);
%     plot([t(k) t(k)], [0 6],'k--'); hold off;
%     xlim([0 t(end)]);
%     xlabh = xlabel('Time, t','FontSize',18);
%     set(xlabh, 'Position',[120, -0.5, 1])
%     ylabel('Current (B) | Voltage (G)','FontSize',18)
%     set(gca,'FontSize',16)
%     set(gca,'Position',[0.08 0.08 0.9 0.18])
    
%     F = getframe(gcf);
%     writeVideo(vidObj,F);