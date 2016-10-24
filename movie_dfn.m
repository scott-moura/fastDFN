%% Create Movie of Doyle-Fuller-Newman Model Results
%   Created Oct 24, 2016 by Scott Moura
%
%   This m-file creates a movie of the solid and electrolyte
%   concentrations (top and bottom), respectively.

close all;
fs = 16;

createMovie = 1;

%% Create Video Object
if(createMovie)
    v = VideoWriter('vids/Dis_5C.avi');
    v.FrameRate = 8;
    open(v);
end

%% Animation

% Vectors
x_vec = (0:(p.Nx))/(p.Nx);

for k = 1:NT

    fighan = figure(2);
    set(fighan,'Position',[291 85 875 578]);
    clf;
    
    % Create a uicontrol of type "text"
    mTextBox = uicontrol('style','text');
    titleStr = sprintf('Time : %2.0f sec | Voltage : %1.2f V',k,Volt(k));
    set(mTextBox,'String',titleStr,'FontSize',fs,'Position',[330 540 235 30]);
    
    %%% ANODE SOLID PHASE CONCENTRATION
    subplot(2,3,1)
    cla
    plot(1:Nn,c_ss_n(:,k)/p.c_s_n_max,'r-', 1:Nn,c_avg_n(:,k)/p.c_s_n_max,'b-','LineWidth',3);
%     plot(1:Nn,phi_s_n(:,k),'LineWidth',2);
    xlim([1 Nn])
%    ylim([min(min(c_ss_n/p.c_s_n_max)) max(max(c_ss_n/p.c_s_n_max))])
    ylim([0 1]);
%     ylabel('$$c_{s,n}(x,t)$$','interpreter','latex','FontSize',fs)
%     ylabel('$$\phi_{s,n}$$','interpreter','latex')
    ylabel('Solid Conc. of Li [%]','FontSize',fs)
%     xlabel('Space across cell, x [node no.]','interpreter','latex','FontSize',fs)
%    xlabel('Space across cell, x [node no.]','FontSize',fs)
    set(gca,'XTickLabel','');
%     legend({'$$c_{ss,n}(x,t)$$';'$$c_{avg,n}(x,t)$$'},'interpreter','latex','FontSize',fs,0);
    title('\bf ANODE','fontsize',fs)
    set(gca,'FontSize', fs)
    set(gca,'Position',[0.1300 0.4775 0.3100 0.4475]);
    
    
    %%% ANODE SOLID PHASE CONCENTRATION
    subplot(2,3,3)
    cla 
    plot(1:Np,c_ss_p(:,k)/p.c_s_p_max,'r-', 1:Np,c_avg_p(:,k)/p.c_s_p_max,'b-','LineWidth',3);
%     plot(1:Np,phi_s_p(:,k),'LineWidth',2)
    xlim([1 Np])
%     ylim([min(min(c_ss_p/p.c_s_p_max)) max(max(c_ss_p/p.c_s_p_max))])
    ylim([0 1]);
%     ylabel('$$c_{s,p}(x,t)$$','interpreter','latex','FontSize',fs)
%     ylabel('$$\phi_{s,p}$$','interpreter','latex')
    ylabel('Solid Conc. of Li [%]','FontSize',fs)
%     xlabel('Space across cell, x [node no.]','interpreter','latex','FontSize',fs)
%    xlabel('Space across cell, x [node no.]','FontSize',fs)
    set(gca,'XTickLabel','');
%     legend({'$$c_{ss,p}(x,t)$$';'$$c_{avg,p}(x,t)$$'},'interpreter','latex','FontSize',fs,0);
    title('\bf CATHODE','fontsize',fs)
    set(gca,'FontSize', fs)
    set(gca,'Position',[0.5950 0.4775 0.3100 0.4475]);
    
    
    
    %%% SOLID PHASE CONCENTRATION
    subplot(2,3,[4 5 6])
    cla
    plot(x_vec, c_ex(:,k)*1e-3,'LineWidth',3,'color',[0 0.6 0]);
%     plot(1:(p.Nx-3), phi_e(:,k),'LineWidth',2)
%     plot(1:Nx,etaall, 'LineWidth',2)
    xlim([0 1])
%     ylabel('$$c_e(x,t)$$','interpreter','latex','FontSize',fs)
%     ylabel('$$\phi_e$$','interpreter','latex')
    ylim([min(min(c_ex*1e-3)), max(max(c_ex*1e-3))])
%     ylim([min(min(phi_e)), max(max(phi_e))])
    ylabel('Elec Conc. of Li [kmol/m^3]','FontSize', fs)
%     xlabel('Space across cell, x [node no.]','interpreter','latex','FontSize',fs)
    xlabel('Space across cell, x/L','FontSize',fs)
    set(gca,'FontSize', fs)

    pause(0.1);
    
    % Save Frame into video
    if(createMovie)
        frame = getframe(fighan);
        writeVideo(v,frame);
    end

end

%% Close Video
if(createMovie)
    close(v);
end


