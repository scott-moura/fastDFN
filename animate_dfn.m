%% Animate Doyle-Fuller-Newman Model Results
%   Created May 12, 2016 by Scott Moura
%
%   This m-file plots an animation of the solid and electrolyte
%   concentrations (top and bottom), respectively. The middle frame plots
%   a user-selected electrochemical state, given by the flag option below.
%
%   1 : j_n,    molar flux
%   2 : i_e,    electrolyte current
%   3 : eta,    overpotential
%   4 : phi_s,  solid potential
%   5 : phi_e,  electrolyte potential


flag = 1;

close all;
fs = 16;


%% Animation

for k = 1:NT

    figure(2)
    set(gcf,'Position',[103 3 1151 673]);
    
    subplot(3,3,1)
    cla
    plot(1:Nn,c_ss_n(:,k)/p.c_s_n_max,'r-', 1:Nn,c_avg_n(:,k)/p.c_s_n_max,'b-','LineWidth',2);
%     plot(1:Nn,phi_s_n(:,k),'LineWidth',2);
    xlim([1 Nn])
    ylim([min(min(c_ss_n/p.c_s_n_max)) max(max(c_ss_n/p.c_s_n_max))])
%     ylim([min(min(phi_s_n)) max(max(phi_s_n))])
%     ylabel('$$c_{s,n}(x,t)$$','interpreter','latex','FontSize',fs)
%     ylabel('$$\phi_{s,n}$$','interpreter','latex')
    ylabel('Solid Conc. [Stoic, -]','FontSize',fs)
%     xlabel('Space across cell, x [node no.]','interpreter','latex','FontSize',fs)
    xlabel('Space across cell, x [node no.]','FontSize',fs)
%     legend({'$$c_{ss,n}(x,t)$$';'$$c_{avg,n}(x,t)$$'},'interpreter','latex','FontSize',fs,0);
    title('\bf ANODE','fontsize',fs)
    set(gca,'FontSize', fs)
    set(gca,'Position',[0.1300 0.7093 0.31 0.2157]);
    
    subplot(3,3,3)
    cla 
    plot(1:Np,c_ss_p(:,k)/p.c_s_p_max,'r-', 1:Np,c_avg_p(:,k)/p.c_s_p_max,'b-','LineWidth',2);
%     plot(1:Np,phi_s_p(:,k),'LineWidth',2)
    xlim([1 Np])
    ylim([min(min(c_ss_p/p.c_s_p_max)) max(max(c_ss_p/p.c_s_p_max))])
%     ylim([min(min(phi_s_p)) max(max(phi_s_p))])
%     ylabel('$$c_{s,p}(x,t)$$','interpreter','latex','FontSize',fs)
%     ylabel('$$\phi_{s,p}$$','interpreter','latex')
    ylabel('Solid Conc. [Stoic, -]','FontSize',fs)
%     xlabel('Space across cell, x [node no.]','interpreter','latex','FontSize',fs)
    xlabel('Space across cell, x [node no.]','FontSize',fs)
%     legend({'$$c_{ss,p}(x,t)$$';'$$c_{avg,p}(x,t)$$'},'interpreter','latex','FontSize',fs,0);
    title('\bf CATHODE','fontsize',fs)
    set(gca,'FontSize', fs)
    set(gca,'Position',[0.5950 0.7093 0.31 0.2157]);
    
    jall = [jn(:,k); zeros(p.Nxs-1,1); jp(:,k)];
    etaall = [eta_n(:,k); zeros(p.Nxs-1,1); eta_p(:,k)];
    ieall = [i_en(:,k); I(k)*ones(p.Nxs-1,1); i_ep(:,k)];
    phisall = [phi_s_n(:,k); zeros(p.Nxs-1,1); phi_s_p(:,k)];
    phieall = phi_e(:,k);
    set(gca,'FontSize', fs)

    subplot(3,3,[4 5 6])
    cla
    
    switch flag
        case 1
            plot(1:Nx,jall,'LineWidth',2);
            ylabel('Molar Flux [mol/(m^2*s)]');
            ylim([1.05*min(min([jn; jp])), 1.05*max(max([jn; jp]))]);
        case 2
            plot(1:Nx,ieall,'LineWidth',2);
            ylabel('Elec. Current [A/m^2]');
            ylim([min(min([i_en; i_ep])), 1.05*max(max([i_en; i_ep]))]);
        case 3
            plot(1:Nx,etaall,'LineWidth',2);
            ylabel('Overpotential [V]');
            ylim([1.05*min(min([eta_n; eta_p])), 1.05*max(max([eta_n; eta_p]))]);
        case 4
            plot(1:Nx,phisall,'LineWidth',2);
            ylabel('Solid Potential [V]');
            ylim([min(min([phi_s_n; phi_s_p])), 1.05*max(max([phi_s_n; phi_s_p]))]);
        case 5
            plot(1:(Nx+2),phieall,'LineWidth',2);
            ylabel('Elec. Potential [V]');
            ylim([min(min(phi_e)), 1.05*max(max(phi_e))]);
    end
            
    xlim([1 Nx])
    xlabel('Space across cell, x [node no.]','FontSize',fs)
    set(gca,'FontSize', fs)
    
    tit_time = sprintf('\b Time : %3.0f',k);
    title(tit_time,'fontsize',fs);

%     subplot(3,3,6)
%     cla
%     plot(1:Np,jp(:,k),'LineWidth',2);
%     xlim([1 Np])
% %     ylim([min(min(i_ep)), max(max(i_ep))])
%     ylabel('$$i_{e,p}(x,t)$$','interpreter','latex')

    subplot(3,3,[7 8 9])
    cla
    plot(1:(p.Nx+1), c_ex(:,k),'LineWidth',2);
%     plot(1:(p.Nx-3), phi_e(:,k),'LineWidth',2)
%     plot(1:Nx,etaall, 'LineWidth',2)
    xlim([1 p.Nx+1])
%     ylabel('$$c_e(x,t)$$','interpreter','latex','FontSize',fs)
%     ylabel('$$\phi_e$$','interpreter','latex')
    ylim([min(min(c_ex)), max(max(c_ex))])
%     ylim([min(min(phi_e)), max(max(phi_e))])
    ylabel('Elec Conc. [mol/m^3]','FontSize', fs)
%     xlabel('Space across cell, x [node no.]','interpreter','latex','FontSize',fs)
    xlabel('Space across cell, x [node no.]','FontSize',fs)
    set(gca,'FontSize', fs)

    pause(0.1);

end