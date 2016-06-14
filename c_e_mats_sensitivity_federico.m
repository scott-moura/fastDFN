%% Matrices for Li Diffusion in Electrolyte Phase, c_e(x,t)
%   Created July 15, 2011 by Scott Moura
%   Modified April 22, 2014 by Federico Bribiesca

function [A,B,C, A_corr_sens_p, A_corr_sens_n, A_corr_sens_s] = c_e_mats_sensitivity_federico(p,c_ex)

% Leading Coefficients
alpha_n = 1 / (p.L_n^2 * p.delta_x_n^2);
alpha_s = 1 / (p.L_s^2 * p.delta_x_s^2);
alpha_p = 1 / (p.L_p^2 * p.delta_x_p^2);

beta_n = (1 - p.t_plus) / (p.epsilon_e_n * p.Faraday * 2 * p.L_n * p.delta_x_n);
beta_s = (1 - p.t_plus) / (p.epsilon_e_s * p.Faraday * 2 * p.L_s * p.delta_x_s);
beta_p = (1 - p.t_plus) / (p.epsilon_e_p * p.Faraday * 2 * p.L_p * p.delta_x_p);

%% Electrolyte Diffusion
De = electrolyteDe(c_ex);

De0 = De(1);
De_n = De(2:p.Nxn);
De_ns = De(p.Nxn+1);
De_s = De(p.Nxn+2 : p.Nxn+2+p.Nxs-2);
De_sp = De(p.Nxn+2+p.Nxs-1);
De_p = De(p.Nxn+2+p.Nxs : end-1);
DeN = De(end);

%% Block Matrices
M1n = zeros(p.Nxn-1);
for idx = 1:p.Nxn-1
    
    if(idx == 1)
        M1n(idx,idx) = alpha_n * -(De0 + De_n(idx+1));
        M1n(idx,idx+1) = alpha_n * De_n(idx+1);
    elseif(idx == p.Nxn-1)
        M1n(idx,idx-1) = alpha_n * De_n(idx-1);
        M1n(idx,idx) = alpha_n * -(De_n(idx-1) + De_ns);
    else
        M1n(idx,idx-1) = alpha_n * De_n(idx-1);
        M1n(idx,idx) = alpha_n * -(De_n(idx-1) + De_n(idx+1));
        M1n(idx,idx+1) = alpha_n * De_n(idx+1);
    end
    
end

M1s = zeros(p.Nxs-1);
for idx = 1:p.Nxs-1
    
    if(idx == 1)
        M1s(idx,idx) = alpha_s * -(De_ns + De_s(idx+1));
        M1s(idx,idx+1) = alpha_s * De_s(idx+1);
    elseif(idx == p.Nxs-1)
        M1s(idx,idx-1) = alpha_s * De_s(idx-1);
        M1s(idx,idx) = alpha_s * -(De_s(idx-1) + De_sp);
    else
        M1s(idx,idx-1) = alpha_s * De_s(idx-1);
        M1s(idx,idx) = alpha_s * -(De_s(idx-1) + De_s(idx+1));
        M1s(idx,idx+1) = alpha_s * De_s(idx+1);
    end
    
end

M1p = zeros(p.Nxp-1);
for idx = 1:p.Nxp-1
    
    if(idx == 1)
        M1p(idx,idx) = alpha_p * -(De_sp + De_p(idx+1));
        M1p(idx,idx+1) = alpha_p * De_p(idx+1);
    elseif(idx == p.Nxp-1)
        M1p(idx,idx-1) = alpha_p * De_p(idx-1);
        M1p(idx,idx) = alpha_p * -(De_p(idx-1) + DeN);
    else
        M1p(idx,idx-1) = alpha_p * De_p(idx-1);
        M1p(idx,idx) = alpha_p * -(De_p(idx-1) + De_p(idx+1));
        M1p(idx,idx+1) = alpha_p * De_p(idx+1);
    end
    
end

rs = [p.Nxn-1; p.Nxs-1; p.Nxp-1];
cs = rs';
M1 = sparse(blkdiagFast(rs,cs,M1n,M1s,M1p));

% M2 : c_e z
M2n = zeros(p.Nxn-1,2);
M2n(1,1) = alpha_n * De0;
M2n(end,end) = alpha_n * De_ns;

M2s = zeros(p.Nxs-1,2);
M2s(1,1) = alpha_s * De_ns;
M2s(end,end) = alpha_s * De_sp;

M2p = zeros(p.Nxp-1,2);
M2p(1,1) = alpha_p * De_sp;
M2p(end,end) = alpha_p * DeN;

M2 = [M2n, zeros(p.Nxn-1,2); ...
      zeros(p.Nxs-1,1), M2s, zeros(p.Nxs-1,1);...
      zeros(p.Nxp-1,2), M2p];
M2 = sparse(M2);  
  
% M3 : i_e
M3 = zeros(p.Nx-3,p.Nx+1);
for idx = 1:p.Nx-3
    
    if(idx <= p.Nxn-1)
        M3(idx,idx) = -beta_n;
        M3(idx,idx+2) = beta_n;
        
    elseif(idx <= p.Nxn+p.Nxs-2)
        M3(idx,idx+1) = -beta_s;
        M3(idx,idx+3) = beta_s;
        
    else
        M3(idx,idx+2) = -beta_p;
        M3(idx,idx+4) = beta_p;
        
    end
    
end

%% Boundary Conditions
N1 = zeros(4,p.Nx-3);
N2 = zeros(4);

N1_BC_n=zeros(4,p.Nx-3); %The normalized impact on the boundary condition is the same for De and \varepsilon (since they enter at the same point)
N1_BC_p=zeros(4,p.Nx-3);
N1_BC_s=zeros(4,p.Nx-3);


% BC1
N1(1,1) = 1/(p.L_n * p.delta_x_n);
N2(1,1) = -1/(p.L_n * p.delta_x_n);


% BC2
N1(2,p.Nxn-1) = p.epsilon_e_n/(p.L_n * p.delta_x_n);
N1(2,p.Nxn) = p.epsilon_e_s/(p.L_s * p.delta_x_s);
N2(2,2) = -p.epsilon_e_n/(p.L_n * p.delta_x_n) - p.epsilon_e_s/(p.L_s * p.delta_x_s);
N1_BC_n(2,p.Nxn-2)=-p.epsilon_e_n/(p.L_n * p.delta_x_n);
N1_BC_n(2,p.Nxn-1)=p.epsilon_e_n/(p.L_n * p.delta_x_n);
N1_BC_s(2,p.Nxn)=p.epsilon_e_n/(p.L_n * p.delta_x_n);
N1_BC_s(2,p.Nxn+1)=-p.epsilon_e_n/(p.L_n * p.delta_x_n);


% BC3
N1(3,p.Nxn+p.Nxs-2) = p.epsilon_e_s/(p.L_s * p.delta_x_s);
N1(3,p.Nxn+p.Nxs-1) = p.epsilon_e_p/(p.L_p * p.delta_x_p);
N2(3,3) = -p.epsilon_e_s/(p.L_s * p.delta_x_s) - p.epsilon_e_p/(p.L_p * p.delta_x_p);
N1_BC_p(3,p.Nxn+p.Nxs-1)=-p.epsilon_e_p/(p.L_p * p.delta_x_p);
N1_BC_p(3,p.Nxn+p.Nxs)=p.epsilon_e_p/(p.L_p * p.delta_x_p);
N1_BC_s(3,p.Nxn+p.Nxs-3)=p.epsilon_e_p/(p.L_p * p.delta_x_p);
N1_BC_s(3,p.Nxn+p.Nxs-2)=-p.epsilon_e_p/(p.L_p * p.delta_x_p);


% BC4
N1(4,end) = -1/(p.L_p * p.delta_x_p);
N2(4,4) = 1/(p.L_p * p.delta_x_p);

%% A,B Matrics
A = sparse(M1 - M2*(N2\N1));
B = sparse(M3);
C = sparse(-N2\N1);
A_corr_sens_p = sparse(-M2*(N2\N1_BC_p));
A_corr_sens_n = sparse(-M2*(N2\N1_BC_n));
A_corr_sens_s = sparse(-M2*(N2\N1_BC_s));


