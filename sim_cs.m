%% Simulate distribution of Li concentration in solid using Pade Approx.
%   Created Oct 2, 2012 by Scott Moura

function [c_sr_n, c_sr_p] = sim_cs(p,t,jn,jp,csn0,csp0)

% Params
Nn = p.Nxn - 1;
Np = p.Nxp - 1;
Nr = 101;
NT = length(t);
p.Nr = Nr;
p.delta_r_n = p.R_s_n/(Nr-1);
p.delta_r_p = p.R_s_p/(Nr-1);

% Load Pade Approximations
% [A_n,A_p,B_n,B_p,C_n,C_p,D_n,D_p] = spm_plant_norm_mats(p);
% A_n = p.D_s_n / p.R_s_n^2 * A_n;
% B_n = -p.D_s_n / p.R_s_n^2 * B_n * p.a_s_n * p.Faraday * p.L_n * p.Area;
% 
% A_p = p.D_s_n / p.R_s_n^2 * A_p;
% B_p = -p.D_s_n / p.R_s_n^2 * B_p * p.a_s_p * p.Faraday * p.L_p * p.Area;

%% Block Matrices for Solid Diffusion
% Electrochemical Model Parameters
alpha_n = (p.D_s_n / 1) / (p.delta_r_n^2);
alpha_p = (p.D_s_p / 1) / (p.delta_r_p^2);

% Block matrices
M1_n = zeros(p.Nr-1);
M1_p = zeros(p.Nr-1);

for idx = 1:(p.Nr-1)
    
    % Upper diagonal
    if(idx ~= 1)
        M1_n(idx,idx-1) = (idx-1)/idx * alpha_n;
        M1_p(idx,idx-1) = (idx-1)/idx * alpha_p;  
    end
    
    % Main diagonal
    M1_n(idx,idx) = -2*alpha_n;
    M1_p(idx,idx) = -2*alpha_p;
    
    % Lower diagonal
    if(idx ~= (p.Nr-1))
        M1_n(idx,idx+1) = (idx+1)/idx * alpha_n;
        M1_p(idx,idx+1) = (idx+1)/idx * alpha_p;
    end
end

M2_n = zeros(p.Nr-1,2);
M2_p = zeros(p.Nr-1,2);

M2_n(end,end) = p.Nr/(p.Nr-1) * alpha_n;
M2_p(end,end) = p.Nr/(p.Nr-1) * alpha_p;

N1 = zeros(2,p.Nr-1);
N1(1,1) = 1;
N1(end,end) = -1;

N2 = diag([-1,1]);

N3_n = [0; -p.delta_r_n/p.D_s_n];
N3_p = [0; -p.delta_r_p/p.D_s_p];

% A,B matrices for each electrode
A_n = M1_n - M2_n*(N2\N1);
A_p = M1_p - M2_p*(N2\N1);

B_n = M2_n*(N2\N3_n);
B_p = M2_p*(N2\N3_p);

% C,D matrices for each electrode
C_n = zeros(1,p.Nr-1);
C_p = zeros(1,p.Nr-1);

C_n(end) = 1;
C_p(end) = 1;

D_n = -p.delta_r_n/p.D_s_n;
D_p = -p.delta_r_p/p.D_s_p;

%%
% Preallocate arrays
c_sr_n = zeros(Nn,Nr,NT);
c_sr_p = zeros(Np,Nr,NT);

% Simulate Particle at x = ii
for ii = 1:Nn
    
    sys = ss(A_n,B_n,C_n,D_n);
    [y,~,x] = lsim(sys,jn(ii,:),t);
    c_sr_n(ii,:,:) = ([x, y])' + csn0;
    
end

% Simulate Particle at x = ii
for ii = 1:Np
    
    sys = ss(A_p,B_p,C_p,D_p);
    [y,~,x] = lsim(sys,jp(ii,:),t);
    c_sr_p(ii,:,:) = ([x, y])' + csp0;
    
end

% 
% for k = 1:NT
%     
%     figure(1)
%     clf
%     contourf(c_sr_n(:,:,k)',51);
%     pause(0.1)
%     
% end


% 
% for k = 1:NT
%     
%     figure(2)
%     plot(1:Nr, c_sr_n(50,:,k))
%     
% end