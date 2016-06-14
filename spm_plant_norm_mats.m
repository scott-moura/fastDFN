%% Matrices for Li Diffusion in Single Particle Model, c_s_p(r,t) & c_s_n(r,t)
%   Created July 21, 2011 by Scott Moura

function [A_n,A_p,B_n,B_p,C_n,C_p,D_n,D_p,varargout] = spm_plant_norm_mats(p)

% Electrochemical Model Parameters
alpha_n = 1 / (p.delta_r_n^2);
alpha_p = (p.D_s_p / p.R_s_p^2) * (p.R_s_n^2 / p.D_s_n) / (p.delta_r_p^2);

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

N3_n = [0; (p.delta_r_n * p.R_s_n)/(p.D_s_n * p.Faraday * p.a_s_n * p.Area * p.L_n)];
N3_p = [0; -(p.delta_r_p * p.R_s_p)/(p.D_s_p * p.Faraday * p.a_s_p * p.Area * p.L_p)];

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

D_n = -(p.delta_r_n * p.R_s_n)/(p.D_s_n * p.Faraday * p.a_s_n * p.Area * p.L_n);
D_p = (p.delta_r_p * p.R_s_p)/(p.D_s_p * p.Faraday * p.a_s_p * p.Area * p.L_p);

% M1,M2,N1,N2,N3 matrices as additional output
if(nargout == 9)
    
    MN_mats.M1_n = M1_n;
    MN_mats.M1_p = M1_p;
    
    MN_mats.M2_n = M2_n;
    MN_mats.M2_p = M2_p;
    
    MN_mats.N1 = N1;
    MN_mats.N2 = N2;
    MN_mats.N3_n = N3_n;
    MN_mats.N3_p = N3_p;
    
    varargout{1} = MN_mats;
end