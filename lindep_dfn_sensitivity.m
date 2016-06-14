%% Linear Dependence of Parameter Sensitivity for DFN Model
%   Created Apr 30, 2014 by Scott Moura
%
%   This code analyzes linear dependence in the DFN model parameter
%   sensitivity equations. The goal is to find which subset of parameters
%   are uniquely identifiable (i.e. linearly dependent)
%
%   UNCERTAIN PARAMETERS, theta
%   1  : D_s_n
%   2  : D_s_p
%   3  : D_e_n
%   4  : D_e_s
%   5  : D_e_p
%   6  : (1-t_plus)
%   7  : 1/sig_n
%   8  : 1/sig_p
%   9  : 1/kappa
%   10 : (1 + d ln f_ca / d ln c_e)
%   11 : k_n
%   12 : k_p
%   13 : R_f_n
%   14 : R_f_p
%   15 : epsilon_e_n
%   16 : epsilon_e_s
%   17 : epsilon_e_p
%   18 : c_s_n_max
%   19 : c_s_p_max
%   20 : h
%   21 : 1/(rho_avg * C_p)
%
%   OUTPUTS OF INTEREST, y
%   1  : Volt
%   2  : SOC
%   3  : T

clc;
clear;

%% Load Sensitivities
fn = 'data/sensitivity/zero_sensitivity.mat';
load(fn);
disp(['Loaded Sensitivity data file:  ' fn]); 

% Parse output sensitivity data (S3)
dfn_fn = out.fn;
S3 = out.S3;
clear out;

% Vector Sizes
Nt = 21;
NT = size(S3,3);

% Parse out sensitivity for each output
S_volt = squeeze(S3(1,:,:))';
S_soc = squeeze(S3(2,:,:))';
S_temp = squeeze(S3(3,:,:))';
clear S3;

%% Parameter Ranking by successive orthogonalization

% Use QR decomposition in economy mode
[Q,R,E] = qr(S_volt,0);

% Extract Diagonal
D = diag(R);