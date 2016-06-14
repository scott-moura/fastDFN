%% Output Eqn for Constrained Vars in Doyle-Fuller-Newman Model
%   Created Aug 30, 2012 by Scott Moura
%   Constraints are output in negative null form

function [C1,C2,D,E] = cons_dfn(p)

%% Constraint Parameters
Imin = -700;
Imax = 700;

theta_n_min = 0.01;
theta_n_max = 0.9; %0.9
theta_p_min = 0.5;
theta_p_max = 0.99;

c_e_min = 150;
c_e_max = 3000;

T_min = 273;
T_max = 273 + 50;

Us = 0;

%% Parse out states
Ncsn = p.PadeOrder * (p.Nxn-1);
Ncsp = p.PadeOrder * (p.Nxp-1);
Nce = p.Nx - 3;
Nc = Ncsn+Ncsp+Nce;
Nn = p.Nxn - 1;
Np = p.Nxp - 1;
Nnp = Nn+Np;
Nx = p.Nx - 3;
Ns = Nx - Nn - Np;
Nz = 3*Nnp + Nx;

%% Preallocate Matrices
Ny = 2 + 2*Nn+2*Np + 2*(p.Nx+1) + 2 + 1;

C1 = zeros(Ny,Nc+1);
C2 = zeros(Ny,Nz);
D = zeros(Ny,1);

%% Indices for constraints
indy_I = 1:2;
indy_cssn = 3:2+(2*Nn);
indy_cssp = indy_cssn(end) + (1:(2*Np));
indy_ce = indy_cssp(end) + (1:(2*(p.Nx+1)));
indy_T = indy_ce(end) + (1:2);
indy_etas = indy_T(end) + 1;

%% Indices for states
ind_csn = 1:Ncsn;
ind_csp = Ncsn+1:Ncsn+Ncsp;

ind_cs = 1:Ncsn+Ncsp;
ind_ce = Ncsn+Ncsp+1:Nc;

ind_phi_s_n = 1:Nn;
ind_phi_s_p = Nn+1:Nnp;

ind_ien = Nnp+1:Nnp+Nn;
ind_iep = Nnp+Nn+1:2*Nnp;

ind_phi_e = 2*Nnp+1 : 2*Nnp+Nx;

ind_jn = 2*Nnp+Nx+1 : 2*Nnp+Nx+Nn;
ind_jp = 2*Nnp+Nx+Nn+1 : Nz;

%% Current Limits
% Minimum Current
D_Imin = -1;
E_Imin = Imin;

% Maximum Current
D_Imax = 1;
E_Imax = -Imax;

% Compile D matrix
D(indy_I,:) = [D_Imin; D_Imax];

%% Depletion/Saturation of Li concentration in solid
Cell_Ccsn = cell(Nn,1);
Cell_Ccsp = cell(Np,1);

% Loop through each "comb tooth" in anode
for idx = 1:Nn
    Cell_Ccsn{idx} = p.C_csn(1,:);
end
C1_cssn = blkdiag(Cell_Ccsn{:});

% Loop through each "comb tooth" in anode
for idx = 1:Np
    Cell_Ccsp{idx} = p.C_csp(1,:);
end
C1_cssp = blkdiag(Cell_Ccsp{:});

% Neg electrode: Min c_ss 
C1_cssnmin = -C1_cssn / p.c_s_n_max;
E_cssnmin = theta_n_min * ones(Nn,1);

% Neg electrode: Max c_ss
C1_cssnmax = C1_cssn / p.c_s_n_max;
E_cssnmax = -theta_n_max * ones(Nn,1);

% Pos electrode: Min c_ss
C1_csspmin = -C1_cssp / p.c_s_p_max;
E_csspmin = theta_p_min * ones(Np,1);

% Pos electrode: Max c_ss
C1_csspmax = C1_cssp / p.c_s_p_max;
E_csspmax = -theta_p_max * ones(Np,1);

% Compile C1 matrix
C1(indy_cssn,ind_csn) = [C1_cssnmin; C1_cssnmax];
C1(indy_cssp,ind_csp) = [C1_csspmin; C1_csspmax];

%% Stress from spatial Li concentration gradient in solid

%% Depletion/Saturation of Li concentration in solid
C1_ce1 = p.C_ce(1,:);
C1_ce2 = [speye(Nn), zeros(Nn, Nx-Nn)];
C1_ce3 = p.C_ce(2,:);
C1_ce4 = [zeros(Ns,Nn), speye(Ns), zeros(Ns,Np)];
C1_ce5 = p.C_ce(3,:);
C1_ce6 = [zeros(Np, Nx-Np), speye(Np)];
C1_ce7 = p.C_ce(4,:);

C1_ce = [C1_ce1; C1_ce2; C1_ce3; C1_ce4; C1_ce5; C1_ce6; C1_ce7];

% Minimum ce
C1_cemin = -C1_ce;
E_cemin = c_e_min * ones(Nx+4,1);

% Maximum ce
C1_cemax = C1_ce;
E_cemax = -c_e_max * ones(Nx+4,1);

% Compile C1 matrix
C1(indy_ce,ind_ce) = [C1_cemin; C1_cemax];

%% Temperature Limits
% Minimum Current
C1_Tmin = -1;
E_Tmin = T_min;

% Maximum Current
C1_Tmax = 1;
E_Tmax = -T_max;

% Compile C1 matrix
C1(indy_T,end) = [C1_Tmin; C1_Tmax];

%% Side rxn overpotential @ x = (L^-) - \Delta x 
C2_etas_psn = [zeros(1,Nn-1), -1];
C2_etas_pe = [zeros(1,Nn-1), 1, zeros(1,Nx-Nn)];
E_etas = Us;

% Compile C2 matrix
C2(indy_etas, ind_phi_s_n) = C2_etas_psn;
C2(indy_etas, ind_phi_e) = C2_etas_pe;

%% Compile Matrices
E = [E_Imin; E_Imax; E_cssnmin; E_cssnmax; E_csspmin; E_csspmax;...
     E_cemin; E_cemax; E_Tmin; E_Tmax; E_etas];
 
%% Create sparse matrices
C1 = sparse(C1);
C2 = sparse(C2);
D = sparse(D);
E = sparse(E);
