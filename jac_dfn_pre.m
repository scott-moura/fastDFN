%% Jacobian for Doyle-Fuller-Newman Model
%   Created May 30, 2012 by Scott Moura
%   Preprocessed elements, i.e. non-state dependent

function [f_x, f_z, g_x, g_z] = jac_dfn_pre(p)


%% Parse out states
Ncsn = p.PadeOrder * (p.Nxn-1);
Ncsp = p.PadeOrder * (p.Nxp-1);
Nce = p.Nx - 3;
Nc = Ncsn+Ncsp+Nce;
Nn = p.Nxn - 1;
Np = p.Nxp - 1;
Nnp = Nn+Np;
Nx = p.Nx - 3;
Nz = 3*Nnp + Nx + 2;

ind_csn = 1:Ncsn;
ind_csp = Ncsn+1:Ncsn+Ncsp;

ind_cs = 1:Ncsn+Ncsp;
ind_ce = Ncsn+Ncsp+1:Nc;

ind_phi_s_n = 1:Nn;
ind_phi_s_p = Nn+1:Nnp;

ind_ien = Nnp+1:Nnp+Nn;
ind_iep = Nnp+Nn+1:2*Nnp;

ind_phi_e = 2*Nnp+1 : 2*Nnp+Nx+2;

ind_jn = 2*Nnp+Nx+3 : 2*Nnp+Nx+2+Nn;
ind_jp = 2*Nnp+Nx+3+Nn : Nz;

%% Preallocate Jacobian
f_x = zeros(Nc+1);
f_z = zeros(Nc+1,Nz);
g_x = zeros(Nz,Nc+1);
g_z = zeros(Nz);

%% Li Diffusion in Solid Phase: c_s(x,r,t)
% Preallocate
Cell_Acsn = cell(Nn,1);
Cell_Bcsn = cell(Nn,1);
Cell_Acsp = cell(Np,1);
Cell_Bcsp = cell(Np,1);

% Loop through each "comb tooth" in anode
for idx = 1:Nn
    Cell_Acsn{idx} = p.A_csn;
    Cell_Bcsn{idx} = p.B_csn;
end
f_x(ind_csn,ind_csn) = blkdiag(Cell_Acsn{:});
f_z(ind_csn,ind_jn) = blkdiag(Cell_Bcsn{:});

% Loop through each "comb tooth" in cathode
for idx = 1:Np
    Cell_Acsp{idx} = p.A_csp;
    Cell_Bcsp{idx} = p.B_csp;
end
f_x(ind_csp,ind_csp) = blkdiag(Cell_Acsp{:});
f_z(ind_csp,ind_jp) = blkdiag(Cell_Bcsp{:});

%% Potential in Solid Phase: phi_s(x,t)
g_z(ind_phi_s_n,ind_phi_s_n) = p.F1_psn;
g_z(ind_phi_s_p,ind_phi_s_p) = p.F1_psp;

g_z(ind_phi_s_n,ind_ien) = p.F2_psn(:,2:end-1);
g_z(ind_phi_s_p,ind_iep) = p.F2_psp(:,2:end-1);

%% Electrolyte Current: i_e(x,t)
g_z(ind_ien,ind_ien) = p.F1_ien;
g_z(ind_iep,ind_iep) = p.F1_iep;

g_z(ind_ien,ind_jn) = p.F2_ien;
g_z(ind_iep,ind_jp) = p.F2_iep;

%%
f_x = sparse(f_x);
f_z = sparse(f_z);
g_x = sparse(g_x);
g_z = sparse(g_z);

% spy(Jac)
% pause;
