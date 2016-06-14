%% Solve algebraic eqns for algebraic vars, using Newton's Method
%   Created May 22, 2012 by Scott Moura

function out = newtonAlgVars(x,z0,Cur,p);

% Discretization parameters
Nn = p.Nxn - 1;
Np = p.Nxp - 1;
Nnp = Nn+Np;
Nx = p.Nx - 3;
Nz = 3*Nnp + Nx;

% Parse out algebraic vars
phi_s_n = z(1:Nn);
phi_s_p = z(Nn+1 : Nnp);
i_en = z(Nnp+1 : Nnp+Nn);
i_ep = z(Nnp+Nn+1 : 2*Nnp);
phi_e = z(2*Nnp+1 : 2*Nnp+Nx);
jn = z(2*Nnp+Nx+1 : 2*Nnp+Nx+Nn);
jp = z(2*Nnp+Nx+Nn+1 : end);

%%% ALGEBRAIC EQUATIONS
% Solid Potential: phi_s(x,t)
[F1_psn,F1_psp,F2_psn,F2_psp,G_psn,G_psp] = phi_s_mats(p);

% Add BCs to i_0
i_enn = [0; i_en; Cur];
i_epp = [Cur; i_ep; 0];

g_phi_s_n = F1_psn*phi_s_n + F2_psn*i_enn + G_psn*Cur;
g_phi_s_p = F1_psp*phi_s_p + F2_psp*i_epp + G_psp*Cur;

% Electrolyte Current: i_e(x,t)
[F1_ien,F1_iep,F2_ien,F2_iep,F3_ien,F3_iep] = i_e_mats(p);

g_i_en = F1_ien*i_en + F2_ien*jn + F3_ien*Cur;
g_i_ep = F1_iep*i_ep + F2_iep*jp + F3_iep*Cur;

% Potential in Electrolyte Phase: phi_e(x,t)
% System matrices
[F1_pe,F2_pe,F3_pe] = phi_e_mats_new(p,c_ex);

% Algebraic eqns (semi-explicit form)
phi_e_dot = F1_pe*phi_e + F2_pe*i_ex + F3_pe*log(c_ex);

% Aggregate g's
g = [g_phi_s_n; g_phi_s_p; g_i_en; g_i_ep; g_phi_e; g_jn; g_jp];

% Jacobian
[~,Jac_z] = jac_dfn(t,x,p);

% Newton iteration
Delta_z = -Jac_z\g;
z_nxt = z + Delta_z;