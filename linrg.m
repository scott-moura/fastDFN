%% Explicit Computation of \beta for Linear DFN Reference Governor
%   Created September 11, 2012 by Scott Moura

function [beta] = linrg(p, x0,z0,I0, x,z,Ir)

%% Parameters
% Simulation Horizon
Ts = 5;

%% Jacobian matrices, evaluated at previous states values
[A11, A12, A21, A22, B1, B2] = jac_dfn_federico(x0,z0,I0, p.f_x,p.f_z,p.g_x,p.g_z, p);%changed to federico

% (A,B) matrices for linear system
A = A11 - A12*(A22\A21);
B = B1 - A12*(A22\B2);

%% Form "G,F" matrices (Discrete Time Method)
% Initial xtilde0
xtilde0 = x - x0;

% Convolution term
tau = 0:p.delta_t:Ts;
NT = length(tau);

% Block Matrices
MM = eye(size(A)) - p.delta_t/2 * A;
AA = eye(size(A)) + p.delta_t/2 * A;

AMinv = AA/MM;

AAA = MM \ (AMinv)^(NT-1) * AA;
BB = zeros([size(AMinv), NT]);
for ii = 1:NT
    BB(:,:,ii) = AMinv^(NT-ii);
end
BBB = MM \ (sum(BB,3)) * B;

% Constraint output equation matrices
[C1,C2,D,E] = cons_dfn(p);

F = (C1*BBB - C2*(A22\(A21*BBB + B2)) + D)*Ir;
G = -C1*(x0 + AAA*xtilde0 - BBB*I0) ...
    -C2*(z0 - A22 \ (A21*AAA*xtilde0 - (A21*BBB+B2)*I0)) - E;

%% Form "G,F" matrices (Exponential Matrix Method)
% % Initial xtilde0
% xtilde0 = x - x0;
% 
% % Convolution term
% tau = 0:p.delta_t:Ts;
% 
% Phi_k = zeros([size(A), length(tau)]);
% for k = 1:(length(tau))
%     Phi_k(:,:,k) = expm(A*(Ts - tau(k)));
% end
% L = sum(Phi_k(:,:,2:end) + Phi_k(:,:,1:end-1) ,3) * p.delta_t/2*B;
% 
% % State transition matrix
% Phi = Phi_k(:,:,1);
% 
% % Constraint output equation matrices
% [C1,C2,D,E] = cons_dfn(p);
% 
% G = -C1*(x0 + Phi*xtilde0 - L*I0) ...
%     - C2*(z0 - A22 \ ( A21*(Phi*xtilde0 - L*I0) - B2*I0 )) - E;
% F = (C1*L - C2*(A22\(A21*L+B2)) + D)*Ir;


%% Compute G by F
% xtildeTs = Phi*xtilde0 + L*(Ir-I0);
% xTs = x0 + xtildeTs;
% zTs = z0 - A22 \ (A21*xtildeTs + B2*(Ir-I0));
% yTs = C1*xTs + C2*zTs + D*Ir + E;

xtildeTs = AAA*xtilde0 + BBB*(Ir-I0);
xTs = x0 + xtildeTs;
zTs = z0 - A22 \ (A21*xtildeTs + B2*(Ir-I0));
yTs = C1*xTs + C2*zTs + D*Ir + E;

GbyF = G./F;
GbyF(F < 0) = -G(F<0)./F(F<0);

%% Solve optimization
beta = min([1; GbyF]);
fprintf(1,'beta = %1.2f\n',beta);

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

%% Indices for constraints
indy_I = 1:2;
indy_cssn = 3:2+(2*Nn);
indy_cssp = indy_cssn(end) + (1:(2*Np));
indy_ce = indy_cssp(end) + (1:(2*(p.Nx+1)));
indy_T = indy_ce(end) + (1:2);
indy_etas = indy_T(end) + 1;

% if(beta < 1)
%     ind = find(GbyF < 1)
%     [indy_ce(1) indy_ce(end)]
%     GbyF(ind)
%     x(Ncsn+Ncsp+1 : Ncsn+Ncsp+Nce)
%     yTs(indy_ce)
%     pause;
%     
% end

