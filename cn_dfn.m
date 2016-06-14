%% Crank-Nicolson Eqns for Doyle-Fuller-Newman Model
%   Created May 22, 2012 by Scott Moura

function [x_nxtf, z_nxtf, varargout] = cn_dfn(x,z,Cur_vec,p)

Cur_prv = Cur_vec(1);
Cur = Cur_vec(2);
Cur_nxt = Cur_vec(3);

%% Parameters

% Newton params
maxIters = 100;
tol = 1e-5;

% % Sparse Linear System Solver Params
% linEqnTol = 1e-6; % Default is 1e-6
% linEqnMaxIter = 100;  % Default is 20

% tol_cs = 1e-4 * 1e-3 * 0.01;
% tol_ce = 1e-3 * 0.01;
% tol_ps = 0.01;
% tol_ie = 0.01;
% tol_pe = 0.01;
% tol_jn = 1e-6;

Nx = length(x);
Nz = length(z);

% Preallocate
x_nxt = zeros(Nx,maxIters);
z_nxt = zeros(Nz,maxIters);

relres = zeros(maxIters,1);
relres(1) = 1;

%% Solve for consistent ICs at k
if(Cur ~= Cur_prv)
    
    disp('Solving for consistent ICs');
    
    % Preallocate
    z_cons = zeros(Nz,maxIters);
    z_cons(:,1) = z;
    
    for idx = 1:(maxIters-1)
        
        % DAE eqns for current time-step
        [~,g] = dae_dfn(x,z_cons(:,idx),Cur,p);
        
        % Jacobian of DAE eqns
        [~,~,~,g_z] = jac_dfn(x,z_cons(:,idx),Cur,p.f_x,p.f_z,p.g_x,p.g_z,p);

        % Newton Iteration
        Delta_z = -(g_z\g);
%         [LL,UU] = ilu(g_z,struct('type','ilutp','droptol',1e-6,'udiag',1));
%         [Delta_z,flag] = bicg(g_z,-g,linEqnTol,linEqnMaxIter,LL,UU);
        z_cons(:,idx+1) = z_cons(:,idx) + Delta_z;
        
        % Check stopping criterion
        relres_z = norm(Delta_z,inf) / norm(z,inf);
        if(relres_z < tol)
            break;
        elseif(idx == (maxIters-1))
            fprintf(1,'Warning: Max Newton Iters Reached | RelChange = %3.2f%%\n',relres_z*100);
        end
        
    end
    
    z = z_cons(:,idx+1);
    
end

%% Solve Nonlinear System using Newton's Method
% DAE eqns for current time-step
[f,~] = dae_dfn(x,z,Cur,p);

% Initialize next x,z
x_nxt(:,1) = x;
z_nxt(:,1) = z;

% Iterate Newton's Method
FAC = [];
for idx = 1:(maxIters-1)

    % DAE eqns for next time-step
    [f_nxt, g_nxt] = dae_dfn(x_nxt(:,idx), z_nxt(:,idx), Cur_nxt, p);

    % Nonlinear system of equations
    F1 = x - x_nxt(:,idx) + p.delta_t/2 * (f + f_nxt);
    F2 = g_nxt;
    F = [F1; F2];

    % Jacobian of DAE eqns
    [f_x, f_z, g_x, g_z] = jac_dfn(x_nxt(:,idx),z_nxt(:,idx),Cur_nxt,p.f_x,p.f_z,p.g_x,p.g_z, p);
    
    % Jacobian of implicit time-stepped eqns
    F1_x = -speye(length(x)) + p.delta_t/2 * f_x;
    F1_z = p.delta_t/2 * f_z;
    F2_x = g_x;
    F2_z = g_z;
    
%     D = sparse(diag(max(abs(g_z),[],2)));
%     
%     [Delta_x,flx,rrx,itx,rrvx] = bicg(F1_x,-F1 - F1_z*Delta_z,linEqnTol,linEqnMaxIter);
%     [Delta_z,flz,rrz,itz,rrvz] = bicg(F2_z,-F2 - F2_x*Delta_x, linEqnTol,linEqnMaxIter);
%     
%     if(flx ~= 0)
%         disp('X: linear solver did not converge')
%     end
%     if(flz ~= 0)
%         disp('Z: linear solver did not converge')
%     end
%     
    J = [F1_x, F1_z; F2_x, F2_z];
    
    % Check Jacobian against numjac
%     Janalytic = [f_x, f_z; g_x, g_z];
%     
%     y = [x_nxt(:,idx); z_nxt(:,idx)];
%     FTY = dae_dfn_numjac(0,y,Cur_nxt,p);
%     thresh = eps * ones(length(y),1);
%     
%     [Jnum,FAC] = numjac(@(t,y) dae_dfn_numjac(t,y,Cur_nxt,p), 0, y, FTY, thresh, FAC, 0);
%     
%     Jerr = abs(Janalytic - Jnum);
% 
%     Jerr(Jerr < 1e-4) = 0;
%     spy(Jerr)
%     assignin('base','Janalytic',Janalytic);
%     assignin('base','Jnum',Jnum);
%     assignin('base','Jerr',Jerr);
%     
%     pause;
    
    % 1/ maximum across rows in Jacobian
%     D = sparse(diag(1./max(abs(J),[],2)));
%     DJ = D*J;
%     DF = D*F;

    % Newton Iteration
%     [LL,UU] = ilu(DJ,struct('type','ilutp','droptol',1e-6,'udiag',1));
    
    Delta_y = -(J\F); 
    
%     [Delta_y,fl,rr,it,rrv] = bicg(DJ,-DF,linEqnTol,linEqnMaxIter,LL,UU);
%     [Delta_y,fl,rr,it,rrv] = gmres(DJ,-DF,[],1e-12,linEqnMaxIter,LL,UU);
%     [Delta_y,fl(1),rr(1),it(1),rrv{1}] = bicg(J,-F,linEqnTol,linEqnMaxIter,LL,UU);
%     [~,fl(2),rr(2),it(2),rrv{2}] = bicgstab(J,-F,linEqnTol,linEqnMaxIter,LL,UU);
%     [~,fl(3),rr(3),it(3),rrv{3}] = bicgstabl(J,-F,linEqnTol,linEqnMaxIter,LL,UU);
%     [~,fl(4),rr(4),it(4),rrv{4}] = cgs(J,-F,linEqnTol,linEqnMaxIter,LL,UU);
%     [~,fl(5),rr(5),it(5),rrv{5}] = qmr(J,-F,linEqnTol,linEqnMaxIter,LL,UU);
%     [~,fl(6),rr(6),it(6),rrv{6}] = tfqmr(J,-F,linEqnTol,linEqnMaxIter,LL,UU);

%     fl
%     log10(rr)
%     rrv
%     
%     figure(10)
%     clf
%     semilogy(1:length(rrv{1}),rrv{1}./norm(F),'g.');
%     hold on;
%     semilogy(1:length(rrv{2}),rrv{2}./norm(F),'m.');
%     semilogy(1:length(rrv{3}),rrv{3}./norm(F),'c.');
%     semilogy(1:length(rrv{4}),rrv{4}./norm(F),'r.');
%     semilogy(1:length(rrv{5}),rrv{5}./norm(F),'b.');
%     semilogy(1:length(rrv{6}),rrv{6}./norm(F),'k.');
%     
%     legend('bicg','bicgstab','bicgstabl','cgs','qmr','tfqmr')
%     pause;

    x_nxt(:,idx+1) = x_nxt(:,idx) + Delta_y(1:Nx);
    z_nxt(:,idx+1) = z_nxt(:,idx) + Delta_y(Nx+1:end);
    
    % Check stopping criterion #1
    y = [x_nxt(:,idx+1); z_nxt(:,idx+1)];
    relres(idx+1) = norm(Delta_y,inf) / norm(y,inf);
    
    if( (relres(idx+1) < tol) && (norm(F,inf) < tol) )
        break;
    elseif(idx == (maxIters-1))
        fprintf(1,'Warning: Max Newton Iters Reached | RelChange = %1.4e%%\n',relres(end)*100);
    end
    
end

%% Output next x,z
x_nxtf = x_nxt(:,idx+1);
z_nxtf = z_nxt(:,idx+1);

newtonStats.iters = idx;
newtonStats.relres = relres;
% newtonStats.condJac = condest(J);
varargout{1} = newtonStats;

