function [x,flag,relres,iter,resvec] = bicgFast(A,b,tol,maxit,M1,M2,x0,varargin)
%BICG   BiConjugate Gradients Method.
%   X = BICG(A,B) attempts to solve the system of linear equations A*X=B
%   for X.  The N-by-N coefficient matrix A must be square and the right
%   hand side column vector B must have length N.
%
%   X = BICG(AFUN,B) accepts a function handle AFUN instead of the matrix A.
%   AFUN(X,'notransp') accepts a vector input X and returns the
%   matrix-vector product A*X while AFUN(X,'transp') returns A'*X. In all
%   of the following syntaxes, you can replace A by AFUN.
%
%   X = BICG(A,B,TOL) specifies the tolerance of the method.  If TOL is []
%   then BICG uses the default, 1e-6.
%
%   X = BICG(A,B,TOL,MAXIT) specifies the maximum number of iterations.  If
%   MAXIT is [] then BICG uses the default, min(N,20).
%
%   X = BICG(A,B,TOL,MAXIT,M) and X = BICG(A,B,TOL,MAXIT,M1,M2) use the
%   preconditioner M or M=M1*M2 and effectively solve the system
%   inv(M)*A*X = inv(M)*B for X. If M is [] then a preconditioner is not
%   applied. M may be a function handle MFUN such that MFUN(X,'notransp')
%   returns M\X and MFUN(X,'transp') returns M'\X.
%
%   X = BICG(A,B,TOL,MAXIT,M1,M2,X0) specifies the initial guess. If X0 is
%   [] then BICG uses the default, an all zero vector.
%
%   [X,FLAG] = BICG(A,B,...) also returns a convergence FLAG:
%    0 BICG converged to the desired tolerance TOL within MAXIT iterations
%    1 BICG iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 BICG stagnated (two consecutive iterates were the same).
%    4 one of the scalar quantities calculated during BICG became
%      too small or too large to continue computing.
%
%   [X,FLAG,RELRES] = BICG(A,B,...) also returns the relative residual
%   NORM(B-A*X)/NORM(B). If FLAG is 0, then RELRES <= TOL.
%
%   [X,FLAG,RELRES,ITER] = BICG(A,B,...) also returns the iteration number
%   at which X was computed: 0 <= ITER <= MAXIT.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = BICG(A,B,...) also returns a vector of
%   the residual norms at each iteration including NORM(B-A*X0).
%
%   Example:
%      n = 100; on = ones(n,1); A = spdiags([-2*on 4*on -on],-1:1,n,n);
%      b = sum(A,2); tol = 1e-8; maxit = 15;
%      M1 = spdiags([on/(-2) on],-1:0,n,n);
%      M2 = spdiags([4*on -on],0:1,n,n);
%      x = bicg(A,b,tol,maxit,M1,M2);
%   Or, use this matrix-vector product function
%      %---------------------------------------------%
%      function y = afun(x,n,transp_flag)
%      if strcmp(transp_flag,'transp')      % y = A'*x
%         y = 4 * x;
%         y(1:n-1) = y(1:n-1) - 2 * x(2:n);
%         y(2:n) = y(2:n) - x(1:n-1);
%      elseif strcmp(transp_flag,'notransp') % y = A*x
%         y = 4 * x;
%         y(2:n) = y(2:n) - 2 * x(1:n-1);
%         y(1:n-1) = y(1:n-1) - x(2:n);
%      end
%      %---------------------------------------------%
%   as input to BICG:
%      x1 = bicg(@(x,tflag)afun(x,n,tflag),b,tol,maxit,M1,M2);
%
%   Class support for inputs A,B,M1,M2,X0 and the output of AFUN:
%      float: double
%
%   See also BICGSTAB, BICGSTABL, CGS, GMRES, LSQR, MINRES, PCG, QMR,
%   SYMMLQ, TFQMR, ILU, FUNCTION_HANDLE.

%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 1.19.4.9 $ $Date: 2010/11/17 11:29:40 $

% Check for an acceptable number of input arguments
if nargin < 2
    error(message('MATLAB:bicg:NotEnoughInputs'));
end

% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = iterchk(A);
if strcmp(atype,'matrix')
    % Check matrix and right hand side vector inputs have appropriate sizes
    [m,n] = size(A);
    if (m ~= n)
        error(message('MATLAB:bicg:NonSquareMatrix'));
    end
    if ~isequal(size(b),[m,1])
        error(message('MATLAB:bicg:RSHsizeMismatchCoeffMatrix', m));
    end
else
    m = size(b,1);
    n = m;
    if ~iscolumn(b)
        error(message('MATLAB:bicg:RSHnotColumn'));
    end
end

% Assign default values to unspecified parameters
if nargin < 3 || isempty(tol)
    tol = 1e-6;
end
warned = 0;
if tol < eps
    warning('MATLAB:bicg:tooSmallTolerance', ...
       'Input tol is smaller than eps and may not be achieved by BICG\n         Try to use a bigger tolerance');
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning(message('MATLAB:bicg:tooBigTolerance'));
    warned = 1;
    tol = 1-eps;
end
if nargin < 4 || isempty(maxit)
    maxit = min(n,20);
end

% Check for all zero right hand side vector => all zero solution
n2b = norm(b);                     % Norm of rhs vector, b
if (n2b == 0)                      % if    rhs vector is all zeros
    x = zeros(n,1);                % then  solution is all zeros
    flag = 0;                      % a valid solution has been obtained
    relres = 0;                    % the relative residual is actually 0/0
    iter = 0;                      % no iterations need be performed
    resvec = 0;                    % resvec(1) = norm(b-A*x) = norm(0)
    if (nargout < 2)
        itermsg('bicg',tol,maxit,0,flag,iter,NaN);
    end
    return
end

if ((nargin >= 5) && ~isempty(M1))
    existM1 = 1;
    [m1type,m1fun,m1fcnstr] = iterchk(M1);
    if strcmp(m1type,'matrix')
        if ~isequal(size(M1),[m,m])
            error(message('MATLAB:bicg:WrongPrecondSize', m));
        end
    end
else
    existM1 = 0;
    m1type = 'matrix';
end

if ((nargin >= 6) && ~isempty(M2))
    existM2 = 1;
    [m2type,m2fun,m2fcnstr] = iterchk(M2);
    if strcmp(m2type,'matrix')
        if ~isequal(size(M2),[m,m])
            error(message('MATLAB:bicg:WrongPrecondSize', m));
        end
    end
else
    existM2 = 0;
    m2type = 'matrix';
end

if ((nargin >= 7) && ~isempty(x0))
    if ~isequal(size(x0),[n,1])
        error(message('MATLAB:bicg:WrongInitGuessSize', n));
    else
        x = x0;
    end
else
    x = zeros(n,1);
end

if ((nargin > 7) && strcmp(atype,'matrix') && ...
        strcmp(m1type,'matrix') && strcmp(m2type,'matrix'))
    error(message('MATLAB:bicg:TooManyInputs'));
end

% Set up for the method
flag = 1;
xmin = x;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
tolb = tol * n2b;                  % Relative tolerance
r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:},'notransp');
normr = norm(r);                   % Norm of residual
normr_act = normr;

if (normr <= tolb)                 % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2b;
    iter = 0;
    resvec = normr;
    if (nargout < 2)
        itermsg('bicg',tol,maxit,0,flag,iter,relres);
    end
    return
end

rt = r;                            % Shadow residual
resvec = zeros(maxit+1,1);         % Preallocate vector for norms of residuals
resvec(1) = normr;                 % resvec(1) = norm(b-A*x0)
normrmin = normr;                  % Norm of residual from xmin
rho = 1;
stag = 0;                          % stagnation of the method
moresteps = 0;
maxmsteps = min([floor(n/50),5,n-maxit]);
maxstagsteps = 3;

% loop over maxit iterations (unless convergence or failure)

for ii = 1 : maxit
    if existM1
        y = iterapp('mldivide',m1fun,m1type,m1fcnstr,r,varargin{:},'notransp');
        if any(~isfinite(y))
            flag = 2;
            break
        end
    else % no preconditioner
        y = r;
    end
    if existM2
        z = iterapp('mldivide',m2fun,m2type,m2fcnstr,y,varargin{:},'notransp');
        yt = iterapp('mldivide',m2fun,m2type,m2fcnstr,rt,varargin{:},'transp');
        if any(~isfinite(z)) || any(~isfinite(yt))
            flag = 2;
            break
        end
    else
        z = y;
        yt = rt;
    end
    if existM1
        zt = iterapp('mldivide',m1fun,m1type,m1fcnstr,yt,varargin{:},'transp');
        if any(~isfinite(zt))
            flag = 2;
            break
        end
    else
        zt = yt;
    end
    
    rho1 = rho;
    rho = rt' * z;
    if rho == 0 || isinf(rho)
        flag = 4;
        break
    end
    if ii == 1
        p = z;
        pt = zt;
    else
        beta = rho / rho1;
        if beta == 0 || isinf(beta)
            flag = 4;
            break
        end
        p = z + beta * p;
        pt = zt + conj(beta) * pt;
    end
    
    q = iterapp('mtimes',afun,atype,afcnstr,p,varargin{:},'notransp');
    if strcmp(atype,'matrix')
        qt = A' * pt;
    else
        qt = iterapp('mtimes',afun,atype,afcnstr,pt,varargin{:},'transp');
    end
    ptq = pt' * q;
    if ptq == 0
        flag = 4;
        break
    else
        alpha = rho / ptq;
    end
    if isinf(alpha)
        flag = 4;
        break
    end
    
    % Check for stagnation of the method
    if abs(alpha)*norm(p) < eps*norm(x)
        stag = stag + 1;
    else
        stag = 0;
    end
    
    x = x + alpha * p;             % form the new iterate
    r = r - alpha * q;
    rt = rt - conj(alpha) * qt;
    normr = norm(r);
    normr_act = normr;
    resvec(ii+1) = normr;
    
    % check for convergence
    if (normr <= tolb || stag >= maxstagsteps || moresteps)
        r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:},'notransp');
        normr_act = norm(r);
        resvec(ii+1,1) = normr_act;
        if (normr_act <= tolb)
            flag = 0;
            iter = ii;
            break
        else
            if stag >= maxstagsteps && moresteps == 0
                stag = 0;
            end
            moresteps = moresteps + 1;
            if moresteps >= maxmsteps
                if ~warned
                    warning('MATLAB:bicg:tooSmallTolerance', ...
                        'Input tol may be smaller than eps*cond(A) and might not be achieved by BICG\n         Try to use a bigger tolerance');
                end
                flag = 3;
                iter = ii;
                break;
            end
        end
    end
    
    if normr_act < normrmin        % update minimal norm quantities
        normrmin = normr_act;
        xmin = x;
        imin = ii;
    end
    
    if stag >= maxstagsteps
        flag = 3;
        break
    end
    
end                                % for ii = 1 : maxit

% returned solution is first with minimal residual
if flag == 0
    relres = normr_act / n2b;
else
    r = b - iterapp('mtimes',afun,atype,afcnstr,xmin,varargin{:},'notransp');
    if flag == 1
        r_comp = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:},'notransp');
        normr_act = norm(r_comp);
    end
    if norm(r) <= normr_act
        x = xmin;
        iter = imin;
        relres = norm(r) / n2b;
    else
        iter = ii;
        relres = normr_act / n2b;
    end
end

% truncate the zeros from resvec
if flag <= 1 || flag == 3
    resvec = resvec(1:ii+1);
else
    resvec = resvec(1:ii);
end

% only display a message if the output flag is not used
if nargout < 2
    itermsg('bicg',tol,maxit,ii,flag,iter,relres);
end
