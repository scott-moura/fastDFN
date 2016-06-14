%% Symbolic computation of Matrices for Li Diffusion in Solid Phase, c_s(r,t)
% Created by Federico B. Apr. 29, 2014


% Load Pade Approximations
N = 3;
% load(['ws/pade_cs' num2str(N) '.mat']);

syms R_s_n D_s_n R_s_p D_s_p;
p.R_s_n=R_s_n;
p.D_s_n=D_s_n;
p.R_s_p=R_s_p;
p.D_s_p=D_s_p;

% Replace symbolic expressions with numerical values
b_vecn = symsubsnum(p.R_s_n, p.D_s_n); %New Jan 10, 2014
a_vecn = symsubsden(p.R_s_n, p.D_s_n); %New Jan 10, 2014

b_vecp = symsubsnum(p.R_s_p, p.D_s_p); %New Jan 10, 2014
a_vecp = symsubsden(p.R_s_p, p.D_s_p); %New Jan 10, 2014

% b_vecn = subs(pade_cs.num(:),{'R';'D'},[p.R_s_n; p.D_s_n]); %Changed Jan 10, 2014, was pade_cs.num
% a_vecn = subs(pade_cs.den(:),{'R';'D'},[p.R_s_n; p.D_s_n]); %Changed Jan 10, 2014, was pade_cs.den

% b_vecp = subs(pade_cs.num(:),{'R';'D'},[p.R_s_p; p.D_s_p]); %Changed Jan 10, 2014, was pade_cs.num
% a_vecp = subs(pade_cs.den(:),{'R';'D'},[p.R_s_p; p.D_s_p]); %Changed Jan 10, 2014, was pade_cs.den

bp_vecn = b_vecn / a_vecn(end);
ap_vecn = a_vecn / a_vecn(end);

bp_vecp = b_vecp / a_vecp(end);
ap_vecp = a_vecp / a_vecp(end);

% Create state-space representation
An = sym(diag(ones(N-1,1),1));
Ap = sym(diag(ones(N-1,1),1));
for idx = 1:N
    An(N,idx) = -ap_vecn(idx);
    Ap(N,idx) = -ap_vecp(idx);
end

Bn = zeros(N,1);
Bp = zeros(N,1);
Bn(N) = 1;
Bp(N) = 1;

C1n = bp_vecn';
C2n = b_vecn(1) * ap_vecn(2:end)';
Cn = [C1n; C2n];

C1p = bp_vecp';
C2p = b_vecp(1) * ap_vecp(2:end)';
Cp = [C1p; C2p];

%% Convert to Jordan-Form   NEW Apr 24, 2014
% Similarity Matrix, computed symbolically
Vn = [-(p.R_s_n^4*(9*sqrt(2429)-457))/(381150*p.D_s_n^2), (p.R_s_n^4*(9*sqrt(2429)+457))/(381150*p.D_s_n^2), 1;...
       p.R_s_n^2*(sqrt(2429)-63)/(2310*p.D_s_n), -p.R_s_n^2*(sqrt(2429)+63)/(2310*p.D_s_n), 0;...
       1, 1, 0];
Vp = [-(p.R_s_p^4*(9*sqrt(2429)-457))/(381150*p.D_s_p^2), (p.R_s_p^4*(9*sqrt(2429)+457))/(381150*p.D_s_p^2), 1;...
       p.R_s_p^2*(sqrt(2429)-63)/(2310*p.D_s_p), -p.R_s_p^2*(sqrt(2429)+63)/(2310*p.D_s_p), 0;...
       1, 1, 0];

An1 = Vn\An*Vn;
Ap1 = Vp\Ap*Vp;

Bn1 = Vn\Bn;
Bp1 = Vp\Bp;

Cn1 = Cn*Vn;
Cp1 = Cp*Vp;

% Perform additional transformation that scales third state, such that it's
% exactly \bar{c}_s^\pm
Vn2 = diag([1 1 1/Cn1(2,3)]);
Vp2 = diag([1 1 1/Cp1(2,3)]);

An2 = Vn2\An1*Vn2;
Ap2 = Vp2\Ap1*Vp2;

Bn2 = Vn2\Bn1;
Bp2 = Vp2\Bp1;

Cn2 = Cn1*Vn2;
Cp2 = Cp1*Vp2;

An = An2;
An_normalized=An/p.D_s_n; %added by Federico for sensitivity analysis (equivalent to having D_s_n=1)
Ap = Ap2;
Ap_normalized=Ap/p.D_s_p; %added by Federico for sensitivity analysis (equivalent to having D_s_p=1)
Bn = Bn2;
Bp = Bp2;
Cn = Cn2;
Cp = Cp2;
   
