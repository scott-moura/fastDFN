function [F,s] = expm_new(A,schur_fact)
%EXPM_NEW  Matrix exponential.
%   EXPM_NEW(A) is the matrix exponential of A computed using
%   an improved scaling and squaring algorithm with a Pade approximation.
%   It exploits triangularity (if any) of A.
%   EXPM_NEW(A,1) uses an initial transformation to complex Schur form.
%   EXPM_NEW(A,2) uses an initial transformation to real Schur form if A
%   is real.

%   Reference: A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring
%   Algorithm for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 31(3):
%   970-989, 2009.

%   Awad H. Al-Mohy and Nicholas J. Higham, April 20, 2010.

upper_triang = isequal(A,triu(A));
lower_triang = isequal(A,tril(A));

if nargin < 2 || isempty(schur_fact) || upper_triang || lower_triang
   schur_fact = false;
end

if schur_fact == 1
    [Q,A] = schur(A,'complex');
elseif schur_fact == 2
    [Q,A] = schur(A);  % Complex Schur form will result if A is complex.
end

% Coefficients of leading terms in the backward error functions h_{2m+1}.
Coeff = [1/100800, 1/10059033600, 1/4487938430976000,...
         1/5914384781877411840000, 1/113250775606021113483283660800000000];

u = eps/2;
n = length(A);
have_A4 = 0;   % prior evaluation of A4
have_A6 = 0;   % prior evaluation of A6

                m_vals = [3 5 7 9 13];
                % theta_m for m=1:13.
                theta = [%3.650024139523051e-008
                         %5.317232856892575e-004
                          1.495585217958292e-002   % m_vals = 3
                         %8.536352760102745e-002
                          2.539398330063230e-001   % m_vals = 5
                         %5.414660951208968e-001
                          9.504178996162932e-001   % m_vals = 7
                         %1.473163964234804e+000
                          2.097847961257068e+000   % m_vals = 9
                         %2.811644121620263e+000
                         %3.602330066265032e+000
                         %4.458935413036850e+000
                          4.250000000000000e+000]; % m_vals = 13

s = 0;


A2 = A*A;
eta1 = max( normAm(A2,2)^(1/4), normAm(A2,3)^(1/6) );
t = eval_alpha(A,1);

    if eta1 <= theta(1) && t == 0
        F = PadeApproximantOfDegree(m_vals(1));
        schur_backtrans
        return
    end

A4 = A2*A2; have_A4 = 1;
eta2 = max( norm(A4,1)^(1/4), normAm(A2,3)^(1/6) );
t = eval_alpha(A,2);

    if eta2 <= theta(2) && t == 0
        F = PadeApproximantOfDegree(m_vals(2));
        schur_backtrans
        return
    end

A6 = A2*A4; have_A6 = 1;
eta3 = max(norm(A6,1)^(1/6), normAm(A4,2)^(1/8));
h = zeros(4,1);

h(3) = eval_alpha(A,3);
h(4) = eval_alpha(A,4);

for i = 3:4
    if eta3 <= theta(i) && h(i) == 0
        F = PadeApproximantOfDegree(m_vals(i));
        schur_backtrans
        return
    end
end

eta4 = max(normAm(A4,2)^(1/8), normAm(A2,5)^(1/10));
eta5 = min(eta3,eta4);
s = max( ceil(log2(eta5/theta(end))), 0 ); % Zero must be here
t = eval_alpha(A/2^s, 5);
s = s + t;
A = A/2^s;  A2 = A2/2^(2*s); A4 = A4/2^(4*s); A6 = A6/2^(6*s); % Scaling

F = PadeApproximantOfDegree(m_vals(end));
if lower_triang, A = A'; F = F'; end

if upper_triang || lower_triang || schur_fact
    F = expm_sqtri(A,F,s);
    if lower_triang, F = F'; end
else
    for k= 1:s
        F = F*F;
    end
end
schur_backtrans

%%%%Nested Functions%%%%

    function t = eval_alpha(A,k)
        alpha = Coeff(k)*normAm(abs(A),2*m_vals(k)+1)/norm(A,1);
        t = max( ceil(log2(alpha/u)/(2*m_vals(k))),0);
    end

    function F = PadeApproximantOfDegree(m)
        %PADEAPPROXIMANTOFDEGREE  Pade approximant to exponential.
        %   F = PADEAPPROXIMANTOFDEGREE(M) is the degree M diagonal
        %   Pade approximant to EXP(A), where M = 3, 5, 7, 9 or 13.
        %   Series are evaluated in decreasing order of powers, which is
        %   in approx. increasing order of maximum norms of the terms.

        c = getPadeCoefficients;

        % Evaluate Pade approximant.
        switch m

            case {3, 5, 7, 9}

                Apowers = cell(ceil((m+1)/2),1);
                Apowers{1} = eye(n);
                Apowers{2} = A2; start = 3;
                if have_A4
                    Apowers{3} = A4;  start = 4;
                end
                if have_A6
                    Apowers{4} = A6;  start = 5;
                end
                for j = start:ceil((m+1)/2)
                    Apowers{j} = Apowers{j-1}*Apowers{2};
                end

                U = zeros(n); V = zeros(n);

                for j = m+1:-2:2
                    U = U + c(j)*Apowers{j/2};
                end
                U = A*U;
                for j = m:-2:1
                    V = V + c(j)*Apowers{(j+1)/2};
                end

            case 13

                % For optimal evaluation need different formula for m >= 12.
                U = A * (A6*(c(14)*A6 + c(12)*A4 + c(10)*A2) ...
                         + c(8)*A6 + c(6)*A4 + c(4)*A2 + c(2)*eye(n) );

                V = A6*(c(13)*A6 + c(11)*A4 + c(9)*A2) ...
                    + c(7)*A6 + c(5)*A4 + c(3)*A2 + c(1)*eye(n);

        end
        F = (-U+V)\(U+V);
        function c = getPadeCoefficients
            % GETPADECOEFFICIENTS Coefficients of numerator P of Pade approximant
            %    C = GETPADECOEFFICIENTS returns coefficients of numerator
            %    of [M/M] Pade approximant, where M = 3,5,7,9,13.
            switch m
                case 3
                    c = [120, 60, 12, 1];
                case 5
                    c = [30240, 15120, 3360, 420, 30, 1];
                case 7
                    c = [17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1];
                case 9
                    c = [17643225600, 8821612800, 2075673600, 302702400, 30270240, ...
                         2162160, 110880, 3960, 90, 1];
                case 13
                    c = [64764752532480000, 32382376266240000, 7771770303897600, ...
                         1187353796428800,  129060195264000,   10559470521600, ...
                         670442572800,      33522128640,       1323241920,...
                         40840800,          960960,            16380,  182,  1];
            end
        end
    end

    function schur_backtrans
    if schur_fact, F = Q*F*Q'; end
    end

%%%%Nested Functions%%%%
end
