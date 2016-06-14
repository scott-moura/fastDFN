%% Pade Solver for Diffusion Eqn
%   Created on May 21, 2012 by Scott Moura
clear

% Pade Order
N = 2;

%% Take Maclaurin Series of original trancendental TF
load('ws/maclaurin_cs_tf.mat')

if(length(c) < (2*(N+1)+1))

    syms R D real positive
    syms s
    
    % This is actually G_star, i.e. G_star = s*G(s)
    G = s * R/D * tanh(R*sqrt(s/D)) / (tanh(R*sqrt(s/D)) - R*sqrt(s/D));

    % Take Maclaurin Series
    c = sym(zeros(2*(N+1)+1,1));
    for k = 0:(2*(N+1))

        % Take derivatives
        dG = simplify(diff(G,'s',k));
        
        % Take limit at s = 0
        c(k+1) = limit(dG,'s',0) / factorial(k);
        pretty(c(k+1));

    end
    save('ws/maclaurin_cs_tf.mat','c');
    
end

% % Numerator and Denominator
% num = sym(zeros(N+1,1));
% den = sym(zeros(N+1,1));
% for k = 0:N
%    
%     num(k+1) = sym(['b' num2str(k)]);
%     if(k == 0)
%         den(1) = 1;
%     else
%         den(k+1) = sym(['a' num2str(k)]);
%     end
%     
% end
% 
% % Construct polynomial equation
% R = sum(den .* s_vec) * (s.^(0:(2*(N+1))) * c) - sum(num.*s_vec);
% RR = collect(R,'s');

% Construct Matrices
A11 = sym(eye(N+1));

A12 = sym(zeros(N+1,N));
for ii = 2:N+1
    for jj = 1:N
        if(ii > jj)
            A12(ii,jj) = -c(ii-1 - (jj-1));
        end
    end
end

A21 = sym(zeros(N,N+1));

A22 = sym(zeros(N));
for ii = 1:N
    for jj = 1:N
        
        A22(ii,jj) = -c(N +(ii-1) - (jj-1) + 1);
        
    end
end

A = [A11, A12; A21, A22];
B = c(1:(2*N+1));
pade_coeffs = A\B;

for k = 0:N
    
    fprintf(1,'b%1.0f = ',k);
    disp(pade_coeffs(k+1))
end

for k = 0:N
    fprintf(1,'a%1.0f = ',k);
    if(k == 0)
        disp('1');
    else
        disp(pade_coeffs(N+k+1));
    end

end

% Form Transfer Function object
num = pade_coeffs(1:N+1);
den = [0; 1; pade_coeffs(N+2:end)];


pade_cs.order = N+1;
pade_cs.num = num;
pade_cs.den = den;
save(['ws/pade_cs' num2str(N+1) '.mat'],'pade_cs');

