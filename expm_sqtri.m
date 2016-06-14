function X = expm_sqtri(T,F,s)
% EXPM_SQTRI   Squaring phase of scaling and squaring method.
%   X = EXPM_SQTRI(T/2^s,F,s) carries out the squaring phase
%   of the scaling and squaring method for an upper quasitriangular T,
%   given T/2^s and a Pade approximant F to e^{T/2^s}.
%   It corrects the diagonal blocks blocks at each step.

%   This M-file exploits Code Fragment 2.1 and Code Fragment 2.2 of the
%   reference below.

%   Reference: A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring
%   Algorithm for the Matrix Exponential,SIAM J. Matrix Anal. Appl. 31(3):
%   970-989, 2009.

%   Awad H. Al-Mohy and Nicholas J. Higham, April 19, 2010.

n = length(T);
k = 1;
% To turn off exact superdiagonal computation force "istriangular = 0".
istriangular = isequal(T,triu(T));

if n > 1
   c = abs(diag(T,-1)) > 0;    % sum(c) = number of 2by2 full blocks
   % NUMBLK blocks with i'th block in rows/cols INDX{i}.
   numblk = n - sum(c);         % The number of blocks
   indx = cell(numblk,1);
   if c(end) == 0
       indx{end} = n; c = [c ; 0];
   end
   for j = 1:numblk
       if c(k)
           indx{j} = k:k+1; k = k+2;
       else
           indx{j} = k; k = k+1;
       end
   end
end

for i = 0:s
    if i > 0, F = F*F; end
    if istriangular
       % Compute diagonal and first superdiagonal.
       for j = 1:2:n
           if j < n
              F(j:j+1,j:j+1) = expmT2by2( 2^i * T(j:j+1,j:j+1) );
           else
              F(n,n) = exp(2^i * T(n,n));
           end
       end
    else
       % Quasitriangular case: compute (block) diagonal only.
       for j = 1:numblk
           F(indx{j},indx{j}) = expm2_by_2( 2^i * T(indx{j},indx{j}) );
       end
    end
end

X = F;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = expm2_by_2(A)
% EXPM2_BY_2  Exponential for a general 2-by-2 matrix A.

if length(A) == 1
    X = exp(A);
else
    X = expm2by2full(A);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = expm2by2full(A)

% EXPM2BY2FULL   Exponential of 2-by-2 full matrix.

a = A(1,1);
b = A(1,2);
c = A(2,1);
d = A(2,2);

delta = sqrt((a-d)^2 + 4*b*c);

X = exp((a+d)/2)  * ...
      [ cosh(delta/2) + (a-d)/2*sinch(delta/2),  b*sinch(delta/2)
        c*sinch(delta/2),  cosh(delta/2) + (d-a)/2*sinch(delta/2) ];

%%%%%%%%%%%%%%%%%%%%%
function y = sinch(x)
    if x == 0
       y = 1;
    else
        y = sinh(x)/x;
    end

%%%%%%%%%%%%%%%%%%%%%%%%
function X = expmT2by2(A)
%EXPMT2BY2    Exponential of 2-by-2 upper triangular matrix.
%   EXPMT2BY2(A) is the exponential of the 2-by-2 upper triangular matrix A.

% Modified from FUNM (EXPM2by2).

a1 = A(1,1);
a2 = A(2,2);

ave = (a1+a2)/2; df  = abs(a1-a2)/2;

if max(ave,df) < log(realmax)
   % Formula fine unless it overflows.
   x12 = A(1,2)*exp( (a1+a2)/2 ) * sinch( (a2-a1)/2 );
else
   % Formula can suffer cancellation.
   x12 = A(1,2)*(exp(a2)-exp(a1))/(a2-a1);
end

X = [exp(a1)  x12
       0      exp(a2)];
