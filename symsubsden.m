function n = symsubsden (R,D)
%   Direct Encoding of pade_cs.den from pade_cs3.mat
%   Created Jan 10, 2014 by Hector Perez 

n = [0 ;
     1 ;
 (3*R^2)/(55*D);
 R^4/(3465*D^2)];