function n = symsubsnum (R,D)
%   Direct Encoding of pade_cs.num from pade_cs3.mat
%   Created Jan 10, 2014 by Hector Perez 

n = [-3/R ;
  -(4*R)/(11*D) ;
 -R^3/(165*D^2)];
