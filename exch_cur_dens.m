%% Exchange Current Density function i_0
%   Created July 12, 2011 by Scott Moura

function [i_0n,i_0p] = exch_cur_dens(p,c_ss_n,c_ss_p,c_e)

% Parse out concentrations in anode and cathode
c_e_n = c_e(1:p.Nxn-1);
c_e_p = c_e(p.Nxn+p.Nxs-1:end);

% Compute exchange current density
i_0n = p.k_n * ((p.c_s_n_max - c_ss_n) .* c_ss_n .* c_e_n).^p.alph;
i_0p = p.k_p * ((p.c_s_p_max - c_ss_p) .* c_ss_p .* c_e_p).^p.alph;