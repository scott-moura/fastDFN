%% Save Doyle-Fuller-Newman Model Results
%   Created June 13, 2012 by Scott Moura

out.date = date;
out.time = t;
out.cur = I;
out.volt = Volt;
out.soc = SOC;
out.temp = T;
out.c_s_n = c_s_n;
out.c_s_p = c_s_p;
out.c_ss_n = c_ss_n;
out.c_ss_p = c_ss_p;
out.c_e = c_e;
% out.c_avg_n = c_avg_n;
% out.c_avg_p = c_avg_p;
% out.eta_n = eta_n;
% out.eta_p = eta_p;
% out.eta_s_Ln = eta_s_Ln;
% out.ce0p = c_ex(end,:);

fileName = input('Save filename? ','s');
save([fileName '.mat'],'out');


