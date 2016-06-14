%% Save Sensitivity Analysis Results of DFN
%   Created February 23, 2014 by Scott Moura

% Meta Data
out.date = date;

% DFN Data Filename
out.fn = fn;

% Sensitivity matrices
out.S1 = S1;
out.S2 = S2;
out.S3 = S3;

fileName = input('Save filename? ','s');
save([fileName '.mat'],'out');