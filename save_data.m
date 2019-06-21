function data=save_data(data,root_path,dir,info)

%function data=save_data(data,root_path,info)
%
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 01/25/2017 last modified 01/26/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%save path
save_path=[root_path(1:end-9) 'Data\' dir];

%datafile
%for Joy
%datafile=['FRLFP_' data{1}.info.align '_' data{1}.info.datafile];
%datafile=['FRLFPtrials_' data{1}.info.align '_' data{1}.info.datafile];
datafile=['RAWtrials_' data.info.align '_' data.info.datafile];


%for Sanjeev
%datafile=[info.datafile(1:end-4) '_sanjeev'];


%save data
display(['saving data in ' save_path datafile]);
save([save_path datafile], 'data')

