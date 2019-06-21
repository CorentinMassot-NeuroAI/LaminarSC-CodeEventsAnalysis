function meanci=get_ci(data,ci_level)

%function meanci=get_ci(data,ci_level)
%   compute mean and ci from data for mean operation!!
%
%WARNING: change data_nnan by mean computed from original sample because
%otherwise use mean(bootstrapped mean values) as an estimate of mean from
%original sample
%
% ci_level: % of confidence interval (e.g. 0.95 for 95%)
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 12/28/2017 last modified 12/28/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%remove nan values
data_nnan=data(find(~isnan(data)));

%alpha
alpha=1-ci_level;   

%T-score for level
ts = tinv([1-alpha/2],length(data_nnan)-1);

%mean data_nnan
%WARNING: change data_nnan by mean computed from original sample
mean_data_nnan=mean(data_nnan);

%Standard Error
sem = std(data_nnan)/sqrt(length(data_nnan));               


%mean and confidence intervals
meanci = [mean_data_nnan ts*sem]; 
