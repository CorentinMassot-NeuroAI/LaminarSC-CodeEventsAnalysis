function [trials_avgn]=get_trials_avg_normalized(trials_avg,trials_bsl_avg,signal,info)

%[trials_avgn]=get_trials_normalized(trials_avg,trials_bsl_avg,info)
%   normalizes trial average by baseline
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 01/22/2017 last modified 01/22/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%trials_avg normalized
trials_avgn=zeros(size(trials_avg));

%compute baseline
baseline=mean(trials_bsl_avg,2);

%normalize across channels
for ch=1:size(trials_avg,1)
    trials_avgn(ch,:)=trials_avg(ch,:)-baseline(ch);
end
display([signal info.align ' normalized to baseline'])