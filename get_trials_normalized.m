function [trials_n]=get_trials_normalized(trials,trials_bsl,signal,info)

%[trials_avgn]=get_trials_normalized(trials_avg,trials_bsl_avg,signal,info)
%   normalizes trials by baseline
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 01/22/2017 last modified 08/30/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%trials normalized
trials_n=zeros(size(trials));

%normalize trials across channels
for ch=1:size(trials,1)
    for t=1:size(trials,2)
        trials_n(ch,t,:)=(squeeze(trials(ch,t,:))-mean(squeeze(trials_bsl(ch,t,:))))';
    end
end
display([signal info.align ' trials normalized to baseline'])