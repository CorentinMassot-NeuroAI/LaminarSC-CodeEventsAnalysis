function [trials_avg trials_var]=get_trials_avg(trials)

%function get_trials_avg(trials)
%   compute and plot trial average of trials
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 07/06/2016 last modified 07/06/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nchannels ntrials triallen]=size(trials);
trials_avg=zeros(nchannels,triallen);
trials_var=zeros(nchannels,triallen);
for ch=1:nchannels
    %take into account NaN values
    trials_avg(ch,:)=nanmean(squeeze(trials(ch,:,:)));
    ntrials_ch=sum(~isnan(squeeze(trials(ch,:,:))),1);
    
    %std
    %trials_var(ch,:)=sqrt(nanvar(squeeze(trials(ch,:,:)),0,1));
    %ste
    trials_var(ch,:)=sqrt(nanvar(squeeze(trials(ch,:,:)),0,1))./ntrials_ch;
end

