function [H_all]=get_bsignif(trials,trials_a,bsls,bsls_a,alpha)

%function [allH]=get_bsignif(trials,trials_a,bsls,bsls_a,alpha)
%   compute significance of activity in trials compare to trials_a, usually
%   inRf compare to anti-RF
%
%  bsls,bsls_a; baseline activity for normalization is required
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 11/01/2016 last modified 01/22/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nchannels ntrials triallen]=size(trials);
[nchannels ntrials_a triallen_a]=size(trials_a);
H_all=zeros(1,nchannels);
P_all=zeros(1,nchannels);

%  histo=nan(max(ntrials,ntrials_a),2);
%  figure;

for ch=nchannels:-1:1
    %avg values and correction for baseline across trials
    trials_s=squeeze(trials(ch,:,:));
    trials_a_s=squeeze(trials_a(ch,:,:));
    
%         %mean
%         trials_avg=nanmean(trials_s,2);
%         trials_a_avg=nanmean(trials_a_s,2);
    
    %peak (of abs for lfp)
    trials_avg=nanmax(abs(trials_s),[],2);
    trials_a_avg=nanmax(abs(trials_a_s),[],2);

    if ~isempty(bsls)
    bsls_avg=nanmean(squeeze(bsls(ch,:,:)),2);
    bsls_a_avg=nanmean(squeeze(bsls_a(ch,:,:)),2);
    
    %WARNING take abs for LFP but should be good for spk too
    trials_avg=abs(trials_avg - bsls_avg);
    trials_a_avg=abs(trials_a_avg - bsls_a_avg);
    end
    
%         %histo
%         histo(1:ntrials,1)=trials_avg';
%         histo(1:ntrials_a,2)=trials_a_avg';
%         edges=[-50:10:250];
%         %edges=[0:10:400];
%         hist=histc(histo,edges);
%         bar(edges,hist,'histc')
    
    %ranksum test of significance   
    [P_all(ch) H_all(ch)]=ranksum(trials_avg,trials_a_avg,'alpha',alpha,'tail','right');
    
    %dprime
    %H_all(ch)=compute_dprime(trials_avg,trials_a_avg,1);

%      display(['H of ch ' num2str(ch) ' :' num2str(H_all(ch))]);
%      pause
end

%for ranksum test
%P_all
H_all=P_all<alpha;
