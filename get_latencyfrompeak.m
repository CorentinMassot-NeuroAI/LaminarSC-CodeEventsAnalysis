function [platencies]=get_latencyfrompeak(trials,peaks,type,info,alpha)

%function [latencies]=get_latencyfrompeak(trials,peaks,type,alpha)
%   get latency of burst of activity in trials starting from peak activity
%
% type: 'spk' or 'lfp'
%
% see also get_latency get_peakfromlatency
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 11/09/2016 last modified 07/24/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

debug=0;
%time ttest stays positive
pos_time=1;

%dprime values
switch info.align
    case 'targ'
        dprime_spk=1;
        dprime_lfp=1;
        winsearch=60
                   
    case 'sacc'
        dprime_spk=1;
        dprime_lfp=1;
        winsearch=30;
end     

[nchannels ~]=size(trials);


H=zeros(nchannels,1);
platencies=zeros(nchannels,3);
platencies(:,2)=1;

if debug,
    histo=nan(winsearch,2);
end

%trials_w=zeros(nchannels,winsearch);
for ch=1:nchannels
    %check peak
    %if peaks(ch,1)>winsearch/2+1+1
    if peaks(ch,1)>3*winsearch/4+1+1
        %peak window
        %peak_w=trials(ch,peaks(ch,1)-winsearch/2+1:peaks(ch,1)+winsearch/2);
        peak_w=trials(ch,peaks(ch,1)-floor(3*winsearch/4)+1:peaks(ch,1)+winsearch/3);
        
        for w=peaks(ch,1):-1:winsearch+1
            trials_w=trials(ch,w-winsearch+1:w);
            
            if platencies(ch,1)~=pos_time,%channel did not reach significance for pos_time
                if debug
                    %plot
                    hdlfig=figure;
                    subplot(1,2,1);hold on;
                    plot(1:winsearch,trials_w(:),'b');
                    plot(1:winsearch,peak_w(:),'r');
                    
                    %histo
                    subplot(1,2,2);hold on;
                    histo(:,1)=trials_w(:)';
                    histo(:,2)=peak_w(:)';
                    edges=[-50:10:250];
                    hist=histc(histo,edges);
                    bar(edges,hist,'histc')
                end
                %one tailed ttest
                switch type
                    case 'spk'
                        %H(ch)=ttest2(trials_w(:),peak_w(:),'alpha',alpha,'tail','left');
                        %[P(ch) H(ch)]=ranksum(trials_w(ch,:),peaks(ch,:),'alpha',alpha);
                        H(ch)=compute_dprime(trials_w(:),peak_w(:),dprime_spk);
                        
                    case 'lfp'
                        %H(ch)=ttest2(trials_w(:),peak_w(:),'alpha',alpha,'tail','both');
                        %[P(ch) H(ch)]=ranksum(trials_w(ch,:),peaks(ch,:),'alpha',alpha);
                        H(ch)=compute_dprime(trials_w(:),peak_w(:),dprime_lfp);
                        
                end
                
                if H(ch)==1
                    %count
                    platencies(ch,1)=platencies(ch,1)+1;
                    
                    if platencies(ch,1)==1
                        
                        %latency
                        platencies(ch,2)=w-round(winsearch/2);
                        
                        %magnitude at latency
                        platencies(ch,3)=mean(trials_w(:));
                    end
                else
                    platencies(ch,1)=0;
                    platencies(ch,2)=1;
                    platencies(ch,3)=0;
                end
                
                %             display(['H of ch ' num2str(ch) ' :' num2str(H(ch))]);
                if debug
                    ch
                    platencies(ch,:)
                    pause
                    %close(hdlfig)
                end
            end
            
        end
        
        
        %check if all latencies have been found
        if sum(platencies(:,1)~=pos_time)==0,
            %         w
            return;
        end
    end
end


