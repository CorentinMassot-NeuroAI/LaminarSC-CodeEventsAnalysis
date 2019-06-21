function [peaks]=get_peakfromlatency(trials,latencies,type)

%function [peaks]=get_peakfromlatency(trials,type,alpha)
%   get peaks of burst of activity in trials after having found first
%   estimation of latencies using get_latency
%
% type: 'spk' or 'lfp'
%
% see also get_latency and get_latencyfrompeak
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 11/08/2016 last modified 06/21/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

debug=0;

[nchannels triallen]=size(trials);

peaks=zeros(nchannels,2);
trials_r=[];trials_rs=[];dtrials_rs=[];
for ch=1:nchannels
    if latencies(ch)~=0
        %reducing trials with latencies
        trials_r=trials(ch,latencies(ch):end);
        
        %smoothing
        trials_rs=smooth(trials_r,20,'rloess');
        
        %diff
        dtrials_rs=diff(trials_rs',1,2);
        dtrials_rss=smooth(dtrials_rs,8,'moving');
        
        if debug
            %plot
            figure;hold on;
            plot(trials_r,'b');
            plot(trials_rs,'r');
            plot(dtrials_rs*10,'g');
            plot(dtrials_rss*10,'m');
        end
        
        %peak (zeros crossing of dtrial_rs)
        switch type
            case 'fr'
                start=30;delta=8;maxwin=15;
                found=0;
                for i=start:length(dtrials_rss)
                    if ~found
                        if sum(dtrials_rss(i-1-delta:i-1)>0)>delta-1 & sum(dtrials_rss(i+1:i+1+delta)<0)>delta-1
                            %find real max (after smoothing true max may be shifted)
                            [peaks(ch,2) peaks(ch,1)]=max(trials_r(i-maxwin:i+maxwin));
                            peaks(ch,1)=peaks(ch,1)+i-maxwin-1;
                            found=1;
                        end
                    end
                end
                
                
                
            case 'lfp'
                found=0;
                
                %polarity
                [maxval_lfp maxind_lfp]=max(trials_r);
                [minval_lfp minind_lfp]=min(trials_r);
                if maxval_lfp>-minval_lfp
                    polarity=1;
                else
                    polarity=-1;
                end
                
                if polarity==1
                    start=50;delta=12;maxwin=15;
                    for i=start:length(dtrials_rss)
                        if ~found
                            dtrials_rss(find(trials_r<0))=nan;
                            if (nansum(dtrials_rss(i-1-delta:i-1)>0)>delta-1 & nansum(dtrials_rss(i+1:i+1+delta)<0)>delta-1)
                                %find real max (after smoothing true max may be shifted)
                                [peaks(ch,2) peaks(ch,1)]=max(trials_r(i-maxwin:i+maxwin));
                                peaks(ch,1)=peaks(ch,1)+i-maxwin-1;
                                found=1;
                            end
                        end
                    end
                elseif polarity==-1
                    start=minind_lfp+20;delta=12;maxwin=15;
                    for i=start:-1:1+1+delta
                        if found~=2
                            dtrials_rss(find(trials_r>0))=nan;
                            if (nansum(dtrials_rss(i-1-delta:i-1)<0)>delta-1 & nansum(dtrials_rss(i+1:i+1+delta)>0)>delta-1)
                                %find real min (after smoothing true max may be shifted)
                                [peaks(ch,2) peaks(ch,1)]=min(trials_r(i-maxwin:i+maxwin));
                                peaks(ch,1)=peaks(ch,1)+i-maxwin-1;
                                found=found+1;
                            end
                        end
                    end
                end
        end

        if debug
            line([peaks(ch,1) peaks(ch,1)],[peaks(ch,2)-10 peaks(ch,2)+10])
            grid
            pause
        end
        
        %timing correction
        peaks(ch,1)=peaks(ch,1)+latencies(ch)-1;
    
    end
end



