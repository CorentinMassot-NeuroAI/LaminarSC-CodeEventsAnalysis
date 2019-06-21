function [seltrials]=get_seltrials_features(data,srtlims,slopemax,repcond)

%[seltrials]=get_seltrials_features(data,info,targslist)
%   select trials according trial features from data recorded with a laminar probe (LMA)
%
% repcon: 1 exclude repeats -1 include only repeats
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 05/01/2018 last modified 05/01/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initializations
seltrials=[];
ntrials_cond=zeros(1,5);

for t=1:numel(data)
    %t
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get features
    %srt
    srt_t=data(t).offline.statstrial.srt;
    
    %peak saccade velocity
    
    %final fixation position
    
    %trend
    slope_t=data(t).offline.statstrial.trend.fitparams(:,2);%slope
    
    %repeat
    rep_t=data(t).offline.repeat;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %conditions
    sel=[];
    
    %srt
    if ~isempty(srtlims)
        if srt_t>srtlims(1) & srt_t<srtlims(2)
            ntrials_cond(1)=ntrials_cond(1)+1;
            sel(1)=1;
        else
            sel(1)=0;
        end
    else
        sel(1)=1;
    end
    
    %trend
    if ~isempty(slopemax)
        if sum(slope_t>slopemax)==0
            ntrials_cond(2)=ntrials_cond(2)+1;
            sel(2)=1;
        else
            sel(2)=0;
        end
    else
        sel(2)=1;
    end
    
    %repeat
    if repcond==1 %exclude repeats
        if rep_t==0
            ntrials_cond(3)=ntrials_cond(3)+1;
            sel(3)=1;
        else
            sel(3)=0;
        end
    elseif repcond==-1 %include only repeats
        if rep_t==1
            ntrials_cond(3)=ntrials_cond(3)+1;
            sel(3)=1;
        else
            sel(3)=0;
        end
    else
        sel(3)=1;
    end
    
     %seltrials
    if sel(1)==1 & sel(2)==1 & sel(3)==1
        seltrials=[seltrials t];
    end
end



%info
display(['Trials features:']);
display(['# of trials within srt limits: ' num2str(ntrials_cond(1)) '/' num2str(numel(data))]);
display(['# of trials within slope limits: ' num2str(ntrials_cond(2)) '/' num2str(numel(data))]);
if repcond==1,display(['# of trials that are not repeats: ' num2str(ntrials_cond(3)) '/' num2str(numel(data))]);
elseif repcond==-1,display(['# of trials that are only repeats: ' num2str(ntrials_cond(3)) '/' num2str(numel(data))]);
end 