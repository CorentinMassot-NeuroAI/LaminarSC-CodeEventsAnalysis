function [trials_c index_c] = clean_trials(trials,signal)

%function [trials_listc index_c] = clean_trials(trials_list,signal)
%  clean list of trials by removing bad channels
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 05/25/2016 last modified 05/25/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch signal
    case 'fr'
        minamp=2;
    case 'lfp'
        minamp=1%4;%;%8;%10;
    case 'raw'
        minamp=0;
        
    case 'no'
        minamp=0;
end;

i=0;
trials_c=[];index_c=[];
for ch=1:size(trials,1)
    if (max(abs(trials(ch,:))-min(abs(trials(ch,:))))>minamp)
        i=i+1;
        trials_c(i,:)=trials(ch,:);
        index_c(i)=ch;
    end
end


