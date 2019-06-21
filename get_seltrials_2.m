function [seltrials]=get_seltrials(data,option)

%[seltrials]=get_seltrials(data,info,targslist)
%   select trials according to success and other criteria from data recorded with a laminar probe (LMA)
%
%  option: 'rpt' if use .rpt file (Eyemove)
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 01/18/2016 last modified 01/18/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%initializations
seltrials=[];
ntrials_cond=zeros(1,5);

for t=1:numel(data)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get events
    %sacc onset
    if strcmp(option,'rpt')
        srt=get_event(data(t),'sacc','rpt');
    end
    
    %target onset
    targ=get_event(data(t),'targCode');
    
    %go onset
    go=get_event(data(t),'goCode');
    
    %fixation target
    fixtarget=get_event(data(t),'fixtarget');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %conditions
    %analyzing only successful trials defined by MonkeyHost
    if 1%data(t).success==1,
        ntrials_cond(1)=ntrials_cond(1)+1;
        
        %srt~=[]
        if ~isempty(srt)
            ntrials_cond(2)=ntrials_cond(2)+1;
            
            %go<srt<fixtarget
            if srt>go
                ntrials_cond(3)=ntrials_cond(3)+1;
                %ntrials_cond(3)=ntrials_cond(2);
                
                %& srt<fixtarget %NOTE: condition on fixtarget: use behavioral data instead
                
                %saccade duration if <100ms
                sacc_dur=fixtarget-srt;
                if sacc_dur<100
                    ntrials_cond(4)=ntrials_cond(4)+1;
                    
                    
                    
%                     %saccade reaction time
%                     rt=srt-go;
%                     if rt>100
%                         ntrials_cond(5)=ntrials_cond(5)+1;
                        
                        seltrials=[seltrials t];
                        
                        
%                     end
                end
            end
        end
    end
    
end



%info
display(['# of selected trials: ' num2str(length(seltrials)) '/' num2str(numel(data))]);
display(['# of successful trials: ' num2str(ntrials_cond(1)) '/' num2str(numel(data))]);
display(['# of trials without missing srt: ' num2str(ntrials_cond(2)) '/' num2str(ntrials_cond(1))])
display(['# of trials with go<srt<fixtarget: ' num2str(ntrials_cond(3)) '/' num2str(ntrials_cond(2))]);
display(['# of trials with sacc duration <100ms: ' num2str(ntrials_cond(4)) '/' num2str(ntrials_cond(3))]);
%display(['# of trials with sacc reaction time >200ms: ' num2str(ntrials_cond(5)) '/' num2str(ntrials_cond(4))]);


