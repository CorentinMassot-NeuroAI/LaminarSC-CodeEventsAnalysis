function allstats=get_allbehavstats(data,seltrials,targslist,option)

%function event=get_behavstats(trial,code,option)
%   get statistics from behvioral data
%
% seltrials: see select_trials
% option:'rpt' use .rpt file (EyeMove)
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 01/23/2017 last modified 01/23/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%initializations
ntargs=size(targslist,1);
allstats=cell(ntargs,1);
ntrials=zeros(ntargs,1);

%get data only for selected_trials
for t=seltrials
    
    %get target index for the trial
    pos=data(t).offline.targpos;
    it=find(pos(1,1)==targslist(:,1) & pos(1,2)==targslist(:,2));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get events
    %sacc onset
    if strcmp(option,'rpt')
        srt=get_event(data(t),'sacc','rpt');
    end
    
    %go onset
    go=get_event(data(t),'goCode');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %saccade reaction time relative to go signal
    stats=srt-go;
    
    
    %creating list of stats
    if ~isempty(allstats{it}),allstatsaux=allstats{it};else allstatsaux=[];end;
    ntrials(it)=ntrials(it)+1;
    allstatsaux(ntrials(it))=stats;
    allstats{it}=allstatsaux;
    
end

%allstatsaux




