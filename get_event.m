function event=get_event(trial,eventname,option)

%function event=get_event(trial,code,option)
%   get event time based on code in statTransitions in trial
%
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 01/09/2017 last modified 01/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

event=[];

if strcmp(eventname,'sacc'),
    %saccade onset based on rpt file (EyeMove)
    if strcmp(option,'rpt')
        meas=get_event(trial,'measCode');
        srt=trial.behavrpt.srts+6; %WARNING:6ms added to correct for eyemove
        if isempty(srt)
            %display('Warning Saccade reaction time is empty (from .rpt file)!')
        else
            event=srt+meas;
        end
        
        %alternative
        %event=trial.behavrpt.saccTime
        
    end
    
elseif strcmp(eventname,'peak'),
    %saccade peak velocity computed offline (see compute_gazedata)
    event=trial.offline.gazepos_events.peak(1,1);
    
elseif strcmp(eventname,'targ_pburst'),
    code=eval(['trial.params.' 'targCode']);
    codeind=find(trial.stateTransitions(1,:) == code);
    event_targ=double(trial.stateTransitions(2,codeind));
    
    %visual burst computed offline (see compute_pburst)
    event=event_targ+trial.offline.targ_pburst_trial.b_begin;
    %remove channel without detected burst (0 by default)
    %TO DO: replace by nan in compute_pburst
    indnan=find(trial.offline.targ_pburst_trial.b_begin==0 & trial.offline.targ_pburst_trial.b_surprise==0);
    event(indnan)=0;
    
elseif strcmp(eventname,'targ_pburst_ch'),
    code=eval(['trial.params.' 'targCode']);
    codeind=find(trial.stateTransitions(1,:) == code);
    event_targ=double(trial.stateTransitions(2,codeind));
    
    %ch_align of visual burst computed offline (see compute_pburst)
    ch_align=trial.offline.targ_pburst_ch_align.ch;
    
    if trial.offline.targ_pburst_trial.b_begin(ch_align)==0 & trial.offline.targ_pburst_trial.b_surprise(ch_align)==0;
        %remove trial without detected burst on ch_align (0 by default)
        %TO DO: replace by nan in compute_pburst
        event=0;
    else
        event=event_targ+trial.offline.targ_pburst_trial.b_begin(ch_align);
    end
    
    
elseif strcmp(eventname,'targ_rsburst'),
    code=eval(['trial.params.' 'targCode']);
    codeind=find(trial.stateTransitions(1,:) == code);
    event_targ=double(trial.stateTransitions(2,codeind));
    
    %visual burst computed offline (see compute_pburst)
    event=event_targ+trial.offline.targ_rsburst_trial.b_begin;
    %remove channel without detected burst (0 by default)
    %TO DO: replace by nan in compute_pburst
    indnan=find(trial.offline.targ_rsburst_trial.b_begin==0 & trial.offline.targ_rsburst_trial.b_surprise==0);
    event(indnan)=0;
    
elseif strcmp(eventname,'targ_rsburst_ch'),
    code=eval(['trial.params.' 'targCode']);
    codeind=find(trial.stateTransitions(1,:) == code);
    event_targ=double(trial.stateTransitions(2,codeind));
    
    %ch_align of visual burst computed offline (see compute_pburst)
    ch_align=trial.offline.targ_rsburst_ch_align.ch;
    if trial.offline.targ_rsburst_trial.b_begin(ch_align)==0 & trial.offline.targ_rsburst_trial.b_surprise(ch_align)==0;
        %remove trial without detected burst on ch_align (0 by default)
        event=0;
    else
        event=event_targ+trial.offline.targ_rsburst_trial.b_begin(ch_align);
    end
    
    
    
elseif strcmp(eventname,'sacc_pburst'),
    %saccade onset based on rpt file (EyeMove)
    if strcmp(option,'rpt')
        meas=get_event(trial,'measCode');
        srt=trial.behavrpt.srts+6; %WARNING:6ms added to correct fro eyemove
        if isempty(srt)
            %display('Warning Saccade reaction time is empty (from .rpt file)!')
        else
            event_sacc=srt+meas;
        end
        
        %alternative
        %event=trial.behavrpt.saccTime
        
    end
    
    %visual burst computed offline (see compute_pburst)
    event=event_sacc+trial.offline.targ_pburst_trial.b_begin;
    
elseif strcmp(eventname,'sacc_rsburst'),
    %saccade onset based on rpt file (EyeMove)
    if strcmp(option,'rpt')
        meas=get_event(trial,'measCode');
        srt=trial.behavrpt.srts+6; %WARNING:6ms added to correct fro eyemove
        if isempty(srt)
            %display('Warning Saccade reaction time is empty (from .rpt file)!')
        else
            event_sacc=srt+meas;
        end
        
        %alternative
        %event=trial.behavrpt.saccTime
        
    end
    
    %visual burst computed offline (see compute_pburst)
    event=event_sacc+trial.offline.targ_rsburst_trial.b_begin;
    
    
    
else %if ~strcmp(eventname,'sacc')
    code=eval(['trial.params.' eventname]);
    codeind=find(trial.stateTransitions(1,:) == code);
    event=double(trial.stateTransitions(2,codeind));
    
end