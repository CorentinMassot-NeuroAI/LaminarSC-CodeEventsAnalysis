%function compute_gazedata

%function compute_gazedata
%   Compute trial-by-trial gaze data recorded with a laminar probe (LMA)
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 08/09/2017 last modified 08/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set paths
[root_path data_path save_path]=set_paths;

%screen size
scrsz = get(groot,'ScreenSize');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
%print figures and save data
savefigs=0;
figtype='epsc2';%'png';%'epsc2';

%alignement
%alignlist={'no' 'targ' 'go' 'sacc'};
alignlist={'sacc'};

%window of analysis
%wind=[];%all
%wind=[-100 400];%targ align
wind=[-200 200];%sacc align

%sigma FR
sigma_FR=6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get data
datalist=load_data_gandhilab(data_path);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analyzing data
dlist=get_dlist;

data=[];
info=[];
for d=dlist(10:end)
    
    %get data and info
    clear('data');
    info.datafile=datalist{d};
    load ([data_path info.datafile]);
    display(info.datafile)
    
    %getting channel mapping and discard selected bad channels
    [info.chmap info.nchannels info.depths]=get_chmap(data(1).info.electrode{2},[]);
    %getting trial type
    info.trialtype=data(1).sequence(1);
    %getting list of targets
    targslist=data(1).offline.targslist;
    %targets index
    targs_ind=get_targsindex(targslist,info);
    
    %target tuning (after compute_tuning)
    info.targ_tuning=data(1).offline.targ_tuning;
    
    %select trials
    seltrials=get_seltrials(data,'rpt');
    
    
    %loop across all alignements
    for al=1:numel(alignlist)
        info.align=alignlist{al};
        
        %get behavioral data with specific alignement independently of target
        [gazepos,events]=get_alldatagaze_align_notarg(data,seltrials,info);
        
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %loop on trials
        info.ntrials=numel(seltrials)
        
        for t=seltrials,
            
            figtrials=figure('Position',[1 100 scrsz(3)-100 scrsz(4)-200]);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %display all targets
            hdlfig=subplot(1,2,1);hold on;
            display_alltargets(targslist,info,hdlfig);
            
            
            %%
            %%%%%%%%%%%%%%%%%%
            %gaze data
            gazepos_t=gazepos{t};
            events_t=events{t};
            event_align=get_eventalign(events_t,info.align);
            
            if ~isempty(event_align)
                
                hdlfig1=subplot(1,2,1);hold on;
                hdlfig2=subplot(1,2,2);hold on;
                
                %%
                %%%%%%%%%%%%%%%%%%
                %get new gaze position events
                display_gazepos=0;
                events_t_new=get_gazepos_events(gazepos_t,events_t,event_align,wind,info,display_gazepos,hdlfig1,hdlfig2);
                
                %plot gaze data
                %plot_gazedata(gazepos_t,events_t,event_align,wind,info,hdlfig1,'',hdlfig2,'XY Eye Traces');
                %subplot(1,2,1);
                %ax=axis;
                %axis([ax(1) ax(2) -4 4]);
                
            end
                       
            %update data
            %WARNING: error: should use lut_trials_targ to ensure right trial number
            
            display('ERROR in SAVING NOT SAVED!')
            %display('NOT SAVED!')
            %velocity peak
            %eval(['data(t).offline.gazepos_events.peak=events_t_new.peak;']);
           
            
            %pause
            close(figtrials)

        end
    end
    
    %save updated data
    update_data(1,0,0,data,data_path,info.datafile,[],[]);
    
end
