%analysis_pburst

%function analysis_pburst
%   analyze bursting activity based on burst detection (pburst)
%   mostly to check if data saved by compute_pburst is ok
%
% see also compute_pburst 
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 10/29/2017 last modified 10/29/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set paths
[root_path data_path save_path]=set_paths;

%screen size
scrsz = get(groot,'ScreenSize');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
%print figures and save data
savedata=0;


%alignement
%alignlist={'no' 'targ' 'go' 'sacc'};
%alignlist={'targ' 'sacc'};
%alignlist={'targ_pburst_ch' 'sacc'};
%alignlist={'targ_rsburst_ch' 'sacc'};
alignlist={'sacc'};

%windows of analysis
%wind_no=[-200 1500];
wind_targ=[0 600];
wind_sacc=[-300 200];%[-200 200];

%algo
algolist={'pburst' 'rsburst'}

%window of burst
windb_targ=[30 150];
windb_sacc=[-100 0];

%vshift
vshift_fr=150;

%sigma FR
sigma_FR=6;

%burst trheshold
thresh=10;%threshold spk/s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get data
datalist=load_data_gandhilab(data_path);

%colorlist
colorlist=get_colorlist;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analyzing data
dlist=get_dlist

%hdlfigallbsignifs=figure;hold on;
data=[];
info=[];
allbob_avg=[];allbobmode_avg=[];allbob_avgn=[];allbobmode_avgn=[];
dd=0;
for d=dlist(1:end)%1:numel(datalist)
    %counter
    dd=dd+1;
    
    
    data=[];
    %get data and info
    info.datafile=datalist{d};
    load ([data_path info.datafile]);
    display(info.datafile)
    
    %getting channel mapping and discard selected bad channels
    [info.chmap info.nchannels info.depths]=get_chmap(data(1).info.electrode{2},[]);
    
    %getting trial type
    info.trialtype=data(1).sequence(1);
    %getting list of targets
    targslist=data(1).offline.targslist;
    
    %target tuning (after compute_tuning)
    targ_tuning=data(1).offline.targ_tuning;
    
    %select trials
    seltrials=get_seltrials(data,'rpt');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Data aligned on target and saccade onset
    for al=1:numel(alignlist)
        info.align=alignlist{al};
        
        switch info.align
            case 'targ'
                wind=wind_targ;
            case 'targ_pburst_ch'
                wind=wind_targ;
            case 'sacc'
                wind=wind_sacc;
        end
        
        [alltrials_spk info.aligntime lut_trials]=get_alltrials_align(data,seltrials,wind,'spk',info,targslist,sigma_FR,1);
        [alltrials_fr ~]=get_alltrials_align(data,seltrials,wind,'fr',info,targslist,sigma_FR,1);
        [alltrials_lfp ~]=get_alltrials_align(data,seltrials,wind,'lfp',info,targslist,sigma_FR,1);
        [allgazepos,allevents]=get_alldatagaze_align(data,seltrials,info,targslist);
        
        %target
        info.targ=targ_tuning;
        
        %selection of tuning data only
        trials_fr=alltrials_fr{info.targ};
        trials_spk=alltrials_spk{info.targ};
        trials_lfp=alltrials_lfp{info.targ};
        gazepos=allgazepos{info.targ};
        events=allevents{info.targ};
        lut_trials_targ=lut_trials{info.targ};
        
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Plot each trial
        
        [info.nchannels info.ntrials info.triallen]=size(trials_fr);
        titlestr={info.datafile ; ['FR ' info.align ' t' num2str(info.targ) ' #trials:' num2str(info.ntrials)]};
        for t=1:info.ntrials,
            
            figtrials=figure('Position',scrsz);hdlfig=subplot(1,1,1);hold on;
            
            trials_fr_t=squeeze(trials_fr(:,t,:));
            trials_spk_t=squeeze(trials_spk(:,t,:));
            events_t=events{t};
            event_align=get_eventalign(events_t,info.align);
            
            %plot
            [range ~]=plot_trials(trials_fr_t,[],[1:info.nchannels],vshift_fr,events_t,event_align,info,hdlfig,titlestr,'-',1);
            
            %Pburst (see compute_pburst)
            for alg=1:2
                algo=algolist{alg};
                b_begin_plot=zeros(info.nchannels,2);b_end_plot=zeros(info.nchannels,2);
                switch algo
                    case 'pburst'
                        eval(['b_begin=data(' num2str(lut_trials_targ(t)) ').offline.' info.align '_pburst_trial.b_begin;']);
                        eval(['b_end=data(' num2str(lut_trials_targ(t)) ').offline.' info.align '_pburst_trial.b_end;']);
                        b_begin_plot(:,1)=b_begin;%-info.aligntime;
                        b_end_plot(:,1)=b_end;%-info.aligntime;
                        plot_events_ch(b_begin_plot,[],vshift_fr,range,info,hdlfig,'n','-',3);
                        plot_events_ch(b_end_plot,[],vshift_fr,range,info,hdlfig,'n','-',3);
                    case 'rsburst'
                        eval(['b_begin=data(' num2str(lut_trials_targ(t)) ').offline.' info.align '_rsburst_trial.b_begin;']);
                        eval(['b_end=data(' num2str(lut_trials_targ(t)) ').offline.' info.align '_rsburst_trial.b_end;']);
                        b_begin_plot(:,1)=b_begin;%-info.aligntime;
                        b_end_plot(:,1)=b_end;%-info.aligntime;
                        plot_events_ch(b_begin_plot,[],vshift_fr,range,info,hdlfig,'n',':',5);
                        plot_events_ch(b_end_plot,[],vshift_fr,range,info,hdlfig,'n',':',5);
                end
            end
            pause
            close(figtrials)
        end
    end
end
        
        