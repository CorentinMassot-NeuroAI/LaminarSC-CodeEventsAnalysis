%compute_peaks_pburst

%function compute_peaks_pburst
%   analyze peaks of spiking activity and LFP based on burst detection (pburst)
%
% see also compute_pburst compute_bsignif
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 11/07/2017 last modified 11/07/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set paths
[root_path data_path save_path]=set_paths;

%screen size
scrsz = get(groot,'ScreenSize');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
%print figures and save data
savedata=0;

%display figures
disp=1;

%alignement
%alignlist={'targ_pburst'};
%alignlist={'targ_pburst_ch'};
%alignlist={'targ_rsburst_ch'};
alignlist={'sacc'};

%windows of analysis
wind_targ_pburst=[-50 150];%[-200 300];
wind_sacc=[-100 100];

%windows baseline
wt=100;
wind_targ_pburst_bsl=[-50-wt -50 ];
ws=50;
wind_sacc_bsl_go=[-ws 0];

%vshift
vshift_fr=150;
vshift_lfp=30;

%sigma FR
sigma_FR=6;

%shift ripple temporal correction
shift_ripple=4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get data
datalist=load_data_gandhilab(data_path);

%colorlist
colorlist=get_colorlist;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analyzing data
dlist=get_dlist

dd=0;
for d=dlist(1:end)%1:numel(datalist)

    %counter
    dd=dd+1;
    
    %get data and info
    data=[];info=[];
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
            case 'targ_pburst_ch'
                wind=wind_targ_pburst;
                wind_bsl=wind_targ_pburst_bsl;
                
                %trials
                [alltrials_fr info.aligntime lut_trials]=get_alltrials_align(data,seltrials,wind,'fr',info,targslist,sigma_FR,1);
                %correct for shift introduced by ripple filtering
                [alltrials_lfp ~]=get_alltrials_align(data,seltrials,wind+shift_ripple,'lfp',info,targslist,sigma_FR,1);
                %baseline
                [alltrials_fr_bsl aligntime_bsl]=get_alltrials_align(data,seltrials,wind_bsl,'fr',info,targslist,sigma_FR,1);
                [alltrials_lfp_bsl aligntime_bsl]=get_alltrials_align(data,seltrials,wind_bsl,'lfp',info,targslist,sigma_FR,1);
                
            case 'sacc'
                wind=wind_sacc;
                wind_bsl=wind_sacc_bsl_go;
                
                %trials
                [alltrials_fr info.aligntime lut_trials]=get_alltrials_align(data,seltrials,wind,'fr',info,targslist,sigma_FR,1);
                %correct for shift introduced by ripple filtering
                [alltrials_lfp ~]=get_alltrials_align(data,seltrials,wind+shift_ripple,'lfp',info,targslist,sigma_FR,1);
                %baseline
                alignaux=info.align;
                info.align='go';
                [alltrials_fr_bsl aligntime_bsl]=get_alltrials_align(data,seltrials,wind_bsl,'fr',info,targslist,sigma_FR,1);
                [alltrials_lfp_bsl aligntime_bsl]=get_alltrials_align(data,seltrials,wind_bsl,'lfp',info,targslist,sigma_FR,1);
                info.align=alignaux;
                
        end
        
        %target
        info.targ=targ_tuning;
        
        %selection of tuning data only
        trials_fr=alltrials_fr{info.targ};
        trials_lfp=alltrials_lfp{info.targ};
        trials_fr_bsl=alltrials_fr_bsl{info.targ};
        trials_lfp_bsl=alltrials_lfp_bsl{info.targ};
        
        [allgazepos,allevents]=get_alldatagaze_align(data,seltrials,info,targslist);
        events=allevents{info.targ};
        lut_trials_targ=lut_trials{info.targ};
        
        [info.nchannels info.ntrials info.triallen]=size(trials_fr);
        
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %normalization per trial
        trials_fr_n=get_trials_normalized(trials_fr,trials_fr_bsl,'FR',info);
        trials_lfp_n=get_trials_normalized(trials_lfp,trials_lfp_bsl,'LFP',info);
        %compute average normalized trials
        [trials_fr_n_avg trials_fr_n_var]=get_trials_avg(trials_fr_n);
        [trials_lfp_n_avg trials_lfp_n_var]=get_trials_avg(trials_lfp_n);
        
        d
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %plot spk lfp aligned on burst per channel
        info.nchs=info.nchannels;
        for t=1:info.ntrials
          
            latencies=zeros(1,info.nchs);
            for ch=1:info.nchs
                
                %trial
                trials_fr_tch=squeeze(trials_fr_n(ch,t,:))';
                trials_lfp_tch=squeeze(trials_lfp_n(ch,t,:))';
                
                if disp
                    figtrial=figure('Position',[scrsz(1) scrsz(4)*0.5-100 scrsz(3) scrsz(4)*0.5]);
                    hdlfigtrial=subplot(1,1,1);hold on;
                    
                    %events_t=events{t};
                    %event_align=get_eventalign(events_t,info.align);
                    
                    
                    %fr
                    titlestr={info.datafile ; ['FR ' info.align ' t' num2str(info.targ) ' ch:' num2str(ch) ' trial:' num2str(t) '/' num2str(info.ntrials)]};
                    info.color=colorlist(1,:);info.nchannels=1;
                    %plot_trials(trials_fr_tch,[],[],[],events_t,event_align,info,hdlfigtrial,titlestr,[],[]);
                    [rangefr ~]=plot_trials(trials_fr_tch/max(trials_fr_tch),[],[],[],[],[],info,hdlfigtrial,titlestr,[],[]);
                    
                    
                    %lfp
                    info.color=colorlist(4,:);info.nchannels=1;
                    %plot_trials(trials_lfp_tch,[],[],[],[],[],info,hdlfigtrial,titlestr,[],[]);
                    [rangelfp ~]=plot_trials(trials_lfp_tch/(-min(trials_lfp_tch)),[],[],[],[],[],info,hdlfigtrial,titlestr,[],[]);
                    
                    %axis
                    ylabel(['Channel ' num2str(ch)]);
                    set(gca,'ytick',[],'yticklabel',[]);
                end
                
                %pause
                %lfp latency trials
                shoulder_lfp=[];lat_lfp=[];lfp_minmax=[];
                if ~isempty(trials_lfp_tch) & ~isnan(trials_lfp_tch)
                    winmin=[-10 30];winsize=20;alpha=0.01;H1_count=10;
                    [shoulder_lfp lat_lfp lfp_minmax]=get_latency_trialch_backspline(trials_lfp_tch,winmin,winsize,'lfp',info,alpha,H1_count);
                    
                    %correction for timing
                    shoulder_lfp(1)=shoulder_lfp(1)-info.aligntime;
                    lat_lfp(1)=lat_lfp(1)-info.aligntime;
                    lfp_minmax(1)=lfp_minmax(1)-info.aligntime;
                    
                    if disp
                        rangeplot=[rangefr(1) rangefr(2) min([rangefr(3) rangelfp(3)]) max([rangefr(4) rangelfp(4)])];
                        
                        shoulder_plot=zeros(1,2);shoulder_plot(1,1)=shoulder_lfp(1);
                        plot_events_ch(shoulder_plot,[],[],rangeplot,info,hdlfigtrial,'n','-',3,colorlist(4,:));
                        lat_plot=zeros(1,2);lat_plot(1,1)=lat_lfp(1,1);
                        plot_events_ch(lat_plot,[],[],rangeplot,info,hdlfigtrial,'n','--',2,'k');
                        lfpmin_plot=zeros(1,2);lfpmin_plot(1,1)=lfp_minmax(1);
                        plot_events_ch(lfpmin_plot,[],1,rangeplot,info,hdlfigtrial,'n','--',2,'k');
                    end
                    
                end
                
                %display pburst
                if disp
                    if ~isempty(trials_fr_tch) & ~isnan(trials_lfp_tch)
                        b_begin=data(lut_trials_targ(t)).offline.targ_pburst_trial.b_begin;
                        ch_align=data(1).offline.targ_pburst_ch_align.ch;
                        pburst_plot=zeros(1,2);pburst_plot(1,1)=b_begin(ch)-b_begin(ch_align);
                        plot_events_ch(pburst_plot,[],[],rangeplot,info,hdlfigtrial,'n','-',2,colorlist(1,:));
                    end
                end
                
                
                %latencies list
                if ~isempty(shoulder_lfp)
                    latencies(ch)=shoulder_lfp(1);
                end
               
            end
            
            %%
            %Update data
            switch info.align
                case 'targ_pburst_ch'
                    field=['targ_pburst_lfplatency'];
                    eval(['data(' num2str(lut_trials_targ(t)) ').offline.' field '_trial=latencies;']);
                    %eval(['data(' num2str(lut_trials_targ(t)) ').offline.' field '_trial']);
                    %eval(['data(' num2str(lut_trials_targ(t)) ').offline']);
                case 'targ_rsburst_ch'
                    field=['targ_rsburst_lfplatency'];
                    eval(['data(' num2str(lut_trials_targ(t)) ').offline.' field '_trial=latencies;']);
            end
            
            %pause
            
            if disp
            pause;
            close(figtrial)
            end
        end
    end
    
    
    %save updated data
    display('NOT SAVED!')
    %update_data(1,0,0,data,data_path,info.datafile,[],[]);
    
    %pause
    
end



