%function compute_pburst

%function compute_pburst
%   Compute P value of burst from laminar data
%   based on Hanes et al. 1995 (P_Burst.m)
%
% see also compute_bsignif
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 09/27/2017 last modified 09/27/2017
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
%alignlist={'no' 'targ' 'go' 'sacc'};
alignlist={'targ'};
%alignlist={'sacc'};

%windows of analysis
%wind_no=[-200 1500];
wind_targ=[-100 300];
wind_sacc=[-300 100];%[-200 200];

%window of burst
windb_targ=[30 150];
windb_sacc=[-100 0];

%vshift
vshift_fr=50;

%sigma FR
sigma_FR=6;

%burst threshold
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
    %targets index
    targs_ind=get_targsindex(targslist,info);
    
    %target tuning (after compute_tuning)
    targ_tuning=data(1).offline.targ_tuning;
    
    %select trials
    seltrials=get_seltrials(data,'rpt');
    %seltrials=[1:numel(data)];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Data aligned on target and saccade onset
    for al=1:numel(alignlist)
        info.align=alignlist{al};
        
        switch info.align
            case 'targ'
                wind=wind_targ;
                windb=windb_targ;
            case 'sacc'
                wind=wind_sacc;
                windb=windb_sacc;
        end
        
        [alltrials_fr info.aligntime ~]=get_alltrials_align(data,seltrials,wind,'fr',info,targslist,sigma_FR,1);
        [alltrials_spk info.aligntime lut_trials]=get_alltrials_align(data,seltrials,wind,'spk',info,targslist,sigma_FR,1);
        [allgazepos,allevents]=get_alldatagaze_align(data,seltrials,info,targslist);
        

        %analysis of trials for each target
        ratio=zeros(2,length(targs_ind),info.nchannels);msurprises=zeros(2,length(targs_ind),info.nchannels);
        ch_align_pburst={};
        for tg=targ_tuning%targs_ind%
            
            info.targ=tg;
            trials_fr=alltrials_fr{info.targ};
            trials_spk=alltrials_spk{info.targ};
            
            gazepos=allgazepos{info.targ};
            events=allevents{info.targ};
            
            lut_trials_targ=lut_trials{info.targ};
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Burst analysis and plot of each trial
            [info.nchannels info.ntrials info.triallen]=size(trials_fr);
            
            %plot each trial
            titlestr={info.datafile ; ['FR ' info.align ' t' num2str(info.targ) ' #trials:' num2str(info.ntrials)]};
            allb_begin=zeros(2,info.ntrials,info.nchannels);allb_end=zeros(2,info.ntrials,info.nchannels);allb_surprise=zeros(2,info.ntrials,info.nchannels);
            for t=1:info.ntrials,
                
                %t
                %pause
                %display(['Processed:' num2str(round(100*t/info.ntrials)) '%'])
                
                if disp,figtrial=figure('Position',scrsz);hdlfig=subplot(1,1,1);hold on;end;
                
                trials_fr_t=squeeze(trials_fr(:,t,:));
                trials_spk_t=squeeze(trials_spk(:,t,:));
                
                events_t=events{t};
                event_align=get_eventalign(events_t,info.align);
                %plot
                if disp,
                    [range ~]=plot_trials(trials_fr_t,[],[1:info.nchannels],vshift_fr,events_t,event_align,info,hdlfig,titlestr,'-',1);
                    %[range ~]=plot_trials(trials_spk_t*vshift_fr,[],[1:info.nchannels],vshift_fr,events_t,event_align,info,hdlfig,titlestr,'-',1);
                end
   
                %%
                %measure Poisson surprise
                b_begin=zeros(1,16);b_end=zeros(1,16);b_surprise=zeros(1,16);
                for algo=1:2
                    for ch=1:info.nchannels
                        %timestamps
                        trials_ts=find(trials_spk_t(ch,:)>=1);%-info.aligntime;
                        
                        switch algo
                            case 1
                                %p_burst (Hanes et al.)
                                [bob,eob,sob]=p_burst(trials_ts,min(trials_ts),max(trials_ts));
                            
                            case 2
                                %rs_burst (Gourevich et al.)
                                if ~isempty(trials_ts)
                                    [sob,lob,bob]=rs_burst(trials_ts);
                                    eob=bob+lob-1;
                                else
                                    bob=[];eob=[];sob=[];
                                end
                        end
                        
                        if ~isempty(bob)
                            bobs=trials_ts(bob)-info.aligntime;
                            eobs=trials_ts(eob)-info.aligntime;
                            
                            %selection of burts in windb window with max surprise
                            %ind=[1:numel(bobs)]
                            ind=find(bobs>windb(1) & bobs<windb(2));
                            
                            %in case 2 candidate bursts
                            [v ind_max]=max(sob(ind));
                            ind=ind(ind_max);
                          
                            if ~isempty(ind)
                                %[val ind]=max(sobs);
                                bob_sel=bobs(ind);
                                eob_sel=eobs(ind);
                                sob_sel=sob(ind);
                                
                                %list
                                b_begin(ch)=bob_sel;
                                b_end(ch)=eob_sel;
                                b_surprise(ch)=sob_sel;
                                
                            else
                                %TO DO: replace by NaN
                                b_begin(ch)=0;
                                b_end(ch)=0;
                                b_surprise(ch)=0;
                                
                            end
                            
                        else
                            %TO DO: replace by NaN
                            b_begin(ch)=0;
                            b_end(ch)=0;
                            b_surprise(ch)=0;
                            
                        end
                    end
                    
                    
                    %%
                    %lists 
                    allb_begin(algo,t,:)=b_begin(:);
                    allb_surprise(algo,t,:)=b_surprise(:);
                    allb_end(algo,t,:)=b_end(:);
                    
                    %%
                    if disp
                        %plot burst events
                        figure(figtrial);
                        b_begin_plot=zeros(info.nchannels,2);
                        b_end_plot=zeros(info.nchannels,2);
                        b_begin_plot(:,1)=b_begin;
                        b_end_plot(:,1)=b_end;
                        if algo==1,
                            plot_events_ch(b_begin_plot,[],vshift_fr,range,info,hdlfig,'n','-',3,[]);
                            plot_events_ch(b_end_plot,[],vshift_fr,range,info,hdlfig,'n','-',3,[]);
                        elseif algo==2
                            plot_events_ch(b_begin_plot,[],vshift_fr,range,info,hdlfig,'n',':',5,[]);
                            plot_events_ch(b_end_plot,[],vshift_fr,range,info,hdlfig,'n',':',5,[]);
                        end
                    end
                    
                    
                    %%
                    %Update data
                    b_trial=struct;
                    b_trial.b_begin=b_begin;
                    b_trial.b_end=b_end;
                    b_trial.b_surprise=b_surprise;
                    if algo==1,
                        eval(['data(' num2str(lut_trials_targ(t)) ').offline.' info.align '_pburst_trial=b_trial;']);
                    elseif algo==2
                        eval(['data(' num2str(lut_trials_targ(t)) ').offline.' info.align '_rsburst_trial=b_trial;']);
                    end
                   
                    
                end
                %pause
                if disp
                    pause
                    close(figtrial)
                end
            end
             
             
            %%
            %compute ratio msurpises and ch_align_pburst
            begins=[];surprises=[];
            for algo=1:2
                for ch=1:info.nchannels
                    begins=allb_begin(algo,:,ch);
                    surprises=allb_surprise(algo,:,ch);
                    
                    %ratio of detected bursts
                    ratio(algo,tg,ch)=sum(abs(begins)>0)/length(begins);
                    
                    %msurprises
                    msurprises(algo,tg,ch)=nanmean(surprises(find(abs(begins)>0)));
                end
            end
        end
        %%
        %ch_align_pburst
        for algo=1:2
            ch_align=struct;
            [ch_align.ratio ch_align.ch]=max(ratio(algo,targ_tuning,:));
            ch_align.surprise=msurprises(algo,targ_tuning,ch_align.ch);
            ch_align_pburst{algo}=ch_align;
        end
     
        %squeeze(ratio(2,targ_tuning,:))'
        %squeeze(msurprises(2,targ_tuning,:))'
        %ch_align_pburst
  
        %%
        %update data
        for algo=1:2,
            if algo==1,
                data=update_data(0,1,0,data,data_path,info.datafile,[info.align '_pburst_ratio'],squeeze(ratio(1,:,:)));
                data=update_data(0,1,0,data,data_path,info.datafile,[info.align '_pburst_msurprises'],squeeze(msurprises(1,:,:)));
                data=update_data(0,1,0,data,data_path,info.datafile,[info.align '_pburst_ch_align'],ch_align_pburst{1});
           
            elseif algo==2,
                data=update_data(0,1,0,data,data_path,info.datafile,[info.align '_rsburst_ratio'],squeeze(ratio(2,:,:)));
                data=update_data(0,1,0,data,data_path,info.datafile,[info.align '_rsburst_msurprises'],squeeze(msurprises(2,:,:)));
                data=update_data(0,1,0,data,data_path,info.datafile,[info.align '_rsburst_ch_align'],ch_align_pburst{2});
              end
        end
        
        %save updated data
        display('NOT SAVED!')
        %update_data(1,0,0,data,data_path,info.datafile,[],[]);

        pause
    end
end

