%function savedata_events_trials

%function savedata_events_trials
%   save data trial-by-trial recorded with a laminar probe (LMA)
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 02/10/2017 last modified 02/10/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set paths
[root_path data_path save_path]=set_paths;

%screen size
scrsz = get(groot,'ScreenSize');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
%print figures and save data
savedata=0;
savefigs=0;
figtype='epsc2';%'png';%'epsc2';

%alignement
%alignlist={'no' 'targ' 'go' 'sacc'};
alignlist={'no'};

%window of analysis
wind=[];%all
%wind=[-100 600];%targ align
%wind=[-600 250];%sacc align

%sigma FR
sigma_FR=6;

%data_save directory
data_save_dir='Data_SC_Sanjeev\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get data
datalist=load_data_gandhilab(data_path);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analyzing data
dlist=get_dlist;

data=[];
info=[];
for d=dlist
    
    %get data and info
    clear('data');clear('data_save');
    
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
    
    %     %target, target index and target in anti-RF
    %     targs_ind_flip=fliplr(targs_ind);
    %     info.targ_tuning_a=targs_ind_flip(find(targs_ind==info.targ_tuning));
    
    
    %select trials
    seltrials=get_seltrials(data,'rpt');
    
    %new data
    data_aux=data(seltrials);
    
    
    
    %%%%%%%%%%%%%%%%%%
    %new structure
    %Uday's structure for Sanjeev: There are 9 datasets, each in the variable 'data', with the following fields:
    data_save={};
    
    for t=1:numel(data_aux)
        % trialSpikeTimestamps: cell(ntrials,nchannels), with spike times on each trial per channel, referenced to beginning of the trial
        for ch=1:info.nchannels
            data_save.trialSpikeTimestamps{t,ch}=data_aux(t).spikeTimestamps{ch};
        end
        
        % trialTargetLocations: double(ntrials,2), with x,y target locations on each trial
        data_save.trialTargetLocations(t,:)=data_aux(t).targets.window(2,1:2);
        
        % trialEventTimes: double(ntrials,8), with corresponding event times for each trial
        data_save.trialEventTimes(t,:)=data_aux(t).stateTransitions(2,:);
        
        % trialGazePosition: cell(ntrials,1), with H and V gaze positions for each trial
        data_save.trialGazePosition{t,1}=data_aux(t).gazePosition;
        
        % trialSaccTime: double(ntrials,1), corresponding to sacc time for each trial
        data_save.trialSaccTime(t,1)=data_aux(t).behavrpt.saccTime;
        
    end
    
    % targetIN: double(1,2), with the x,y of the receptive field target. You can use this to sort the other data based on IN or OUT (ipsi/contra) condition
    data_save.targetIN(1,:)=targslist(info.targ_tuning,:);
    %data_save.targetOUT(1,:)=targslist(info.targ_tuning_a,:);
    
    
    % EventNames: cell(1,8), with labels indicating relevant events during the trial
    data_save.EventNames(1,:)=data_aux(1).params;
    
    % neuronType: cell(1,nchannels), indicating which channels I think there are meaningful neurons in according to my sorting procedure (combination of automatic firing rate criteria and manual detection), and their type (V/M/VM) classifications. Feel free to use your own procedures if you think consistency is needed.
    for ch=1:info.nchannels
        targ_b=data_aux(1).offline.targ_bsignif(ch);
        sacc_b=data_aux(1).offline.sacc_bsignif(ch);
        
        if targ_b==1 & sacc_b==1
            data_save.neuronType{1,ch}='VM';
        elseif targ_b==1
            data_save.neuronType{1,ch}='V';
        elseif sacc_b==1
            data_save.neuronType{1,ch}='M';
        else
            data_save.neuronType{1,ch}='';
        end            
        
        targ_t=data_aux(1).offline.targ_bthresh(ch);
        sacc_t=data_aux(1).offline.sacc_bthresh(ch);
        
        if (targ_b==1 & targ_t>=1) | (sacc_b==1 & sacc_t>=1)
            data_save.goodChannel(1,ch)=1;
        else
            data_save.goodChannel(1,ch)=0;
        end            
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %save
    save_data(data_save,root_path,data_save_dir,info);
    %pause
end







%     %MISC
%     %loop across all alignements
%     for al=1:numel(alignlist)
%         info.align=alignlist{al};
%
%         %get all neural and behavioral data with specific alignement
%         [alltrials_spk_tuning,info.aligntime]=get_alltrials_align(data,seltrials,[],'fr',info,targslist,sigma_FR,1);
%         [alltrials_spk,info.aligntime]=get_alltrials_align(data,seltrials,wind,'fr',info,targslist,sigma_FR,1);
%         [alltrials_lfp,info.aligntime]=get_alltrials_align(data,seltrials,wind,'lfp',info,targslist,sigma_FR,1);
%         [allgazepos,allevents]=get_alldatagaze_align(data,seltrials,info,targslist);
%
%
%
%         %%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %analysis of trials for each target
%         for tg=info.targ_tuning;%targs_ind,
%
%             %neural and behavioral signals for target tg
%             trials_spk=alltrials_spk{tg};
%             [info.nchannels,info.ntrials,info.triallen]=size(trials_spk);
%             trials_lfp=alltrials_lfp{tg};
%             [info.nchannels,info.ntrials,info.triallen]=size(trials_lfp);
%             gazepos=allgazepos{tg};
%             events=allevents{tg};
%
%             %loop on trials
%             for t=1:info.ntrials,
%
%                 figtrials=figure('Position',[1 100 scrsz(3)-100 scrsz(4)-200]);
%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %display all targets
%                 hdlfig=subplot(2,3,1);hold on;
%                 display_alltargets(targslist,info,hdlfig);
%
%                 %compute target tuning
%                 hdlfig=subplot(2,3,3);hold on;
%                 plot_targtuning(alltrials_spk_tuning,targs_ind,info,hdlfig,'Target tuning');
%
%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %target index
%                 info.targ=tg;
%
%
%                 %%
%                 %%%%%%%%%%%%%%%%%%
%                 %gaze data
%                 gazepos_t=gazepos{t};
%                 events_t=events{t};
%                 event_align=get_eventalign(events_t,info.align);
%
%                 if ~isempty(event_align)
%
%                     hdlfig1=subplot(2,3,1);hold on;
%                     hdlfig2=subplot(2,3,4);hold on;
%                     events_t.peak=plot_gazedata(gazepos_t,events_t,event_align,wind,info,hdlfig1,'',hdlfig2,'XY Eye Traces');
%
%                     %%
%                     %%%%%%%%%%%%%%%%%%
%                     %spk
%                     trials_spk_t=squeeze(trials_spk(:,t,:));
%                     %remove trials with amplitude that is too small
%                     [trials_spk_tc index_spk_tc]=clean_trials(trials_spk_t,'fr');
%                     hdlfig=subplot(2,3,2);hold on;
%                     titlestr={info.datafile ; ['FR ' info.align ' t' num2str(info.targ) ' trial:' num2str(t) '/' num2str(info.ntrials)]};
%                     plot_trials(trials_spk_tc,[],index_spk_tc,[],events_t,event_align,info,hdlfig,titlestr);
%                     %plot_trials(trials_spk_t,[],[],[],events_t,event_align,info,hdlfig,titlestr);
%
%
%                     %%%%%%%%%%%%%%%%%%
%                     %lfp
%                     trials_lfp_t=squeeze(trials_lfp(:,t,:));
%                     %remove trials with amplitude that is too small
%                     [trials_lfp_tc index_lfp_tc]=clean_trials(trials_lfp_t,'lfp');
%                     hdlfig=subplot(2,3,5);hold on;
%                     titlestr='LFP';
%                     plot_trials(trials_lfp_tc,[],index_lfp_tc,[],events_t,event_align,info,hdlfig,titlestr);
%

%                        pause
%                 end
%
%                 close(figtrials)
%             end
%             pause
%
%             end
%             end
