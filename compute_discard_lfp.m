%function compute_discard_lfp

%function compute_discard_lfp
%   create list of discarded channels with bad lfp
%
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 08/12/2018 last modified 08/12/2018
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
%alignlist={'targ'};
alignlist={'targ_pburst_ch' };
%alignlist={'sacc' };


%window of analysis
wind_targ=[-10 340];
wind_targ_pburst=[-50 300];%targ_pburst_ch align
%wind_sacc=[-200 0];%%pre saccadic
%wind_sacc=[-50 100];%[-50 10];%movement burst

%bsl windows (see also compute_vmi)
%targ
wt=100;%50
wind_targ_bsl=[50-wt 50];
wind_targ_pburst_bsl=[-50-wt -50 ];%[30-wt 30];
%sacc
ws=50;
wind_sacc_bsl_go=[-ws 0];%[50-ws 50];%



%sigma FR
sigma_FR=6;

%vshift
vshift=30;

%shift ripple temporal correction
shift_ripple=4


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get data
datalist=load_data_gandhilab(data_path);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analyzing data
dlist=get_dlist

data=[];info=[];
allcsd={};allcsdbsl={};
dd=0;
for d=dlist(1:end)
    %counter
    dd=dd+1;
    %get data and info
    info.datafile=datalist{d};
    load ([data_path info.datafile]);
    display(info.datafile)
    
    %getting channel mapping and discard selected bad channels
    discard=[];%11;%[]
    [info.chmap info.nchannels info.depths]=get_chmap(data(1).info.electrode{2},discard);
    %getting trial type
    info.trialtype=data(1).sequence(1);
    %getting list of targets
    targslist=data(1).offline.targslist;
    %targets index
    targs_ind=get_targsindex(targslist,info);
    
    %target tuning (after compute_tuning)
    targ_tuning=data(1).offline.targ_tuning;
    
    
    %select trials
    %seltrials=get_seltrials(data,'rpt');
    seltrials_d=get_seltrials(data,'rpt');
    
    %select trials according to features: srt trend repeat
    seltrials_f=get_seltrials_features(data(seltrials_d),[200 400],[],0);
    
    seltrials=seltrials_d(seltrials_f);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Data aligned on target and saccade onset
    for al=1:numel(alignlist)
        info.align=alignlist{al};
        
        switch info.align
            case 'targ'
                [alltrials_lfp info.aligntime_w ~]=get_alltrials_align(data,seltrials,wind_targ+shift_ripple,'lfp',info,targslist,sigma_FR,1);
                [alltrials_lfp_bsl info.aligntime_bsl ~]=get_alltrials_align(data,seltrials,wind_targ_bsl+shift_ripple,'lfp',info,targslist,sigma_FR,0);
                
            case 'targ_pburst_ch'
                [alltrials_lfp info.aligntime_w ~]=get_alltrials_align(data,seltrials,wind_targ_pburst+shift_ripple,'lfp',info,targslist,sigma_FR,1);
                [alltrials_lfp_bsl info.aligntime_bsl ~]=get_alltrials_align(data,seltrials,wind_targ_pburst_bsl+shift_ripple,'lfp',info,targslist,sigma_FR,0);
                
            case 'sacc'
                [alltrials_lfp info.aligntime_w ~]=get_alltrials_align(data,seltrials,wind_sacc+shift_ripple,'lfp',info,targslist,sigma_FR,1);
                
                alignaux=info.align;
                info.align='go';
                [alltrials_lfp_bsl info.aligntime_bsl ~]=get_alltrials_align(data,seltrials,wind_sacc_bsl_go+shift_ripple,'lfp',info,targslist,sigma_FR,0);
                info.align=alignaux;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %analysis of trials for each target
        tg=targ_tuning;%targs_ind
        
        %target index
        info.targ=tg;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %plot lfp
        %figures
        figcsd=figure('Position',[1 100 scrsz(3)-100 scrsz(4)-200]);
        
        trials_lfp=alltrials_lfp{tg};
        trials_lfp_bsl=alltrials_lfp_bsl{tg};
        info.aligntime=info.aligntime_w;
        maxplot=[];
        
        [info.nchannels info.ntrials info.triallen]=size(trials_lfp);
        
        
        %compute average trials
        trials_lfp_n=get_trials_normalized(trials_lfp,trials_lfp_bsl,'lfp',info);
        [trials_lfp_avg trials_spk_var]=get_trials_avg(trials_lfp_n);
        
        figure(figcsd);hdlfig=subplot(1,3,1);hold on;
        titlestr={info.datafile ; ['LFP' info.align ' t' num2str(info.targ) ' #t' num2str(info.ntrials)]};
        range=plot_trials(trials_lfp_avg,[],[],vshift,[],[],info,hdlfig,titlestr,[],[]);
        
        
        %%
        %%%%%%%%%%%%%%%%%%
        %select discarded channels
        discardlist=[];
        ok=0;
        while ~ok
            discardlist = str2num(input('List of discarded LFP channels:','s'));
            
            %%
            %%%%%%%%%%%%%%%%%%
            %plot LFP and CSD with discarded channels
            aux=[1:info.nchannels];
            aux(discardlist)=nan;
            selectlist=aux(find(~isnan(aux)))
            
            figcsd2=figure('Position',[1 100 scrsz(3)-100 scrsz(4)-200]);
            hdlfig=subplot(1,3,2);hold on;
            titlestr='LFP with discarded channels';
            range=plot_trials(trials_lfp_avg(selectlist,:),[],selectlist,vshift,[],[],info,hdlfig,titlestr,[],[]);
            
            hdlfig=subplot(1,3,3);hold on;
            titlestr='CSD with discarded channels';
            [csd zs maxplot]=plot_csdtrials('lfp',trials_lfp_avg(selectlist,:),[],selectlist,[],[],maxplot,info,hdlfig,titlestr);
            
            %ok?
            display('list of discarded LFP channels: ', num2str(discardlist));
            ok = str2num(input('Are you ok with list of discarded LFP channels? 1/0:','s'));
            close(figcsd2)
            
        end
        
        %%
        %%%%%%%%%%%%%%%%%%
        %Update data
        %data=update_data(0,1,0,data,data_path,info.datafile,['discardlfp' info.align],csdfeat);
        data=update_data(0,1,0,data,data_path,info.datafile,['discardlfp'],discardlist);
        
        d
        dd
        %pause
        close all
        
        
    end
    
    %%%%%%%%%%%%%%%%%%
    %saving updated data
    update_data(1,0,0,data,data_path,info.datafile,[],[]);
    %display('NOT SAVED!')
    %pause
end


