%function analysis_events

%function analysis_events
%   Analysis of data averaged over trials recorded with a laminar probe (LMA)
%
%
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 11/03/2015 last modified 01/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set paths
[root_path data_path save_path]=set_paths;

%screen size
scrsz = get(groot,'ScreenSize');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
%print figures
savefigs=0;
figtype='epsc2';%'png';%'epsc2';

%alignement
%alignlist={'no' 'targ' 'go' 'sacc' 'targ_pburst' 'targ_rsburst' 'sacc_pburst' 'sacc_rsburst'};
%alignlist={'targ_pburst_ch' 'targ_rsburst_ch' 'targ' };
%alignlist={'targ'};
alignlist={'targ_pburst_ch'};
%alignlist={'sacc' };


%window of analysis
%paper Functional organization
%wind=[];%all
%wind_targ=[-50 300];%targ align
wind_targ=[-150 260];%targ align
%wind_targ=[-50 100];%targ align
%wind_sacc=[-500 200];%sacc align

%paper LFP/CSD
%wind_targ=[0 200];%targ align


wind_targ_pburst=[-50 350];%[-50 300];%targ_pburst_ch align
%wind_sacc=[-200 0];%[-200 -50];%sacc align buildup
wind_sacc=[-50 100];%[-50 10];%sacc align burst
%wind_sacc=[-10 150];%sacc align transsaccadic


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

%shift
shift_spk=50;%10;%50;%50%80;%100;%78.9386;%[]
shift_lfp=30;%280000000;%30;%28.6944;%[]

%newdata directory
savedata=0;
newdata_dir='Data_SC_Joy\';

%shift ripple temporal correction
shift_ripple=4

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
    clear('data')
    info.datafile=datalist{d};
    load ([data_path info.datafile]);
    display(info.datafile)
    
    %getting channel mapping and discard selected bad channels
    [info.chmap info.nchannels info.depths]=get_chmap(data(1).info.electrode{2},[]);
    %[info.chmap info.nchannels info.depths]=get_chmap([9:16  1:8],[]);
     
     
    %getting trial type
    info.trialtype=data(1).sequence(1);
    %getting list of targets
    targslist=data(1).offline.targslist;
    %targets index
    targs_ind=get_targsindex(targslist,info);
    %%target tuning (after compute_tuning)
    info.targ_tuning=data(1).offline.targ_tuning;
    
    %select trials
    seltrials=get_seltrials(data,'rpt');
    %seltrials_d=get_seltrials(data,'rpt');
    
    %select trials according to features: srt trend repeat
    %seltrials_f=get_seltrials_features(data(seltrials_d),[0 10000],[],0);
    %seltrials_f=get_seltrials_features(data(seltrials_d),[0 10000],[]);
    %seltrials_f=get_seltrials_features(data(seltrials_d),[200 290],[]);
    %seltrials_f=get_seltrials_features(data(seltrials_d),[350 500],[],0);
    %seltrials_f=get_seltrials_features(data(seltrials_d),[300 340],[],1);
    %seltrials_f=get_seltrials_features(data(seltrials_d),[250 300],[],0);
    %seltrials_f=get_seltrials_features(data(seltrials_d),[300 350],[],0);
    %seltrials_f=get_seltrials_features(data(seltrials_d),[200 400],[],0);
  
    
    %seltrials_f=get_seltrials_features(data(seltrials_d),[0 100],[],0);
    %seltrials_f=get_seltrials_features(data(seltrials_d),[40 80],[],0);
    %seltrials_f=get_seltrials_features(data(seltrials_d),[40 50],[],0);
    %seltrials_f=get_seltrials_features(data(seltrials_d),[50 80],[],0);
    
    %seltrials_f=get_seltrials_features(data(seltrials_d),[200 300],[],0);
    %seltrials_f=get_seltrials_features(data(seltrials_d),[200 400],[],0);

    %seltrials=seltrials_d(seltrials_f);
    
    %loop across all alignements
    for al=1:numel(alignlist)
        info.align=alignlist{al};
        
        %get alltrials with specific alignement
        [alltrials_spk_tuning info.aligntime ~]=get_alltrials_align(data,seltrials,[],'fr',info,targslist,sigma_FR,1);
        %[alltrials_spk info.aligntime ~]=get_alltrials_align(data,seltrials,wind,'fr',info,targslist,sigma_FR,1);
        %[alltrials_lfp ~]=get_alltrials_align(data,seltrials,wind+shift_ripple,'lfp',info,targslist,sigma_FR,1);
        [allstats]=get_allbehavstats(data,seltrials,targslist,'rpt');
        
        
         switch info.align
                case 'targ'
                    [alltrials_spk info.aligntime ~]=get_alltrials_align(data,seltrials,wind_targ,'fr',info,targslist,sigma_FR,1);
                    [alltrials_spk_bsl ~]=get_alltrials_align(data,seltrials,wind_targ_bsl,'fr',info,targslist,sigma_FR,0);
                    
                    [alltrials_lfp info.aligntime ~]=get_alltrials_align(data,seltrials,wind_targ+shift_ripple,'lfp',info,targslist,sigma_FR,1);
                    [alltrials_lfp_bsl ~]=get_alltrials_align(data,seltrials,wind_targ_bsl+shift_ripple,'lfp',info,targslist,sigma_FR,0);
                    
                case 'targ_pburst_ch'
                    [alltrials_spk info.aligntime ~]=get_alltrials_align(data,seltrials,wind_targ_pburst,'fr',info,targslist,sigma_FR,1);
                    [alltrials_spk_bsl ~]=get_alltrials_align(data,seltrials,wind_targ_pburst_bsl,'fr',info,targslist,sigma_FR,0);
                    
                    [alltrials_lfp info.aligntime ~]=get_alltrials_align(data,seltrials,wind_targ_pburst+shift_ripple,'lfp',info,targslist,sigma_FR,1);
                    [alltrials_lfp_bsl ~]=get_alltrials_align(data,seltrials,wind_targ_pburst_bsl+shift_ripple,'lfp',info,targslist,sigma_FR,0);
                    
                case 'sacc'
                    [alltrials_spk info.aligntime ~]=get_alltrials_align(data,seltrials,wind_sacc,'fr',info,targslist,sigma_FR,1);
                    [alltrials_lfp info.aligntime ~]=get_alltrials_align(data,seltrials,wind_sacc+shift_ripple,'lfp',info,targslist,sigma_FR,1);
                    
                    alignaux=info.align;
                    info.align='go';
                    [alltrials_spk_bsl ~]=get_alltrials_align(data,seltrials,wind_sacc_bsl_go,'fr',info,targslist,sigma_FR,0);
                    [alltrials_lfp_bsl ~]=get_alltrials_align(data,seltrials,wind_sacc_bsl_go+shift_ripple,'lfp',info,targslist,sigma_FR,0);
                    info.align=alignaux;
         end
            
        
        %save data
        newdata={};
        
        %% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %analysis of trials for each target
        for tg=info.targ_tuning %targs_ind %[6:10],%info.targ_tuning%
        
            figtrials=figure('Position',[1 100 scrsz(3)-100 scrsz(4)-200]);
            figtrials2=figure('Position',[1 100 scrsz(3)-100 scrsz(4)-200]);
            figtrials3=figure('Position',[1 100 scrsz(3)-100 scrsz(4)-200]);
         
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %display all targets
            figure(figtrials)
            hdlfig=subplot(2,3,1);hold on;
            display_alltargets(targslist,info,hdlfig);
            
%              %compute target tuning
%              figure(figtrials)
%              hdlfig=subplot(2,3,4);hold on;
%              plot_targtuning(alltrials_spk_tuning,targs_ind,info,hdlfig,'Target tuning');
         
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %target index
            info.targ=tg;
               
            %%
            %%%%%%%%%%%%%%%%%%
            %spk
            trials_spk=alltrials_spk{tg};
            trials_spk_bsl=alltrials_spk_bsl{tg};
            [info.nchannels info.ntrials info.triallen]=size(trials_spk);
            
            %             %compute average trials
            %             [trials_spk_avg trials_spk_var]=get_trials_avg(trials_spk);
            %
            %             %remove trials with amplitude that is too small
            %             [trials_spk_avgc index_spk_c]=clean_trials(trials_spk_avg,'fr');
            
            
            %compute baseline-corrected average trials 
            trials_spk_n=get_trials_normalized(trials_spk,trials_spk_bsl,'FR',info);
            [trials_spk_n_avg trials_spk_n_var]=get_trials_avg(trials_spk_n);
            
            
            figure(figtrials)
            hdlfig=subplot(2,3,2);hold on;
            titlestr={info.datafile ; ['FR ' info.align ' t' num2str(info.targ) ' #trials:' num2str(info.ntrials)]};
            plot_trials(trials_spk_n_avg,[],[],shift_spk,[],[],info,hdlfig,titlestr,[],[]);
           
            
            figure(figtrials2);hold on;
            plot_trials(trials_spk_n_avg,[],[],0,[],[],info,hdlfig,titlestr,[],[]);
            
            
            %%
            %%%%%%%%%%%%%%%%%%
            %lfp
            trials_lfp=alltrials_lfp{tg};
            trials_lfp_bsl=alltrials_lfp_bsl{tg};
            [info.nchannels info.ntrials info.triallen]=size(trials_lfp);
            
            
            %             %%
            %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             %HACK to use CSD
            %             %NOTE: compute csd for each trial first and then avg in order to compute latency
            %             %get csd channels in RF
            %             display('USING CSD SIGNALS!!!')
            %             trials_csd=[];
            %             for tcsd=1:info.ntrials
            %                 trials_lfp_aux=squeeze(trials_lfp(:,tcsd,:));
            %                 [csd zs]=get_csdtrials(trials_lfp_aux,[1:length(info.chmap)],info);
            %                 %get csd channels
            %                 [trials_csd_aux depths_ch]=get_csdchannels(csd,info);%TO DO compare with second derivative of LFP
            %                 trials_csd(:,tcsd,:)=trials_csd_aux;
            %              end
            %             %WARNING
            %             trials_lfp=trials_csd;
            %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            
            
            
            %             %compute average trials
            %             [trials_lfp_avg trials_lfp_var]=get_trials_avg(trials_lfp);
            %
            %           
            
            %compute baseline-corrected average trials 
            trials_lfp_n=get_trials_normalized(trials_lfp,trials_lfp_bsl,'lfp',info);
            [trials_lfp_n_avg trials_spk_n_var]=get_trials_avg(trials_lfp_n);
            
            %remove trials with amplitude that is too small
            [trials_lfp_n_avgc index_lfp_c]=clean_trials(trials_lfp_n_avg,'lfp');
            
            
            figure(figtrials)
            hdlfig=subplot(2,3,5);hold on;
            titlestr='LFP';
            plot_trials(trials_lfp_n_avgc,[],index_lfp_c,shift_lfp,[],[],info,hdlfig,titlestr,[],[]);
            
            
            figure(figtrials3);hold on;
            plot_trials(trials_lfp_n_avgc,[],index_lfp_c,0,[],[],info,hdlfig,titlestr,[],[]);
            
            
            
            %%
            %%%%%%%%%%%%%%%%%%
            %plot CSD
            
            figure(figtrials)
            hdlfig=subplot(2,3,6);hold on;
            titlestr='CSD';
            %plot_csdtrials(trials_lfp_n_avg,[1:info.nchannels],[],[],[],info,hdlfig,titlestr);
            plot_csdtrials('lfp',trials_lfp_n_avgc,[],index_lfp_c,[],[],[],info,hdlfig,titlestr);
            
            
            %ch_ref
            %alignment of onset using CSD features (after compute_CSDfeature)
            info.csdfeat_avg_targ=data(1).offline.csdfeat_avg_targ;
            info.zs=data(1).offline.csdzs;
            dref=info.csdfeat_avg_targ(2);
            
            [aux info_r ch_ref dref_conv]=get_data_aligndepth(zeros(info.nchannels,1),dref,info,[]);
            display(dref)
            display(ch_ref)
            display(dref_conv)
            
            %%
            %%%%%%%%%%%%%%%%%%
            %behavioral stats
            stats_t=allstats{tg};
            figure(figtrials)
            hdlfig=subplot(2,3,3);hold on;
            titlestr='Behavioral stats';
            plot_behavstats(stats_t,info,hdlfig,titlestr);
            
            
            %%%%%%%%%%%%%%%%%%
            %save figs
            if savefigs
                saveas(figtrials,[save_path info.datafile '_' info.align '_t' num2str(info.targ) '.' figtype],figtype);
            end;
            
            %%
            %%%%%%%%%%%%%%%%%%
            %save data
            if savedata
                newdata{tg}.lfp=trials_lfp_avgc;
                newdata{tg}.fr=trials_spk_avgc;
                newdata{tg}.info=info;
            else
                d
                pause
                close all
            end
        end
        
        %save data
        if savedata
            save_data(newdata,root_path,newdata_dir,info);
        end
        
    end
    %close all
end



