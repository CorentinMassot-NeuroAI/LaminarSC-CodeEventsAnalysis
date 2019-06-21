%function compute_CSDfeature

%function compute_CSDfeature
%   Compute CSD feature from laminar data
%
% see also compute_tuning
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 10/15/2016 last modified 01/22/2017
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
alignlist={'targ' 'sacc'};

%window of analysis
wind_targ=[-10 340];
wind_sacc=[-100 250];


%sigma FR
sigma_FR=6;

%vshift
vshift=10;%28.6944;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get data
datalist=load_data_gandhilab(data_path);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analyzing data
dlist=get_dlist

hdlfigallvmis=figure;hold on;
data=[];
info=[];
for d=dlist
    %get data and info
    info.datafile=datalist{d};
    load ([data_path info.datafile]);
    display(info.datafile)
    
    %getting channel mapping and discard selected bad channels
    discard=data(1).offline.discardlfp;
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
    seltrials=get_seltrials(data,'rpt');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Data aligned on target and saccade onset
    for al=1:numel(alignlist)
        info.align=alignlist{al};

        switch info.align
            case 'targ'
                [alltrials_lfp_targ aligntime_targ]=get_alltrials_align(data,seltrials,wind_targ,'lfp',info,targslist,sigma_FR,0);
            case 'sacc'
                [alltrials_lfp_sacc aligntime_sacc]=get_alltrials_align(data,seltrials,wind_sacc,'lfp',info,targslist,sigma_FR,0);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get spk data for tuning
    [alltrials_spk_tuning info.aligntime]=get_alltrials_align(data,seltrials,[],'fr',info,targslist,sigma_FR,0);
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %analysis of trials for each target
    for tg=targ_tuning;%targs_ind        
        figtrials=figure('Position',[1 100 scrsz(3)-100 scrsz(4)-200]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %display all targets
        hdlfig=subplot(2,3,1);hold on;
        display_alltargets(targslist,info,hdlfig);
        
        %compute target tuning
        hdlfig=subplot(2,3,4);hold on;
        plot_targtuning(alltrials_spk_tuning,targs_ind,info,hdlfig,'Target tuning');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %target index
        info.targ=tg;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %plot lfp of targ and sacc
        for al=1:numel(alignlist)
            info.align=alignlist{al};
            switch info.align
                case 'targ'
                    trials_lfp=alltrials_lfp_targ{tg};
                    info.aligntime=aligntime_targ;
                case 'sacc'
                    trials_lfp=alltrials_lfp_sacc{tg};
                    info.aligntime=aligntime_sacc;
            end
            
            
            [info.nchannels info.ntrials info.triallen]=size(trials_lfp);
            %compute average trials
            [trials_lfp_avg trials_lfp_var]=get_trials_avg(trials_lfp);
            %remove trials with amplitude that is too small
            [trials_lfp_avgc index_lfp_c]=clean_trials(trials_lfp_avg,'lfp');
            hdlfig=subplot(2,3,al+1);hold on;
            titlestr={info.datafile ; ['LFP ' info.align ' t' num2str(info.targ) ' #trials:' num2str(info.ntrials)]};
            range=plot_trials(trials_lfp_avgc,[],index_lfp_c,vshift,[],[],info,hdlfig,titlestr,'-',1);                  
            %plot_event(wind_vmi,info.aligntime,range,hdlfig);
        
            %%
            %%%%%%%%%%%%%%%%%%
            %plot CSD
            hdlfig=subplot(2,3,al+4);hold on;
            titlestr='CSD';
            
            [csd zs clim]=plot_csdtrials('lfp',trials_lfp_avgc,[],index_lfp_c,[],[],[],info,hdlfig,titlestr);
            %find feature in CSD
            zs=zs*1e-3;
            [csdfeat ~]=findfeatures_CSD(1,1,0,csd,zs,info,hdlfig);     
            figure(figtrials)
            %hdlrect=rectangle('Curvature', [1 1],'Position', [peak.t-1 peak.d-1 2 2]);
            %set(hdlrect,'linewidth',2);
            %line([aligntime aligntime] ,[0 length(zs)]);
        
            %%
            %%%%%%%%%%%%%%%%%%
            %Update data
            %first paper (functional analysis of SC)
            %data=update_data(0,1,0,data,data_path,info.datafile,['csdfeat_avg_' info.align],csdfeat);

            %second paper _p2 (LFP/CSD analysis of SC)
            data=update_data(0,1,0,data,data_path,info.datafile,['csdfeat_avg_' info.align '_p2'],csdfeat);

        end
        
        
        %pause
        %saving updated data
        data=update_data(0,1,0,data,data_path,info.datafile,['csdzs'],zs);
        update_data(1,0,0,data,data_path,info.datafile,[],[]);
        %display('NOT SAVED!')
        close(figtrials)
    
    end
end

