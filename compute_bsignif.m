                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  %function compute_bsignif
%function compute_bsignif
%   Compute burst significance from laminar data
%
% see also compute_tuning
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 06/23/2016 last modified 06/23/2017
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
%alignlist={'targ_pburst_ch' 'sacc' };
%alignlist={'targ_rsburst_ch' 'sacc' };
%alignlist={'targ_pburst_ch' 'targ_rsburst_ch'};
alignlist={'targ' 'sacc' };
%alignlist={'sacc' };

%windows of analysis
%plot
wind_targ=[-50 450];
wind_sacc=[-300 200];%[-200 200];

% %to compute bsignif
% %target 20deg
% wind_targ_bsignif=[100 200];
% wind_targ_bsl=[-50 50];   
% wind_sacc_bsignif=[-25 75];
%
% %target 3deg
% wind_targ_bsignif=[100 150];
% wind_targ_bsl=[0 50];
% wind_sacc_bsignif=[-25 25];
%
% %wind_sacc_bsl=[-150 -100]

%vshift
vshift_spk=100;

%sigma FR
sigma_FR=6;

%burst trheshold
thresh=10;%threshold spk/s

%alpha test statistic
alpha=0.001
        
%thresholds for pburst constraints (just for display purposes, not involved in bsingif and thresh computation)
thresh_ratios=0.15;%0.15
thresh_surprises=4

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
for d=dlist%1:numel(datalist)
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
    targs_ind_flip=fliplr(targs_ind);
    
    %target tuning (after compute_tuning)
    targ_tuning=data(1).offline.targ_tuning;
    
    %select trials
    seltrials=get_seltrials(data,'rpt');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute size of adaptive window for analysis
    %     pos=targslist(targ_tuning,:);
    %     amp=sqrt(pos(1)^2+pos(2)^2);
    %     if amp>5,
    %         wadapt=floor(p(2)+amp*p(1));
    %     else
    %         wadapt=ceil(p(2)+5*p(1));
    %     end
    
    wt=100;%50
    wind_targ_bsignif=[110 110+wt];
    wind_targ_bsl=[50-wt 50];
    
    wind_targ_pburst=[0 wt];
    wind_targ_pburst_bsl=[-50-wt -50 ];%[30-wt 30];
    
    ws=50;
    wind_sacc_bsignif=[25-ws 25];%[-25 -25+wadapt];
    wind_sacc_bsl=[-150-ws -150];%[50-ws 50];%
    wind_sacc_bsl_go=[-ws 0];%[50-ws 50];%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Data aligned on target and saccade onset
    for al=1:numel(alignlist)
        info.align=alignlist{al};
        
        switch info.align
            case 'targ'
                [alltrials_spk_targ aligntime_targ]=get_alltrials_align(data,seltrials,wind_targ,'fr',info,targslist,sigma_FR,1);
                %[alltrials_spk_targ_bsignif aligntime_targ_bsignif]=get_alltrials_align(data,seltrials,wind_targ_bsignif,'fr',info,targslist,sigma_FR,0);
                
                [alltrials_spk_targ_bsignif aligntime_targ_bsignif]=get_alltrials_align(data,seltrials,wind_targ_bsignif,'fr',info,targslist,sigma_FR,0);
                
                [alltrials_spk_targ_bsl aligntime_targ_bsl]=get_alltrials_align(data,seltrials,wind_targ_bsl,'fr',info,targslist,sigma_FR,0);
            
                %bsl normalization
                %[alltrials_spk_sacc_bsl aligntime_sacc_bsl]=get_alltrials_align(data,seltrials,wind_sacc_bsl,'fr',info,targslist,sigma_FR,0);
           
            case 'targ_pburst_ch'
                [alltrials_spk_targ aligntime_targ]=get_alltrials_align(data,seltrials,wind_targ,'fr',info,targslist,sigma_FR,1);
                
                [alltrials_spk_targ_pburst aligntime_targ_pburst]=get_alltrials_align(data,seltrials,wind_targ_pburst,'fr',info,targslist,sigma_FR,0);
                
                %alignaux=info.align;
                %info.align='targ';
                [alltrials_spk_targ_pburst_bsl aligntime_targ_pburst_bsl]=get_alltrials_align(data,seltrials,wind_targ_pburst_bsl,'fr',info,targslist,sigma_FR,0);
                %info.align=alignaux;
                
            case 'targ_rsburst_ch'
                [alltrials_spk_targ aligntime_targ]=get_alltrials_align(data,seltrials,wind_targ,'fr',info,targslist,sigma_FR,1);
                
                [alltrials_spk_targ_rsburst aligntime_targ_rsburst]=get_alltrials_align(data,seltrials,wind_targ_pburst,'fr',info,targslist,sigma_FR,0);
                
                %alignaux=info.align;
                %info.align='targ';
                [alltrials_spk_targ_rsburst_bsl aligntime_targ_rsburst_bsl]=get_alltrials_align(data,seltrials,wind_targ_pburst_bsl,'fr',info,targslist,sigma_FR,0);
                %info.align=alignaux;
                
                
            case 'sacc'
                [alltrials_spk_sacc aligntime_sacc]=get_alltrials_align(data,seltrials,wind_sacc,'fr',info,targslist,sigma_FR,1);
                [alltrials_spk_sacc_bsignif aligntime_sacc_bsignif]=get_alltrials_align(data,seltrials,wind_sacc_bsignif,'fr',info,targslist,sigma_FR,0);
                alignaux=info.align;
                info.align='go';
                [alltrials_spk_sacc_bsl aligntime_sacc_bsl]=get_alltrials_align(data,seltrials,wind_sacc_bsl_go,'fr',info,targslist,sigma_FR,0);
                info.align=alignaux;
        end
    end
    
     [allgazepos,allevents]=get_alldatagaze_align(data,seltrials,info,targslist);
       
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %target, target index and target in anti-RF
    info.targ=targ_tuning;
    info.targ_ind=find(targs_ind==targ_tuning);
    targ_tuning_a=targs_ind_flip(info.targ_ind);
    
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot spk of targ and sacc
    figbsignif=figure('Position',[1 100 scrsz(3)-100 scrsz(4)-200]);
    for al=1:numel(alignlist)
        info.align=alignlist{al};
        switch info.align
            case 'targ'
                info.aligntime=aligntime_targ;
                wind_bsignif=wind_targ_bsignif;
                trials_spk=alltrials_spk_targ{targ_tuning};
                
                wind_bsl=wind_targ_bsl;
                trials_spk_bsl=alltrials_spk_targ_bsl{targ_tuning};
                 
                
                
            case 'targ_pburst_ch'
                info.aligntime=aligntime_targ_pburst;
                wind_bsignif=wind_targ_pburst;
                trials_spk=alltrials_spk_targ_pburst{targ_tuning};
                
                wind_bsl=wind_targ_pburst_bsl;
                trials_spk_bsl=alltrials_spk_targ_pburst_bsl{targ_tuning};
                 
            case 'targ_rsburst_ch'
                info.aligntime=aligntime_targ_rsburst;
                wind_bsignif=wind_targ_pburst;
                trials_spk=alltrials_spk_targ_rsburst{targ_tuning};
                
                wind_bsl=wind_targ_pburst_bsl;
                trials_spk_bsl=alltrials_spk_targ_rsburst_bsl{targ_tuning};
                   
                
            case 'sacc'
                info.aligntime=aligntime_sacc;
                wind_bsignif=wind_sacc_bsignif;
                trials_spk=alltrials_spk_sacc{targ_tuning};
                
                wind_bsl=wind_sacc_bsl;
                trials_spk_bsl=alltrials_spk_sacc_bsl{targ_tuning};

                
        end
        
        gazepos=allgazepos{info.targ};
        events=allevents{info.targ};
        
        %%baseline
        %wind_bsl=wind_targ_bsl;
        %trials_spk_bsl=alltrials_spk_targ_bsl{targ_tuning};
        
        
        [info.nchannels info.ntrials info.triallen]=size(trials_spk);
        %compute average trials
        [trials_spk_avg trials_spk_var]=get_trials_avg(trials_spk);
        %remove trials with amplitude that is too small
        %[trials_spk_avgc index_spk_c]=clean_trials(trials_spk_avg,'fr');
        %compute average trials of baseline
        [trials_spk_bsl_avg trials_spk_bsl_var]=get_trials_avg(trials_spk_bsl);
        %normalize average trials by baseline
        %trials_spk_avgcn=get_trials_normalized(trials_spk_avgc,trials_spk_bsl_avg,'fr',info);
        trials_spk_avgn=get_trials_avg_normalized(trials_spk_avg,trials_spk_bsl_avg,'fr',info);
        
        %plot
        hdlfig=subplot(1,2,al);hold on;
        titlestr={info.datafile ; ['FR ' info.align ' t' num2str(info.targ) ' #trials:' num2str(info.ntrials)]};
        [range ~]=plot_trials(trials_spk_avgn,[],[1:info.nchannels],vshift_spk,[],[],info,hdlfig,titlestr,[],[]);
        %plot wind_bsignif limits
        plot_event(wind_bsignif,info.aligntime,range,1,hdlfig);
        %plot wind_bsl limits
        %if strcmp(info.align,'targ')
        plot_event(wind_bsl,info.aligntime,range,3,hdlfig);
        %end
        %same axis
        if al==1,
            axis tight;
            ax1=axis;
        else
            axis tight;
            ax2=axis;
            axis([ax2(1) ax2(2) ax1(3) ax1(4)]);
        end
        
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         %Update data thresh_al
        %         %thresh=5;%threshold spk/s
        %         thresh_al=mean(trials_spk_avgn,2)>thresh; %select channel when at any time the mean FR crossed the threshold
        %         %pause
        %         data=update_data(0,1,0,data,data_path,info.datafile,field,thresh_al);
    end
    
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %compute threshold of targ and sacc
    %TODO: fix plot of targ signals
    
    figthresh=figure('Position',[1 100 scrsz(3)-100 scrsz(4)-200]);
    for al=1:numel(alignlist)
        info.align=alignlist{al};
        switch info.align
            case 'targ'
                info.aligntime=aligntime_targ_bsignif;
                wind_bsignif=wind_targ_bsignif;
                trials_spk=alltrials_spk_targ_bsignif{targ_tuning};
                
                wind_bsl=wind_targ_bsl;
                trials_spk_bsl=alltrials_spk_targ_bsl{targ_tuning};
                 
                field='targ_bthresh';
                
            case 'targ_pburst_ch'
                info.aligntime=aligntime_targ_pburst;
                wind_bsignif=wind_targ_pburst;
                trials_spk=alltrials_spk_targ_pburst{targ_tuning};
                
                wind_bsl=wind_targ_pburst_bsl;
                trials_spk_bsl=alltrials_spk_targ_pburst_bsl{targ_tuning};
                 
                field='targ_pburstch_bthresh';
                
            case 'targ_rsburst_ch'
                info.aligntime=aligntime_targ_rsburst;
                wind_bsignif=wind_targ_pburst;
                trials_spk=alltrials_spk_targ_rsburst{targ_tuning};
                
                wind_bsl=wind_targ_pburst_bsl;
                trials_spk_bsl=alltrials_spk_targ_rsburst_bsl{targ_tuning};
                 
                field='targ_rsburstch_bthresh';
                                
            case 'sacc'
                info.aligntime=aligntime_sacc_bsignif;
                wind_bsignif=wind_sacc_bsignif;
                trials_spk=alltrials_spk_sacc_bsignif{targ_tuning};
                
                wind_bsl=wind_sacc_bsl;
                trials_spk_bsl=alltrials_spk_sacc_bsl{targ_tuning};

                field='sacc_bthresh';
        end
        
        gazepos=allgazepos{info.targ};
        events=allevents{info.targ};
        
        %%baseline
        %wind_bsl=wind_targ_bsl;
        %trials_spk_bsl=alltrials_spk_targ_bsl{targ_tuning};
        
        
        [info.nchannels info.ntrials info.triallen]=size(trials_spk);
        %normalization of each trial 
        trials_spk_n=get_trials_normalized(trials_spk,trials_spk_bsl,'FR',info);
         
        %compute average trials
        [trials_spk_avg trials_spk_var]=get_trials_avg(trials_spk);
        %remove trials with amplitude that is too small
        %[trials_spk_avgc index_spk_c]=clean_trials(trials_spk_avg,'fr');
        %compute average trials of baseline
        [trials_spk_bsl_avg trials_spk_bsl_var]=get_trials_avg(trials_spk_bsl);
        %normalize average trials by baseline
        %trials_spk_avgcn=get_trials_normalized(trials_spk_avgc,trials_spk_bsl_avg,'fr',info);
        trials_spk_avgn=get_trials_avg_normalized(trials_spk_avg,trials_spk_bsl_avg,'fr',info);
        
        %plot
        hdlfig=subplot(1,2,al);hold on;
        titlestr={info.datafile ; ['FR ' info.align ' t' num2str(info.targ) ' #trials:' num2str(info.ntrials)]};
        [range ~]=plot_trials(trials_spk_avgn,[],[1:info.nchannels],vshift_spk,[],[],info,hdlfig,titlestr,[],[]);
        %plot wind_bsignif limits
        plot_event(wind_bsignif,info.aligntime,range,1,hdlfig);
        %plot wind_bsl limits
        %if strcmp(info.align,'targ')
        plot_event(wind_bsl,info.aligntime,range,3,hdlfig);
        %end
        %same axis
        if al==1,
            axis tight;
            ax1=axis;
        else
            axis tight;
            ax2=axis;
            axis([ax2(1) ax2(2) ax1(3) ax1(4)]);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Update data thresh_al
        %thresh_al_trials
        %thresh_al_trials=nanmedian(mean(trials_spk_n,[],3),2)>thresh;
        thresh_al_trials=nanmedian(mean(trials_spk_n,3),2)>thresh;
        data=update_data(0,1,0,data,data_path,info.datafile,[field '_trials'],thresh_al_trials);

        %thresh_al
        %thresh_al=mean(trials_spk_avgn,[],2)>thresh; %select channel when at any time the mean FR crossed the threshold
        thresh_al=mean(trials_spk_avgn,2)>thresh; %select channel when at any time the mean FR crossed the threshold
        data=update_data(0,1,0,data,data_path,info.datafile,field,thresh_al);
        
        %pause
    end
   
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %compute significance of burst of activity for the different alignements
    figure(figbsignif);
    for al=1:numel(alignlist)
        info.align=alignlist{al};
        switch info.align
            case 'targ'
                field='targ_bsignif';
                wind_bsignif=wind_targ_bsignif;
                info.aligntime=aligntime_targ;
                
                %in RF
                trials_spk_bsignif=alltrials_spk_targ_bsignif{targ_tuning};
                trials_spk_bsl=alltrials_spk_targ_bsl{targ_tuning};
                
                %in anti-RF
                trials_spk_bsignif_a=alltrials_spk_targ_bsignif{targ_tuning_a};
                trials_spk_bsl_a=alltrials_spk_targ_bsl{targ_tuning_a};
                
                %thesh_al
                %thresh_al=data(1).offline.targ_bthresh;
                thresh_al=data(1).offline.targ_bthresh_trials;
                
            case 'targ_pburst_ch'
                field='targ_pburstch_bsignif';
                wind_bsignif=wind_targ_pburst;
                info.aligntime=aligntime_targ_pburst;
                
                %in RF
                trials_spk_bsignif=alltrials_spk_targ_pburst{targ_tuning};
                trials_spk_bsl=alltrials_spk_targ_pburst_bsl{targ_tuning};
                
                %in anti-RF
                trials_spk_bsignif_a=alltrials_spk_targ_pburst{targ_tuning_a};
                trials_spk_bsl_a=alltrials_spk_targ_pburst_bsl{targ_tuning_a};
                
                %thesh_al
                %thresh_al=data(1).offline.targ_pburstch_bthresh;
                thresh_al=data(1).offline.targ_pburstch_bthresh_trials;
                
                
            case 'targ_rsburst_ch'
                field='targ_rsburstch_bsignif';
                wind_bsignif=wind_targ_pburst;
                info.aligntime=aligntime_targ_rsburst;
                
                %in RF
                trials_spk_bsignif=alltrials_spk_targ_rsburst{targ_tuning};
                trials_spk_bsl=alltrials_spk_targ_rsburst_bsl{targ_tuning};
                
                %in anti-RF
                trials_spk_bsignif_a=alltrials_spk_targ_rsburst{targ_tuning_a};
                trials_spk_bsl_a=alltrials_spk_targ_rsburst_bsl{targ_tuning_a};
                
                %thesh_al
                %thresh_al=data(1).offline.targ_rsburstch_bthresh;
                thresh_al=data(1).offline.targ_rsburstch_bthresh_trials;
                
                
            case 'sacc'
                field='sacc_bsignif';
                wind_bsignif=wind_sacc_bsignif;
                info.aligntime=aligntime_sacc;
                
                %in RF
                trials_spk_bsignif=alltrials_spk_sacc_bsignif{targ_tuning};
                trials_spk_bsl=alltrials_spk_sacc_bsl{targ_tuning};
                
                %in anti-RF
                trials_spk_bsignif_a=alltrials_spk_sacc_bsignif{targ_tuning_a};
                trials_spk_bsl_a=alltrials_spk_sacc_bsl{targ_tuning_a};
                
                %thesh_al
                %thresh_al=data(1).offline.sacc_bthresh;
                thresh_al=data(1).offline.sacc_bthresh_trials;
                
        end
        %%baseline
        %trials_spk_bsl={targ_tuning};
        %trials_spk_bsl_a=alltrials_spk_targ_bsl{targ_tuning_a};
        
        
        %Peak significance
        %compare to baseline
        H_al=get_bsignif(trials_spk_bsignif,trials_spk_bsl,[],[],alpha);
        %compare to anti-RF (NOTE not good because sometime clear burst and not selected)
        %H_al=get_bsignif(trials_spk_bsignif,trials_spk_bsignif_a,trials_spk_bsl,trials_spk_bsl_a,0.01);
        
        %plot bsignif
        hdlfig=subplot(1,2,al);hold on;            
        center_bsignif(:,1)=H_al'.*(round(wind_bsignif(2)-wind_bsignif(1))/2+wind_bsignif(1)+1*[1:info.nchannels]');
        center_bsignif(find(center_bsignif(:,1)==0),1)=-info.aligntime;
        center_bsignif(:,2)=0;%vshift_spk;
        plot_events_ch(center_bsignif,[],vshift_spk,range,info,hdlfig,[],'-',2,'');
        
%         %plot thresh_al
%         hdlfig=subplot(1,2,al);hold on;            
%         center_bsignif(:,1)=(thresh_al).*(round(wind_bsignif(2)-wind_bsignif(1))/2+wind_bsignif(1)+5+1*[1:info.nchannels]');
%         center_bsignif(find(center_bsignif(:,1)==0),1)=-info.aligntime;
%         center_bsignif(:,2)=0;%vshift_spk;
%         plot_events_ch(center_bsignif,[],vshift_spk,range,info,hdlfig,[],'-.',2);
        
        
        %plot thresh_al_trials
        hdlfig=subplot(1,2,al);hold on;            
        center_bsignif(:,1)=(thresh_al_trials).*(round(wind_bsignif(2)-wind_bsignif(1))/2+wind_bsignif(1)+5+1*[1:info.nchannels]');
        center_bsignif(find(center_bsignif(:,1)==0),1)=-info.aligntime;
        center_bsignif(:,2)=0;%vshift_spk;
        plot_events_ch(center_bsignif,[],vshift_spk,range,info,hdlfig,[],'-.',2,'');
        
        %%
        %plot pburst_signif
        switch info.align
            case 'targ_pburst_ch'
                ratios=data(1).offline.targ_pburst_ratio(targ_tuning,:)>thresh_ratios;
                surprises=data(1).offline.targ_pburst_msurprises(targ_tuning,:)>thresh_surprises;
                
            case 'targ_rsburst_ch'
                ratios=data(1).offline.targ_rsburst_ratio(targ_tuning,:)>thresh_ratios;
                surprises=data(1).offline.targ_rsburst_msurprises(targ_tuning,:)>thresh_surprises;
                
            case 'sacc'
                ratios=data(1).offline.sacc_pburst_ratio(targ_tuning,:)>thresh_ratios;
                surprises=data(1).offline.sacc_pburst_msurprises(targ_tuning,:)>thresh_surprises;
        end
        pburst_signif=[];
        pburst_signif=(ratios & surprises)';
        hdlfig=subplot(1,2,al);hold on;
        center_bsignif(:,1)=(pburst_signif).*(round(wind_bsignif(2)-wind_bsignif(1))/2+wind_bsignif(1)+15+1*[1:info.nchannels]');
        center_bsignif(find(center_bsignif(:,1)==0),1)=-info.aligntime;
        center_bsignif(:,2)=0;%vshift_spk;
        plot_events_ch(center_bsignif,[],vshift_spk,range,info,hdlfig,[],'-',2,'');
        
        
        %Update data H_al
        data=update_data(0,1,0,data,data_path,info.datafile,field,H_al);
        
    end
    %data(1).offline.targ_bsignif
    %data(1).offline.sacc_bsignif
    
    %%
    %update data
    %display('NOT SAVED!')
    update_data(1,0,0,data,data_path,info.datafile,[],[]);
    
    %pause
    close(figbsignif)
    close(figthresh)
   
    
end

