%function compute_vmi

%function compute_vmi
%   Compute VMI index from laminar data
%
% see also compute_tuning
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 10/14/2016 last modified 01/19/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TO DO: select the same trials in 'sacc' than in 'targ_pburst_ch'


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

%burst classification
classiflist={''};
%classiflist={'vis' 'vm' 'mov'};

%alignment
%alignlist={'no' 'targ' 'go' 'sacc'};
%alignlist={'targ' 'sacc'};
alignlist={'targ_pburst_ch' 'sacc'};


%windows of analysis
%plot
wind_targ=[-150 350];
wind_sacc=[-300 200];%[-200 200];



% %to compute vmi
% %target 20deg
% wind_targ_vmi=[100 200];
% wind_targ_bsl=[-50 50];
% wind_sacc_vmi=[-25 75];
%
% %target 3deg
% wind_targ_vmi=[100 150];
% wind_targ_bsl=[0 50];
% wind_sacc_vmi=[-25 25];
%
% %wind_sacc_bsl=[-150 -100]


%adaptive window
[p,polystats] = polyfit([5 20],[55 100],1);


%vshift
vshift_spk=50;%100

%sigma FR
sigma_FR=6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get data
datalist=load_data_gandhilab(data_path);

%colorlist
colorlist=get_colorlist;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analyzing data
dlist=get_dlist

%hdlfigallvmis=figure;hold on;
data=[];
info=[];
for cl=1:numel(classiflist)
    info.classif=classiflist{cl}
    
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
        wind_targ_vmi=[110 110+wt];
        wind_targ_bsl=[50-wt 50];
        
        wind_targ_pburst_vmi=[0 wt];
        wind_targ_pburst_bsl=[-50-wt -50 ];%[30-wt 30];
        
        
        ws=50;
        wind_sacc_vmi=[25-ws 25];%[-25 -25+wadapt];
        %wind_sacc_bsl=[-150-ws -150];%[50-ws 50];%
        wind_sacc_bsl_go=[-ws 0];%[50-ws 50];%
        
        
        
        %     ws=25;
        %     wind_sacc_vmi=[-ws 0];
        %     wind_sacc_bsl=[-150-ws -150];%[50-ws 50];%
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Data aligned on target and saccade onset
        for al=1:numel(alignlist)
            info.align=alignlist{al};
            
            switch info.align
                case 'targ'
                    [alltrials_spk_targ aligntime_targ]=get_alltrials_align(data,seltrials,wind_targ,'fr',info,targslist,sigma_FR,1);
                    [alltrials_spk_targ_vmi aligntime_targ_vmi]=get_alltrials_align(data,seltrials,wind_targ_vmi,'fr',info,targslist,sigma_FR,0);
                    [alltrials_spk_targ_bsl aligntime_targ_bsl]=get_alltrials_align(data,seltrials,wind_targ_bsl,'fr',info,targslist,sigma_FR,0);
                    
                    %bsl normalization
                    %[alltrials_spk_sacc_bsl aligntime_sacc_bsl]=get_alltrials_align(data,seltrials,wind_sacc_bsl,'fr',info,targslist,sigma_FR,0);
                    
                case 'targ_pburst_ch'
                    [alltrials_spk_targ_pburst aligntime_targ_pburst]=get_alltrials_align(data,seltrials,wind_targ,'fr',info,targslist,sigma_FR,1);
                    [alltrials_spk_targ_pburst_vmi aligntime_targ_pburst_vmi]=get_alltrials_align(data,seltrials,wind_targ_pburst_vmi,'fr',info,targslist,sigma_FR,0);
                    [alltrials_spk_targ_pburst_bsl aligntime_targ_pburst_bsl]=get_alltrials_align(data,seltrials,wind_targ_pburst_bsl,'fr',info,targslist,sigma_FR,0);
                    
                    
                case 'sacc'
                    [alltrials_spk_sacc aligntime_sacc]=get_alltrials_align(data,seltrials,wind_sacc,'fr',info,targslist,sigma_FR,1);
                    [alltrials_spk_sacc_vmi aligntime_sacc_vmi]=get_alltrials_align(data,seltrials,wind_sacc_vmi,'fr',info,targslist,sigma_FR,0);
                    
                    alignaux=info.align;
                    info.align='go';
                    [alltrials_spk_sacc_bsl aligntime_sacc_bsl]=get_alltrials_align(data,seltrials,wind_sacc_bsl_go,'fr',info,targslist,sigma_FR,0);
                    info.align=alignaux;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %get spk data for tuning
        [alltrials_spk_tuning info.aligntime]=get_alltrials_align(data,seltrials,[],'fr',info,targslist,sigma_FR,0);
        
         
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %analysis of trials for each target
        figtrials=figure('Position',[1 100 scrsz(3)-100 scrsz(4)-200]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %display all targets
        hdlfig=subplot(2,3,1);hold on;
        display_alltargets(targslist,info,hdlfig);
        
        %     %compute target tuning
        %     hdlfig=subplot(2,3,4);hold on;
        %     plot_targtuning(alltrials_spk_tuning,targs_ind,info,hdlfig,'Target tuning');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %target, target index and target in anti-RF
        info.targ=targ_tuning;
        info.targ_ind=find(targs_ind==targ_tuning);
        targ_tuning_a=targs_ind_flip(info.targ_ind);
        
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %plot spk of targ and sacc
        for al=1:numel(alignlist)
            info.align=alignlist{al};
            switch info.align
                case 'targ'
                    info.aligntime=aligntime_targ;
                    wind_vmi=wind_targ_vmi;
                    trials_spk=alltrials_spk_targ{targ_tuning};
                    
                    wind_bsl=wind_targ_bsl;
                    trials_spk_bsl=alltrials_spk_targ_bsl{targ_tuning};
                    
                    %field='targ_bthresh';
                    
                case 'targ_pburst_ch'
                    info.aligntime=aligntime_targ_pburst;
                    wind_vmi=wind_targ_pburst_vmi;
                    trials_spk=alltrials_spk_targ_pburst{targ_tuning};
                    
                    wind_bsl=wind_targ_pburst_bsl;
                    trials_spk_bsl=alltrials_spk_targ_pburst_bsl{targ_tuning};
                    
                case 'sacc'
                    info.aligntime=aligntime_sacc;
                    wind_vmi=wind_sacc_vmi;
                    trials_spk=alltrials_spk_sacc{targ_tuning};
                    
                    %wind_bsl=wind_sacc_bsl;
                    trials_spk_bsl=alltrials_spk_sacc_bsl{targ_tuning};
                    
                    %field='sacc_bthresh';
            end
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
            %trials_spk_avgcn=get_trials_avg_normalized(trials_spk_avgc,trials_spk_bsl_avg,'fr',info);
            trials_spk_avgn=get_trials_avg_normalized(trials_spk_avg,trials_spk_bsl_avg,'fr',info);
            
            
            
            
            %plot
            hdlfig=subplot(2,3,al+1);hold on;
            titlestr={info.datafile ; ['FR ' info.align ' t' num2str(info.targ) ' #trials:' num2str(info.ntrials)]};
            [range ~]=plot_trials(trials_spk_avgn,[],[1:info.nchannels],vshift_spk,[],[],info,hdlfig,titlestr,[],[]);
            %plot wind_vmi limits
            plot_event(wind_vmi,info.aligntime,range,1,hdlfig);
            %plot wind_bsl limits
            if strcmp(info.align,'targ')
                plot_event(wind_bsl,info.aligntime,range,3,hdlfig);
            end
            %same axis
            if al==1,
                axis tight;
                ax1=axis;
            else
                axis tight;
                ax2=axis;
                axis([ax2(1) ax2(2) ax1(3) ax1(4)]);
            end
            
            
            %         %NOTE: done in compute_bsignif
            %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %         %Update data thresh_al
            %         thresh=10;%threshold spk/s
            %         thresh_al=sum(trials_spk_avgn>thresh,2)'; %select channel when at any time the mean FR crossed the threshold
            %         data=update_data(0,1,0,data,data_path,info.datafile,field,thresh_al);
        end
        
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %index of not-significant activity: select channel if at least one the burst has a significant activity
        %and classification selection
        %NOTE: here find indexes of not-significant bursts
        
        switch alignlist{1}
            case 'targ'
                targ_bsignif=data(1).offline.targ_bsignif;sacc_bsignif=data(1).offline.sacc_bsignif;
                targ_bthresh=data(1).offline.targ_bthresh';sacc_bthresh=data(1).offline.sacc_bthresh';
                %ind_bsignif=find(targ_bsignif==0 & sacc_bsignif==0 )
                ind_targ_bsignif=find(targ_bsignif==0 | targ_bthresh==0);
                ind_sacc_bsignif=find(sacc_bsignif==0 | sacc_bthresh==0);
                ind_bsignif=find((targ_bsignif==0 | targ_bthresh==0) & (sacc_bsignif==0 | sacc_bthresh==0));
                
                
            case 'targ_pburst_ch'
                %no classification
                if strcmp(info.classif,'')
                    %thresholds
                    thresh_ratios=0.15;%0.15
                    thresh_surprises=4;
                    %targ_pburst
                    ratios_targ=data(1).offline.targ_pburst_ratio(targ_tuning,:)>thresh_ratios;
                    surprises_targ=data(1).offline.targ_pburst_msurprises(targ_tuning,:)>thresh_surprises;
                    bsignif_targ=data(1).offline.targ_pburstch_bsignif;
                    bthresh_targ=data(1).offline.targ_pburstch_bthresh_trials';
                    %sacc
                    ratios_sacc=data(1).offline.sacc_pburst_ratio(targ_tuning,:)>thresh_ratios;
                    surprises_sacc=data(1).offline.sacc_pburst_msurprises(targ_tuning,:)>thresh_surprises;
                    bsignif_sacc=data(1).offline.sacc_bsignif;
                    bthresh_sacc=data(1).offline.sacc_bthresh_trials';
                    
                    %targ_bsignif=(ratios_targ & surprises_targ & bsignif_targ & bthresh_targ);
                    %sacc_bsignif=(ratios_sacc & surprises_sacc & bsignif_sacc & bthresh_sacc);
                    targ_bsignif=(bsignif_targ & bthresh_targ);
                    sacc_bsignif=(bsignif_sacc & bthresh_sacc);
                    
                    %indexes not significant
                    ind_targ_pburst_bsignif=find(targ_bsignif==0);
                    ind_sacc_bsignif=find(sacc_bsignif==0);
                    ind_bsignif=find(targ_bsignif==0 & sacc_bsignif==0);
                    
                else
                    %%%%%%%%%%%%%%%%%
                    %NOTE: all before could be removed and replaced by classif lists
                    classif_vis=data(1).offline.classif_pburst_vis;
                    classif_vm=data(1).offline.classif_pburst_vm;
                    classif_mov=data(1).offline.classif_pburst_mov;
                    
                    ind_notvis_bsignif=find(classif_vis==0);
                    ind_notvm_bsignif=find(classif_vm==0);
                    ind_notmov_bsignif=find(classif_mov==0);
                    
                    switch info.classif
                        case 'vis'
                            ind_targ_pburst_bsignif=ind_notvis_bsignif;
                            ind_sacc_bsignif=ind_notvis_bsignif;
                            ind_bsignif=ind_notvis_bsignif
                        case 'vm'
                            ind_targ_pburst_bsignif=ind_notvm_bsignif;
                            ind_sacc_bsignif=ind_notvm_bsignif;
                            ind_bsignif=ind_notvm_bsignif
                        case 'mov'
                            ind_targ_pburst_bsignif=ind_notmov_bsignif;
                            ind_sacc_bsignif=ind_notmov_bsignif;
                            ind_bsignif=ind_notmov_bsignif
                    end
                end
                %pause
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % vmi index
        %loop on vmi index type
        allspk_vmi_avg=zeros(2,length(targslist),info.nchannels);
        allspk_vmi_avg_bsignif=zeros(2,length(targslist),info.nchannels);
        allspk_vmi_var=zeros(2,length(targslist),info.nchannels);
        for v=2%[2 4] %1:4
            %get spk and baseline data to compute vmi
            for al=1:numel(alignlist)
                info.align=alignlist{al};
                switch info.align
                    case 'targ'
                        info.aligntime=aligntime_targ_vmi;
                        
                        %in RF
                        trials_spk_vmi=alltrials_spk_targ_vmi{targ_tuning};
                        trials_spk_bsl=alltrials_spk_targ_bsl{targ_tuning};
                        
                        field_targ=[info.classif ''];
                        
                    case 'targ_pburst_ch'
                        info.aligntime=aligntime_targ_pburst_vmi;
                        
                        %in RF
                        trials_spk_vmi=alltrials_spk_targ_pburst_vmi{targ_tuning};
                        trials_spk_bsl=alltrials_spk_targ_pburst_bsl{targ_tuning};
                        
                        field_targ=[info.classif 'pburst_'];
                        
                    case 'sacc'
                        info.aligntime=aligntime_sacc_vmi;
                        
                        %in RF
                        trials_spk_vmi=alltrials_spk_sacc_vmi{targ_tuning};
                        trials_spk_bsl=alltrials_spk_sacc_bsl{targ_tuning};
                        
                end
                %%baseline
                %trials_spk_bsl=alltrials_spk_targ_bsl{targ_tuning};
                
                %spk
                [info.nchannels info.ntrials info.triallen]=size(trials_spk_vmi);
                %compute average trials
                [trials_spk_vmi_avg trials_spk_vmi_var]=get_trials_avg(trials_spk_vmi);
                %remove trials with amplitude that is too small
                %[trials_spk_avgc index_spk_c]=clean_trials(trials_spk_avg,'spk');
                
                %compute average trials of baseline
                [trials_spk_bsl_avg trials_spk_bsl_var]=get_trials_avg(trials_spk_bsl);
                
                
                %%%%%%%%%%%%%%%%%%%%
                %different vmis
                switch v
                    case 1
                        %mean
                        allspk_vmi_avg(al,targ_tuning,:)=mean(trials_spk_vmi_avg,2);
                        field='vmis_mean';
                    case 2
                        %mean-baseline
                        trials_spk_vmi_avg_aux=trials_spk_vmi_avg-mean(trials_spk_bsl_avg,2);
                        
                        allspk_vmi_avg(al,targ_tuning,:)=mean(trials_spk_vmi_avg_aux,2);
                        %allspk_vmi_avg(al,targ_tuning,:)=abs(mean(trials_spk_vmi_avg_aux,2));
                        
                        allspk_vmi_var(al,targ_tuning,:)=var(trials_spk_vmi_avg_aux,[],2);
                        
                        %field='vmis_mean_bsl';
                        field='vmis_mean_bslbefore';
                        %field='vmis_mean_bslbefore_2';
                        %field_dp='dprimes_mean_bsl';%'dprimes_mean_bslbefore';
                    case 3
                        %peak
                        allspk_vmi_avg(al,targ_tuning,:)=max(trials_spk_vmi_avg,[],2);
                        field='vmis_peak';
                    case 4
                        %peak-baseline
                        allspk_vmi_avg(al,targ_tuning,:)=max(trials_spk_vmi_avg,[],2)-mean(trials_spk_bsl_avg,2);
                        allspk_vmi_var(al,targ_tuning,:)=var(trials_spk_vmi_avg,[],2);
                        %field='vmis_peak_bsl';
                        field='vmis_peak_bslbefore';%
                        %field_dp='dprimes_peak_bsl';%'dprimes_peak_bslbefore';
                end
                
                %%%%%%%%%%%%%%%%%%%%
                %Update data allspk_vmi_avg_bsignif
                %correction using significance of burst activity
                allspk_vmi_avg_bsignif(al,targ_tuning,:)=allspk_vmi_avg(al,targ_tuning,:);
                switch info.align
                    case 'targ'
                        allspk_vmi_avg_bsignif(al,targ_tuning,ind_targ_bsignif)=nan;
                    case 'targ_pburst_ch'
                        allspk_vmi_avg_bsignif(al,targ_tuning,ind_targ_pburst_bsignif)=nan;
                    case 'sacc'
                        allspk_vmi_avg_bsignif(al,targ_tuning,ind_sacc_bsignif)=nan;
                end
                
                %%%%%%%%%%%%%%%%%%%%
                figure(figtrials);hdlfig=subplot(2,3,5);hold on;
                plot(squeeze(allspk_vmi_avg(al,targ_tuning,:)),[1:info.nchannels],'linewidth',2,'color',colorlist(al+v,:));
                
            end
            figure(figtrials);hdlfig=subplot(2,3,5);hold on;
            axis tight;xlabel('Avg FR (spk/s)');ylabel('Channels');
            
            %update data allspk_vmi_avg_bsignif
            field2=['vmi.' field_targ field '_' num2str(wt) '_' num2str(ws)];
            data=update_data(0,1,0,data,data_path,info.datafile,[field2 '_bursts'],allspk_vmi_avg_bsignif);
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %compute vmi
            %NOTE: make a subfunction
            vmis=zeros(length(targslist),info.nchannels);
            vmis_amp=zeros(length(targslist),info.nchannels);
            vmis_c=zeros(length(targslist),info.nchannels);
            vmis_bsignif=zeros(length(targslist),info.nchannels);
            dprimes=zeros(length(targslist),info.nchannels);
            dprimes_bsignif=zeros(length(targslist),info.nchannels);
            
            for tg=targ_tuning,%targs_ind,
                
                %Visuo-Motor Index formula
                vmis(tg,:)=squeeze((allspk_vmi_avg(2,tg,:)-allspk_vmi_avg(1,tg,:))./(allspk_vmi_avg(2,tg,:)+allspk_vmi_avg(1,tg,:)));
                
                %             %correction using threshold on activity
                %             vmis_amp(tg,:)=max(squeeze(allspk_vmi_avg(2,tg,:)),squeeze(allspk_vmi_avg(1,tg,:)));
                %             vmis_amp(tg,:)=vmis_amp(tg,:)/max(vmis_amp(tg,:),[],2);
                %             ind=find(vmis_amp(tg,:)<0.2);%threshold
                %             vmis_c(tg,:)=vmis(tg,:);
                %             vmis_c(tg,ind)=nan;
                
                
                
                %correction using significance of burst activity
                %NOTE could have avoid this step by using allspk_vmi_avg_bsignif
                vmis_bsignif(tg,:)=vmis(tg,:);
                vmis_bsignif(tg,ind_bsignif)=nan;
                
                %if value beyond 1 or -1 because of normalization force them to
                %be 1 or -1
                vmis_bsignif(tg,find(vmis_bsignif(tg,:)>1))=1;
                vmis_bsignif(tg,find(vmis_bsignif(tg,:)<-1))=-1;
                
                %dprime
                dprimes(tg,:)=squeeze((allspk_vmi_avg(2,tg,:)-allspk_vmi_avg(1,tg,:))./sqrt(0.5*(allspk_vmi_var(2,tg,:)+allspk_vmi_var(1,tg,:))));
                dprimes_bsignif(tg,:)=dprimes(tg,:);
                dprimes_bsignif(tg,ind_bsignif)=nan;
                
            end
            figure(figtrials);hdlfig=subplot(2,3,6);hold on;
            %plot(vmis(targ_tuning,:),[1:info.nchannels],'linewidth',5,'color',colorlist(1,:));
            %plot(vmis_amp(targ_tuning,:),[1:info.nchannels],'linewidth',2,'color',colorlist(3,:));
            %plot(vmis_c(targ_tuning,:),[1:info.nchannels],'linewidth',2,'color',colorlist(5,:));
            
            %plot vmis
            plot_vmis(vmis_bsignif,targ_tuning,'-',1,3,info,hdlfig,[]);
            
            %         figure(hdlfigallvmis)
            %         %plot(vmis_c(targ_tuning,:),[1:info.nchannels],'linewidth',3,'color',colorlist(1,:));
            %
            %         plot_vmis(vmis_bsignif,targ_tuning,'-',1,3,info,hdlfig,[]);
            %         axis([-1 1 1 info.nchannels]);
            %         xlabel('VMI');ylabel('Channels');
            
            %plot dprime
            figure(figtrials);hdlfig=subplot(2,3,4);hold on;
            plot_dprimes(dprimes_bsignif,targ_tuning,'-',1,3,info,hdlfig,[]);
            
            %%
            %%%%%%%%%%%%%%%%%%
            
            %update data vmis
            field2=['vmi.' field_targ field '_' num2str(wt) '_' num2str(ws)];
            data=update_data(0,1,0,data,data_path,info.datafile,field2,vmis);
            %%data=update_data(0,1,0,data,data_path,info.datafile,[field2 '_c'],vmis_c);
            data=update_data(0,1,0,data,data_path,info.datafile,[field2 '_bsignif'],vmis_bsignif);
            
            %%update data dprime
            %field3=['vmi.' field_dp '_' num2str(wt) '_' num2str(ws)];
            %data=update_data(0,1,0,data,data_path,info.datafile,[field3 '_bsignif'],dprimes_bsignif);
            
            %pause
            
        end
        
        display('NOT SAVED!')
        %update_data(1,0,0,data,data_path,info.datafile,[],[]);
        pause
        %close(figtrials)
        
        
    end
end

