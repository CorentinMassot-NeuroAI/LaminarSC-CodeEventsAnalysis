%function compute_onset_buildup_avg

%function compute_onset_buildup_avg
%   Compute onset of build-up and burst activity of average spiking activity recorded with a
%   laminar probe (LMA) 
%
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 12/15/2016 last modified 01/23/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set paths
[root_path data_path save_path]=set_paths;

save_path=[root_path 'Results\Results_SC_onsetbuildup\'];


%screen size
scrsz = get(groot,'ScreenSize');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
%print figures and save data
savedata=0;
savefigs=0;
figtype='png';%'epsc2';

%display plot
disp_plot=0;

%alignement
%alignlist={'no' 'targ' 'go' 'sacc'};
%alignlist={'targ'};
alignlist={'sacc'};


%windows of analysis (do not change)
%wind_targ=[-50 250];%[-50 150];
%wind_sacc=[-400 100];%[-100 200]
wind_sacc=[-600 100];%[-100 200]

%windows baseline
wt=100;
%wind_targ_pburst_bsl=[-50-wt -50 ];
wind_bsl_sacc=[-100 0];

% %starting time for latency estimation
% starting_targ=50;
% starting_sacc=-149;

% %gaussian window for latency
% gw_width=2.5;

% %alpha of ttest
% alpha=0.01;%0.01;

%vshift
vshift_spk=100;
vshift_lfp=30;%29;


%sigma FR
sigma_FR=6;

% %shift ripple temporal correction
% shift_ripple=4


%start .pptx file
savepptx=0;
if savepptx
    isOpen  = exportToPPTX();
    if ~isempty(isOpen),
        % If PowerPoint already started, then close first and then open a new one
        exportToPPTX('close');
    end
    
    exportToPPTX('new','Dimensions',[12 6], ...
        'Title','SC onset buildup', ...
        'Author','Corentin', ...
        'Subject','Automatically generated PPTX file from output of analysis_onset_buildup_avg.m', ...
        'Comments',' ');

    %tmp filename
    file=[save_path 'onset_buildup_tmp' '.' figtype];
    
    %pptx filename
    datenow=datestr(now);
    datenow=[datenow(1:11) '-' datenow(13:14) 'h' datenow(16:17) 'm' datenow(19:20) 's']
    filepptx=[save_path 'SC_onsetbuildup-' datenow]

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get data
datalist=load_data_gandhilab(data_path);


%colorlist
colorlist=get_colorlist;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analyzing data
dlist=get_dlist


data=[];
info=[];
dd=0;
for d=dlist([1:end])%dlist(20:end)
    %counter
    dd=dd+1;
    
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
    amplist=sqrt(targslist(:,1).^2+targslist(:,2).^2);
    
    %targets index
    targs_ind=get_targsindex(targslist,info);
    targs_ind_flip=fliplr(targs_ind);
    
    %target tuning (after compute_tuning)
    info.targ_tuning=data(1).offline.targ_tuning;
    
    %select trials
    seltrials=get_seltrials(data,'rpt');
    
    %bursts significance
    thresh_ratios=0.15;%0.15
    thresh_surprises=4;
    %     %targ_pburst
    %     ratios_targ=data(1).offline.targ_pburst_ratio(info.targ_tuning,:)>thresh_ratios;
    %     surprises_targ=data(1).offline.targ_pburst_msurprises(info.targ_tuning,:)>thresh_surprises;
    %     bsignif_targ=data(1).offline.targ_pburstch_bsignif;
    %     bthresh_targ=data(1).offline.targ_pburstch_bthresh_trials';
    %sacc
    ratios_sacc=data(1).offline.sacc_pburst_ratio(info.targ_tuning,:)>thresh_ratios;
    surprises_sacc=data(1).offline.sacc_pburst_msurprises(info.targ_tuning,:)>thresh_surprises;
    bsignif_sacc=data(1).offline.sacc_bsignif;
    bthresh_sacc=data(1).offline.sacc_bthresh_trials';
    
    %     targ_bsignif=(ratios_targ & surprises_targ & bsignif_targ & bthresh_targ)
    %sacc_bsignif=(ratios_sacc & surprises_sacc & bsignif_sacc & bthresh_sacc);
    sacc_bsignif=(bsignif_sacc & bthresh_sacc);
    
    
    
    %loop across all alignements
    aux_spk=[];aux_lfp=[];auxp_spk=[];auxp_lfp=[];auxv_spk=[];auxv_lfp=[];auxvbsl_spk=[];auxvbsl_lfp=[];
    for al=1%1:numel(alignlist)
        info.align=alignlist{al};
        
        %get all neural and behavioral data with specific alignement
        switch info.align
            case 'sacc'
                wind=wind_sacc;
                wind_bsl=wind_bsl_sacc;
                burst_bsignif=sacc_bsignif;
                %starting=starting_sacc;
        end
        
        %signals
        [alltrials_spk,aligntime_spk]=get_alltrials_align(data,seltrials,wind,'fr',info,targslist,sigma_FR,1);
        %         [alltrials_lfp,aligntime_lfp]=get_alltrials_align(data,seltrials,wind+shift_ripple,'lfp',info,targslist,sigma_FR,1);
        
        %baseline
        alignaux=info.align;
        info.align='go';
        [alltrials_spk_bsl aligntime_bsl]=get_alltrials_align(data,seltrials,wind_bsl,'fr',info,targslist,sigma_FR,1);
        %         [alltrials_lfp_bsl aligntime_bsl]=get_alltrials_align(data,seltrials,wind_bsl,'lfp',info,targslist,sigma_FR,1);
        info.align=alignaux;
        
        
        %burst_bsignif
        burst_bsignif=double(burst_bsignif);
        burst_bsignif(find(burst_bsignif==0))=nan;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %analysis of trials for each target
        for tg=info.targ_tuning;%targs_ind,
            %target, target index and target in anti-RF
            info.targ=tg;
            info.targ_ind=find(targs_ind==info.targ_tuning);
            tga=targs_ind_flip(info.targ_ind);
            info.targ_tuning_a=tga;
            
            
            %neural and behavioral signals for target tg
            trials_spk=alltrials_spk{tg};
            trials_spk_bsl=alltrials_spk_bsl{tg};
            
            %             trials_lfp=alltrials_lfp{tg};
            %             trials_lfp_bsl=alltrials_lfp_bsl{tg};
            
            
            %figure
            figtrials=figure('Position',[scrsz(3)/3 100 scrsz(3)/2 scrsz(4)-200]);
            
            
            
            %%%%%%%%%%%%%%%%%%
            %spk
            [info.nchannels info.ntrials info.triallen]=size(trials_spk);
            info.aligntime=aligntime_spk;
            
            %compute average of normalized trials in RF
            trials_spk_n=get_trials_normalized(trials_spk,trials_spk_bsl,'FR',info);
            [trials_spk_n_avg trials_spk_n_var]=get_trials_avg(trials_spk_n);

            %compute average trials in RF
            %[trials_spk_avg trials_spk_var]=get_trials_avg(trials_spk);
            %[trials_spk_bsl_avg trials_spk_bsl_var]=get_trials_avg(trials_spk_bsl);
            
            
            %plot avg and ste
            figure(figtrials)
            hdlfig=subplot(1,1,1);hold on;
            [range_spk vshift_spk]=plot_trials(trials_spk_n_avg,[],[],[],[],[],info,hdlfig,[],[],[]);
            plot_trials(trials_spk_n_avg+trials_spk_n_var,[],[],[],[],[],info,hdlfig,[],'-',1);
            plot_trials(trials_spk_n_avg-trials_spk_n_var,[],[],[],[],[],info,hdlfig,[],'-',1);
            
            
            %%%%%%%%%%%%%%%%%%
            %bootstrap on trials
            ntrials=size(trials_spk,2);
            trials_boot=[];fitparamsd_boot=[];peak_boot=[];elb1_boot=[];elb2_boot=[];   
            inflect1_boot=[];inflect2_boot=[];inflect3_boot=[];
            resnorm1_boot={};resnorm2_boot={};resnorm2_boot={};
            fitparams1_boot=[];fitparams2_boot=[];fitparams3_boot=[];
            resnormfit1_boot=[];resnormfit2_boot=[];resnormfit3_boot=[];
                
            nboot=100;%1000%20;%200;
            for b=1:nboot
                d
                dd
                b
                
                %random selection with replacement
                samples=randi(ntrials,1,ntrials);
                %samples=[1:ntrials];display('No bootstrapping!');
                
                
                trials_spk_boot=trials_spk(:,samples,:);
                trials_spk_bsl_boot=trials_spk_bsl(:,samples,:);
                
                %%%%%
                %compute average of normalized trials in RF
                [info.nchannels info.ntrials info.triallen]=size(trials_spk_boot);
                info.aligntime=aligntime_spk;
                trials_spk_n_boot=get_trials_normalized(trials_spk_boot,trials_spk_bsl_boot,'FR',info);
                [trials_spk_n_avg_boot trials_spk_n_var_boot]=get_trials_avg(trials_spk_n_boot);
                
                
                %%%%%
                %detrend activity (in case of anticipatory activity after go cue)
                %linear function
                fun = @(x,xdata)x(1)+x(2)*xdata;
                options=optimset('Display','off');
                
                snip_b=-500+info.aligntime;%-300+info.aligntime;
                snip_e=-300+info.aligntime;%-150+info.aligntime;
                trial_snip=trials_spk_n_avg_boot(:,snip_b:snip_e);
                trial_snipd=nan(info.nchannels,size(trials_spk_n_avg_boot,2));
                fitparams_d=[];
                for ch=1:info.nchannels
                    xp1=[1:length(trial_snip(ch,:))];
                    x0p1=[xp1(1) trial_snip(ch,1)];
                    [fitparams1 resnorm1] = lsqcurvefit(fun,x0p1,xp1,trial_snip(ch,:),[],[],options);
                    %display linear trend
                    plot(xp1+snip_b,fun(fitparams1,xp1),'b-')
                    trial_snipd(ch,xp1+snip_b:xp1+snip_b+length(trial_snip(ch,:))-1)=fun(fitparams1,xp1);
                    fitparams_d(ch,:)=fitparams1;%detrending parameters
                end
                %plot_trials(trial_snipd,[],[],vshift_spk,[],[],info,hdlfig,[],'--',2);

                %detrend trials
                trials_spk_n_avg_boot_d=trials_spk_n_avg_boot;
                trials_spk_boot_d=trials_spk_boot;
                xd=[1:length(trials_spk_n_avg_boot_d(1,snip_b:end))];
                for ch=1:info.nchannels
                    trials_spk_n_avg_boot_d(ch,snip_b:end)=trials_spk_n_avg_boot_d(ch,snip_b:end)-fun(fitparams_d(ch,:),xd);
                    for t=1:size(trials_spk_boot_d,2)
                        trials_spk_boot_d(ch,t,snip_b:end)=squeeze(trials_spk_boot_d(ch,t,snip_b:end))'-fun(fitparams_d(ch,:),xd);
                    end
                end
                %plot_trials(trials_spk_n_avg_boot_d,[],[],vshift_spk,[],[],info,hdlfig,[],'--',2);
                
                
                
                %%%%%
                %peak activity
                [peak_vals peak_times]=max(trials_spk_n_avg_boot(:,1:info.aligntime+30),[],2);
                peak(:,1)=peak_times;
                peak(:,2)=peak_vals;
                peak(:,1)=peak(:,1)-info.aligntime;%correction for timing
                
                %plot peak
                peak_plot=[];
                peak_plot(:,1)=peak(:,1);
                peak_plot(:,2)=peak(:,2);
                plot_events_ch(peak_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','--',1,'k');
                %plot_events_ch(peak_plot,[],vshift_spk,range_spk,info,hdlfig,'n','-',1,'');
                
                %%
                %%%%%
                %elb_1 (spk onsets based on significance of activity on detrended trials)
                %forward approach
                alpha=0.01;
                search_wind(1)=info.aligntime-300;%info.aligntime-250;
                latcount_max=100;%duration significance is maintained
                latcount_min=latcount_max;%min duration not to reset latcount
                elb_1=[];
                for ch=1:info.nchannels
                    search_wind(2)=peak(ch,1)+info.aligntime;
                    %onset=get_latency_trials(squeeze(trials_spk_boot(ch,:,:)),squeeze(trials_spk_bsl_boot(ch,:,:)),search_wind,1,'fr',info,alpha,latcount_max,latcount_min);
                    onset=get_latency_trials(squeeze(trials_spk_boot_d(ch,:,:)),squeeze(trials_spk_bsl_boot(ch,:,:)),search_wind,1,'fr',info,alpha,latcount_max,latcount_min);
                    
                    if ~isnan(onset(1))
                        elb_1(ch,1)=onset(1)-info.aligntime;%correction for timing
                        elb_1(ch,2)=trials_spk_n_avg_boot(ch,onset(1));
                    else
                        elb_1(ch,1:2)=nan;
                    end
                end
                
                %plot elb_1
                if disp_plot
                    figure(figtrials)
                    onset_plot=[];
                    onset_plot(:,1)=elb_1(:,1);
                    onset_plot(:,2)=elb_1(:,2);
                    %plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','--',1,'');
                    plot_events_ch(onset_plot,[],vshift_spk,range_spk,info,hdlfig,'n','--',1,'');
                end
                
%                 %%%%%
%                 %elb_2 (backward approach)
%                 %winsize=10;
%                 alpha=0.01;
%                 latcount_max=100;%duration significance is maintained
%                 latcount_min=round(latcount_max);%min duration not to reset latcount
%                 elb_2=[];
%                 for ch=1:info.nchannels
%                     peak_ch=peak(ch,1)+info.aligntime;
%                     onset=get_latency_trials_back(squeeze(trials_spk_boot(ch,:,:)),peak_ch,'fr',info,alpha,latcount_max,latcount_min);
%                     if ~isnan(onset(1))
%                         elb_2(ch,1)=onset(1)-info.aligntime;%correction for timing
%                         elb_2(ch,2)=trials_spk_n_avg_boot(ch,onset(1));
%                     else
%                         elb_2(ch,1:2)=nan;
%                     end
%                 end
%                 
%                 %plot elb_2
%                 if disp_plot
%                     figure(figtrials)
%                     onset_plot=[];
%                     onset_plot(:,1)=elb_2(:,1);
%                     onset_plot(:,2)=elb_2(:,2);
%                     %plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','-',1,'');
%                     plot_events_ch(onset_plot,[],vshift_spk,range_spk,info,hdlfig,'n','-',1,'');
%                 end
                
                %%
                %%%%%
                %inflection 1 point between elb_1 and Peak 
                anchor_b=elb_1(:,1)+info.aligntime;
                anchor_e=peak(:,1)+info.aligntime;
                min_bins=10;
                inflect_1=nan(info.nchannels,3);
                resnorm_1=cell(1,info.nchannels);
                fitparams_1=nan(info.nchannels,4);
                resnormfit_1=nan(info.nchannels,2);
                %figresnorm=figure;hold on;
                for ch=1:info.nchannels
                    if ~isnan(anchor_b(ch)) & ~isnan(anchor_e(ch)) & (anchor_e(ch)-anchor_b(ch)+1)>(2*min_bins)
                        %ch
                        [inflect,resnorm,resnorm_fit,fitparams]=get_inflection_2pwlr(trials_spk_n_avg_boot(ch,:),anchor_b(ch),anchor_e(ch),min_bins);
                        inflect_1(ch,1)=inflect(1)-info.aligntime;%correction for timing
                        inflect_1(ch,2:3)=inflect(2:3);
                        resnorm_1{ch}=resnorm;
                        fitparams_1(ch,:)=fitparams;
                        resnormfit_1(ch,:)=resnorm_fit;
                        
                        %plot resnorm
%                         if ~isempty(resnorm_1{ch}) %trial_sel was long enough
%                             resnorm_plot=nan(1,info.aligntime+length(resnorm_1{ch})-1);
%                             resnorm_plot(1,info.aligntime:end)=resnorm_1{ch};
%                             resnorm_plot=resnorm_plot/max(resnorm_plot);
%                             figure(figresnorm);subplot(1,3,1);hold on;
%                             plot(resnorm_plot,'color',colorlist(ch,:),'linewidth',1);
%                             title('Residuals of 2pwlr for inflect 1');xlabel('Time (ms)');ylabel('R2');
%                         end
                    end
                end
                
                %channel without inflection point
                %ind_ninflect=find(fitparamslist(:,2)<0.09*0.001)
                
                %plot inflection
                if disp_plot
                    figure(figtrials)
                    inflect_plot=[];
                    inflect_plot(:,1)=inflect_1(:,1);
                    inflect_plot(:,2)=inflect_1(:,2);
                    %inflect_spk_plot(ind_ninflect,:)=nan;
                    %plot_events_ch(inflect_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','-',2,'');
                    plot_events_ch(inflect_plot,[],vshift_spk,range_spk,info,hdlfig,'n','-',2,'');
                end
                
%                 %%
%                 %%%%%
%                 %inflection 2 point between elb_2 and elb_2 - 250 
%                 inflect2_wind=2*abs(elb_2(:,1)-elb_1(:,1));%250;%150;
%                 anchor_b=nanmax(-250,elb_2(:,1)-inflect2_wind)+info.aligntime;
%                 anchor_e=elb_2(:,1)+info.aligntime;
%                 min_bins=10;
%                 inflect_2=nan(info.nchannels,3);
%                 resnorm_2=cell(1,info.nchannels);
%                 fitparams_2=nan(info.nchannels,4);
%                 resnormfit_2=nan(info.nchannels,2);
%                 for ch=1:info.nchannels
%                     if ~isnan(anchor_b(ch)) & ~isnan(anchor_e(ch)) & (anchor_e(ch)-anchor_b(ch)+1)>(2*min_bins)
%                         %ch
%                         [inflect,resnorm,resnorm_fit,fitparams]=get_inflection_2pwlr(trials_spk_n_avg_boot(ch,:),anchor_b(ch),anchor_e(ch),min_bins);
%                         inflect_2(ch,1)=inflect(1)-info.aligntime;%correction for timing
%                         inflect_2(ch,2:3)=inflect(2:3);
%                         resnorm_2{ch}=resnorm;
%                         fitparams_2(ch,:)=fitparams;
%                         resnormfit_2(ch,:)=resnorm_fit;
%                         
%                         %plot resnorm
% %                         if ~isempty(resnorm_2{ch}) %trial_sel was long enough
% %                             resnorm_plot=nan(1,info.aligntime+length(resnorm_2{ch})-1);
% %                             resnorm_plot(1,info.aligntime:end)=resnorm_2{ch};
% %                             resnorm_plot=resnorm_plot/max(resnorm_plot);
% %                             figure(figresnorm);subplot(1,3,2);hold on;
% %                             plot(resnorm_plot,'color',colorlist(ch,:),'linewidth',1);
% %                             title('Residuals of 2pwlr for inflect 2');xlabel('Time (ms)');ylabel('R2');
% %                         end
%                     end
%                 end
%                 
%                 %channel without inflection point
%                 %ind_ninflect=find(fitparamslist(:,2)<0.09*0.001)
%                 
%                 %plot inflection
%                 if disp_plot
%                     figure(figtrials)
%                     inflect_plot=[];
%                     inflect_plot(:,1)=inflect_2(:,1);
%                     inflect_plot(:,2)=inflect_2(:,2);
%                     %inflect_spk_plot(ind_ninflect,:)=nan;
%                     %plot_events_ch(inflect_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','--',2,'');
%                     plot_events_ch(inflect_plot,[],vshift_spk,range_spk,info,hdlfig,'n','--',2,'');
%                 end
                
                %%
                %%%%%
                %inflection 3 point between elb_1 and peak 
                anchor_b=max(-150,elb_1(:,1)-150)+info.aligntime;
                anchor_e=peak(:,1)+info.aligntime;
                min_bins=10;
                inflect_3=nan(info.nchannels,3);
                resnorm_3=cell(1,info.nchannels);
                fitparams_3=nan(info.nchannels,4);
                resnormfit_3=nan(info.nchannels,2);
                for ch=1:info.nchannels
                    if ~isnan(anchor_b(ch)) & ~isnan(anchor_e(ch)) & (anchor_e(ch)-anchor_b(ch)+1)>(2*min_bins)
                        %ch
                        [inflect,resnorm,resnorm_fit,fitparams]=get_inflection_2pwlr(trials_spk_n_avg_boot(ch,:),anchor_b(ch),anchor_e(ch),min_bins);
                        inflect_3(ch,1)=inflect(1)-info.aligntime;%correction for timing
                        inflect_3(ch,2:3)=inflect(2:3);
                        resnorm_3{ch}=resnorm;
                        fitparams_3(ch,:)=fitparams;
                        resnormfit_3(ch,:)=resnorm_fit;
                        
                        %plot resnorm
%                         if ~isempty(resnorm_3{ch}) %trial_sel was long enough
%                             resnorm_plot=nan(1,info.aligntime+length(resnorm_3{ch})-1);
%                             resnorm_plot(1,info.aligntime:end)=resnorm_3{ch};
%                             resnorm_plot=resnorm_plot/max(resnorm_plot);
%                             figure(figresnorm);subplot(1,3,3);hold on;
%                             plot(resnorm_plot,'color',colorlist(ch,:),'linewidth',1);
%                             title('Residuals of 2pwlr for inflect 2');xlabel('Time (ms)');ylabel('R2');
%                         end
                    end
                end
                
                %channel without inflection point
                %ind_ninflect=find(fitparamslist(:,2)<0.09*0.001)
                
                %plot inflection
                if disp_plot
                    figure(figtrials)
                    inflect_plot=[];
                    inflect_plot(:,1)=inflect_3(:,1);
                    inflect_plot(:,2)=inflect_3(:,2);
                    %inflect_spk_plot(ind_ninflect,:)=nan;
                    %plot_events_ch(inflect_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','--',2,'');
                    plot_events_ch(inflect_plot,[],vshift_spk,range_spk,info,hdlfig,'n','-',3,'');
                end
                
                %%%%%
                %keep bootstrapped estimations
                trials_boot(b,:,:)=trials_spk_n_avg_boot;
                fitparamsd_boot(b,:,:)=fitparams_d;
                peak_boot(b,:,:)=peak;
                elb1_boot(b,:,:)=elb_1;
                elb2_boot(b,:,:)=nan;%elb_2;
                inflect1_boot(b,:,:)=inflect_1;
                inflect2_boot(b,:,:)=nan;%inflect_2;
                inflect3_boot(b,:,:)=inflect_3;
                resnorm1_boot{b}=resnorm_1;
                resnorm2_boot{b}=nan;%resnorm_2;
                resnorm3_boot{b}=resnorm_3;
                fitparams1_boot(b,:,:)=fitparams_1;
                fitparams2_boot(b,:,:)=nan;%fitparams_2;
                fitparams3_boot(b,:,:)=fitparams_3;
                resnormfit1_boot(b,:,:)=resnormfit_1;
                resnormfit2_boot(b,:,:)=nan;%resnormfit_2;
                resnormfit3_boot(b,:,:)=resnormfit_3;
                
                if disp_plot
                    %pause
                end
                
            end%boot
            %pause
            %save results
            if ~disp_plot
                suffixe='_detrend_500_300_2';%'_detrend_500_300';'_newdfix';%'_newd';%'_new';%'_inflect2_150';%''
                namesave=[save_path 'results_onsetbuildup_' info.datafile(1:end-4) '_nboot' num2str(nboot) suffixe];
                save(namesave, 'trials_boot','fitparamsd_boot','peak_boot','elb1_boot','elb2_boot',...
                    'inflect1_boot','inflect2_boot','inflect3_boot',...
                    'resnorm1_boot','resnorm2_boot','resnorm3_boot',...
                    'fitparams1_boot','fitparams2_boot','fitparams3_boot',...
                    'resnormfit1_boot','resnormfit2_boot','resnormfit3_boot');
                
                
                close all
            end
            
%             %%
%             %%%%%%%%%%%%%%%%%%
%             %plot avg boot
%             trials_boot_avg=squeeze(mean(trials_boot,1));
%             trials_boot_var=squeeze(std(trials_boot,1)/size(trials_boot,1));
%             
%             %plot avg and ste
%             figtrialsbootall=figure
%             hdlfig=subplot(1,1,1);hold on;
%             for b=1:nboot
%             [range_spk vshift_spk]=plot_trials(squeeze(trials_boot(b,:,:)),[],[],[],[],[],info,hdlfig,[],[],[]);
%             end
%             
%             %%
%             %plot avg and ste
%             figtrialsboot=figure
%             hdlfig=subplot(1,1,1);hold on;
%             [range_spk vshift_spk]=plot_trials(trials_boot_avg,[],[],[],[],[],info,hdlfig,[],[],[]);
%             plot_trials(trials_boot_avg+trials_boot_var,[],[],[],[],[],info,hdlfig,[],'-',1);
%             plot_trials(trials_boot_avg-trials_boot_var,[],[],[],[],[],info,hdlfig,[],'-',1);
%             
%             %plot peak
%             figure(figtrialsboot)
%             onset_plot=[];
%             onset_plot(:,1)=squeeze(mean(peak_boot(:,:,1),1));
%             onset_plot(:,2)=squeeze(mean(peak_boot(:,:,2),1))
%             plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','--',1,'k');
% 
%             %plot elb_1
%             figure(figtrialsboot)
%             onset_plot=[];
%             onset_plot(:,1)=squeeze(median(elb1_boot(:,:,1),1));
%             onset_plot(:,2)=squeeze(median(elb1_boot(:,:,2),1))
%             plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','--',1,'');
%             
%             %plot elb_2
%             figure(figtrialsboot)
%             onset_plot=[];
%             onset_plot(:,1)=squeeze(median(elb2_boot(:,:,1),1));
%             onset_plot(:,2)=squeeze(median(elb2_boot(:,:,2),1))
%             plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','-',1,'');
%             
%             
%             %%
%             %plot elbows
%             figelbs=figure('Position',[1 100 scrsz(3) scrsz(4)-500]);
%             hold on;
%             
%             %plot onset elb_1
%             elb1_boot_avg=squeeze(mean(elb1_boot(:,:,1),1))
%             elb1_boot_var=squeeze(std(elb1_boot(:,:,1),1));%/size(elb1_boot,1));
%             subplot(1,3,1);hold on;
%             errorbar(1:info.nchannels,elb1_boot_avg,elb1_boot_var);
%             subplot(1,3,2);hold on;
%             plot(1:info.nchannels,elb1_boot_avg);
%             boxplot(squeeze(elb1_boot(:,:,1)));
%             title('Elbow 1')
%             
%             %plot onset elb_2
%             elb2_boot_avg=squeeze(mean(elb2_boot(:,:,1),1))
%             elb2_boot_var=squeeze(std(elb2_boot(:,:,1),1));%/size(elb2_boot,1));
%             subplot(1,3,1);hold on;
%             errorbar(1:info.nchannels,elb2_boot_avg,elb2_boot_var);
%             subplot(1,3,3);hold on;
%             plot(1:info.nchannels,elb2_boot_avg);
%             boxplot(squeeze(elb2_boot(:,:,1)));
%             title('Elbows 2')
%             
%             
%             subplot(1,3,1)
%             axis([1 info.nchannels -280 20]);
%             xlabel('Channel');ylabel('Time from saccade onset (ms)')
%             title('Elbows 1 and 2')
%             legend({'Elbow 1' 'Elbow 2'},'Location','Southeast')
%             
%             %subplot(1,2,2)
%             %axis([1 info.nchannels -280 20]);            
%             %xlabel('Channel');ylabel('Time from saccade onset (ms)')
%             %title('Elbows 1 and 2')
%             %legend({'Elbow 1' 'Elbow 2'},'Location','Southeast')
%             
%             %%
%             %plot distribution of Elbows' onset
%             %distr elb_1
%             figdistr1=figure;
%             figdistr2=figure;
%             
%             for ch=1:info.nchannels
%                 figure(figdistr1)
%                 subplot(4,4,ch);hold on;
%                 histo=reshape(elb1_boot(:,ch,1),1,numel(elb1_boot(:,ch,1)));
%                 edges=[-300:1:50];
%                 hist=histc(histo,edges);
%                 bar(edges,hist,'histc')
%                 axis([-300 50 0 10])
%                 axis tight
%                 if ch==1
%                     title(['Elbow 1 ' num2str(ch)])
%                 else
%                     title(num2str(ch))
%                 end
%                 
%                 %distr elb_2
%                 figure(figdistr2)
%                 subplot(4,4,ch);hold on;
%                 histo=reshape(elb2_boot(:,ch,1),1,numel(elb2_boot(:,ch,1)));
%                 edges=[-350:1:50];
%                 hist=histc(histo,edges);
%                 bar(edges,hist,'histc')
%                 axis([-300 50 0 10])
%                 axis tight 
%                 if ch==1
%                     title(['Elbow 2 ' num2str(ch)])
%                 else
%                     title(num2str(ch))
%                 end
%                
%             end
%             
%             %%
%             %plot amplitude differences
%             figamps=figure('Position',[1 100 scrsz(3) scrsz(4)-500]);
%             hold on;
%             
%             %plot amplitude diff baseline elb_1
%             elb1_boot_avg=squeeze(mean(elb1_boot(:,:,2),1))
%             elb1_boot_var=squeeze(std(elb1_boot(:,:,2),1));%/size(elb1_boot,1));
%             subplot(1,3,1);hold on;
%             errorbar(1:info.nchannels,elb1_boot_avg,elb1_boot_var);
%             axis([1 info.nchannels 0 max(elb1_boot_avg)+max(elb1_boot_var)]);      
%             xlabel('Channel');ylabel('FR (spk/s)')
%             title('Elbow1-Baseline')
%             
%             
%             %plot amplitude diff elb_1 elb_2
%             diffamp=elb2_boot(:,:,2)-elb1_boot(:,:,2);
%             diffamp_avg=squeeze(mean(diffamp,1))
%             diffamp_var=squeeze(std(diffamp,1));%/size(diffamp,1));
%             subplot(1,3,2);hold on;
%             errorbar(1:info.nchannels,diffamp_avg,diffamp_var);
%             axis([1 info.nchannels 0 (max(diffamp_avg)+max(diffamp_var))]);      
%             xlabel('Channel');ylabel('FR (spk/s)')
%             title('Elbow2-Elbow1')
%             
%             %plot amplitude diff elb_2 peak
%             diffamp=peak_boot(:,:,2)-elb1_boot(:,:,2);
%             diffamp_avg=squeeze(mean(diffamp,1))
%             diffamp_var=squeeze(std(diffamp,1));%/size(diffamp,1));
%             subplot(1,3,3);hold on;
%             errorbar(1:info.nchannels,diffamp_avg,diffamp_var);
%             axis([1 info.nchannels 0 (max(diffamp_avg)+max(diffamp_var))]);      
%             xlabel('Channel');ylabel('FR (spk/s)')
%             title('Peak-Elbow2')
%             
%             
% 
%             %%
%             %save in powerpoint
%             %             %start .pptx file
%             %             savepptx=1;
%             %             if savepptx
%             %                 isOpen  = exportToPPTX();
%             %                 if ~isempty(isOpen),
%             %                     % If PowerPoint already started, then close first and then open a new one
%             %                     exportToPPTX('close');
%             %                 end
%             %
%             %                 exportToPPTX('new','Dimensions',[12 6], ...
%             %                     'Title','SC onset buildup', ...
%             %                     'Author','Corentin', ...
%             %                     'Subject','Automatically generated PPTX file from output of analysis_onset_buildup_avg.m', ...
%             %                     'Comments',' ');
%             %
%             %                 %tmp filename
%             %                 file=[save_path 'onset_buildup_tmp' '.' figtype];
%             %
%             %                 %pptx filename
%             %                 file_pptx=[save_path 'SC_onsetbuildup'];
%             %
%             %             end
% 
%             
%             %%
%             if savepptx,
%                 savetopptx(figtrialsboot,file,figtype,{info.datafile ;' trials bootstrapped'});
%                 savetopptx(figelbs,file,figtype,{info.datafile ;' elbows'});
%                 savetopptx(figdistr1,file,figtype,{info.datafile ; ' distribution of elbow 1'});
%                 savetopptx(figdistr2,file,figtype,{info.datafile ;' distribution of elbow 2'});
%                 savetopptx(figamps,file,figtype,{info.datafile ;' amplitude difference of elbows'});
%             end
%             
%             %             %%
%             %             if savepptx
%             %                 %close .pptx
%             %                 newFile = exportToPPTX('saveandclose',filepptx)
%             %             end
% 
%             close all
%             
%             
%            
%          
%             

% %     pause
        end
    end
end

%%
if savepptx
    %close .pptx
    newFile = exportToPPTX('saveandclose',filepptx)
end


