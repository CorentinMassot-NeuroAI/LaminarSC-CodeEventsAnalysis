%function analysis_onset_buildup_avg

%function analysis_onset_buildup_avg
%   Analysis of onset of build-up and burst activity of average spiking activity recorded with a
%   laminar probe (LMA)
%
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 11/29/2016 last modified 01/23/2018
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



%alignement
%alignlist={'no' 'targ' 'go' 'sacc'};
alignlist={'sacc'};


%windows of analysis (do not change)
%for first results
%wind_sacc=[-400 250];%[-400 100]

%for detrend_500_300
wind_sacc=[-600 100];

%windows baseline
%wt=100;
%wind_targ_pburst_bsl=[-50-wt -50 ];
%wind_bsl_sacc=[-200 -100];%[-200 150];


%vshift
vshift_spk=100;
vshift_lfp=30;%29;


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

allpeak={};allelb1={};allelb2={};allinflect1={};allinflect2={};allinflect3={};allaccum={};allburst={};
classif=struct;allscat={};
data=[];
info=[];
dd=0;
for d=dlist(1:end)
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
    
    %targets index
    targs_ind=get_targsindex(targslist,info);
    
    %target tuning (after compute_tuning)
    info.targ_tuning=data(1).offline.targ_tuning;
    info.targ=info.targ_tuning;
    
    %align 'sacc'
    info.align=alignlist{1};
    
    %aligntime
    info.aligntime=abs(min(wind_sacc));
    
    %sacc bursts significance
    thresh_ratios=0.15;%0.15
    thresh_surprises=4;
    ratios_sacc=data(1).offline.sacc_pburst_ratio(info.targ_tuning,:)>thresh_ratios;
    surprises_sacc=data(1).offline.sacc_pburst_msurprises(info.targ_tuning,:)>thresh_surprises;
    bsignif_sacc=data(1).offline.sacc_bsignif;
    bthresh_sacc=data(1).offline.sacc_bthresh_trials';
    %burst_bsignif=(ratios_sacc & surprises_sacc & bsignif_sacc & bthresh_sacc);
    burst_bsignif=(bsignif_sacc & bthresh_sacc);
    
    burst_bsignif=double(burst_bsignif);
    burst_bsignif(find(burst_bsignif==0))=nan;
    
    %ntrials
    info.ntrials=0;
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %load results onsetbuildup for average trials
    %     nboot=1;
    %     nameload=[save_path 'results_onsetbuildup_' info.datafile(1:end-4) '_nboot' num2str(nboot)];
    %     load(nameload)
    %     %'trials_boot','peak_boot','elb1_boot','elb2_boot','inflect1_boot','inflect2_boot','inflect3_boot','resnorm1_boot','resnorm2_boot','resnorm3_boot');
    %
    %     trials_avg=squeeze(trials_boot(1,:,:));
    %     peak_avg=squeeze(peak_boot(1,:,:));
    %     elb1_avg=squeeze(elb1_boot(1,:,:));
    %     elb2_avg=squeeze(elb2_boot(1,:,:));
    %     inflect1_avg=squeeze(inflect1_boot(1,:,1:2));
    %     inflect2_avg=squeeze(inflect2_boot(1,:,1:2));
    %     inflect3_avg=squeeze(inflect3_boot(1,:,1:2));
    %     %resnorm1_avg=resnorm1_boot{1};
    %     %resnorm2_avg=resnorm2_boot{1};
    %     %resnorm3_avg=resnorm3_boot{1};
    %
    %%
    %%%%%%%%%%%%%%%%%%
    %     %plot avg and ste
    %     figtrials=figure('Position',[scrsz(3)/6 150 scrsz(3)/1.5 scrsz(4)-300]);
    %     hdlfig=subplot(1,1,1);hold on;
    %     [range_spk vshift_spk]=plot_trials(trials_avg,[],[],[],[],[],info,hdlfig,[],[],[]);
    %     %plot_trials(trials_avg+trials_var,[],[],[],[],[],info,hdlfig,[],'-',1);
    %     %plot_trials(trials_avg-trials_var,[],[],[],[],[],info,hdlfig,[],'-',1);
    %
    %     %plot peak
    %     plot_events_ch(peak_avg.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','--',1,'k');
    %
    %     %plot elb_1
    %     plot_events_ch(elb1_avg.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','--',1,'');
    %
    %     %plot elb_2
    %     plot_events_ch(elb2_avg.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','-',1,'');
    %
    %     %plot inflect 1
    %     plot_events_ch(inflect1_avg.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','-',2,'');
    %
    %     %plot inflect 2
    %     plot_events_ch(inflect2_avg.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','--',2,'');
    %
    %     %plot inflect 3
    %     plot_events_ch(inflect3_avg.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','-',2,'k');
    
    
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load results onsetbuildup for bootstrapped trials
    nboot=100;%551;%100;%20;
    %suffixe='_newd';%'_newdfix';%fix point for 2pwlr    %'_new';%'_inflect2_150';%''
    suffixe='_detrend_500_300';%'_detrend_500_300_2';%change detrending and elb1 windows
    
    nameload=[save_path 'results_onsetbuildup_' info.datafile(1:end-4) '_nboot' num2str(nboot) suffixe]
    load(nameload)
    
    %'trials_boot','peak_boot','elb1_boot','elb2_boot',...
    %'inflect1_boot','inflect2_boot','inflect3_boot',...
    %'resnorm1_boot','resnorm2_boot','resnorm3_boot',...
    %'fitparams1_boot','fitparams2_boot','fitparams3_boot',...
    %'resnormfit1_boot','resnormfit2_boot','resnormfit3_boot');
    
    %%
    %%%%%%%%%%%%%
    %plot all bootstrapped trials
    %     %plot all trials boot
    %     figtrialsbootall=figure
    %     hdlfig=subplot(1,1,1);hold on;
    %     for b=1:nboot
    %         [range_spk vshift_spk]=plot_trials(squeeze(trials_boot(b,:,:)),[],[],[],[],[],info,hdlfig,[],[],[]);
    %     end
    
    %%
    %%%%%%%%%%%%%
    %plot trials avg and ste
    trials_boot_avg=squeeze(mean(trials_boot,1));
    trials_boot_var=squeeze(std(trials_boot,1)/size(trials_boot,1));
    figtrialsboot=figure;
    hdlfig=subplot(1,1,1);hold on;
    [range_spk vshift_spk]=plot_trials(trials_boot_avg,[],[],[],[],[],info,hdlfig,[],[],[]);
    plot_trials(trials_boot_avg+trials_boot_var,[],[],[],[],[],info,hdlfig,[],'-',1);
    plot_trials(trials_boot_avg-trials_boot_var,[],[],[],[],[],info,hdlfig,[],'-',1);
    figtrialsboot2=figure;
    hdlfig=subplot(1,1,1);hold on;
    [range_spk vshift_spk]=plot_trials(trials_boot_avg,[],[],[],[],[],info,hdlfig,[],[],[]);
    plot_trials(trials_boot_avg+trials_boot_var,[],[],[],[],[],info,hdlfig,[],'-',1);
    plot_trials(trials_boot_avg-trials_boot_var,[],[],[],[],[],info,hdlfig,[],'-',1);
    figtrialsboot3=figure;
    hdlfig=subplot(1,1,1);hold on;
    [range_spk vshift_spk]=plot_trials(trials_boot_avg,[],[],[],[],[],info,hdlfig,[],[],[]);
    plot_trials(trials_boot_avg+trials_boot_var,[],[],[],[],[],info,hdlfig,[],'-',1);
    plot_trials(trials_boot_avg-trials_boot_var,[],[],[],[],[],info,hdlfig,[],'-',1);
    
    %     %%%%%%%%%%%%%
    %     %plot all bootstrapped onsets
    %     for b=1:nboot
    %         %plot peak
    %         figure(figtrialsboot)
    %         onset_plot=[];onset_plot(:,1)=squeeze(peak_boot(b,:,1));onset_plot(:,2)=squeeze(peak_boot(b,:,2));
    %         plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','--',1,'k');
    %         %plot elb_1
    %         onset_plot=[];onset_plot(:,1)=squeeze(elb1_boot(b,:,1));onset_plot(:,2)=squeeze(elb1_boot(b,:,2));
    %         plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','--',1,'');
    %         %plot elb_2
    %         onset_plot=[];onset_plot(:,1)=squeeze(elb2_boot(b,:,1));onset_plot(:,2)=squeeze(elb2_boot(b,:,2));
    %         plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','-',1,'');
    %
    %         %plot inflect 1
    %         figure(figtrialsboot2)
    %         onset_plot=[];onset_plot(:,1)=squeeze(inflect1_boot(b,:,1));onset_plot(:,2)=squeeze(inflect1_boot(b,:,2));
    %         plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','-',2,'');
    %         %plot inflect 2
    %         onset_plot=[];onset_plot(:,1)=squeeze(inflect2_boot(b,:,1));onset_plot(:,2)=squeeze(inflect2_boot(b,:,2));
    %         plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','--',2,'');
    %         %plot inflect 3
    %         onset_plot=[];onset_plot(:,1)=squeeze(inflect3_boot(b,:,1));onset_plot(:,2)=squeeze(inflect3_boot(b,:,2));
    %         plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','-',3,'');
    %
    %     end
    
    
    %%
    %%%%%%%%%%%%%
    %compute CI
    ci_level=0.95;%0.8;%95%
    p_level=0.01;
    peak_ci=[];elb1_ci=[];elb2_ci=[];inflect1_ci=[];inflect2_ci=[];inflect3_ci=[];
    fitslopes1_boot=[];fitslopes1_1_ci=[];fitslopes1_2_ci=[];fitslopes1_diff_ci=[];tsignif_fitslopes1=[];
    %fitslopes2_boot=[];fitslopes2_1_ci=[];fitslopes2_2_ci=[];
    fitslopes3_boot=[];fitslopes3_1_ci=[];fitslopes3_2_ci=[];fitslopes3_diff_ci=[];tsignif_fitslopes3=[];
    tsignif_elb1_inflect1=[];%tsignif_elb1_inflect3=[];
    
    elb2_boot=nan(size(elb1_boot));
    inflect2_boot=nan(size(elb1_boot));
    
    
    for ch=1:info.nchannels
        peak_ci(ch,1:3)=get_ci2(peak_boot(:,ch,1),ci_level);
        elb1_ci(ch,1:3)=get_ci2(elb1_boot(:,ch,1),ci_level);
        %elb1 ci normalization search window of elb1 is 250ms (see compute_onset_buildup_avg)
        %elb1_ci(ch,4:5)=(elb1_ci(ch,2:3)-elb1_ci(ch,1))./250+elb1_ci(ch,1);
        %elb1 ci normalization search window of elb1 is 300ms (see compute_onset_buildup_avg)
        elb1_ci(ch,4:5)=(elb1_ci(ch,2:3)-elb1_ci(ch,1))./300+elb1_ci(ch,1);
      
        elb2_ci(ch,1:3)=get_ci2(elb2_boot(:,ch,1),ci_level);
        inflect1_ci(ch,1:3)=get_ci2(inflect1_boot(:,ch,1),ci_level);
        %inflect1 ci normalization
        inflect1_ci(ch,4:5)=(inflect1_ci(ch,2:3)-inflect1_ci(ch,1))./mean(abs(elb1_boot(:,ch,1)-peak_boot(:,ch,1)))+inflect1_ci(ch,1);
        inflect2_ci(ch,1:3)=get_ci2(inflect2_boot(:,ch,1),ci_level);
        inflect3_ci(ch,1:3)=get_ci2(inflect3_boot(:,ch,1),ci_level);
        %inflect3 ci normalization search window of inflect3 is max(-150,elb_1(:,1)-150) to peak (see compute_onset_buildup_avg)
        inflect3_ci(ch,4:5)=(inflect3_ci(ch,2:3)-inflect3_ci(ch,1))./mean(abs(max(-150,elb1_boot(:,ch,1)-150)-peak_boot(:,ch,1)))+inflect3_ci(ch,1);
        
        %linear fit slopes
        %fitslopes1_boot(:,ch,1)=get_fitangles(squeeze(fitparams1_boot(:,ch,1:2)));
        %fitslopes1_boot(:,ch,2)=get_fitangles(squeeze(fitparams1_boot(:,ch,3:4)));
        fitslopes1_boot(:,ch,1)=squeeze(fitparams1_boot(:,ch,2));
        fitslopes1_boot(:,ch,2)=squeeze(fitparams1_boot(:,ch,4));
        fitslopes1_1_ci(ch,:)=get_ci2(fitslopes1_boot(:,ch,1),ci_level);
        fitslopes1_2_ci(ch,:)=get_ci2(fitslopes1_boot(:,ch,2),ci_level);
        fitslopes1_diff_ci(ch,:)=get_ci2(fitslopes1_boot(:,ch,2)-fitslopes1_boot(:,ch,1),ci_level);
        %[H1 p1]=get_testsignif(fitslopes1_boot(:,ch,1),fitslopes1_boot(:,ch,2),p_level);
        tsignif_fitslopes1(ch,1)=get_testsignifci(fitslopes1_1_ci(ch,:),fitslopes1_2_ci(ch,:),0);
        %for example [28] in paper
        %tsignif_fitslopes1(ch,1)=get_testsignifci(fitslopes1_1_ci(ch,:),fitslopes1_2_ci(ch,:),-0.5);

        %         fitslopes2_boot(:,ch,1)=get_fitangles(squeeze(fitparams2_boot(:,ch,1:2)));
        %         fitslopes2_boot(:,ch,2)=get_fitangles(squeeze(fitparams2_boot(:,ch,3:4)));
        %         fitslopes2_1_ci(ch,:)=get_ci2(fitslopes2_boot(:,ch,1),ci_level);
        %         fitslopes2_2_ci(ch,:)=get_ci2(fitslopes2_boot(:,ch,2),ci_level);
        
        %fitslopes3_boot(:,ch,1)=get_fitangles(squeeze(fitparams3_boot(:,ch,1:2)));
        %fitslopes3_boot(:,ch,2)=get_fitangles(squeeze(fitparams3_boot(:,ch,3:4)));
        fitslopes3_boot(:,ch,1)=squeeze(fitparams3_boot(:,ch,2));
        fitslopes3_boot(:,ch,2)=squeeze(fitparams3_boot(:,ch,4));
        fitslopes3_1_ci(ch,:)=get_ci2(fitslopes3_boot(:,ch,1),ci_level);
        fitslopes3_2_ci(ch,:)=get_ci2(fitslopes3_boot(:,ch,2),ci_level);
        fitslopes3_diff_ci(ch,:)=get_ci2(fitslopes3_boot(:,ch,2)-fitslopes3_boot(:,ch,1),ci_level);
        %tsignif_fitslopes3(ch,1:2)=get_testsignif(fitslopes3_boot(:,ch,1),fitslopes3_boot(:,ch,2),p_level);
        tsignif_fitslopes3(ch,1)=get_testsignifci(fitslopes3_1_ci(ch,:),fitslopes3_2_ci(ch,:),0);
        %for example [28] in paper
        %tsignif_fitslopes3(ch,1)=get_testsignifci(fitslopes3_1_ci(ch,:),fitslopes3_2_ci(ch,:),-0.5);
        
        
        %ranksum tests between elb1 and inflect1
        %tsignif_elb1_inflect1(ch,:)=get_testsignif(elb1_boot(:,ch,1),inflect1_boot(:,ch,1),p_level);
        tsignif_elb1_inflect1(ch,:)=get_testsignifci(elb1_ci(ch,[1 4 5]),inflect1_ci(ch,[1 4 5]),0);
        %to reproduce old results in power pointSC_onsetbuildup-VG-all-22-Feb-2018.pptx
        %tsignif_elb1_inflect1(ch,:)=get_testsignifci(elb1_ci(ch,:),inflect1_ci(ch,:),0);
        %for example [28] in paper
        %tsignif_elb1_inflect1(ch,:)=get_testsignifci(elb1_ci(ch,1:3),inflect1_ci(ch,1:3),10);
       
        %tsignif_elb1_inflect3(ch,:)=get_testsignif(elb1_boot(:,ch,1),inflect3_boot(:,ch,1),p_level);
        
    end
    
    %%
        %%%%%%%%%%%%%
        %plot onsets mean and ci
        linetype_ci={'-' '--' '--','--','--'};
        linewidth_ci=[4 2 2 2 2];
    
        %peak mean and ci
       figure(figtrialsboot)
       color_ci={'k' 'k' 'k'};
        for c=1:3
            onset_plot=[];onset_plot(:,1)=peak_ci(:,c);
            for ch=1:info.nchannels
                if ~isnan(onset_plot(ch,1)) & round(onset_plot(ch,1))+info.aligntime <= size(trials_boot_avg,2)
                    onset_plot(ch,2)=trials_boot_avg(ch,round(onset_plot(ch,1))+info.aligntime);%peak_avg(:,2);
                end
            end
            plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n',linetype_ci{c},linewidth_ci(c),color_ci{c});
        end
        %elb1 mean and ci
        figure(figtrialsboot)
        color_ci={'b' 'b' 'b'};
        for c=1:3
            onset_plot=[];onset_plot(:,1)=elb1_ci(:,c);
            for ch=1:info.nchannels
                if ~isnan(onset_plot(ch,1))
                    onset_plot(ch,2)=trials_boot_avg(ch,round(onset_plot(ch,1))+info.aligntime);%elb1_avg(:,2);
                end
            end
            plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n',linetype_ci{c},linewidth_ci(c),color_ci{c});
        end
    %     %elb2 mean and ci
    %     figure(figtrialsboot)
    %     color_ci={'g' 'g' 'g'};
    %     for c=1:3
    %         onset_plot=[];onset_plot(:,1)=elb2_ci(:,c);
    %         for ch=1:info.nchannels
    %             if ~isnan(onset_plot(ch,1))
    %                 onset_plot(ch,2)=trials_boot_avg(ch,round(onset_plot(ch,1))+info.aligntime);%elb2_avg(:,2);
    %             end
    %         end
    %         plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n',linetype_ci{c},linewidth_ci(c),color_ci{c});
    %     end
    
        %inflect1 mean and ci
       figure(figtrialsboot)
        color_ci={'r' 'r' 'r' 'r' 'r'};
        for c=1:3
            onset_plot=[];onset_plot(:,1)=inflect1_ci(:,c);
            for ch=1:info.nchannels
                if ~isnan(onset_plot(ch,1))
                    onset_plot(ch,2)=trials_boot_avg(ch,round(onset_plot(ch,1))+info.aligntime);%inflect1_avg(:,2);
                end
            end
            plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n',linetype_ci{c},linewidth_ci(c),color_ci{c});
        end
    
%         %inflect1 mean and ci
%         figure(figtrialsboot2)
%         color_ci={'r' 'r' 'r' 'r' 'r'};
%         for c=1:3
%             onset_plot=[];onset_plot(:,1)=inflect1_ci(:,c);
%             for ch=1:info.nchannels
%                 if ~isnan(onset_plot(ch,1))
%                     onset_plot(ch,2)=trials_boot_avg(ch,round(onset_plot(ch,1))+info.aligntime);%inflect1_avg(:,2);
%                 end
%             end
%             plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n',linetype_ci{c},linewidth_ci(c),color_ci{c});
%         end
    
    
%         %inflect2 mean and ci
%         figure(figtrialsboot2)
%         color_ci={'m' 'm' 'm'};
%         for c=1:3
%             onset_plot=[];onset_plot(:,1)=inflect2_ci(:,c);
%             for ch=1:info.nchannels
%                 if ~isnan(onset_plot(ch,1))
%                     onset_plot(ch,2)=trials_boot_avg(ch,round(onset_plot(ch,1))+info.aligntime);%inflect2_avg(:,2);
%                 end
%             end
%             plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n',linetype_ci{c},linewidth_ci(c),color_ci{c});
%         end
    
%inflect3 mean and ci
figure(figtrialsboot3)
color_ci={'c' 'c' 'c'};
for c=1:3
    onset_plot=[];onset_plot(:,1)=inflect3_ci(:,c);
    for ch=1:info.nchannels
        if ~isnan(onset_plot(ch,1))
            onset_plot(ch,2)=trials_boot_avg(ch,round(onset_plot(ch,1))+info.aligntime);%inflect3_avg(:,2);
        end
    end
    plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n',linetype_ci{c},linewidth_ci(c),color_ci{c});
end
    
    %     %%
    %     %%%%%%%%%%%%%
    %     %plot stats for session
    %     %NOTE: here theshold on CI on both sides
    %     %     %ci_level=0.8
    %     %     peak_t=20;
    %     %     elb1_t=30;
    %     %     elb2_t=30;
    %     %     inflect1_t=20;
    %     %     inflect1n_t=1.2;%1.3
    %     %     inflect2_t=50;%30;
    %     %     inflect3_t=20;
    %
    %     %ci_level=0.95
    %     peak_t=20;
    %     elb1_t=40;
    %     elb1n_t=0.15;
    %     elb2_t=40;
    %     inflect1_t=20;
    %     inflect1n_t=0.15;%1.3
    %     inflect2_t=50;%30;
    %     inflect3_t=20;
    %     inflect3n_t=0.15;
    %
    %     %     peak_t=inf;
    %     %     elb1_t=inf;
    %     %     elb2_t=inf;
    %     %     inflect1_t=inf;
    %     %     inflect1n_t=inf;
    %     %     inflect2_t=inf;%30;
    %     %     inflect3_t=inf;
    %
    %     %     peak_t=0;
    %     %     elb1_t=0;
    %     %     elb2_t=0;
    %     %     inflect1_t=0;
    %     %     inflect1n_t=0;
    %     %     inflect2_t=0;
    %     %     inflect3_t=0;
    %
    %
    %     %plot ci of onsets
    %     figelbs=figure('Position',[1 100 scrsz(3) scrsz(4)-500]);
    %     hold on;
    %     %peak
    %     subplot(1,7,1);hold on;
    %     errorbar(1:info.nchannels,zeros(16,1),abs(peak_ci(:,2)-peak_ci(:,1)),abs(peak_ci(:,3)-peak_ci(:,1)));
    %     h=line([0 17],[peak_t peak_t]);set(h,'color','r');
    %     h=line([0 17],[-peak_t -peak_t]);set(h,'color','r');
    %     grid minor;axis([0 17 -100 100]);
    %     xlabel('Channel');ylabel('CI (ms)')
    %     title('peak')
    %     %elb1
    %     subplot(1,7,2);hold on;
    %     errorbar(1:info.nchannels,zeros(16,1),abs(elb1_ci(:,2)-elb1_ci(:,1)),abs(elb1_ci(:,3)-elb1_ci(:,1)));
    %     h=line([0 17],[elb1_t elb1_t]);set(h,'color','r');
    %     h=line([0 17],[-elb1_t -elb1_t]);set(h,'color','r');
    %     grid minor;axis([0 17 -100 100]);
    %     xlabel('Channel');ylabel('CI (ms)')
    %     title('elb1')
    %     %elb1 normalized
    %     subplot(1,7,3);hold on;
    %     hdl_e=errorbar(1:info.nchannels,zeros(16,1),abs(elb1_ci(:,4)-elb1_ci(:,1)),abs(elb1_ci(:,5)-elb1_ci(:,1)));
    %     set(hdl_e,'color','m')
    %     h=line([0 17],[elb1n_t elb1n_t]);set(h,'color','r');
    %     h=line([0 17],[-elb1n_t -elb1n_t]);set(h,'color','r');
    %     grid minor;axis([0 17 -0.5 0.5]);
    %     xlabel('Channel');ylabel('CI (ms)')
    %     title('elb1 normalized')
    %     %     %elb2
    %     %     subplot(1,6,4);hold on;
    %     %     errorbar(1:info.nchannels,zeros(16,1),abs(elb2_ci(:,2)-elb2_ci(:,1)),abs(elb2_ci(:,3)-elb2_ci(:,1)));
    %     %     h=line([0 17],[elb2_t elb2_t]);set(h,'color','r');
    %     %     h=line([0 17],[-elb2_t -elb2_t]);set(h,'color','r');
    %     %     grid minor;axis([0 17 -100 100]);
    %     %     xlabel('Channel');ylabel('CI (ms)')
    %     %     title('elb2')
    %     %inflect1
    %     subplot(1,7,4);hold on;
    %     errorbar(1:info.nchannels,zeros(16,1),abs(inflect1_ci(:,2)-inflect1_ci(:,1)),abs(inflect1_ci(:,3)-inflect1_ci(:,1)));
    %     h=line([0 17],[inflect1_t inflect1_t]);set(h,'color','r');
    %     h=line([0 17],[-inflect1_t -inflect1_t]);set(h,'color','r');
    %     grid minor;axis([0 17 -100 100]);
    %     xlabel('Channel');ylabel('CI (ms)')
    %     title('inflect1')
    %     %inflect1 normalized
    %     subplot(1,7,5);hold on;
    %     hdl_e=errorbar(1:info.nchannels,zeros(16,1),abs(inflect1_ci(:,4)-inflect1_ci(:,1)),abs(inflect1_ci(:,5)-inflect1_ci(:,1)));
    %     set(hdl_e,'color','m')
    %     h=line([0 17],[inflect1n_t inflect1n_t]);set(h,'color','r');
    %     h=line([0 17],[-inflect1n_t -inflect1n_t]);set(h,'color','r');
    %     grid minor;axis([0 17 -0.5 0.5]);
    %     xlabel('Channel');ylabel('CI (ms)')
    %     title('inflect1 normalized')
    %     %     %inflect2
    %     %     subplot(1,6,7);hold on;
    %     %     errorbar(1:info.nchannels,zeros(16,1),abs(inflect2_ci(:,2)-inflect2_ci(:,1)),abs(inflect2_ci(:,3)-inflect2_ci(:,1)));
    %     %     h=line([0 17],[inflect2_t inflect2_t]);set(h,'color','r');
    %     %     h=line([0 17],[-inflect2_t -inflect2_t]);set(h,'color','r');
    %     %     grid minor;axis([0 17 -100 100]);
    %     %     xlabel('Channel');ylabel('CI (ms)')
    %     %     title('inflect2')
    %     %inflect3
    %     subplot(1,7,6);hold on;
    %     errorbar(1:info.nchannels,zeros(16,1),abs(inflect3_ci(:,2)-inflect3_ci(:,1)),abs(inflect3_ci(:,3)-inflect3_ci(:,1)));
    %     h=line([0 17],[inflect3_t inflect3_t]);set(h,'color','r');
    %     h=line([0 17],[-inflect3_t -inflect3_t]);set(h,'color','r');
    %     grid minor;axis([0 17 -100 100]);
    %     xlabel('Channel');ylabel('CI (ms)')
    %     title('inflect3')
    %     %inflect3 normalized
    %     subplot(1,7,7);hold on;
    %     hdl_e=errorbar(1:info.nchannels,zeros(16,1),abs(inflect3_ci(:,4)-inflect3_ci(:,1)),abs(inflect3_ci(:,5)-inflect3_ci(:,1)));
    %     set(hdl_e,'color','m')
    %     h=line([0 17],[inflect3n_t inflect3n_t]);set(h,'color','r');
    %     h=line([0 17],[-inflect3n_t -inflect3n_t]);set(h,'color','r');
    %     grid minor;axis([0 17 -0.5 0.5]);
    %     xlabel('Channel');ylabel('CI (ms)')
    %     title('inflect3 normalized')
    %
    %     %pause
    
    %%
    %%%%%%%%%%%%%
    %plot stats for session
    %NOTE: here threshold on range of ci (sum of both sides)
    %ci_level=0.95
    peak_t=40;
    elb1_t=80;
    %elb1n_t=0.4;%for results 
    elb1n_t=0.6;%0.4;%0.6;%0.3;0.45;
    %for example [28 ]in paper
    %elb1n_t=0.4 %0.5;%0.6;%0.3;0.45;
    
    elb2_t=80;
    inflect1_t=40;
    inflect1n_t=0.5;%0.6;%0.3;0.45;%1.3
    inflect1adiff_t=25;
    inflect2_t=100;%30;
    inflect3_t=40;
    inflect3n_t=0.5;%0.6;%0.3;0.45;
    inflect3adiff_t=25;
    
%     peak_t=inf;
%     elb1_t=inf;
%     elb1n_t=inf;
%     elb2_t=inf;
%     inflect1_t=inf;
%     inflect1adiff_t=-inf;
%     inflect1n_t=inf;
%     inflect2_t=inf;%30;
%     inflect3_t=inf;
%     inflect3n_t=inf;
%     inflect3adiff_t=-inf;
    
    %     peak_t=0;
    %     elb1_t=0;
    %     elb2_t=0;
    %     inflect1_t=0;
    %     inflect1n_t=0;
    %     inflect2_t=0;
    %     inflect3_t=0;
    
    
    %plot ci of onsets
    figelbs=figure('Position',[1 100 scrsz(3) scrsz(4)-500]);
    hold on;
    %peak
    subplot(1,6,1);hold on;
    errorbar(1:info.nchannels,zeros(16,1),zeros(info.nchannels,1),abs(peak_ci(:,2)-peak_ci(:,3)));
    h=line([0 17],[peak_t peak_t]);set(h,'color','r');
    grid minor;axis([0 17 0 200]);
    xlabel('Channel');ylabel('CI (ms)')
    title('peak')
%     %elb1
%     subplot(1,7,2);hold on;
%     errorbar(1:info.nchannels,zeros(16,1),zeros(info.nchannels,1),abs(elb1_ci(:,2)-elb1_ci(:,3)));
%     h=line([0 17],[elb1_t elb1_t]);set(h,'color','r');
%     grid minor;axis([0 17 0 200]);
%     xlabel('Channel');ylabel('CI (ms)')
%     title('elb1')
    %elb1 normalized
    subplot(1,6,2);hold on;
    hdl_e=errorbar(1:info.nchannels,zeros(16,1),zeros(info.nchannels,1),abs(elb1_ci(:,4)-elb1_ci(:,5)));
    set(hdl_e,'color','m')
    h=line([0 17],[elb1n_t elb1n_t]);set(h,'color','r');
    grid minor;axis([0 17 0 1]);
    xlabel('Channel');ylabel('CI (ms)')
    title('elb1 normalized')
   
    %     %elb2
    %     subplot(1,6,4);hold on;
    %     errorbar(1:info.nchannels,zeros(16,1),zeros(info.nchannels,1),abs(elb2_ci(:,2)-elb2_ci(:,3)));
    %     h=line([0 17],[elb2_t elb2_t]);set(h,'color','r');
    %     grid minor;axis([0 17 0 200]);
    %     xlabel('Channel');ylabel('CI (ms)')
    %     title('elb2')
   
%     %inflect1
%     subplot(1,7,4);hold on;
%     errorbar(1:info.nchannels,zeros(16,1),zeros(info.nchannels,1),abs(inflect1_ci(:,2)-inflect1_ci(:,3)));
%     h=line([0 17],[inflect1_t inflect1_t]);set(h,'color','r');
%     grid minor;axis([0 17 0 200]);
%     xlabel('Channel');ylabel('CI (ms)')
%     title('inflect1')
    %inflect1 normalized
    subplot(1,6,3);hold on;
    hdl_e=errorbar(1:info.nchannels,zeros(16,1),zeros(info.nchannels,1),abs(inflect1_ci(:,4)-inflect1_ci(:,5)));
    set(hdl_e,'color','m')
    h=line([0 17],[inflect1n_t inflect1n_t]);set(h,'color','r');
    grid minor;axis([0 17 0 1]);
    xlabel('Channel');ylabel('CI (ms)')
    title('inflect1 normalized')
   
    %     %inflect2
    %     subplot(1,6,7);hold on;
    %     errorbar(1:info.nchannels,zeros(16,1),zeros(info.nchannels,1),abs(inflect2_ci(:,2)-inflect2_ci(:,3)));
    %     h=line([0 17],[inflect2_t inflect2_t]);set(h,'color','r');
    %     grid minor;axis([0 17 0 200]);
    %     xlabel('Channel');ylabel('CI (ms)')
    %     title('inflect2')
    
    %     %inflect3
    %     subplot(1,7,6);hold on;
    %     errorbar(1:info.nchannels,zeros(16,1),zeros(info.nchannels,1),abs(inflect3_ci(:,2)-inflect3_ci(:,3)));
    %     h=line([0 17],[inflect3_t inflect3_t]);set(h,'color','r');
    %     grid minor;axis([0 17 0 200]);
    %     xlabel('Channel');ylabel('CI (ms)')
    %     title('inflect3')
    %inflect3 normalized
    subplot(1,6,5);hold on;
    hdl_e=errorbar(1:info.nchannels,zeros(16,1),zeros(info.nchannels,1),abs(inflect3_ci(:,4)-inflect3_ci(:,5)));
    set(hdl_e,'color','m')
    h=line([0 17],[inflect3n_t inflect3n_t]);set(h,'color','r');
    grid minor;axis([0 17 0 1]);
    xlabel('Channel');ylabel('CI (ms)')
    title('inflect3 normalized')
    
    
    %plot slopes of inflections
    %figslopes=figure('Position',[1 100 scrsz(3) scrsz(4)-500]);
    %hold on;
    %inflect1 slopes
    subplot(1,6,4);hold on;
    hdl_e=errorbar(1:info.nchannels,fitslopes1_1_ci(:,1),fitslopes1_1_ci(:,3)-fitslopes1_1_ci(:,1),fitslopes1_1_ci(:,2)-fitslopes1_1_ci(:,1));
    set(hdl_e,'color','b')
    hdl_e=errorbar(1:info.nchannels,fitslopes1_2_ci(:,1),fitslopes1_2_ci(:,3)-fitslopes1_2_ci(:,1),fitslopes1_2_ci(:,2)-fitslopes1_2_ci(:,1));
    set(hdl_e,'color','r')
    h=line([0 17],[0 0]);set(h,'color','k');
    grid minor;axis([0 17 -1 10]);
    xlabel('Channel');ylabel('slopes+CI (spk/s^2)')
    title('inflect1 slopes')
    %     %inflect1 diff slopes
    %     subplot(1,4,2);hold on;
    %     hdl_e=errorbar(1:info.nchannels,fitslopes1_diff_ci(:,1),fitslopes1_diff_ci(:,2),fitslopes1_diff_ci(:,3));
    %     set(hdl_e,'color','b')
    %     h=line([0 17],[0 0]);set(h,'color','k');
    %     h=line([0 17],[inflect1adiff_t inflect1adiff_t]);set(h,'color','r');
    %     grid minor;axis([0 17 -10 10]);
    %     xlabel('Channel');ylabel('slopes+CI (spk/s^2)')
    %     title('inflect1 diff slopes')
    %inflect3 slopes
    subplot(1,6,6);hold on;
    hdl_e=errorbar(1:info.nchannels,fitslopes3_1_ci(:,1),fitslopes3_1_ci(:,3)-fitslopes3_1_ci(:,1),fitslopes3_1_ci(:,2)-fitslopes3_1_ci(:,1));
    set(hdl_e,'color','b')
    hdl_e=errorbar(1:info.nchannels,fitslopes3_2_ci(:,1),fitslopes3_2_ci(:,3)-fitslopes3_2_ci(:,1),fitslopes3_2_ci(:,2)-fitslopes3_2_ci(:,1));
    set(hdl_e,'color','r')
    h=line([0 17],[0 0]);set(h,'color','k');
    grid minor;axis([0 17 -1 10]);
    xlabel('Channel');ylabel('slopes+CI (spk/s^2)')
    title('inflect3 slopes')
    %     %inflect3 diff slopes
    %     subplot(1,4,4);hold on;
    %     hdl_e=errorbar(1:info.nchannels,fitslopes3_diff_ci(:,1),fitslopes3_diff_ci(:,2),fitslopes3_diff_ci(:,3));
    %     set(hdl_e,'color','b')
    %     h=line([0 17],[0 0]);set(h,'color','k');
    %     h=line([0 17],[inflect3adiff_t inflect3adiff_t]);set(h,'color','r');
    %     grid minor;axis([0 17 -10 10]);
    %     xlabel('Channel');ylabel('slopes+CI (spk/s^2)')
    %     title('inflect3 diff slopes')
 
  
    %plot regression between CI and slopes diff
    figciadiff=figure;
    subplot(1,2,1);hold on;
    %regression 1
    xr=fitslopes1_diff_ci(:,1);
    yr=abs(inflect1_ci(:,4)-inflect1_ci(:,5));
    x_thresh=inflect1adiff_t;
    y_thresh=inflect1n_t;
    nnan=find(~isnan(xr) & ~isnan(yr));
    xr=xr(nnan);yr=yr(nnan);
    [p_fit rsq_fit yr_fit]=get_regressioncoefs(xr,yr,1);
    plot(xr,yr,'ok')
    plot(xr,yr_fit,'r-')
    %h=line([x_thresh x_thresh],[0 max(yr)+0.1]);set(h,'color','k','linestyle','--');
    h=line([0 max(xr)+1],[y_thresh y_thresh]);set(h,'color','k','linestyle','--');  
    axis([0 max(xr)+1 0 max(yr)+0.1])
    text(max(xr)-1,max(yr),['r2=' num2str(round(rsq_fit,2))])
    ylabel('CI normalized (ms)');xlabel('slopes diff (spk/s^2)');
    title('inflect 1')
    subplot(1,2,2);hold on;
    %regression 3
    xr=fitslopes3_diff_ci(:,1);
    yr=abs(inflect3_ci(:,4)-inflect3_ci(:,5));
    x_thresh=inflect3adiff_t;
    y_thresh=inflect3n_t;
    nnan=find(~isnan(xr) & ~isnan(yr));
    xr=xr(nnan);yr=yr(nnan);
    [p_fit rsq_fit yr_fit]=get_regressioncoefs(xr,yr,1);
    plot(xr,yr,'ok')
    plot(xr,yr_fit,'r-')
    %h=line([x_thresh x_thresh],[0 max(yr)+0.1]);set(h,'color','k','linestyle','--');
    h=line([0 max(xr)+1],[y_thresh y_thresh]);set(h,'color','k','linestyle','--');  
    axis([0 max(xr)+1 0 max(yr)+0.1])
    text(max(xr)-1,max(yr),['r2=' num2str(round(rsq_fit,2))])
    ylabel('CI normalized (ms)');xlabel('slopes diff (spk/s^2)');
    title('inflect 3')
    
    tsignif_fitslopes1'
    tsignif_fitslopes3'
    %pause
    
    %     %%
    %     %%%%%%%%%%%%%
    %     %save in powerpoint
    %     if savepptx,
    %         savetopptx(figtrialsboot,file,figtype,{info.datafile ;' Elbow 1 and 2'});
    %         savetopptx(figtrialsboot2,file,figtype,{info.datafile ;' Inflection 1 and 2'});
    %         %savetopptx(figtrialsboot3,file,figtype,{info.datafile ;' Inflection 3'});
    %         savetopptx(figelbs,file,figtype,{info.datafile ;' CI of onsets and slopes'});
    %         savetopptx(figciadiff,file,figtype,{info.datafile ;' Regression CI vs. slopes diff'});
    %
    %     end
    
    
    
    %%
    %     %%%%%%%%%%%%%
%     %onsets selection by threshold on both side of ci
%     %peak
%     peak_aux=(sum(abs(peak_ci(:,2:3)-peak_ci(:,1))<peak_t,2)==2);
%     peak_sel=peak_ci(:,1:3);
%     peak_sel(find(peak_aux==0),:)=nan;
%     %     %elb1
%     %     elb1_aux=(sum(abs(elb1_ci(:,2:3)-elb1_ci(:,1))<elb1_t,2)==2);
%     %     elb1_sel=elb1_ci(:,1:3);
%     %     elb1_sel(find(elb1_aux==0),:)=nan;
%     %elb1 normalized
%     elb1_aux=(sum(abs(elb1_ci(:,4:5)-elb1_ci(:,1))<elb1n_t,2)==2);
%     elb1_sel=elb1_ci(:,[1 4 5]);
%     elb1_sel(find(elb1_aux==0),:)=nan;
%     %elb2
%     elb2_aux=(sum(abs(elb2_ci(:,2:3)-elb2_ci(:,1))<elb2_t,2)==2);
%     elb2_sel=elb2_ci(:,1:3);
%     elb2_sel(find(elb2_aux==0),:)=nan;
%     %     %inflect1
%     %     inflect1_aux=(sum(abs(inflect1_ci(:,2:3)-inflect1_ci(:,1))<inflect1_t,2)==2);
%     %     inflect1_sel=inflect1_ci(:,1:3);
%     %     inflect1_sel(find(inflect1_aux==0),:)=nan;
%     %inflect1 normalized
%     inflect1_aux=(sum(abs(inflect1_ci(:,4:5)-inflect1_ci(:,1))<inflect1n_t,2)==2);
%     inflect1_sel=inflect1_ci(:,[1 4 5]);
%     inflect1_sel(find(inflect1_aux==0),:)=nan;
%     %inflect2
%     inflect2_aux=(sum(abs(inflect2_ci(:,2:3)-inflect2_ci(:,1))<inflect2_t,2)==2);
%     inflect2_sel=inflect2_ci(:,1:3);
%     inflect2_sel(find(inflect2_aux==0),:)=nan;
%     %     %inflect3
%     %     inflect3_aux=(sum(abs(inflect3_ci(:,2:3)-inflect3_ci(:,1))<inflect3_t,2)==2);
%     %     inflect3_sel=inflect3_ci(:,1:3);
%     %     inflect3_sel(find(inflect3_aux==0),:)=nan;
%     %inflect3 normalized
%     inflect3_aux=(sum(abs(inflect3_ci(:,4:5)-inflect3_ci(:,1))<inflect3n_t,2)==2);
%     inflect3_sel=inflect3_ci(:,[1 4 5]);
%     inflect3_sel(find(inflect3_aux==0),:)=nan;
    
    
    %onsets selection by threshold on range of ci
    %peak
    peak_aux=abs(peak_ci(:,2)-peak_ci(:,3))<peak_t;
    peak_sel=peak_ci(:,1:3);
    peak_sel(find(peak_aux==0),:)=nan;
    
    %elb1
    elb1_thresh=abs(elb1_ci(:,4)-elb1_ci(:,5))<elb1n_t;
    elb1_sel=elb1_ci(:,[1 4 5]);
    elb1_sel(find(elb1_thresh==0),:)=nan;
    
    %elb2
    elb2_sel=nan(16,3);
    
    %inflect1
    %inflect1_thresh=abs(inflect1_ci(:,4)-inflect1_ci(:,5))<inflect1n_t;
    %inflect1_thresh=abs(inflect1_ci(:,4)-inflect1_ci(:,5))<inflect1n_t & fitslopes1_diff_ci(:,1)>inflect1adiff_t;
    inflect1_thresh=tsignif_fitslopes1;
    inflect1_sel=inflect1_ci(:,[1 4 5]);
    inflect1_sel(find(inflect1_thresh==0),:)=nan;
    
    %elb2
    inflect2_sel=nan(16,3);
    
    %inflect3
    %inflect3_thresh=abs(inflect3_ci(:,4)-inflect3_ci(:,5))<inflect3n_t; 
    %inflect3_thresh=abs(inflect3_ci(:,4)-inflect3_ci(:,5))<inflect3n_t & fitslopes3_diff_ci(:,1)>inflect3adiff_t;
    inflect3_thresh=tsignif_fitslopes3;
    inflect3_sel=inflect3_ci(:,[1 4 5]);
    inflect3_sel(find(inflect3_thresh==0),:)=nan;
    
    
    %     %plot selected onsets
    %     linetype_ci={'-' '--' '--'};
    %     linewidth_ci=[4 2 2];
    % %     color_ci={'k' 'k' 'k'};
    % %     %peak mean and ci
    % %     figure(figtrialsboot)
    % %     for c=1:3
    % %         onset_plot=[];onset_plot(:,1)=peak_sel(:,c);
    % %         for ch=1:info.nchannels
    % %             if ~isnan(onset_plot(ch,1))
    % %                 onset_plot(ch,2)=trials_boot_avg(ch,round(onset_plot(ch,1))+info.aligntime);%peak_avg(:,2);
    % %             end
    % %         end
    % %         plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n',linetype_ci{c},linewidth_ci(c),color_ci{c});
    % %     end
    %     %elb1 mean and ci
    %     figure(figtrialsboot)
    %     color_ci={'k' 'k' 'k'};
    %     for c=1:3
    %         onset_plot=[];onset_plot(:,1)=elb1_sel(:,c);
    %         for ch=1:info.nchannels
    %             if ~isnan(onset_plot(ch,1))
    %                 onset_plot(ch,2)=trials_boot_avg(ch,round(onset_plot(ch,1))+info.aligntime);%elb1_avg(:,2);
    %             end
    %         end
    %         plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n',linetype_ci{c},linewidth_ci(c),color_ci{c});
    %     end
    % %     %elb2 mean and ci
    % %     figure(figtrialsboot)
    % %     for c=1:3
    % %         onset_plot=[];onset_plot(:,1)=elb2_sel(:,c);
    % %         for ch=1:info.nchannels
    % %             if ~isnan(onset_plot(ch,1))
    % %                 onset_plot(ch,2)=trials_boot_avg(ch,round(onset_plot(ch,1))+info.aligntime);%elb2_avg(:,2);
    % %             end
    % %         end
    % %         plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n',linetype_ci{c},linewidth_ci(c),color_ci{c});
    % %     end
    %     %inflect1 mean and ci
    %     figure(figtrialsboot)
    %     for c=1:3
    %         onset_plot=[];onset_plot(:,1)=inflect1_sel(:,c);
    %         for ch=1:info.nchannels
    %             if ~isnan(onset_plot(ch,1))
    %                 onset_plot(ch,2)=trials_boot_avg(ch,round(onset_plot(ch,1))+info.aligntime);%inflect1_avg(:,2);
    %             end
    %         end
    %         plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n',linetype_ci{c},linewidth_ci(c),color_ci{c});
    %     end
    %     %inflect1 mean and ci
    %     figure(figtrialsboot2)
    %     for c=1:3
    %         onset_plot=[];onset_plot(:,1)=inflect1_sel(:,c);
    %         for ch=1:info.nchannels
    %             if ~isnan(onset_plot(ch,1))
    %                 onset_plot(ch,2)=trials_boot_avg(ch,round(onset_plot(ch,1))+info.aligntime);%inflect1_avg(:,2);
    %             end
    %         end
    %         plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n',linetype_ci{c},linewidth_ci(c),color_ci{c});
    %     end
    %     %inflect2 mean and ci
    %     figure(figtrialsboot2)
    %     for c=1:3
    %         onset_plot=[];onset_plot(:,1)=inflect2_sel(:,c);
    %         for ch=1:info.nchannels
    %             if ~isnan(onset_plot(ch,1))
    %                 onset_plot(ch,2)=trials_boot_avg(ch,round(onset_plot(ch,1))+info.aligntime);%inflect2_avg(:,2);
    %             end
    %         end
    %         plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n',linetype_ci{c},linewidth_ci(c),color_ci{c});
    %     end
    %     %inflect3 mean and ci
    %     figure(figtrialsboot3)
    %     for c=1:3
    %         onset_plot=[];onset_plot(:,1)=inflect3_sel(:,c);
    %         for ch=1:info.nchannels
    %             if ~isnan(onset_plot(ch,1))
    %                 onset_plot(ch,2)=trials_boot_avg(ch,round(onset_plot(ch,1))+info.aligntime);%inflect3_avg(:,2);
    %             end
    %         end
    %         plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n',linetype_ci{c},linewidth_ci(c),color_ci{c});
    %     end
    %
    
    %%
    %%%%%%%%%%%%%
    %classification accum/burst
    
    %selection between elb1/inflect1/inflect3
    %     %threshold each sides separately
    %     %elb1_thresh=(sum(abs(elb1_ci(:,2:3)-elb1_ci(:,1))<elb1_t,2))>=1;
    %     %inflect1_thresh=(sum(abs(inflect1_ci(:,2:3)-inflect1_ci(:,1))<inflect1_t,2))>=1;
    %     %inflect3_thresh=(sum(abs(inflect3_ci(:,2:3)-inflect3_ci(:,1))<inflect3_t,2))>=1;
    %     elb1_thresh=(sum(abs(elb1_ci(:,4:5)-elb1_ci(:,1))<elb1n_t,2))>=1;
    %     inflect1_thresh=(sum(abs(inflect1_ci(:,4:5)-inflect1_ci(:,1))<inflect1n_t,2))>=1;
    %     inflect3_thresh=(sum(abs(inflect3_ci(:,4:5)-inflect3_ci(:,1))<inflect3n_t,2))>=1;
    
    %     %threshold on range of CI
    %     elb1_thresh=abs(elb1_ci(:,4)-elb1_ci(:,5))<elb1n_t;
    %     inflect1_thresh=abs(inflect1_ci(:,4)-inflect1_ci(:,5))<inflect1n_t;
    %     inflect3_thresh=abs(inflect3_ci(:,4)-inflect3_ci(:,5))<inflect3n_t;
    
    classif_t=-50
    accum_onset=nan(info.nchannels,2);burst_onset=nan(info.nchannels,2);
    inflect31_onset=nan(info.nchannels,2);inflect32_onset=nan(info.nchannels,2);
    scat_01=nan(info.nchannels,2);
    scat_10=nan(info.nchannels,2);
    scat_10_3=nan(info.nchannels,2);
    scat_11=nan(info.nchannels,2);
    scat_11_3=nan(info.nchannels,2);
 
    for ch=1:info.nchannels
        if elb1_thresh(ch)~=1 & inflect1_thresh(ch)==1
            %classification of inflect1
            if inflect1_ci(ch,1)>classif_t
                burst_onset(ch,1)=inflect1_ci(ch,1);
                burst_onset(ch,2)=trials_boot_avg(ch,round(burst_onset(ch,1))+info.aligntime);
            else
                accum_onset(ch,1)=inflect1_ci(ch,1);
                accum_onset(ch,2)=trials_boot_avg(ch,round(accum_onset(ch,1))+info.aligntime);
            end
            
            %for scatter plot
            scat_01(ch,1:2)=[0 inflect1_ci(ch,1)];
                     
            display(['ch=' num2str(ch) ' elb1~=1 inflect1=1'])
            display('This condition happens sometimes!')
            %pause

            
            
        elseif elb1_thresh(ch)==1 & inflect1_thresh(ch)~=1 
            %NOTE: optimally replace classif_t boundary by progressive increase
            %of the significance of the slopes difference for inflect3
            %classification of elb1
            if elb1_ci(ch,1)<classif_t
                accum_onset(ch,1)=elb1_ci(ch,1);
                accum_onset(ch,2)=trials_boot_avg(ch,round(accum_onset(ch,1))+info.aligntime);
 
                %for scatter plot
                scat_10(ch,1:2)=[elb1_ci(ch,1) 0];
                
                display(['ch=' num2str(ch) ' elb1=1 inflect1~=1  elb1<classif_t'])
                display('elb1 (<classif_t):')
                elb1_ci(ch,1)
                elb1_ci(ch,1)<classif_t               
                %pause
 
            elseif elb1_ci(ch,1)>classif_t & inflect3_thresh(ch)==1
                if inflect3_ci(ch,1)>classif_t
                    burst_onset(ch,1)=inflect3_ci(ch,1);
                    burst_onset(ch,2)=trials_boot_avg(ch,round(burst_onset(ch,1))+info.aligntime);
                else
                    accum_onset(ch,1)=inflect3_ci(ch,1);
                    accum_onset(ch,2)=trials_boot_avg(ch,round(accum_onset(ch,1))+info.aligntime);
                end
            
                %for scatter plot
                scat_10(ch,1:2)=[elb1_ci(ch,1) 0];
                scat_10_3(ch,1:2)=[inflect3_ci(ch,1) 0];
                
                %for displaying inflect3
                inflect31_onset(ch,1)=inflect3_ci(ch,1);
                inflect31_onset(ch,2)=trials_boot_avg(ch,round(inflect31_onset(ch,1))+info.aligntime);
                
                display(['ch=' num2str(ch) ' elb1=1 inflect1~=1 elb1>classif_t inflect3=1'])
                display('inflect3 (elb1>classif_t):')
                inflect3_ci(ch,1)
                inflect3_ci(ch,1)>classif_t                
                %pause
            
            else
                display(['ch=' num2str(ch) ' elb1=1 inflect1~=1 elb1>classif_t inflect3~=1'])
                display('discarded estimations')
                display('elb1 (>classif_t):')
                elb1_ci(ch,1)
                elb1_ci(ch,1)>classif_t               
                %pause
            end 
            
            
        elseif elb1_thresh(ch)==1 & inflect1_thresh(ch)==1
            if tsignif_elb1_inflect1(ch)==1
                accum_onset(ch,1)=elb1_ci(ch,1);
                accum_onset(ch,2)=trials_boot_avg(ch,round(accum_onset(ch,1))+info.aligntime);
                burst_onset(ch,1)=inflect1_ci(ch,1);
                burst_onset(ch,2)=trials_boot_avg(ch,round(burst_onset(ch,1))+info.aligntime);
                display(['ch=' num2str(ch) ' elb1=1 inflect1=1'])
                display('best scenario')
                %pause
                
                %for scatter plot
                scat_11(ch,1:2)=[elb1_ci(ch,1) inflect1_ci(ch,1)];
                
                
            elseif tsignif_elb1_inflect1(ch)~=1 & inflect3_thresh(ch)==1
                %classification of inflect3
                if inflect3_ci(ch,1)>classif_t
                    burst_onset(ch,1)=inflect3_ci(ch,1);
                    burst_onset(ch,2)=trials_boot_avg(ch,round(burst_onset(ch,1))+info.aligntime);
                else
                    accum_onset(ch,1)=inflect3_ci(ch,1);
                    accum_onset(ch,2)=trials_boot_avg(ch,round(accum_onset(ch,1))+info.aligntime);
                end
                
                %for displaying inflect3
                inflect32_onset(ch,1)=inflect3_ci(ch,1);
                inflect32_onset(ch,2)=trials_boot_avg(ch,round(inflect32_onset(ch,1))+info.aligntime);
                
                %for scatter plot
                scat_11_3(ch,1:2)=[inflect3_ci(ch,1)  inflect3_ci(ch,1)];
                      
                display(['ch=' num2str(ch) ' elb1=1 inflect1=1 tsignif_elb1_inflect1~=1 inflect3=1'])
                display('This condition happens sometimes!')
                display('inflect3 (>classif_t):')
                inflect3_ci(ch,1)
                inflect3_ci(ch,1)>classif_t
                %pause
            
            elseif tsignif_elb1_inflect1(ch)~=1 & inflect3_thresh(ch)~=1
                display(['ch=' num2str(ch) ' elb1=1 inflect1=1 tsignif_elb1_inflect1~=1 inflect3~=1'])
                display('This condition never happens!')
                pause
                
            else
                display(['ch=' num2str(ch) ' elb1=1 inflect1=1 else condition'])
                display('This condition never happens!')
                elb1_thresh(ch)
                inflect1_thresh(ch)
                tsignif_elb1_inflect1(ch)
                inflect3_thresh(ch)
                pause
            
            end
            
        else
            display(['ch=' num2str(ch) ' else condition'])
            elb1_thresh(ch)
            inflect1_thresh(ch)
            %pause
        end
    end
    
    %plot inflect3
    figure(figtrialsboot3)
    plot_events_ch(inflect31_onset.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','--',5,'g');
    plot_events_ch(inflect32_onset.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','--',5,'k');
    
    
    %plot peak
    figure(figtrialsboot2)
    for c=1
        onset_plot=[];onset_plot(:,1)=peak_sel(:,c);
        for ch=1:info.nchannels
            if ~isnan(onset_plot(ch,1)) & round(onset_plot(ch,1))+info.aligntime <= size(trials_boot_avg,2)
                onset_plot(ch,2)=trials_boot_avg(ch,round(onset_plot(ch,1))+info.aligntime);%peak_avg(:,2);
            end
        end
        plot_events_ch(onset_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','-',5,'k');
    end
    
    %plot accum
    figure(figtrialsboot2)
    plot_events_ch(accum_onset.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','-',5,'c');
    %plot burst
    figure(figtrialsboot2)
    plot_events_ch(burst_onset.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','--',5,'m');
    
    
    
    %%
    %%%%%%%%%%%%%
    %save in powerpoint
    if savepptx,
        savetopptx(figtrialsboot,file,figtype,{info.datafile ;' elb1 and inflect1'});
        savetopptx(figtrialsboot3,file,figtype,{info.datafile ;' inflect3'});
        savetopptx(figelbs,file,figtype,{info.datafile ;' CI of onsets and slopes'});
        savetopptx(figtrialsboot2,file,figtype,{info.datafile ;' accum and burst (classification)'});
        %savetopptx(figciadiff,file,figtype,{info.datafile ;' Regression CI vs. slopes diff'});
    end
    
    
    
    %%
    %%%%%%%%%%%%%%%%
    %alignment of onset using CSD features (after compute_CSDfeature)
    info.csdfeat_avg_targ=data(1).offline.csdfeat_avg_targ;
    info.zs=data(1).offline.csdzs;
    dref=info.csdfeat_avg_targ(2);
    
    %peak
    onset_aux=[];
    onset_aux(:,[1 3 4])=peak_sel;
    for ch=1:info.nchannels
        if ~isnan(onset_aux(ch,1))
            onset_aux(ch,2)=trials_boot_avg(ch,round(onset_aux(ch,1))+info.aligntime);
        else
            onset_aux(ch,2)=nan;
        end
    end
    [peak_r info_r ch_ref ~]=get_data_aligndepth(onset_aux,dref,info);
    %elb1
    onset_aux=[];
    onset_aux(:,[1 3 4])=elb1_sel;
    for ch=1:info.nchannels
        if ~isnan(onset_aux(ch,1))
            onset_aux(ch,2)=trials_boot_avg(ch,round(onset_aux(ch,1))+info.aligntime);
        else
            onset_aux(ch,2)=nan;
        end
    end
    [elb1_r info_r ch_ref ~]=get_data_aligndepth(onset_aux,dref,info);
    %elb2
    onset_aux=[];
    onset_aux(:,[1 3 4])=elb2_sel;
    for ch=1:info.nchannels
        if ~isnan(onset_aux(ch,1))
            onset_aux(ch,2)=trials_boot_avg(ch,round(onset_aux(ch,1))+info.aligntime);
        else
            onset_aux(ch,2)=nan;
        end
    end
    [elb2_r info_r ch_ref ~]=get_data_aligndepth(onset_aux,dref,info);
    %inflect1
    onset_aux=[];
    onset_aux(:,[1 3 4])=inflect1_sel(:,1:3);
    for ch=1:info.nchannels
        if ~isnan(onset_aux(ch,1))
            onset_aux(ch,2)=trials_boot_avg(ch,round(onset_aux(ch,1))+info.aligntime);
        else
            onset_aux(ch,2)=nan;
        end
    end
    [inflect1_r info_r ch_ref ~]=get_data_aligndepth(onset_aux,dref,info);
    %inflect2
    onset_aux=[];
    onset_aux(:,[1 3 4])=inflect2_sel;
    for ch=1:info.nchannels
        if ~isnan(onset_aux(ch,1))
            onset_aux(ch,2)=trials_boot_avg(ch,round(onset_aux(ch,1))+info.aligntime);
        else
            onset_aux(ch,2)=nan;
        end
    end
    [inflect2_r info_r ch_ref ~]=get_data_aligndepth(onset_aux,dref,info);
    %inflect3
    onset_aux=[];
    onset_aux(:,[1 3 4])=inflect3_sel;
    for ch=1:info.nchannels
        if ~isnan(onset_aux(ch,1))
            onset_aux(ch,2)=trials_boot_avg(ch,round(onset_aux(ch,1))+info.aligntime);
        else
            onset_aux(ch,2)=nan;
        end
    end
    [inflect3_r info_r ch_ref ~]=get_data_aligndepth(onset_aux,dref,info);
    
    %accum_onset
    onset_aux=[];
    onset_aux(:,[1 2])=accum_onset;
    [accum_r info_r ch_ref dref_conv]=get_data_aligndepth(onset_aux,dref,info);
    %burst_onset
    onset_aux=[];
    onset_aux(:,[1 2])=burst_onset;
    [burst_r info_r ch_ref ~]=get_data_aligndepth(onset_aux,dref,info);
    
    %scatter plot
    [scat_01_r ~]=get_data_aligndepth(scat_01,dref,info);
    [scat_10_r ~]=get_data_aligndepth(scat_10,dref,info);
    [scat_10_3_r ~]=get_data_aligndepth(scat_10_3,dref,info);
    [scat_11_r ~]=get_data_aligndepth(scat_11,dref,info);
    [scat_11_3_r ~]=get_data_aligndepth(scat_11_3,dref,info);
    
    
    %remove non significant burst
    [burst_bsignif_r info_r ch_ref ~]=get_vmis_aligndepth(burst_bsignif,dref,info);
    burst_bsignif_mat=[burst_bsignif_r ; burst_bsignif_r ; burst_bsignif_r ; burst_bsignif_r]';
    peak_r=peak_r.*burst_bsignif_mat;
    elb1_r=elb1_r.*burst_bsignif_mat;
    elb2_r=elb2_r.*burst_bsignif_mat;
    inflect1_r=inflect1_r.*burst_bsignif_mat;
    inflect2_r=inflect2_r.*burst_bsignif_mat;
    inflect3_r=inflect3_r.*burst_bsignif_mat;
    accum_r=accum_r.*[burst_bsignif_r ; burst_bsignif_r]';
    burst_r=burst_r.*[burst_bsignif_r ; burst_bsignif_r]';
    scat_01_r=scat_01_r.*[burst_bsignif_r ; burst_bsignif_r]';
    scat_10_r=scat_10_r.*[burst_bsignif_r ; burst_bsignif_r]';
    scat_10_3_r=scat_10_3_r.*[burst_bsignif_r ; burst_bsignif_r]';
    scat_11_r=scat_11_r.*[burst_bsignif_r ; burst_bsignif_r]';
    scat_11_3_r=scat_11_3_r.*[burst_bsignif_r ; burst_bsignif_r]';
    
    %list of onsets for all sessions
    allpeak{dd}=peak_r;
    allelb1{dd}=elb1_r;
    allelb2{dd}=elb2_r;
    allinflect1{dd}=inflect1_r;
    allinflect2{dd}=inflect2_r;
    allinflect3{dd}=inflect3_r;
    allaccum{dd}=accum_r;
    allburst{dd}=burst_r;
    
    %histos
    classif.accum_only{dd} = ~isnan(accum_r(:,1)) & isnan(burst_r(:,1));
    classif.burst_only{dd} = isnan(accum_r(:,1)) & ~isnan(burst_r(:,1));
    classif.accum_burst{dd} = ~isnan(accum_r(:,1)) & ~isnan(burst_r(:,1));
    classif.nothing{dd} = isnan(accum_r(:,1)) & isnan(burst_r(:,1));
    classif.accum_only_accum_burst_only{dd} = (~isnan(accum_r(:,1)) & isnan(burst_r(:,1))) | (~isnan(accum_r(:,1)) & ~isnan(burst_r(:,1)));
    
    %for scatter plot
    allscat.scat_01{dd}=scat_01_r;
    allscat.scat_10{dd}=scat_10_r;
    allscat.scat_10_3{dd}=scat_10_3_r;
    allscat.scat_11{dd}=scat_11_r;
    allscat.scat_11_3{dd}=scat_11_3_r;
    
    
    
    %     %display in original channel coordinates for checking
    %     display(['accum_only: ' num2str((find(classif.accum_only{dd})-16-(16/2-dref_conv))')])
    %     display(['burst_only: ' num2str((find(classif.burst_only{dd})-16-(16/2-dref_conv))')])
    %     display(['accum_burst: ' num2str((find(classif.accum_burst{dd})-16-(16/2-dref_conv))')])
    %     display(['no classif: ' num2str((find(classif.nothing{dd})-16-(16/2-dref_conv))')])
    
    
    d
    dd
    pause
    close all
    
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%peak, elb1 and inflect1
%plot latencies
fig_l=figure('Position',[scrsz(3)/4 100 scrsz(3)/2 scrsz(4)-200]);hold on;
fig_lfit=figure('Position',[scrsz(3)/4 100 scrsz(3)/2 scrsz(4)-200]);hold on;
[fig_allpeak_l ~]=plot_stats_depths_v('latency',allpeak,[ ],1,'k',[-110 20],fig_l,info,datalist,dlist,fig_l,[])
[fig_allelb1_l ~]=plot_stats_depths_v('latency',allelb1,[ ],1,'b',[-110 20],fig_l,info,datalist,dlist,fig_l,[])
[fig_allinflect1_l ~]=plot_stats_depths_v('latency',allinflect1,[ ],1,'r',[-110 20],fig_l,info,datalist,dlist,fig_l,[])

%plot_stats_depths('latency',allelb2,[-50 0],1,'g',[-100 0],fig_l,info,datalist,dlist)
%plot_stats_depths('latency',allinflect2,[-150 0],1,'m',[-100 0],fig_l,info,datalist,dlist)
%plot_stats_depths('latency',allinflect3,[-50 0],1,'y',[-100 0],fig_l,info,datalist,dlist)
grid;
legend('Peak','','fit','onset elb1','','fit','onset inflect1','','fit','Location','Northwest')


%%
%plot magnitudes
fig_m=figure('Position',[scrsz(3)/4 100 scrsz(3)/2 scrsz(4)-200]);hold on;
fig_mfit=figure('Position',[scrsz(3)/4 100 scrsz(3)/2 scrsz(4)-200]);hold on;
[fig_allpeak_m max_peak_m]=plot_stats_depths_v('magnitude',allpeak,[ ],1,'k',[0 200],fig_m,info,datalist,dlist,fig_m,[])
[fig_allelb1_m ~]=plot_stats_depths_v('magnitude',allelb1,[ ],1,'b',[0 200],fig_m,info,datalist,dlist,fig_m,[])
[fig_allinflect1_m max_inflect1_m]=plot_stats_depths_v('magnitude',allinflect1,[ ],1,'r',[0 200],fig_m,info,datalist,dlist,fig_m,[])
%plot_stats_depths('magnitude',allelb2,[0 500],1,'g',[0 170],fig_m,info,datalist,dlist)
%plot_stats_depths('magnitude',allinflect2,[0 500],1,'m',[0 170],fig_m,info,datalist,dlist)
%plot_stats_depths('magnitude',allinflect3,[0 500],1,'y',[0 170],fig_m,info,datalist,dlist)
grid;
%coef=2.7;%VG
%coef=2;%MG
%coef=max_peak_m/max_inflect1_m;
%display(['coef multiplicator=' num2str(coef)])
%plot_stats_depths_v('magnitude',allinflect1,[ ],coef,'g',[0 180],fig_m,info,datalist,dlist,fig_m,[])
%legend('Peak','','fit','onset elb1','','fit','onset inflect1','','fit',[num2str(coef) 'X onset inflect1'],'','fit','Location','NorthEast')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%peak, accum and burst (classification)
%plot latencies
fig_class_l=figure('Position',[scrsz(3)/4 100 scrsz(3)/2 scrsz(4)-200]);hold on;
fig_class_lfit=figure('Position',[scrsz(3)/4 100 scrsz(3)/2 scrsz(4)-200]);hold on;
[fig_allpeak_l ~]=plot_stats_depths_v('latency',allpeak,[ ],1,'k',[-100 20],fig_class_l,info,datalist,dlist,fig_class_l,[]);
[fig_allaccum_l ~]=plot_stats_depths_v('latency',allaccum,[ ],1,'b',[-110 20],fig_class_l,info,datalist,dlist,fig_class_l,[]);%corr [30]
[fig_allburst_l ~]=plot_stats_depths_v('latency',allburst,[ ],1,'r',[-110 20],fig_class_l,info,datalist,dlist,fig_class_l,[]);
figure(fig_class_l)
grid;
legend('Peak','','fit','onset accum','','fit','onset burst','','fit','Location','NorthWest')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot magnitudes
fig_class_m=figure('Position',[scrsz(3)/4 100 scrsz(3)/2 scrsz(4)-200]);hold on;
fig_class_mfit=figure('Position',[scrsz(3)/4 100 scrsz(3)/2 scrsz(4)-200]);hold on;
[fig_allpeak_m max_peak_m]=plot_stats_depths_v('magnitude',allpeak,[ ],1,'k',[0 200],fig_class_m,info,datalist,dlist,fig_class_m,[])
[fig_allaccum_m ~]=plot_stats_depths_v('magnitude',allaccum,[ ],1,'b',[0 200],fig_class_m,info,datalist,dlist,fig_class_m,[])
[fig_allburst_m max_burst_m]=plot_stats_depths_v('magnitude',allburst,[ ],1,'r',[0 200],fig_class_m,info,datalist,dlist,fig_class_m,[])
%coef=3.2;%VG
%coef=3.5;%MG
coef=max_peak_m/max_burst_m;
display(['coef multiplicator=' num2str(coef)])
plot_stats_depths_v('magnitude',allburst,[ ],coef,'g',[0 200],fig_class_m,info,datalist,dlist,fig_class_m,[])
grid;
legend('Peak','','fit','onset accum','','fit','onset burst','','fit',[num2str(coef) 'X onset burst'],'','fit','Location','NorthEast')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot histos classification
fig_class_h=figure('Position',[scrsz(3)/4 100 scrsz(3)/2 scrsz(4)-200]);hold on;
max_val=length(dlist);
%plot total
ltotal=max([size(classif.accum_only,2) size(classif.burst_only,2) size(classif.accum_burst,2)]);
classif.total={};
for dtot=[1:ltotal],
    classif.total{dtot}=classif.accum_only{dtot}+classif.burst_only{dtot}+classif.accum_burst{dtot};
end
total(1,:)=plot_histos_depths_v(classif.total,classif.total,[ ],'k',[0 max_val],fig_class_h);

total(2,:)=plot_histos_depths_v(classif.accum_only,classif.total,[ ],'g',[0 max_val],fig_class_h);
total(3,:)=plot_histos_depths_v(classif.burst_only,classif.total,[ ],'m',[0 max_val],fig_class_h);
total(4,:)=plot_histos_depths_v(classif.accum_burst,classif.total,[ ],'b',[0 max_val],fig_class_h);
%plot_histos_depths_v(classif.nothing,classif.total,[ ],'k',[0 max_val],fig_class_h)

total(5,:)=plot_histos_depths_v(classif.accum_only_accum_burst_only,classif.total,[ ],'r',[0 max_val],fig_class_h);

total

subplot(1,3,1);
text(10,8,num2str(round(total,2)))
legend('Significant bursts','Build-up only','Burst only','Build-up+Burst','Build-up only and Build-up+Burst','Location','NorthEast')

%proportion of each type of activity


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%scatter plot of Elb1/inflect1/inflect3
figscat=figure;hold on;
for dd=1:size(allscat.scat_01,2)
    aux=allscat.scat_11{dd};
    for ch=1:size(aux,1)
        plot(aux(ch,1),aux(ch,2),'mo')
    end
    
    aux=allscat.scat_11_3{dd};
    for ch=1:size(aux,1)
        plot(aux(ch,1),aux(ch,2),'ko')
    end
 
    aux=allscat.scat_10{dd};
    for ch=1:size(aux,1)
        plot(aux(ch,1),aux(ch,2),'bo')
    end
    
    aux=allscat.scat_10_3{dd};
    for ch=1:size(aux,1)
        plot(aux(ch,1),aux(ch,2),'go')
    end
    
    aux=allscat.scat_01{dd};
    for ch=1:size(aux,1)
        plot(aux(ch,1),aux(ch,2),'ro')
    end
    
    
end
%%
minplot=-150;
hdl=line([0 minplot],[0 minplot]);
set(hdl,'color','k','linestyle','--')
hdl=line([minplot 0],[classif_t classif_t]);
set(hdl,'color','r','linestyle','--')
hdl=line([classif_t classif_t],[minplot 0 ]);
set(hdl,'color','r','linestyle','--')
xlabel('Elb1');ylabel('Inflect1')
title('all events onsets')
axis([minplot 0 minplot 0])


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%histos plot of Elb1/inflect1/inflect3
fighistos=figure;hold on;
hist10=[];hist10_3=[];hist01=[];hist11_3=[];hist11_1=[];hist11_2=[];
minch=14;maxch=32;
for dd=1:size(allscat.scat_01,2)
    aux=allscat.scat_11{dd};
    hist11_1=[hist11_1 aux(minch:maxch,1)'];
    hist11_2=[hist11_2 aux(minch:maxch,2)'];
    
    aux=allscat.scat_11_3{dd};
    hist11_3=[hist11_3 aux(minch:maxch,1)'];
    
    aux=allscat.scat_10{dd};
    hist10=[hist10 aux(minch:maxch,1)'];
    
    aux=allscat.scat_10_3{dd};
    hist10_3=[hist10_3 aux(minch:maxch,1)'];
    
    aux=allscat.scat_01{dd};
    hist01=[hist01 aux(minch:maxch,2)'];
end    
histo=[];
edges=[minplot:2:0];
maxhist=11;


subplot(1,6,1);hold on;
histo(:,1)=hist11_1';
hist=histc(histo,edges);
hdl=bar(edges,hist,'histc');
set(hdl,'Facecolor','b','EdgeColor','b')
hdl=line([classif_t classif_t],[0 maxhist+1]);
set(hdl,'color','k','linestyle','--')
axis([minplot 0 0 max(hist(:,1))+1])
xlabel('Elb1');ylabel('Number of MUA')
axis([-150 0 0 maxhist+1])

subplot(1,6,2);hold on;
histo(:,1)=hist11_2';
hist=histc(histo,edges);
hdl=bar(edges,hist,'histc');
set(hdl,'Facecolor','r','EdgeColor','r')
hdl=line([classif_t classif_t],[0 maxhist+1]);
set(hdl,'color','k','linestyle','--')
axis([minplot 0 0 max(hist(:,1))+1])
xlabel('Inflect1');ylabel('Number of MUA')
axis([-150 0 0 maxhist+1])

subplot(1,6,3);
histo(:,1)=hist11_3';
hist=histc(histo,edges);
hdl=bar(edges,hist,'histc');
set(hdl,'Facecolor','b','EdgeColor','b')
hdl=line([classif_t classif_t],[0 maxhist+1]);
set(hdl,'color','k','linestyle','--')
axis([minplot 0 0 max(hist(:,1))+1])
xlabel('Elb1 Inflect1 Inflect3');ylabel('Number of MUA')
axis([-150 0 0 maxhist+1])

subplot(1,6,4);
histo(:,1)=hist10';
hist=histc(histo,edges);
hdl=bar(edges,hist,'histc');
set(hdl,'Facecolor','b','EdgeColor','b')
hdl=line([classif_t classif_t],[0 maxhist+1]);
set(hdl,'color','k','linestyle','--')
axis([minplot 0 0 maxhist+1])
xlabel('Elb1 not Inflect1');ylabel('Number of MUA')
axis([-150 0 0 maxhist+1])

subplot(1,6,5);
histo(:,1)=hist10_3';
hist=histc(histo,edges);
hdl=bar(edges,hist,'histc');
set(hdl,'Facecolor','g','EdgeColor','g')
hdl=line([classif_t classif_t],[0 maxhist+1]);
set(hdl,'color','k','linestyle','--')
axis([minplot 0 0 max(hist(:,1))+1])
xlabel('Elb1 not Inflect1 Inflect3');ylabel('Number of MUA')
axis([-150 0 0 maxhist+1])

subplot(1,6,6);
histo(:,1)=hist01';
hist=histc(histo,edges);
hdl=bar(edges,hist,'histc');
set(hdl,'Facecolor','r','EdgeColor','r')
hdl=line([classif_t classif_t],[0 maxhist+1]);
set(hdl,'color','k','linestyle','--')
axis([minplot 0 0 max(hist(:,1))+1])
xlabel('not Elb1 Inflect1');ylabel('Number of MUA')
axis([-150 0 0 maxhist+1])


   


%%
%%%%%%%%%%%%%
%save in powerpoint
if savepptx,
    savetopptx(fig_l,file,figtype,{info.datafile ;' Average latencies of elb1/inflect1'});
    savetopptx(fig_m,file,figtype,{info.datafile ;' Average firing rates of elb1/inflect1'});
    savetopptx(fig_class_l,file,figtype,{info.datafile ;' Average latencies of accum/burst'});
    savetopptx(fig_class_m,file,figtype,{info.datafile ;' Average firing rates of accum/burst'});
    %     savetopptx(fig_allpeak_l,file,figtype,{info.datafile ;' All sessions peak latencies'});
    %     savetopptx(fig_allelb1_l,file,figtype,{info.datafile ;' All sessions elb1 latencies'});
    %     savetopptx(fig_allinflect1_l,file,figtype,{info.datafile ;' All sessions inflect1 latencies '});
    %     savetopptx(fig_allpeak_m,file,figtype,{info.datafile ;' All sessions peak firing rates'});
    %     savetopptx(fig_allelb1_m,file,figtype,{info.datafile ;' All sessions elb1 firing rates'});
    %     savetopptx(fig_allinflect1_m,file,figtype,{info.datafile ;' All sessions inflect1 firing rates'});
    
end

%%
if savepptx
    %close .pptx
    newFile = exportToPPTX('saveandclose',filepptx)
end









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MISC

%%
%old classification rules
% %             %classification of inflect3
% %             %(use histo of inflect3 compare with histo of elb1 and inflect1 wihtout inflect3 cases)
% %             if inflect3_ci(ch,1)>-50
% %                 burst_onset(ch,1)=inflect3_ci(ch,1);
% %                 burst_onset(ch,2)=trials_boot_avg(ch,round(burst_onset(ch,1))+info.aligntime);
% %             else
% %                 accum_onset(ch,1)=inflect3_ci(ch,1);
% %                 accum_onset(ch,2)=trials_boot_avg(ch,round(accum_onset(ch,1))+info.aligntime);
% %             end
% %             inflect31_onset(ch,1)=inflect3_ci(ch,1);
% %             inflect31_onset(ch,2)=trials_boot_avg(ch,round(inflect31_onset(ch,1))+info.aligntime);


%         elseif elb1_thresh(ch)==1 & inflect1_thresh(ch)~=1 & inflect3_thresh(ch)==1
%             if tsignif_elb1_inflect3(ch)==1
% %                 accum_onset(ch,1)=elb1_ci(ch,1);
% %                 accum_onset(ch,2)=trials_boot_avg(ch,round(accum_onset(ch,1))+info.aligntime);
% %                 %ch
% %                 %pause
% %
%             elseif tsignif_elb1_inflect3(ch)~=1
% %                     %classification of inflect3
% %                     %(use histo of inflect3 compare with histo of elb1 and inflect1 wihtout inflect3 cases)
% %                     if inflect3_ci(ch,1)>-50
% %                         burst_onset(ch,1)=inflect3_ci(ch,1);
% %                         burst_onset(ch,2)=trials_boot_avg(ch,round(burst_onset(ch,1))+info.aligntime);
% %                     else
% %                         accum_onset(ch,1)=inflect3_ci(ch,1);
% %                         accum_onset(ch,2)=trials_boot_avg(ch,round(accum_onset(ch,1))+info.aligntime);
% %                     end
% %                     inflect33_onset(ch,1)=inflect3_ci(ch,1);
% %                     inflect33_onset(ch,2)=trials_boot_avg(ch,round(inflect33_onset(ch,1))+info.aligntime);
%             end


%%
%%%%%%%%%%%%%
%other critaria than threshold
%     %plot resnorm 1 and 2
%
%     %                         %plot resnorm
%     %                         if ~isempty(resnorm_3{ch}) %trial_sel was long enough
%     %                             resnorm_plot=nan(1,info.aligntime+length(resnorm_3{ch})-1);
%     %                             resnorm_plot(1,info.aligntime:end)=resnorm_3{ch};
%     %                             resnorm_plot=resnorm_plot/max(resnorm_plot);
%     %                             figure(figresnorm);subplot(1,3,3);hold on;
%     %                             plot(resnorm_plot,'color',colorlist(ch,:),'linewidth',1);
%     %                             title('Residuals of 2pwlr for inflect 2');xlabel('Time (ms)');ylabel('R2');
%     %                         end
%
%
%
%     %%
%     %plot amplitude differences
%     figamps=figure('Position',[1 100 scrsz(3) scrsz(4)-500]);
%     hold on;
%
%     %plot amplitude diff baseline elb_1
%     elb1_boot_avg=squeeze(mean(elb1_boot(:,:,2),1))
%     elb1_boot_var=squeeze(std(elb1_boot(:,:,2),1));%/size(elb1_boot,1));
%     subplot(1,3,1);hold on;
%     errorbar(1:info.nchannels,elb1_boot_avg,elb1_boot_var);
%     axis([1 info.nchannels 0 max(elb1_boot_avg)+max(elb1_boot_var)]);
%     xlabel('Channel');ylabel('FR (spk/s)')
%     title('Elbow1-Baseline')
%
%
%     %plot amplitude diff elb_1 elb_2
%     diffamp=elb2_boot(:,:,2)-elb1_boot(:,:,2);
%     diffamp_avg=squeeze(mean(diffamp,1))
%     diffamp_var=squeeze(std(diffamp,1));%/size(diffamp,1));
%     subplot(1,3,2);hold on;
%     errorbar(1:info.nchannels,diffamp_avg,diffamp_var);
%     axis([1 info.nchannels 0 (max(diffamp_avg)+max(diffamp_var))]);
%     xlabel('Channel');ylabel('FR (spk/s)')
%     title('Elbow2-Elbow1')
%
%     %plot amplitude diff elb_2 peak
%     diffamp=peak_boot(:,:,2)-elb1_boot(:,:,2);
%     diffamp_avg=squeeze(mean(diffamp,1))
%     diffamp_var=squeeze(std(diffamp,1));%/size(diffamp,1));
%     subplot(1,3,3);hold on;
%     errorbar(1:info.nchannels,diffamp_avg,diffamp_var);
%     axis([1 info.nchannels 0 (max(diffamp_avg)+max(diffamp_var))]);
%     xlabel('Channel');ylabel('FR (spk/s)')
%     title('Peak-Elbow2')
%




% %linear fit of each dataset
% %fit linear function
% npoly=1;%polynomial number for polyfit
%
% dispfit=1;
% if dispfit, figfit=figure;hold on;end;
% fun = @(x,xdata)x(1)+x(2)*xdata;
% options=optimset('Display','on');
% al=1;%useless
% lims=[16 33];
% dspk=0;dlfp=0;
% fitparams_spk=[];fitparams_lfp=[];lat_spk_all=[];lat_lfp_all=[];lat_spk_index=[];lat_lfp_index=[];
% regressparams_spk=[];regressparams_lfp=[];
% for dd=1:size(allonset_spk,2),
%
%     aux_spk=allonset_spk{dd};
%     lat_spk=squeeze(aux_spk(al,lims(1):lims(2),1));
%     xp1=find(~isnan(lat_spk));
%     yp1=lat_spk(xp1);
%     [yp1_s yp1_i]=sort(yp1);
%     xp1_s=xp1(yp1_i);
%
%     yfit_spk=[];yfit_lfp=[];
%     if ~isempty(xp1_s)
%         dspk=dspk+1;
%         x0p1=[yp1_s(1) xp1_s(1)];
%         [fitparams1 resnorm1] = lsqcurvefit(fun,x0p1,yp1_s,xp1_s,[],[],options);
%         fitparams_spk(dspk,:)=[resnorm1 fitparams1];
%
%         %lists
%         lat_spk_all=[lat_spk_all , lat_spk];
%         lat_spk_index=[lat_spk_index , 1:length(lat_spk)];
%
%         %regression
%         [p rsq yfit_spk]=get_regressioncoefs(yp1_s,xp1_s,npoly)
%         regressparams_spk(dspk,1)=rsq;
%         regressparams_spk(dspk,2:2+npoly)=p;
%         regressparams_spk(dspk,2+npoly+1)=numel(xp1);
%
%         if dispfit
%             figure(figfit)
%             subplot(1,2,1);hold on;
%             plot(yp1_s,xp1_s,'o','color',colorlist(dd,:));
%             %plot(yp1_s,fun(fitparams1,yp1_s),'b-');
%             %plot(yp1_s,xp1_s,'ko',yp1_s,yfit_spk,'b-')
%             plot(yp1_s,yfit_spk,'b-')
%
%             %axis
%             axis([-20 20 1 18])
%             set(gca,'Ytick',[1:2:18],'Yticklabel',[-8:2:8])
%             xlabel('Latency (ms)');ylabel('Channel')
%         end
%
%     end
%
%
%     %     aux_lfp=allonset_lfp{dd};
%     %     lat_lfp=squeeze(aux_lfp(al,lims(1):lims(2),1));
%     %     xp2=find(~isnan(lat_lfp));
%     %     yp2=lat_lfp(xp2);
%     %     [yp2_s yp2_i]=sort(yp2);
%     %     xp2_s=xp2(yp2_i);
%     %     if ~isempty(xp2_s)
%     %         dlfp=dlfp+1;
%     %         x0p2=[yp2(1) xp2(1)];
%     %         [fitparams2 resnorm2] = lsqcurvefit(fun,x0p2,yp2_s,xp2_s,[],[],options);
%     %         fitparams_lfp(dlfp,:)=[resnorm2 fitparams2];
%     %
%     %         %lists
%     %         lat_lfp_all=[lat_lfp_all , lat_lfp];
%     %         lat_lfp_index=[lat_lfp_index , 1:length(lat_lfp)];
%     %
%     %         %regression
%     %         [p rsq yfit_lfp]=get_regressioncoefs(yp2_s,xp2_s,npoly)
%     %         regressparams_lfp(dlfp,1)=rsq;
%     %         regressparams_lfp(dlfp,2:2+npoly)=p;
%     %         regressparams_lfp(dlfp,2+npoly+1)=numel(xp2);
%     %
%     %         if dispfit
%     %             figure(figfit)
%     %             subplot(1,2,2);hold on;
%     %             plot(yp2_s,xp2_s,'o','color',colorlist(dd,:));
%     %             %plot(yp2_s,fun(fitparams2,yp2_s),'r-');
%     %             %plot(yp2_s,xp2_s,'ko',yp2_s,yfit_lfp,'r-')
%     %             %pause
%     %
%     %             %axis
%     %             axis([-20 20 1 18])
%     %             set(gca,'Ytick',[1:2:18],'Yticklabel',[-8:2:8])
%     %             xlabel('Latency (ms)');ylabel('Channel')
%     %
%     %         end
%     %
%     %     end
%
%
% end


% %regression
% ind=find(~isnan(lat_spk_all));
% x=lat_spk_index(ind);
% y=lat_spk_all(ind);
%
% [p_spk rsq_spk]=get_regressioncoefs(x,y)
%
% ind=find(~isnan(lat_lfp_all));
% x=lat_lfp_index(ind);
% y=lat_lfp_all(ind);
%
% [p_lfp rsq_lfp]=get_regressioncoefs(x,y)







%%inflection
%              pause
%             %inflection point from spk onset and peak activity
%             anchor_b=onset_spk(:,1)+info.aligntime-20;
%             anchor_e=peak+info.aligntime;
%             inflect_spk=nan(info.nchannels,3);
%             figresnorm=figure;hold on;
%             fitparamslist=[]
%             for ch=1:info.nchannels
%                 if ~isnan(anchor_b(ch)) & ~isnan(anchor_e(ch))
%                     [inflect,resnorm]=get_inflection_2pwlr(trials_spk_n_avg(ch,:),anchor_b(ch),anchor_e(ch));
%                     inflect_spk(ch,1)=inflect(1)-info.aligntime;%correction for timing
%                     inflect_spk(ch,2:3)=inflect(2:3);
%
%                     %plot resnorm
%                     if ~isempty(resnorm) %trial_sel was long enough
%                         resnormplot=nan(1,info.aligntime+length(resnorm)-1);
%                         resnormplot(1,info.aligntime:end)=resnorm;
%                         resnormplot=resnormplot/max(resnormplot);
%                         figure(figresnorm);
%                         subplot(2,1,1);hold on;
%                         plot(resnormplot,'color',colorlist(ch,:),'linewidth',2);
%                         %plot(inflect_spk(1),inflect_spk(3),'ko');
%                         title('Residuals of 2pwlr');xlabel('Time (ms)');ylabel('R2');
%                         subplot(2,1,2);hold on;
%
%                         resnormplot_s=sort(resnormplot);%sorting
%                         resnormplot_sn=resnormplot_s-resnormplot_s(1);%normalization
%                         plot(resnormplot_sn,'color',colorlist(ch,:),'linewidth',2);
%                         %plot(inflect_spk(1),inflect_spk(3),'ko');
%                         title('Sorted residuals of 2pwlr');xlabel('Index');ylabel('R2');
%                         %pause
%
%                         %%
%                         %quadratic function
%                         resnormplot_r2=resnormplot_sn(find(~isnan(resnormplot_sn)));
%                         %fun = @(x,xdata)x(1)+x(2).*xdata+x(3).*(xdata.*xdata);
%                         fun = @(x,xdata)x(1)+x(2).*(xdata.*xdata);
%                         options=optimset('Display','off');
%                         xp=[1:length(resnormplot_r2)];
%                         x0p=[xp(1) resnormplot_r2(1)];
%                         [fitparams resnorm_fit2] = lsqcurvefit(fun,x0p,xp,resnormplot_r2,[],[],options);
%                         plot(xp,fun(fitparams,xp),'--','color',colorlist(ch,:),'linewidth',2)
%                         fitparamslist(ch,:)=fitparams;
%                         %pause
%
%                     end
%                 end
%             end
%
%             %channel without inflection point
%             ind_ninflect=find(fitparamslist(:,2)<0.09*0.001)
%
%             %plot inflection
%             figure(figtrials)
%             inflect_spk_plot=[];
%             inflect_spk_plot(:,1)=inflect_spk(:,1);
%             inflect_spk_plot(:,2)=inflect_spk(:,2);
%             %inflect_spk_plot(ind_ninflect,:)=nan;
%             plot_events_ch(inflect_spk_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','-',2,'');
%             %plot_events_ch(onset_spk_plot,[],vshift_spk,range_spk,info,hdlfig,'n','-',1,'');
%
%
%

% %%%%%%%%%%
% onset estimation within a trial by comparing within window distribution
% winsize=10;search_wind=[info.aligntime-200 info.aligntime];alpha=0.001;
% onset_spk2=zeros(info.nchannels,2);
% for ch=1:info.nchannels,
%     trials_ch=trials_spk_n_avg(ch,:);
%     onset_spk_ch=get_latency_trialch(trials_ch,winsize,search_wind,'fr',info,alpha);
%     onset_spk2(ch,:)=onset_spk_ch;
% end
% onset_spk2(:,1)=onset_spk2(:,1)-info.aligntime;%correction for timing
%
% %plot onsets
% onset_spk2_plot=[];
% onset_spk2_plot(:,1)=onset_spk2(:,1);
% onset_spk2_plot(:,2)=0;
% plot_events_ch(onset_spk2_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','--',1,'');


%             %%
%             %%%%%%%%%%%%%%%%
%             %lfp
%             [info.nchannels info.ntrials info.triallen]=size(trials_lfp);
%
%             %normalization
%             trials_lfp_n=get_trials_normalized(trials_lfp,trials_lfp_bsl,'LFP',info);
%             %compute average trials in RF and aRF
%             [trials_lfp_n_avg trials_lfp_n_var]=get_trials_avg(trials_lfp_n);
%
%             info.aligntime=aligntime_lfp;
%             hdlfig2=subplot(1,3,2);hold on;
%             [range_lfp vshift_lfp]=plot_trials(trials_lfp_n_avg,[],[],[],[],[],info,hdlfig2,[],[],[]);
%
%             %get lfp latency
%             winmin=[-10 30];winsize=20;alpha=0.01;H1_count=10;
%             lat_lfp=zeros(info.nchannels,2);
%             for ch=1:info.nchannels,
%                 trials_ch=trials_lfp_n_avg(ch,:);
%                 [lat_ch latback minmax]=get_latency_trialch_backspline(trials_ch,winmin,winsize,'fr',info,alpha,H1_count);
%                 lat_lfp(ch,:)=lat_ch;
%             end
%
%             lat_lfp(:,1)=lat_lfp(:,1)-info.aligntime;%correction for timing
%
%             lat_lfp_plot=[];
%             lat_lfp_plot(:,1)=lat_lfp(:,1);
%             lat_lfp_plot(:,2)=0;%lat_spk(:,3);
%
%             %plot_events_ch(lat_lfp_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_lfp,range_lfp,info,hdlfig2,'n','-',1,'');
%             plot_events_ch(lat_lfp_plot,[],vshift_lfp,range_lfp,info,hdlfig2,'n','-',1,'');
%
% %             %plot lfp latency on spk plot
% %             lat_lfp_onspk=lat_lfp+shift_ripple;
% %             lat_lfp_onspk(:,2)=0;%lat_spk(:,3);
% %
% %             %hdlfig=subplot(1,2,1);hold on;
% %             %plot_events_ch(lat_lfp_onspk.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','-.',1,'');
% %             plot_events_ch(lat_lfp_onspk,[],vshift_spk,range_spk,info,hdlfig,'n','-.',1,'');
% %
%             %%get lfp snr
%             %var_lfp=var(trials_lfp_avgn,[],2);
%             %var_lfpbsl=var(trials_lfp_latbsl_avgn,[],2);


