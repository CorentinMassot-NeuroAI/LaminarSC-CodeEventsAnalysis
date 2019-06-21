%function analysis_peak_latencyfrompeak

%function analysis_peak_latencyfrompeak
%   Analysis of peaks and latencies based on first estimation of peak of trial-averaged activity recorded with a
%   laminar probe (LMA) 
%
% measures P along sliding windows, if necessary use a gaussian to put more
% weight on center of windows. find first inflection point of second
% derivative to measure latency
% from there use derivative of signal to find first peak.
%
% see also analysis_latency analysis_peak
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 11/07/2016 last modified 07/24/2017
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
alignlist={'targ'};
%alignlist={'sacc'};


%windows of analysis (do not change)
wind_targ=[0 400];
wind_sacc=[-50 50];%[-100 200]

%targ baseline
wind_targ_bsl=[-50 50];
wind_sacc_bsl=[-200 -150];%see compute_bsignif

%targ/sacc latency baseline
wind_targ_latbsl=[10 110];%overlap a little with beginning of burst and deflection %[0 100];%[-50 50];
wind_sacc_latbsl=[-100 -50];%[-200 -150];%[-150 -100];
    

%gaussian window for latency
gw_width=2.5;

%alpha of ttest
alpha=0.001;

%vshift
vshift_spk=100;
vshift_lfp=30;%29;


%sigma FR
sigma_FR=6;

%correct for filter phase shift introduced by ripple
pshift=4;%ms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get data
datalist=load_data_gandhilab(data_path);


           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analyzing data
dlist=get_dlist

alllat_spk=[];alllat_lfp=[];
allpeaks_spk=[];allpeaks_lfp=[];
allvar_spk=[];allvar_lfp=[];
allvarbsl_spk=[];allvarbsl_lfp=[];

data=[];
info=[];
dd=0;
for d=dlist
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
    
    %select trials
    seltrials=get_seltrials(data,'rpt');
    
    %bursts significance
    %targ_bsignif=data(1).offline.targ_bsignif;
    %sacc_bsignif=data(1).offline.sacc_bsignif;
    targ_bsignif=data(1).offline.targ_bsignif & data(1).offline.targ_bthresh';
    sacc_bsignif=data(1).offline.sacc_bsignif & data(1).offline.sacc_bthresh';
      
    
    
    %loop across all alignements
    aux_spk=[];aux_lfp=[];auxp_spk=[];auxp_lfp=[];auxv_spk=[];auxv_lfp=[];auxvbsl_spk=[];auxvbsl_lfp=[];
    for al=1%1:numel(alignlist)
        info.align=alignlist{al};
        
        %get all neural and behavioral data with specific alignement
        switch info.align
            case 'targ'
                [alltrials_spk,aligntime_spk]=get_alltrials_align(data,seltrials,wind_targ,'fr',info,targslist,sigma_FR,1);
                [alltrials_lfp,aligntime_lfp]=get_alltrials_align(data,seltrials,wind_targ,'lfp',info,targslist,sigma_FR,1);
                
                wind=wind_targ;
                wind_bsl=wind_targ_bsl;
                wind_latbsl=wind_targ_latbsl;
                burst_bsignif=targ_bsignif;
                
                %baseline latency
                [alltrials_spk_latbsl aligntime_latbsl]=get_alltrials_align(data,seltrials,wind_targ_latbsl,'fr',info,targslist,sigma_FR,1);
                [alltrials_lfp_latbsl aligntime_latbsl]=get_alltrials_align(data,seltrials,wind_targ_latbsl,'lfp',info,targslist,sigma_FR,1);
                
                %baseline
                [alltrials_spk_bsl aligntime_bsl]=get_alltrials_align(data,seltrials,wind_bsl,'fr',info,targslist,sigma_FR,1);
                [alltrials_lfp_bsl aligntime_bsl]=get_alltrials_align(data,seltrials,wind_bsl,'lfp',info,targslist,sigma_FR,1);
                
                
            case 'sacc'
                [alltrials_spk aligntime_spk]=get_alltrials_align(data,seltrials,wind_sacc,'fr',info,targslist,sigma_FR,1);
                [alltrials_lfp aligntime_lfp]=get_alltrials_align(data,seltrials,wind_sacc,'lfp',info,targslist,sigma_FR,1);
                
                wind=wind_sacc;
                wind_bsl=wind_sacc_bsl;
                wind_latbsl=wind_sacc_latbsl;
                burst_bsignif=sacc_bsignif;
                
                %baseline latency
                [alltrials_spk_latbsl aligntime_latbsl]=get_alltrials_align(data,seltrials,wind_sacc_latbsl,'fr',info,targslist,sigma_FR,1);
                [alltrials_lfp_latbsl aligntime_latbsl]=get_alltrials_align(data,seltrials,wind_sacc_latbsl,'lfp',info,targslist,sigma_FR,1);
                
                %baseline
                info.align='targ';
                [alltrials_spk_bsl aligntime_bsl]=get_alltrials_align(data,seltrials,wind_bsl,'fr',info,targslist,sigma_FR,1);
                [alltrials_lfp_bsl aligntime_bsl]=get_alltrials_align(data,seltrials,wind_bsl,'lfp',info,targslist,sigma_FR,1);
                info.align=alignlist{al};
                
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %analysis of trials for each target
        for tg=info.targ_tuning;%targs_ind,
            %target index
            info.targ=tg;
            
            %neural and behavioral signals for target tg
            trials_spk=alltrials_spk{tg};
            trials_spk_latbsl=alltrials_spk_latbsl{tg};
            trials_spk_bsl=alltrials_spk_bsl{tg};
            
            trials_lfp=alltrials_lfp{tg};
            trials_lfp_latbsl=alltrials_lfp_latbsl{tg};
            trials_lfp_bsl=alltrials_lfp_bsl{tg};
            
            %figure
            figtrials=figure('Position',[scrsz(3)/3 100 scrsz(3)/2 scrsz(4)-200]);
            
            
            
            %%%%%%%%%%%%%%%%%%
            %spk
            [info.nchannels info.ntrials info.triallen]=size(trials_spk);
            
            %compute average trials
            [trials_spk_avg trials_spk_var]=get_trials_avg(trials_spk);
            
            %average latency baseline
            [trials_spk_latbsl_avg trials_spk_latbsl_var]=get_trials_avg(trials_spk_latbsl);
            
            %baseline
            [trials_spk_bsl_avg trials_spk_bsl_var]=get_trials_avg(trials_spk_bsl);
            
            %normalization
            trials_spk_avgn=get_trials_normalized(trials_spk_avg,trials_spk_bsl_avg,'FR',info);
            trials_spk_latbsl_avgn=get_trials_normalized(trials_spk_latbsl_avg,trials_spk_bsl_avg,'FR',info);
            
%             %difference
%             trials_spk_avgn(2:16,:)=trials_spk_avgn(2:16,:)-trials_spk_avgn(1:15,:);
%             trials_spk_latbsl_avgn(2:16,:)=trials_spk_latbsl_avgn(2:16,:)-trials_spk_latbsl_avgn(1:15,:);
%                         
%             %trials_spk_avgn=trials_spk_avgn-mean(trials_spk_avgn);
            
            info.aligntime=aligntime_spk;
            hdlfig=subplot(1,2,1);hold on;
            titlestr={info.datafile ; ['FR ' info.align ' t' num2str(info.targ) '/' num2str(info.ntrials)]};
            [range vshift_spk]=plot_trials(trials_spk_avgn,[],[],[],[],[],info,hdlfig,titlestr);
            
            
            %%%%%%%
            %smoothing
            %NOTE:latency results are not significantly different w/o smoothing
            %gw=gausswin(gw_size,gw_width);gw=gw/sum(gw);
            for ch=1:size(trials_spk_avgn,1)
                %gausswin
                %trials_spk_avgn(ch,:)=conv(trials_spk_avgn(ch,:),gw','same');
                %moving avg
                trials_spk_avgn(ch,:)=smooth(trials_spk_avgn(ch,:),5,'moving');
                
            end
            
%             %%%%%%%
%             %get spk latency
%             lat_spk=get_latency(trials_spk_avgn,trials_spk_latbsl_avgn,1,'fr',info,alpha);
%             lat_spk(:,2)=lat_spk(:,2)+wind(1)-1;%correction for timing
%             
%             lat_spk_plot=lat_spk;
%             lat_spk_plot(:,3)=0;%lat_spk(:,3);
%             
%             %%plot_events_ch(lat_spk(:,2:3),[],vshift_spk,range,info,hdlfig,[],'-');
%             %plot_events_ch(lat_spk_plot(:,2:3).*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range,info,hdlfig,[],'-');
%             
%             
            %%%%%%%
            %get spk peaks
            %peaks_spk=get_peakfromlatency(trials_spk_avgn,lat_spk(:,2)-wind(1),'spk');
            %starting=lat_spk(:,2)-(wind(1)-1);
            starting=ones(info.nchannels,2);
            peaks_spk=get_peakfromlatency_minmax(trials_spk_avgn,starting,wind(1),info.align,'spk');
            peaks_spk(:,1)=peaks_spk(:,1)+wind(1)-1;%correction for timing
            
            figure(figtrials)
            hdlfig=subplot(1,2,1);hold on;
            %plot_events_ch(peaks_spk,[],vshift_spk,range,info,hdlfig,[],'--');
            %taking care of bsignif for plotting
            peaks_spk_plot=peaks_spk.*[burst_bsignif ; burst_bsignif]';
            peaks_spk_plot(find(burst_bsignif==0),1)=wind(1);peaks_spk_plot(find(burst_bsignif==0),2)=wind(1);       
            plot_events_ch(peaks_spk_plot,[],vshift_spk,range,info,hdlfig,[],'--');
            grid
            
            %%%%%%%
            %spk latency from peak
            plat_spk=get_latencyfrompeak(trials_spk_avgn,peaks_spk-(wind(1)-1),'spk',info,alpha);
            plat_spk(:,2)=plat_spk(:,2)+wind(1)-1;%correction for timing
            %plot_events_ch(plat_spk(:,2:3),[],vshift_spk,range,info,hdlfig,[],'-');
            %taking care of bsignif for plotting
            plat_spk_plot(:,2:3)=plat_spk(:,2:3).*[burst_bsignif ; burst_bsignif]';
            plot_events_ch(plat_spk_plot(:,2:3),[],vshift_spk,range,info,hdlfig,[],'-');
            
            
            %get spk snr
            var_spk=var(trials_spk_avgn,[],2);
            var_spkbsl=var(trials_spk_latbsl_avgn,[],2);
            
            
            
            %%%%%%%%%%%%%%%%
            %lfp
            [info.nchannels info.ntrials info.triallen]=size(trials_lfp);
            
            %compute average trials
            [trials_lfp_avg trials_lfp_var]=get_trials_avg(trials_lfp);
            
            %average latency baseline
            [trials_lfp_latbsl_avg trials_lfp_latbsl_var]=get_trials_avg(trials_lfp_latbsl);
            
            %baseline
            [trials_lfp_bsl_avg trials_lfp_bsl_var]=get_trials_avg(trials_lfp_bsl);
            
            %             %detrend lfp
            %             trials_lfp_avg=detrend(trials_lfp_avg);
            %             trials_lfp_latbsl_avg=detrend(trials_lfp_latbsl_avg);
            %             trials_lfp_bsl_avg=detrend(trials_lfp_bsl_avg);
            
            %normalization
            trials_lfp_avgn=get_trials_normalized(trials_lfp_avg,trials_lfp_bsl_avg,'LFP',info);
            trials_lfp_latbsl_avgn=get_trials_normalized(trials_lfp_latbsl_avg,trials_lfp_bsl_avg,'LFP',info);
            
%             %difference
%             trials_lfp_avgn(2:16,:)=trials_lfp_avgn(2:16,:)-trials_lfp_avgn(1:15,:);
%             trials_lfp_latbsl_avgn(2:16,:)=trials_lfp_latbsl_avgn(2:16,:)-trials_lfp_latbsl_avgn(1:15,:);
%             
%             trials_lfp_avgn(2:16,:)=trials_lfp_avgn(2:16,:)-trials_lfp_avgn(1:15,:);
%             trials_lfp_latbsl_avgn(2:16,:)=trials_lfp_latbsl_avgn(2:16,:)-trials_lfp_latbsl_avgn(1:15,:);
            
%              trials_lfp_avgn=trials_lfp_avgn-mean(trials_lfp_avgn);
%              trials_lfp_latbsl_avgn=trials_lfp_latbsl_avgn-mean(trials_lfp_latbsl_avgn);

            info.aligntime=aligntime_lfp;
            hdlfig=subplot(1,2,2);hold on;
            titlestr={info.datafile ; ['LFP ' info.align ' t' num2str(info.targ) '/' num2str(info.ntrials)]};
            [range vshift_lfp]=plot_trials(trials_lfp_avgn,[],[],[],[],[],info,hdlfig,titlestr);
            
            %%%%%%%
            %smoothing
            %NOTE:latency results are not significantly different w/o smoothing
            %gw=gausswin(gw_size,gw_width);gw=gw/sum(gw);
            for ch=1:size(trials_lfp_avgn,1)
                %gausswin
                %trials_lfp_avgn(ch,:)=conv(trials_lfp_avgn(ch,:),gw','same');
                %moving avg
                trials_lfp_avgn(ch,:)=smooth(trials_lfp_avgn(ch,:),5,'moving');
                
            end
            
%             %%%%%%%
%             %get lfp latency
%             starting=100-wind_latbsl(2)/2;
%             lat_lfp=get_latency(trials_lfp_avgn,trials_lfp_latbsl_avgn,starting,'lfp',info,alpha);
%             lat_lfp(:,2)=lat_lfp(:,2)+wind(1)-1;%correction for timing
%             
%             %%plot_events_ch(lat_lfp(:,2:3),[],vshift_lfp,range,info,hdlfig,[],'-');
%             %plot_events_ch(lat_lfp(:,2:3).*[burst_bsignif ; burst_bsignif]',[],vshift_lfp,range,info,hdlfig,[],'-');
            
            
            %%%%%%%
            %get lfp peaks
            %starting=ones(info.nchannels,2);
            %peaks_lfp=get_peakfromlatency(trials_lfp_avgn,starting,'lfp');
            %starting=lat_lfp(:,2)-(wind(1)-1);
            starting=ones(info.nchannels,2);
            peaks_lfp=get_peakfromlatency_minmax(trials_lfp_avgn,starting,wind(1),info.align,'lfp');
            peaks_lfp(:,1)=peaks_lfp(:,1)+wind(1)-1;%correction for timing
            
            figure(figtrials)
            hdlfig=subplot(1,2,2);hold on;
            %plot_events_ch(peaks_lfp,[],vshift_lfp,range,info,hdlfig,[],'--');
            %taking care of bsignif for plotting
            peaks_lfp_plot=peaks_lfp.*[burst_bsignif ; burst_bsignif]';
            peaks_lfp_plot(find(burst_bsignif==0),1)=wind(1);peaks_lfp_plot(find(burst_bsignif==0),2)=wind(1);       
            plot_events_ch(peaks_lfp_plot,[],vshift_lfp,range,info,hdlfig,[],'--');

            grid
            
            %%%%%%%
            %lfp latency from peak
            plat_lfp=get_latencyfrompeak(trials_lfp_avgn,peaks_lfp-(wind(1)-1),'lfp',info,alpha);
            plat_lfp(:,2)=plat_lfp(:,2)+wind(1)-1;%correction for timing
            %plot_events_ch(plat_lfp(:,2:3),[],vshift_lfp,range,info,hdlfig,[],'-');
            %taking care of bsignif for plotting
            plat_lfp_plot(:,2:3)=plat_lfp(:,2:3).*[burst_bsignif ; burst_bsignif]';
            plot_events_ch(plat_lfp_plot(:,2:3),[],vshift_lfp,range,info,hdlfig,[],'-');
            
          
            
            %correct for filter phase shift introduced by ripple
            %lat_lfp(:,2)=lat_lfp(:,2)-pshift;%correction for pshift by ripple
            peaks_lfp(:,1)=peaks_lfp(:,1)-pshift;%correction for pshift by ripple
            plat_lfp(:,2)=plat_lfp(:,2)-pshift;%correction for pshift by ripple
          
            
            %plot lfp latency on spk plot 
            %(introduced by ripple)
            plat_lfp_onspk=plat_lfp;
            plat_lfp_onspk(:,3)=plat_spk(:,3);
            
            hdlfig=subplot(1,2,1);hold on;
            plot_events_ch((plat_lfp_onspk(:,2:3)).*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range,info,hdlfig,[],'-.');
            
            
            %get lfp snr
            var_lfp=var(trials_lfp_avgn,[],2);
            var_lfpbsl=var(trials_lfp_latbsl_avgn,[],2);

            
            %%%%%%%%%%%%%%%%
            %alignment of lat using CSD features (after compute_CSDfeature)
            info.csdfeat_avg_targ=data(1).offline.csdfeat_avg_targ;
            info.zs=data(1).offline.csdzs;
            dref=info.csdfeat_avg_targ(2);
%             [lataux_r1 info_r ch_ref ~]=get_vmis_aligndepth(lat_spk(:,2)',dref,info);
%             [lataux_r2 info_r ch_ref ~]=get_vmis_aligndepth(lat_spk(:,3)',dref,info);
%             lat_spk_r(:,2)=lataux_r1;
%             lat_spk_r(:,3)=lataux_r2;
%             [lataux_r1 info_r ch_ref ~]=get_vmis_aligndepth(lat_lfp(:,2)',dref,info);
%             [lataux_r2 info_r ch_ref ~]=get_vmis_aligndepth(lat_lfp(:,3)',dref,info);
%             lat_lfp_r(:,2)=lataux_r1;
%             lat_lfp_r(:,3)=lataux_r2;
            
            
            %%%%%%%%%%%%%%%%
            %alignment of peaks using CSD features (after compute_CSDfeature)
            [peaksaux_r1 info_r ch_ref ~]=get_vmis_aligndepth(peaks_spk(:,1)',dref,info);
            [peaksaux_r2 info_r ch_ref ~]=get_vmis_aligndepth(peaks_spk(:,2)',dref,info);
            peaks_spk_r(:,1)=peaksaux_r1;
            peaks_spk_r(:,2)=peaksaux_r2;
            [peaksaux_r1 info_r ch_ref ~]=get_vmis_aligndepth(peaks_lfp(:,1)',dref,info);
            [peaksaux_r2 info_r ch_ref ~]=get_vmis_aligndepth(peaks_lfp(:,2)',dref,info);
            peaks_lfp_r(:,1)=peaksaux_r1;
            peaks_lfp_r(:,2)=peaksaux_r2;
            
            
            %%%%%%%%%%%%%%%%
            %alignment of plat using CSD features (after compute_CSDfeature) 
            [plataux_r1 info_r ch_ref ~]=get_vmis_aligndepth(plat_spk(:,2)',dref,info);
            [plataux_r2 info_r ch_ref ~]=get_vmis_aligndepth(plat_spk(:,3)',dref,info);
            plat_spk_r(:,2)=plataux_r1;
            plat_spk_r(:,3)=plataux_r2;
            [plataux_r1 info_r ch_ref ~]=get_vmis_aligndepth(plat_lfp(:,2)',dref,info);
            [plataux_r2 info_r ch_ref ~]=get_vmis_aligndepth(plat_lfp(:,3)',dref,info);
            plat_lfp_r(:,2)=plataux_r1;
            plat_lfp_r(:,3)=plataux_r2;
            
            
            %%%%%%%%%%%%%%%%
            %alignment of var and var bsl using CSD features (after compute_CSDfeature)
            [var_spk_r info_r ch_ref ~]=get_vmis_aligndepth(var_spk',dref,info);
            [var_spkbsl_r info_r ch_ref ~]=get_vmis_aligndepth(var_spkbsl',dref,info);
            [var_lfp_r info_r ch_ref ~]=get_vmis_aligndepth(var_lfp',dref,info);
            [var_lfpbsl_r info_r ch_ref ~]=get_vmis_aligndepth(var_lfpbsl',dref,info);
            
            
            %%%%%%%%%%%%%%%%%%
            %list of latencies and peaks of activity
%             aux_spk(al,:,:)=lat_spk_r;
%             aux_lfp(al,:,:)=lat_lfp_r;
            auxp_spk(al,:,:)=peaks_spk_r;
            auxp_lfp(al,:,:)=peaks_lfp_r;
            auxplat_spk(al,:,:)=plat_spk_r;
            auxplat_lfp(al,:,:)=plat_lfp_r;
            auxv_spk(al,:)=var_spk_r';
            auxv_lfp(al,:)=var_lfp_r';
            auxvbsl_spk(al,:)=var_spkbsl_r';
            auxvbsl_lfp(al,:)=var_lfpbsl_r';
            
            %pause
            %close(figtrials)
            
        end
        %test of significant burst
        bsignif=find(burst_bsignif==0);
%         aux_spk(al,bsignif,:)=nan;
%         aux_lfp(al,bsignif,:)=nan;
        auxp_spk(al,bsignif,:)=nan;
        auxp_lfp(al,bsignif,:)=nan;
        auxplat_spk(al,bsignif,:)=nan;
        auxplat_lfp(al,bsignif,:)=nan;
        auxv_spk(al,bsignif)=nan;
        auxv_lfp(al,bsignif)=nan;
        auxvbsl_spk(al,bsignif)=nan;
        auxvbsl_lfp(al,bsignif)=nan;
        
    end
%     alllat_spk{dd}=aux_spk;
%     alllat_lfp{dd}=aux_lfp;
    allpeaks_spk{dd}=auxp_spk;
    allpeaks_lfp{dd}=auxp_lfp;
    allplat_spk{dd}=auxplat_spk;
    allplat_lfp{dd}=auxplat_lfp;
    allvar_spk{dd}=auxv_spk;
    allvar_lfp{dd}=auxv_lfp;
    allvarbsl_spk{dd}=auxvbsl_spk;
    allvarbsl_lfp{dd}=auxvbsl_lfp;
    
    %pause
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot peaks timing and magnitude

%latency or magnitude of peaks
peaksplat=1;% peaks 1 plat 2
timag=2;%peaks timing 1 peaks magnitude 2

colorlist=get_colorlist;
nchannels=48;%16;
channels=[1:nchannels];

r_mag=[];
r_lat=[];
figregress=figure('Position',[scrsz(3)/4 100 scrsz(3)/2 scrsz(4)-200]);
ddlist=[1:length(dlist)];
peaks_spk_all=[];peaks_lfp_all=[];var_spk_all=[];var_lfp_all=[];varbsl_spk_all=[];varbsl_lfp_all=[];
for dd=1:size(allpeaks_spk,2),
    info.datafile=datalist{dlist(dd)};
    
    %targ or sacc
    for al=1%1:2
        info.align=alignlist{al};
          
        
        %%%%%%%%%%%%%%%%%%%%%
        if peaksplat==1
            aux_spk=allpeaks_spk{dd};
            aux_lfp=allpeaks_lfp{dd};
            
            if timag==1
                %peaks timing
                peaks_spk=squeeze(aux_spk(al,:,1));
                peaks_lfp=squeeze(aux_lfp(al,:,1));
                
            elseif timag==2
                %peaks magnitude
                peaks_spk=squeeze(aux_spk(al,:,2));
                peaks_lfp=squeeze(aux_lfp(al,:,2));
            end
            
        elseif peaksplat==2
            aux_spk=allplat_spk{dd};
            aux_lfp=allplat_lfp{dd};
            if timag==1
                %peaks timing
                peaks_spk=squeeze(aux_spk(al,:,2));
                peaks_lfp=squeeze(aux_lfp(al,:,2));
                
            elseif timag==2
                %peaks magnitude
                peaks_spk=squeeze(aux_spk(al,:,3));
                peaks_lfp=squeeze(aux_lfp(al,:,3));
            end
        end
        
       
        %var
        var_spk=allvar_spk{dd};
        var_lfp=allvar_lfp{dd};
        varbsl_spk=allvarbsl_spk{dd};
        varbsl_lfp=allvarbsl_lfp{dd};
        
        %lists
        peaks_spk_all=[peaks_spk_all ; peaks_spk];
        peaks_lfp_all=[peaks_lfp_all ; peaks_lfp];
        
        var_spk_all=[var_spk_all ; var_spk];
        var_lfp_all=[var_lfp_all ; var_lfp];
        varbsl_spk_all=[varbsl_spk_all ; varbsl_spk];
        varbsl_lfp_all=[varbsl_lfp_all ; varbsl_lfp];
        
         
        %plot fr peaks
        hdlfig=subplot(2,2,1+3*(al-1));hold on;
        titlestr={info.datafile ; info.align};
        plot(peaks_spk,1:nchannels,'-o','color',colorlist(ddlist(dd),:));
        
        ylabel('Channel');
        if timag==1
            xlabel('FR latency');
        elseif timag==2
            xlabel('FR magnitude')
        end
        axis tight
        
        
        %plot lfp peaks
        hdlfig=subplot(2,2,2+3*(al-1));hold on;
        %titlestr={info.datafile ; info.align};
        plot(peaks_lfp,1:nchannels,'--s','color',colorlist(ddlist(dd),:));
        ylabel('Channel');
        if timag==1
            xlabel('LFP latency');
        elseif timag==2
            xlabel('LFP magnitude')
        end
        axis tight
           
        %         %%%%%%%%%%%%%%%%%%%
        %         %plot regression fr/lfp on latency
        %         titlestr=info.align;
        %         hdlfig=subplot(2,2,3+3*(al-1));hold on;
        %         %plot_regress(lat_spk([1:nchannels]),lat_lfp([1:nchannels]),'Latency',ddlist(dd),info,hdlfig,titlestr);
        %         plot_regress(lat_spk([1:nchannels]),lat_lfp([1:nchannels]),'Latency',ddlist(dd),info,hdlfig,titlestr);
        %         axis tight

    end
    %pause
    %figure(figregress); title(info.datafile);
    
end


%plot average and ci
color_avgconf=['b' 'r']
for sig=1:2
    switch sig
        case 1
            peaks_all=peaks_spk_all;
        case 2
            peaks_all=peaks_lfp_all;
    end
    
    %remove outliers
    switch info.align
        case 'targ'
            if timag==1
                %timing
                peaks_all(find(peaks_all<60 | peaks_all>150))=nan;
            elseif timag==2
                %magnitude
                %peaks_all(find(peaks_all<60 | peaks_all>200))=nan;
            end
        case 'sacc'
            if timag==1
                %timing
                outliers(sig,:)=find(peaks_all<-30 | peaks_all>50);
                peaks_all(outliers(sig,:))=nan;
            elseif timag==2
                %magnitude
                %NOTE: should actually remove the outiers given by the
                %timing
                %peaks_all(find(peaks_all<-200 | peaks_all>200))=nan;
                peaks_all(outliers(sig,:))=nan;
            end
    end
    
    
    %plot mean and 95% confidence intervals
    peaks_avg=nanmean(peaks_all,1);
    
    
    subplot(2,2,sig);hold on;
    subplot(2,2,1);hold on;
    plot(peaks_avg,1:nchannels,color_avgconf(sig),'Linewidth',3);
    
    
    %find channel range
    chs_r=find(~isnan(peaks_avg));
    [vmiss imiss]=find(chs_r(2:end)-chs_r(1:end-1)>1);
    
    %consider only the last consecutive channels
    %     %if ~isempty(imiss),min_ch=chs_r(max(imiss)+1);else min_ch=chs_r(1);end
    %min_ch=chs_r(1);
    %     if ~isempty(imiss),max_ch=chs_r(imiss-1);else max_ch=max(chs_r);end
    %max_ch=chs_r(imiss(1)-1);
    
    %alt
    if timag==1 & sig==1
        min_ch=chs_r(1);
        %min_ch=chs_r(imiss(1)+1);
        max_ch=max(chs_r);
        %max_ch=chs_r(imiss(1)-1);
        %max_ch=chs_r(imiss(2)-1);
        pause
    elseif timag==1 & sig==2
        min_ch=chs_r(1);
        %min_ch=chs_r(imiss(1)+1);
        max_ch=max(chs_r);
        %max_ch=chs_r(imiss(1)-1);
        
    elseif timag==2 & sig==1
        min_ch=chs_r(1);%chs_r(imiss(1)+1);
        max_ch=max(chs_r);
        %max_ch=chs_r(imiss(1)-1);
    elseif timag==2 & sig==2
        min_ch=chs_r(1);
        max_ch=max(chs_r);
        %max_ch=chs_r(imiss(1)-1);
        
    end
    
    %compute ci
    ind=0;peaks_ci=[];
    for ch=min_ch:max_ch,
        ind=ind+1;
        
        aux=(peaks_all(find(~isnan(peaks_all(:,ch))),ch));
        
        if numel(aux)<=1,
            peaks_ci(ind,:)=[peaks_avg(ch) ; peaks_avg(ch)];
        else
            peaks_ci(ind,:) = bootci(2000,{@mean,aux},'type','per');
        end
    end
    fill([peaks_ci(:,1)' fliplr(peaks_ci(:,2)')],[min_ch:max_ch max_ch:-1:min_ch], 1,'facecolor',color_avgconf(sig),'edgecolor','none','facealpha', 0.3);
    
    lims=[16 33];
    subplot(2,2,sig);hold on;
    switch info.align
        case 'targ'
            if timag==1 & sig==1
                axis([50 150 lims(1) lims(2)])
                set(gca,'Ytick',[lims(1):2:lims(2)+1],'Yticklabel',[-8:2:10])
            elseif timag==1 & sig==2
                axis([50 150 lims(1) lims(2)])
                set(gca,'Ytick',[lims(1):2:lims(2)+1],'Yticklabel',[-8:2:10])
                
            elseif timag==2 & sig==1
                axis([-90 150 lims(1) lims(2)])
                set(gca,'Ytick',[lims(1):2:lims(2)+1],'Yticklabel',[-8:2:10])
            elseif timag==2 & sig==2
                axis([-90 150 lims(1) lims(2)])
                set(gca,'Ytick',[lims(1):2:lims(2)+1],'Yticklabel',[-8:2:10])
            end
        case 'sacc'
            if timag==1 & sig==1
                axis([-30 30 lims(1) lims(2)])
                set(gca,'Ytick',[lims(1):2:lims(2)+1],'Yticklabel',[-8:2:10])
            elseif timag==1 & sig==2
                axis([-30 30 lims(1) lims(2)])
                set(gca,'Ytick',[lims(1):2:lims(2)+1],'Yticklabel',[-8:2:10])
                
            elseif timag==2 & sig==1
                axis([-50 250 lims(1) lims(2)])
                set(gca,'Ytick',[lims(1):2:lims(2)+1],'Yticklabel',[-8:2:10])
            elseif timag==2 & sig==2      
                axis([-50 250 lims(1) lims(2)])
                set(gca,'Ytick',[lims(1):2:lims(2)+1],'Yticklabel',[-8:2:10])
            end
    end
            
    display(['average: ' num2str(mean(peaks_avg(min_ch:max_ch))) ' and variance: ' num2str(var(peaks_avg(min_ch:max_ch)))])
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%stats
%channels of interest
min_ch=16;
max_ch=33;
p=[];h=[];

%%latencies without outliers
%lat_spk_all(find(lat_spk_all<60 | lat_spk_all>150))=nan;
%lat_lfp_all(find(lat_lfp_all<60 | lat_lfp_all>150))=nan;

%normality test (Kolmogorov-Smirnov)
%spk
chi=0;
for ch=min_ch:max_ch %limit test at channels of interest
    chi=chi+1;
    [h(chi),p(chi)] = kstest(peaks_spk_all(:,ch));
end
%output
h
p


%lfp
chi=0;
for ch=min_ch:max_ch
    chi=chi+1;
    [h(chi),p(chi)] = kstest(peaks_lfp_all(:,ch));
end
%output
h
p


%conclusions
%all latencies channels have a normal distribution!!




%%%%%%%%%%%%%%%%%%%%%%%%%
%parametric test

chi=0;
for ch=min_ch:max_ch
    chi=chi+1;
    %[p(chi),h(chi),stats] = ranksum(peaks_spk_all(:,ch),peaks_lfp_all(:,ch));
    [h(chi),p(chi),stats] = ttest2(peaks_spk_all(:,ch),peaks_lfp_all(:,ch));
end
h
p

display(['Significant difference for channels:' num2str(find(p<0.01) + min_ch-1 -23)])

%results:
%Significant difference for channels:14 (29)

%%%%%%%%%%%%%%%%%%%%%%%%%
%parametric pair ttest

chi=0;
for ch=min_ch:max_ch
    chi=chi+1;
    [h(chi),p(chi),stats] = ttest(peaks_spk_all(:,ch),peaks_lfp_all(:,ch));
end
h
p

display(['Significant pair ttest difference for channels:' num2str(find(p<0.01) + min_ch-1 -23)])

%results:
%Significant difference for channels:14 (29)



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot ranked histograms
color_avgconf=['.b' '.r']
[npeaks nchannels]=size(peaks_spk_all);

min_ch=16;max_ch=33;

%%remove outliers
%peaks_spk_all(find(peaks_spk_all<60 | peaks_spk_all>150))=nan;
%peaks_lfp_all(find(peaks_lfp_all<60 | peaks_lfp_all>200))=nan;

%rank order
peaks_spk_rk=[];peaks_lfp_rk=[];
for ch=min_ch:max_ch
    [peaks_spk_rk(:,ch) i_rk]=sort(peaks_spk_all(:,ch));
    peaks_lfp_rk(:,ch)=peaks_lfp_all(i_rk,ch);

    %[peaks_lfp_rk(:,ch) i_rk]=sort(peaks_lfp_all(:,ch));
    %peaks_spk_rk(:,ch)=peaks_spk_all(i_rk,ch);
end

figure;hold on;
for ch=min_ch:max_ch,
    subplot(5,4,ch-min_ch+1);hold on;
    if ch==min_ch,title('Peaks timing / ch');end
    plot(peaks_spk_rk(:,ch),1:npeaks,'o','MarkerSize',5,'MarkerFaceColor','b');
    plot(peaks_lfp_rk(:,ch),1:npeaks,'o','MarkerSize',5,'MarkerFaceColor','r');
    switch info.align
        case 'targ'
            axis([50 150 1 npeaks+1]);axis square;
        case 'sacc'
            axis([-30 30 1 npeaks+1]);axis square;
    end
    xlabel(['ch' num2str(ch-min_ch+1)])
end


%rank order by peaks timing differences
peaks_diff=[];
for ch=min_ch:max_ch
    peaks_diff=abs(peaks_spk_all(:,ch)-peaks_lfp_all(:,ch));
    [val_rk i_rk]=sort(peaks_diff);
    peaks_spk_rk(:,ch)=peaks_spk_all(i_rk,ch);
    peaks_lfp_rk(:,ch)=peaks_lfp_all(i_rk,ch);
end

figure;hold on;
for ch=min_ch:max_ch,
    subplot(5,4,ch-min_ch+1);hold on;
    if ch==min_ch,title('peaks timing (sort by diff) / ch');end
    plot(peaks_spk_rk(:,ch),1:npeaks,'o','MarkerSize',5,'MarkerFaceColor','b');
    plot(peaks_lfp_rk(:,ch),1:npeaks,'o','MarkerSize',5,'MarkerFaceColor','r');
    switch info.align
        case 'targ'
            axis([50 150 1 npeaks+1]);axis square;
        case 'sacc'
            axis([-30 30 1 npeaks+1]);axis square;
    end
    xlabel(['ch' num2str(ch-min_ch+1)])
end


% %%
% %peaks timing differences vs. snr
% lat_diff=[];snr_spk=[];snr_lfp=[];
% lat_diff_rk=[];snr_spk_rk=[];snr_lfp_rk=[];
% varbsl_spk_rk=[];varbsl_lfp_rk=[];
% 
% for ch=min_ch:max_ch
%     
%     %latency difference
%     lat_diff(:,ch)=(lat_lfp_all(:,ch)-lat_spk_all(:,ch));
%     [lat_diff_rk(:,ch) i_rk]=sort(lat_diff(:,ch));
%     
%     %latency
%     %[lat_spk_rk(:,ch) i_rk]=sort(lat_spk_all(:,ch));
%     %lat_lfp_rk(:,ch)=lat_lfp_all(i_rk,ch);
% 
%     %[lat_lfp_rk(:,ch) i_rk]=sort(lat_lfp_all(:,ch));
%     %lat_spk_rk(:,ch)=lat_spk_all(i_rk,ch);
%     
%     %snr
%     snr_spk(:,ch)=var_spk_all(:,ch)./varbsl_spk_all(:,ch);
%     snr_lfp(:,ch)=var_lfp_all(:,ch)./varbsl_lfp_all(:,ch);
%     snr_spk_rk(:,ch)=snr_spk(i_rk,ch);
%     snr_lfp_rk(:,ch)=snr_lfp(i_rk,ch);
%     
%     %var
%     varbsl_spk_rk(:,ch)=varbsl_spk_all(i_rk,ch);
%     varbsl_lfp_rk(:,ch)=varbsl_lfp_all(i_rk,ch);
%     
% end
% 
% 
% figure;hold on;
% for ch=min_ch:max_ch,
%     subplot(5,4,ch-min_ch+1);hold on;
%     if ch==min_ch,title('Lat diff / snr');end
%     plot(lat_diff_rk(:,ch),snr_spk_rk(:,ch),'o','MarkerSize',5,'MarkerFaceColor','b');
%     plot(lat_diff_rk(:,ch),snr_lfp_rk(:,ch),'o','MarkerSize',5,'MarkerFaceColor','r');
%     axis([-20 20 0 50]);axis square;
%     xlabel(['ch' num2str(ch-min_ch+1)])
% end
% 
% figure;hold on;
% for ch=min_ch:max_ch,
%     subplot(5,4,ch-min_ch+1);hold on;
%     if ch==min_ch,title('Lat diff / var');end
%     plot(lat_diff_rk(:,ch),varbsl_spk_rk(:,ch),'o','MarkerSize',5,'MarkerFaceColor','b');
%     plot(lat_diff_rk(:,ch),varbsl_lfp_rk(:,ch),'o','MarkerSize',5,'MarkerFaceColor','r');
%     axis([-20 20 0 50]);axis square;
%     xlabel(['ch' num2str(ch-min_ch+1)])
% end









