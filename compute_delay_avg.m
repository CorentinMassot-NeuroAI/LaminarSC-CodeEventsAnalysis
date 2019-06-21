%function compute_delay_avg

%function compute_delay_avg
%   Compute activity during delay period based on average activity recorded with a
%   laminar probe (LMA)
%
%
% see also analysis_delay_avg
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 02/28/2018 last modified 02/22/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TO DO
% display amplitude of visual burst / visual burst / delay period / snippet
% of additional delay for increaseing go cues timing / motor buildup and
% burst / latency of buildup / amplitude of burst

% decrease of visual burst across depth


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set paths
[root_path data_path save_path]=set_paths;
save_path=[root_path 'Results\Results_SC_delay\'];


%screen size
scrsz = get(groot,'ScreenSize');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
%save data
savedata=0;

%display plot
disp_plot=0;

%alignement
%alignlist={'no' 'targ' 'go' 'sacc'};
%alignlist={'targ'};
alignlist={'targ_pburst_ch' };

%windows of analysis (do not change)
wind_targ=[-50 1500];%[-50 150];
wind_go=[-1000 500];

%windows baseline
wt=100;
%wind_bsl_targ=[-wt 0 ];
wind_bsl_targ=[-50-wt -50 ] %newbsl

%wind_bsl_sacc=[-100 0];

%vshift
vshift_spk=100;
vshift_lfp=30;%29;

%sigma FR
sigma_FR=6;


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
    %targ_pburst
    ratios_targ=data(1).offline.targ_pburst_ratio(info.targ_tuning,:)>thresh_ratios;
    surprises_targ=data(1).offline.targ_pburst_msurprises(info.targ_tuning,:)>thresh_surprises;
    bsignif_targ=data(1).offline.targ_pburstch_bsignif;
    bthresh_targ=data(1).offline.targ_pburstch_bthresh_trials';
    %     %sacc
    %     ratios_sacc=data(1).offline.sacc_pburst_ratio(info.targ_tuning,:)>thresh_ratios;
    %     surprises_sacc=data(1).offline.sacc_pburst_msurprises(info.targ_tuning,:)>thresh_surprises;
    %     bsignif_sacc=data(1).offline.sacc_bsignif;
    %     bthresh_sacc=data(1).offline.sacc_bthresh_trials';
    
    %targ_bsignif=(ratios_targ & surprises_targ & bsignif_targ & bthresh_targ)
    targ_bsignif=(bsignif_targ & bthresh_targ);
    %sacc_bsignif=(ratios_sacc & surprises_sacc & bsignif_sacc & bthresh_sacc);
    %sacc_bsignif=(bsignif_sacc & bthresh_sacc);
    
    
    
    %loop across all alignements
    aux_spk=[];aux_lfp=[];auxp_spk=[];auxp_lfp=[];auxv_spk=[];auxv_lfp=[];auxvbsl_spk=[];auxvbsl_lfp=[];
    for al=1%1:numel(alignlist)
        info.align=alignlist{al};
        
        %get all neural and behavioral data with specific alignement
        switch info.align
            case 'targ'
                wind=wind_targ;
                wind_bsl=wind_bsl_targ;
                burst_bsignif=targ_bsignif;
            case 'targ_pburst_ch'
                wind=wind_targ;
                wind_bsl=wind_bsl_targ;
                burst_bsignif=targ_bsignif;
            case 'go'
                wind=wind_go;
                wind_bsl=wind_bsl_targ;
                burst_bsignif=targ_bsignif;
        end
        
        %signals
        [alltrials_spk,aligntime_spk lut_trials]=get_alltrials_align(data,seltrials,wind,'fr',info,targslist,sigma_FR,1);
        %baseline
        [alltrials_spk_bsl aligntime_bsl]=get_alltrials_align(data,seltrials,wind_bsl,'fr',info,targslist,sigma_FR,1);
        
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
            
            
            %%%%%%%%%%%%%%%%%%
            %figure
            figtrials=figure('Position',[scrsz(3)/3 100 scrsz(3)/2 scrsz(4)-200]);
            
            
            %%%%%%%%%%%%%%%%%%
            %spk
            [info.nchannels info.ntrials info.triallen]=size(trials_spk);
            info.aligntime=aligntime_spk;
            
            %compute average of normalized trials in RF
            trials_spk_n=get_trials_normalized(trials_spk,trials_spk_bsl,'FR',info);
            [trials_spk_n_avg trials_spk_n_var]=get_trials_avg(trials_spk_n);

            %plot avg and ste
            figure(figtrials)
            hdlfig=subplot(1,1,1);hold on;
            [range_spk vshift_spk]=plot_trials(trials_spk_n_avg,[],[],[],[],[],info,hdlfig,[],[],[]);
            plot_trials(trials_spk_n_avg+trials_spk_n_var,[],[],[],[],[],info,hdlfig,[],'-',1);
            plot_trials(trials_spk_n_avg-trials_spk_n_var,[],[],[],[],[],info,hdlfig,[],'-',1);
            
            
            %%
            %%%%%%%%%%%%%%%%%%
            %plot all 'go' cues across trials
            switch info.align
                case {'targ','targ_pburst_ch'}
                    go_onsets=zeros(numel(lut_trials{info.targ_tuning}),info.nchannels,2);
                    for t=lut_trials{info.targ_tuning}
                        go_onset = get_event(data(t),'goCode','')-info.aligntime;
                        targ_onset = get_event(data(t),'targCode','')-info.aligntime;
                        go_onsets(t,:,1)=ones(info.nchannels,1)*(go_onset-targ_onset);
                        
                        plot_events_ch(squeeze(go_onsets(t,:,:)).*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','-',1,'b');
                        %plot_events_ch(onsets,[],vshift_spk,range_spk,info,hdlfig,'n','-',1,'k');
                    end
            end
            grid
            
            %soonest go onset
            aux=sort(unique(go_onsets(:)));
            go_snst=aux(2);
            
            
            %%
            %%%%%%%%%%%%%%%%%%
            %peak activity
            [peak_vals peak_times]=max(trials_spk_n_avg(:,1:info.aligntime+200),[],2);
            peaks(:,1)=peak_times;
            peaks(:,2)=peak_vals;
            peaks(:,1)=peaks(:,1)-info.aligntime;%correction for timing
            
            %plot peak
            peak_plot=[];
            peak_plot(:,1)=peaks(:,1);
            peak_plot(:,2)=peaks(:,2);
            plot_events_ch(peak_plot.*[burst_bsignif ; burst_bsignif]',[],vshift_spk,range_spk,info,hdlfig,'n','-',2,'k');
            %plot_events_ch(peak_plot,[],vshift_spk,range_spk,info,hdlfig,'n','-',1,'');

            %latest peak time
            peak_ltst=max(peaks(:,1));
        
            
            %%
            %%%%%%%%%%%%%%%%%%
            %compute binned activity across channels            
            figdelaybins=figure('Position',[scrsz(3)/3 200 scrsz(3)/1.8 scrsz(4)-400]);hold on;
            %figdelaybinsn=figure;hold on;
            binsize=50;
            nbins=floor((go_snst-peak_ltst+1)/binsize);
            delay_bins=zeros(info.nchannels,nbins);delay_bins_n=zeros(info.nchannels,nbins);
            legend_bins={};
            colorlist2 = colormap(jet(nbins));%hot hsv copper
            for bi=1:nbins
                delay_bins(:,bi)=nanmean(trials_spk_n_avg(:,peak_ltst+(bi-1)*binsize:peak_ltst+bi*binsize-1),2);
            
                %without normalization
                figure(figdelaybins)
                subplot(1,2,1);hold on;
                plot(delay_bins(:,bi),1:info.nchannels,'color',colorlist2(bi,:))
                
                %with peak normalization
                delay_bins_n(:,bi)=delay_bins(:,bi)./peaks(:,2);
                figure(figdelaybins)
                subplot(1,2,2);hold on;
                plot(delay_bins_n(:,bi),1:info.nchannels,'color',colorlist2(bi,:))
            
               legend_bins{bi}=num2str(bi);
            end
            
            figure(figdelaybins)
            subplot(1,2,1);
            axis([-20 max(peaks(:,2)) 1 info.nchannels])
            hdl=line([0 0],[1 info.nchannels]);
            set(hdl,'color','k','linestyle','--')
            grid
            legend(legend_bins)
            ylabel('Channel');xlabel('Firing rate (spk/s)');
            title('Delay activity')
            subplot(1,2,2);
            axis([-0.3 1 1 info.nchannels])
            hdl=line([0 0],[1 info.nchannels]);
            set(hdl,'color','k','linestyle','--')
            grid
            ylabel('Channel');xlabel('Normalized Firing rate (peak)');
            
            %%
            %save results
            if ~disp_plot
                %suffixe=[info.align '_' num2str(binsize)];
                suffixe=[info.align '_' num2str(binsize) '_newbsl'];%new bsl [-50-wt -50]
                namesave=[save_path 'results_delay_' info.datafile(1:end-4) '_' suffixe];
                save(namesave, 'delay_bins' , 'peaks' );
            end
            
        
            %             %%
            %             if savepptx,
            %                 savetopptx(figtrialsboot,file,figtype,{info.datafile ;' trials bootstrapped'});
            %                 savetopptx(figelbs,file,figtype,{info.datafile ;' elbows'});
            %                 savetopptx(figdistr1,file,figtype,{info.datafile ; ' distribution of elbow 1'});
            %                 savetopptx(figdistr2,file,figtype,{info.datafile ;' distribution of elbow 2'});
            %                 savetopptx(figamps,file,figtype,{info.datafile ;' amplitude difference of elbows'});
            %             end

            %pause
            close all
        
        end
    end
end

%%
if savepptx
    %close .pptx
    newFile = exportToPPTX('saveandclose',filepptx)
end


