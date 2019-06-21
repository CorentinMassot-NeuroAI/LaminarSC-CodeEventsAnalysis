%function analysis_delay_avg

%function analysis_delay    _avg
%   Analyze activity during delay period based on average activity recorded with a
%   laminar probe (LMA)
%
% see also compute_delay_avg
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 11/29/2016 last modified 11/29/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TO DO
% add actvity across depth after go cue and before build-up onset.
% add comparison with magnitdue of burst onset / peak motor activity
% add remaining delay for later go cue



%set paths
[root_path data_path save_path]=set_paths;
save_path=[root_path 'Results\Results_SC_delay\'];


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
%alignlist={'sacc'};
alignlist={'targ_pburst_ch' };


%windows of analysis (do not change)
wind_sacc=[-400 100];%[-100 200]

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

alldelaybins={};alldelaybinsn={};alldelaybinsn2={};
nbins_list=[];
data=[];info=[];dd=0;
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
    
    %     %sacc bursts significance
    %     thresh_ratios=0.15;%0.15
    %     thresh_surprises=4;
    %     ratios_sacc=data(1).offline.sacc_pburst_ratio(info.targ_tuning,:)>thresh_ratios;
    %     surprises_sacc=data(1).offline.sacc_pburst_msurprises(info.targ_tuning,:)>thresh_surprises;
    %     bsignif_sacc=data(1).offline.sacc_bsignif;
    %     bthresh_sacc=data(1).offline.sacc_bthresh_trials';
    %     %burst_bsignif=(ratios_sacc & surprises_sacc & bsignif_sacc & bthresh_sacc);
    %     burst_bsignif=(bsignif_sacc & bthresh_sacc);
    %
    %targ bursts significance
    thresh_ratios=0.15;%0.15
    thresh_surprises=4;
    ratios_targ=data(1).offline.targ_pburst_ratio(info.targ_tuning,:)>thresh_ratios;
    surprises_targ=data(1).offline.targ_pburst_msurprises(info.targ_tuning,:)>thresh_surprises;
    bsignif_targ=data(1).offline.targ_pburstch_bsignif;
    bthresh_targ=data(1).offline.targ_pburstch_bthresh_trials';
    
    %burst_bsignif=(ratios_targ & surprises_targ & bsignif_targ & bthresh_targ)
    burst_bsignif=(bsignif_targ & bthresh_targ);
    
    
    burst_bsignif=double(burst_bsignif);
    burst_bsignif(find(burst_bsignif==0))=nan;
    
    %ntrials
    info.ntrials=0;
    
    
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load results onsetbuildup for bootstrapped trials
    binsize=50;
    %suffixe=['targ_pburst_ch_' num2str(binsize)];
    suffixe=[info.align '_' num2str(binsize) '_newbsl'];%new bsl [-50-wt -50]
    nameload=[save_path 'results_delay_' info.datafile(1:end-4) '_' suffixe]
    load(nameload)
    
    %'delay_bins' , 'peaks'
    
    
    %%
    %%%%%%%%%%%%%
    %normalize delay bins
    delay_bins_n=[];delay_bins_n2=[];
    nbins=size(delay_bins,2);
    nbins_list(dd)=nbins;
    for bi=1:nbins
        delay_bins_n(:,bi)=delay_bins(:,bi)./peaks(:,2);
        delay_bins_n2(:,bi)=(peaks(:,2)-delay_bins(:,bi))./peaks(:,2);
    end
    
    
%     %%
%     %%%%%%%%%%%%%
%     %plot delay bins
%     figdelaybins=figure('Position',[scrsz(3)/3 200 scrsz(3)/1.8 scrsz(4)-400]);hold on;
%     legend_bins={};
%     colorlist2 = colormap(jet(nbins));%hot hsv copper
%     for bi=1:nbins
%         %without normalization
%         figure(figdelaybins)
%         subplot(1,2,1);hold on;
%         plot(delay_bins(:,bi),1:info.nchannels,'color',colorlist2(bi,:))
%         
%         %witht peak normalization
%         figure(figdelaybins)
%         subplot(1,2,2);hold on;
%         plot(delay_bins_n(:,bi),1:info.nchannels,'color',colorlist2(bi,:))
%         
%         legend_bins{bi}=num2str(bi);
%     end
%     
%     
%     figure(figdelaybins)
%     subplot(1,2,1);
%     axis([-20 max(peaks(:,2)) 1 info.nchannels])
%     hdl=line([0 0],[1 info.nchannels]);
%     set(hdl,'color','k','linestyle','--')
%     grid
%     legend(legend_bins)
%     ylabel('Channel');xlabel('Firing rate (spk/s)');
%     title('Delay activity')
%     subplot(1,2,2);
%     axis([-0.3 1 1 info.nchannels])
%     hdl=line([0 0],[1 info.nchannels]);
%     set(hdl,'color','k','linestyle','--')
%     grid
%     ylabel('Channel');xlabel('Normalized Firing rate (peak)');
%     title(info.datafile)
%     
    
    %     %%
    %     %%%%%%%%%%%%%
    %     %save in powerpoint
    %     if savepptx,
    %         savetopptx(figtrialsboot,file,figtype,{info.datafile ;' elb1 and inflect1'});
    %         savetopptx(figtrialsboot3,file,figtype,{info.datafile ;' inflect3'});
    %         savetopptx(figelbs,file,figtype,{info.datafile ;' CI of onsets and slopes'});
    %         savetopptx(figtrialsboot2,file,figtype,{info.datafile ;' accum and burst (classification)'});
    %         %savetopptx(figciadiff,file,figtype,{info.datafile ;' Regression CI vs. slopes diff'});
    %     end
    
    
    
    %%
    %%%%%%%%%%%%%%%%
    %alignment of onset using CSD features (after compute_CSDfeature)
    info.csdfeat_avg_targ=data(1).offline.csdfeat_avg_targ;
    info.zs=data(1).offline.csdzs;
    dref=info.csdfeat_avg_targ(2);
    
    %delay_bins
    [delay_bins_r info_r ch_ref ~]=get_data_aligndepth(delay_bins,dref,info,[]);
    [delay_bins_n_r info_r ch_ref ~]=get_data_aligndepth(delay_bins_n,dref,info,[]);
    [delay_bins_n2_r info_r ch_ref ~]=get_data_aligndepth(delay_bins_n2,dref,info,[]);
    
    
    %remove non significant burst
    [burst_bsignif_r info_r ch_ref ~]=get_vmis_aligndepth(burst_bsignif,dref,info);
    delay_bins_r=delay_bins_r.*burst_bsignif_r';
    delay_bins_n_r=delay_bins_n_r.*burst_bsignif_r';
    delay_bins_n2_r=delay_bins_n2_r.*burst_bsignif_r';
    
    %list of onsets for all sessions
    alldelaybins{dd}=delay_bins_r;
    alldelaybinsn{dd}=delay_bins_n_r;
    alldelaybinsn2{dd}=delay_bins_n2_r;
    
    
    %pause
    close all
    
end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%min nbins
%nbins_min=max(nbins_list);
nbins_min=min(nbins_list);
%nbins_min=6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot delay_bins
fig_d=figure('Position',[scrsz(3)/4 100 scrsz(3)/2 scrsz(4)-200]);hold on;
colorlist2 = colormap(jet(nbins_min));%hot hsv copper
legend_bins={};
ibi=1;
for bi=1:nbins_min
    fig_all=plot_stats_depths_v('delay_magnitude',alldelaybins,[ ],bi,colorlist2(bi,:),[-20 180],fig_d,info,datalist,dlist,fig_d,[])
   
%     legend_bins{bi}=num2str(bi);
    legend_bins{ibi}=num2str(bi);
    legend_bins{ibi+1}=num2str(bi);
    legend_bins{ibi+2}=num2str(bi);
    ibi=ibi+3;
    %pause
end
grid;
axis tight
axis([16 32 -20 180 -1 1])
title('delay')
legend(legend_bins,'Location','SouthEast')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot delay_bins_n (normalization by peaks)
fig_dn=figure('Position',[scrsz(3)/4 100 scrsz(3)/2 scrsz(4)-200]);hold on;
colorlist2 = colormap(jet(nbins_min));%hot hsv copper
legend_bins={};
ibi=1;
for bi=1:nbins_min
    fig_all=plot_stats_depths_v('delay_magnitude',alldelaybinsn,[ ],bi,colorlist2(bi,:),[-0.5 1],fig_dn,info,datalist,dlist,fig_dn,[])
    
%     legend_bins{bi}=num2str(bi);
    legend_bins{ibi}=num2str(bi);
    legend_bins{ibi+1}=num2str(bi);
    legend_bins{ibi+2}=num2str(bi);
    ibi=ibi+3;
end
grid;
axis tight
axis([16 32 -0.1 1 -1 1])
title('Delay normalized by peak')
legend(legend_bins,'Location','SouthEast')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot delay_bins_n2
fig_dn=figure('Position',[scrsz(3)/4 100 scrsz(3)/2 scrsz(4)-200]);hold on;
colorlist2 = colormap(jet(nbins));%hot hsv copper
legend_bins={};
ibi=1;
for bi=1:nbins_min
    fig_all=plot_stats_depths_v('delay_magnitude',alldelaybinsn2,[ ],bi,colorlist2(bi,:),[-0.5 1],fig_dn,info,datalist,dlist,fig_dn,[])
    
%     legend_bins{bi}=num2str(bi);
    legend_bins{ibi}=num2str(bi);
    legend_bins{ibi+1}=num2str(bi);
    legend_bins{ibi+2}=num2str(bi);
    ibi=ibi+3;
end
grid;
axis tight
title('Delay: diff normalized by peak')
legend(legend_bins,'Location','SouthEast')




%%
% %%%%%%%%%%%%%
% %save in powerpoint
% if savepptx,
%     savetopptx(fig_l,file,figtype,{info.datafile ;' Average latencies of elb1/inflect1'});
%     savetopptx(fig_m,file,figtype,{info.datafile ;' Average firing rates of elb1/inflect1'});
%     savetopptx(fig_class_l,file,figtype,{info.datafile ;' Average latencies of accum/burst'});
%     savetopptx(fig_class_m,file,figtype,{info.datafile ;' Average firing rates of accum/burst'});
%     %     savetopptx(fig_allpeak_l,file,figtype,{info.datafile ;' All sessions peak latencies'});
%     %     savetopptx(fig_allelb1_l,file,figtype,{info.datafile ;' All sessions elb1 latencies'});
%     %     savetopptx(fig_allinflect1_l,file,figtype,{info.datafile ;' All sessions inflect1 latencies '});
%     %     savetopptx(fig_allpeak_m,file,figtype,{info.datafile ;' All sessions peak firing rates'});
%     %     savetopptx(fig_allelb1_m,file,figtype,{info.datafile ;' All sessions elb1 firing rates'});
%     %     savetopptx(fig_allinflect1_m,file,figtype,{info.datafile ;' All sessions inflect1 firing rates'});
%
% end

%%
if savepptx
    %close .pptx
    newFile = exportToPPTX('saveandclose',filepptx)
end



