%function analysis_vmi

%function analysis_vmi
%   Analysis of data averaged over trials recorded with a laminar probe (LMA)
%   Compute VMI index and realign using CSDfeature
%
% see also compute_vmi compute_CSDfeature
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 10/14/2016 last modified 01/22/2017
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
figtype='epsc';%'png'
save_path='C:\Users\Corentin\Work\NeuroPITT\Publications\Figures\';


%alignement
%alignlist={'no' 'targ' 'go' 'sacc'};
%alignlist={'targ' 'sacc'};
alignlist={'targ_pburst_ch' 'sacc'};

%window of analysis
wind_targ=[-50 350];%[-350 600];%[-10 340];
wind_sacc=[-150 250];%[-550 400];%[-100 250];

%max_csdplot
max_csdplot=8.5828e+05;

%sigma FR
sigma_FR=6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get data
datalist=load_data_gandhilab(data_path);

%colorlist
colorlist=get_colorlist;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create list of vmis
dlist=get_dlist

data=[];info=[];
vmis_list=[];
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
    %%targets index
    %targs_ind=get_targsindex(targslist,info);
    
    %target tuning (after compute_tuning)
    targ_tuning=data(1).offline.targ_tuning;
    info.targ_tuning=targ_tuning;
    
    %%%%%%%%%%%%
    %VMI spk (after compute_vmi)
    %%vmis=data(1).offline.vmi.vmis_mean;
    %%vmis=data(1).offline.vmi.vmis_mean_bsl;
    %%vmis=data(1).offline.vmi.vmis_peak;
    %%vmis=data(1).offline.vmi.vmis_peak_bsl;
    
    %%vmis=data(1).offline.vmi.vmis_mean_bsignif;
    %vmis=data(1).offline.vmi.vmis_mean_bsl_bsignif;
    %%vmis=data(1).offline.vmi.vmis_peak_bsignif;
    %vmis=data(1).offline.vmi.vmis_peak_bsl_bsignif;
    
    %constant windows 100ms
    %vmis=data(1).offline.vmi.vmis_mean_bsl_100_bsignif;
    %vmis=data(1).offline.vmi.vmis_peak_bsl_100_bsignif;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %wt 100 ws 50 bslbefore
    %SELECTED
    %vmis=data(1).offline.vmi.vmis_mean_bslbefore_100_50_bsignif;    
    %vmis=data(1).offline.vmi.vmis_peak_bslbefore_100_50_bsignif;
    
    %vmis=data(1).offline.vmi.vmis_mean_bslbefore_2_100_50_bsignif;    
    %vmis=data(1).offline.vmi.vmis_mean_bslbefore_2_100_50_bsignif;  
    
    %wt 100 ws 25 bslbefore (to test effective part of the movemment burst)
    %vmis=data(1).offline.vmi.vmis_mean_bslbefore_100_25_bsignif;    
  
    
    %targ aligned on pburst
    vmis=data(1).offline.vmi.pburst_vmis_mean_bslbefore_100_50_bsignif;    
    %vmis=data(1).offline.vmi.pburst_vmis_peak_bslbefore_100_50_bsignif;
    
    
    %targ aligned on pburst + classification
    %vmis=data(1).offline.vmi.vispburst_vmis_mean_bslbefore_100_50_bsignif;    
    %vmis=data(1).offline.vmi.vmpburst_vmis_mean_bslbefore_100_50_bsignif;    
    %vmis=data(1).offline.vmi.movpburst_vmis_mean_bslbefore_100_50_bsignif;    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %wt 50 ws 50 bslbefore
    %vmis=data(1).offline.vmi.vmis_mean_bslbefore_50_50_bsignif;
    %vmis=data(1).offline.vmi.vmis_peak_bslbefore_50_50_bsignif;
    
    %wt 100 ws 50 bsl
    %vmis=data(1).offline.vmi.vmis_mean_bsl_100_50_bsignif;
    %vmis=data(1).offline.vmi.vmis_peak_bsl_100_50_bsignif;
    
    %wt 50 ws 50 bsl
    %vmis=data(1).offline.vmi.vmis_mean_bsl_50_50_bsignif;
    %vmis=data(1).offline.vmi.vmis_peak_bsl_50_50_bsignif;
    
    %Dprimes
    %wt 100 ws 50 bslbefore
    %vmis=data(1).offline.vmi.dprimes_mean_bslbefore_100_50_bsignif;
    %vmis=data(1).offline.vmi.dprimes_peak_bslbefore_100_50_bsignif;
    
    %wt 50 ws 50 bslbefore
    %vmis=data(1).offline.vmi.dprimes_mean_bslbefore_50_50_bsignif;
    %vmis=data(1).offline.vmi.dprimes_peak_bslbefore_50_50_bsignif;
    
    %wt 100 ws 50 bsl
    %vmis=data(1).offline.vmi.dprimes_mean_bsl_100_50_bsignif;
    %vmis=data(1).offline.vmi.dprimes_peak_bsl_100_50_bsignif;
    
    %wt 50 ws 50 bsl
    %vmis=data(1).offline.vmi.dprimes_mean_bsl_50_50_bsignif;
    %vmis=data(1).offline.vmi.dprimes_peak_bsl_50_50_bsignif;
    
    
    %LFP
    %SELECTED
    %wt 100 ws 50 bslbefore
    %vmis=data(1).offline.vmi.vmislfp_mean_bslbefore_100_50;
    %vmis=data(1).offline.vmi.vmislfp_mean_negbursts_bslbefore_100_50;
    %vmis=data(1).offline.vmi.vmislfp_mean_posbursts_bslbefore_100_50;
    %vmis=data(1).offline.vmi.vmispeaklfp_mean_bslbefore_100_50;
    
    %vmis=data(1).offline.vmi.vmislfp_peak_bslbefore_100_50_bsignif;
    %vmis=data(1).offline.vmi.vmislfp_absmean_bslbefore_100_50_bsignif;
    %vmis=data(1).offline.vmi.vmislfp_abspeak_bslbefore_100_50_bsignif;
    %vmis=data(1).offline.vmi.vmislfp_minmaxpeak_bslbefore_100_50_bsignif;
    %vmis=data(1).offline.vmi.vmislfp_min_bslbefore_100_50_bsignif;
    
    %wt 100 ws 25 bslbefore (to test effective part of the movemment burst)
    %vmis=data(1).offline.vmi.vmislfp_mean_bslbefore_100_25;
    
    %%%%%%%%%%%%
    %bursts
    %NOTE plots are not adapted for the bursts use axis tight and run for
    %both alignment separately
    %vmis=squeeze(data(1).offline.vmi.vmis_mean_bsl_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmis_mean_bsl_bursts(2,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmis_peak_bsl_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmis_peak_bsl_bursts(2,:,:));
    
    %constant windows 100ms
    %vmis=squeeze(data(1).offline.vmi.vmis_mean_bsl_100_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmis_mean_bsl_100_bursts(2,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmis_peak_bsl_100_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmis_peak_bsl_100_bursts(2,:,:));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %wt 100 ws 50 bslbefore
    %SELECTED
    %vmis=squeeze(data(1).offline.vmi.vmis_mean_bslbefore_100_50_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmis_mean_bslbefore_100_50_bursts(2,:,:));
    
    %vmis=squeeze(data(1).offline.vmi.vmis_peak_bslbefore_100_50_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmis_peak_bslbefore_100_50_bursts(2,:,:));
    
    %targ aligned on pburst
    %vmis=squeeze(data(1).offline.vmi.pburst_vmis_mean_bslbefore_100_50_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.pburst_vmis_mean_bslbefore_100_50_bursts(2,:,:));

    %vmis=squeeze(data(1).offline.vmi.pburst_vmis_peak_bslbefore_100_50_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.pburst_vmis_peak_bslbefore_100_50_bursts(2,:,:));

    
    %targ aligned on pburst + classification
    %vmis=squeeze(data(1).offline.vmi.vispburst_vmis_mean_bslbefore_100_50_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.vispburst_vmis_mean_bslbefore_100_50_bursts(2,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmpburst_vmis_mean_bslbefore_100_50_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmpburst_vmis_mean_bslbefore_100_50_bursts(2,:,:));
    %vmis=squeeze(data(1).offline.vmi.movpburst_vmis_mean_bslbefore_100_50_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.movpburst_vmis_mean_bslbefore_100_50_bursts(2,:,:));
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %wt 50 ws 50 bslbefore
    %vmis=squeeze(data(1).offline.vmi.vmis_mean_bslbefore_50_50_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmis_mean_bslbefore_50_50_bursts(2,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmis_peak_bslbefore_50_50_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmis_peak_bslbefore_50_50_bursts(2,:,:));
    
    %wt 100 ws 50 bsl
    %vmis=squeeze(data(1).offline.vmi.vmis_mean_bsl_100_50_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmis_mean_bsl_100_50_bursts(2,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmis_peak_bsl_100_50_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmis_peak_bsl_100_50_bursts(2,:,:));
    
    %wt 50 ws 50 bsl
    %vmis=squeeze(data(1).offline.vmi.vmis_mean_bsl_50_50_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmis_mean_bsl_50_50_bursts(2,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmis_peak_bsl_50_50_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmis_peak_bsl_50_50_bursts(2,:,:));
    
    %LFP
    %wt 100 ws 50 bslbefore
    %SELECTED
    %vmis=squeeze(data(1).offline.vmi.vmislfp_mean_bslbefore_100_50_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmislfp_mean_bslbefore_100_50_bursts(2,:,:));
    
    %vmis=squeeze(data(1).offline.vmi.vmislfp_absmean_bslbefore_100_50_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmislfp_absmean_bslbefore_100_50_bursts(2,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmislfp_abspeak_bslbefore_100_50_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmislfp_abspeak_bslbefore_100_50_bursts(2,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmislfp_minmaxpeak_bslbefore_100_50_bursts(1,:,:));
    %vmis=squeeze(data(1).offline.vmi.vmislfp_minmaxpeak_bslbefore_100_50_bursts(2,:,:));
    
    
    %     %%%%%%%%%%%%
    %     %VMI lfp
    %     %vmis=data(1).offline.vmislfp_bsl_sgfc;
    %     vmis=data(1).offline.vmislfp_bsl_noabs_sgfc;
    %     %vmis=data(1).offline.vmislfp_peak_bsl_sgfc;
    
    
    %CSD features (after compute_CSDfeature)
    info.csdfeat_avg_targ=data(1).offline.csdfeat_avg_targ;
    info.csdfeat_avg_sacc=data(1).offline.csdfeat_avg_sacc;
    
    %list of vmis
    vmis_list(dd).vmis=vmis(:,:);
    vmis_list(dd).info=info;
    
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure
%figvmicsd=figure('Position',[1 100 scrsz(3)-100 scrsz(4)-200]);
figvmicsd=figure;

%color
color_vmis=1;color_conf='blue';
%color_vmis=2;color_conf='red';
%color_vmis=3;color_conf='green';


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %plot vmis without re-alignment
% info.align=alignlist{1};
% %hdlfig=subplot(2,3,2);hold on;
% vmis_alltuning=[];
% for dd=1:size(vmis_list,2),
%     vmis=vmis_list(dd).vmis;
%     info=vmis_list(dd).info;
%     
%     %plot_vmis(vmis,info.targ_tuning,'o',dd,1,info,hdlfig,[]);
%     %plot_vmis(vmis,info.targ_tuning,'-',dd,1,info,hdlfig,[]);
%     
%     %vmis_avg
%     if isempty(vmis_alltuning),
%         vmis_alltuning=vmis(info.targ_tuning,:);
%     else
%         vmis_alltuning=[vmis_alltuning ; vmis(info.targ_tuning,:)];
%     end
%     
%     %line at ch_ref
%     %    hl=line([-1 1] ,[ch_ref ch_ref]);
%     %    set(hl,'Color',colorlist(1,:),'LineStyle','--','Linewidth',1);
%     %pause
% end
% vmis_avg=nanmean(vmis_alltuning,1);
% plot_vmis(vmis_avg,1,'-',color_vmis,3,info,hdlfig,[]);
% %plot_dprimes(vmis_avg,1,'-',color_vmis,3,info,hdlfig,[]);
% 
% % %standard deviation
% % vmis_std=nanstd(vmis_alltuning,1);
% % fill([vmis_avg-vmis_std fliplr(vmis_avg+vmis_std)],[1:info.nchannels info.nchannels:-1:1], 1,'facecolor','blue','edgecolor','none','facealpha', 0.3);
% 
% % %%
% % % %95% confidence interval
% % vmis_ci=[];
% % for ch=1:info.nchannels,
% %     aux=(vmis_alltuning(find(~isnan(vmis_alltuning(:,ch))),ch));
% %     if numel(aux)<=1,
% %         vmis_ci(ch,:)=[vmis_avg(ch) ; vmis_avg(ch)];
% %     else
% %         vmis_ci(ch,:) = bootci(2000,{@mean,aux},'type','per');
% %     end
% % end
% % fill([vmis_ci(:,1)' fliplr(vmis_ci(:,2)')],[1:info.nchannels info.nchannels:-1:1], 1,'facecolor',color_conf,'edgecolor','none','facealpha', 0.3);
% 
% %axis
% lims=[info.chmap(1) info.chmap(end)];
% %axis([-1 1 lims(1) lims(2)]);
% axis([-1 1 1 length(info.chmap)]);
% %axis([-5 5 1 length(info.chmap)]);%dprimes
% 
% %variance
% vmis_var=nanvar(vmis_alltuning,1);
% mvar=nanmean(vmis_var(lims(1):lims(2)));
% title({[info.datafile(20:end-4) ' (' num2str(dd) ' vmis)'] ; ['variance=' num2str(mvar)]});
% 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot with re-alignment
%vmis_listr=[];
for realign=2%1:3 %assuming alignlist={'targ' 'sacc'}
    switch realign
        case 1
            al=1;
            title_al='targ (csdfeat 1)';
        case 2
            al=1;
            title_al='targ (csdfeat 2)';
        case 3
            al=2;
            title_al='sacc';
    end
    info.align=alignlist{1};
    
    vmis_r_alltuning=[];
    ch_ref_list=[];
    %hdlfig=subplot(2,3,3+realign);hold on;%figure;hold on;
    hdlfig=subplot(1,1,1);hold on;
    for dd=1:size(vmis_list,2),
        vmis=vmis_list(dd).vmis;
        
        %
        
        info=vmis_list(dd).info;
        %dref
        switch realign
            case 1
                dref=info.csdfeat_avg_targ(1);
            case 2
                dref=info.csdfeat_avg_targ(2);
            case 3
                dref=info.csdfeat_avg_sacc;
        end
        
        
        %realign vmis
        if ~isempty(dref)
            [vmis_r info_r ch_ref dref_conv]=get_vmis_aligndepth(vmis,dref,info);
            
            %%list of all ch_ref
            %ch_ref_list(dd)=dref_conv;
            
            %vmis_r(:,1+info.nchannels:info.nchannels+ch_ref)=nan;
            %vmis_r(:,1:ch_ref-1)=nan;
            %vmis_r(:,-1+ch_ref:end)=nan;
            
            %normalization to plot bursts
            %vmis_r(info_r.targ_tuning,:)=vmis_r(info_r.targ_tuning,:)/max(abs(vmis_r(info_r.targ_tuning,:)));
            
            %%plot realigned vmis
            %plot_vmis(vmis_r,info_r.targ_tuning,'o',dd,1,info_r,hdlfig,[]);
            %plot_vmis(vmis_r,info_r.targ_tuning,'-',dd,1,info_r,hdlfig,[]);
           
            %vmis_listr
            %vmis_r_list(al,dd).vmis=vmis_r;
            %vmis_r_list(al,dd).info=info_r;
            
            %vmis_r_avg
            if isempty(vmis_r_alltuning),
                vmis_r_alltuning=vmis_r(info_r.targ_tuning,:);
            else
                vmis_r_alltuning=[vmis_r_alltuning ; vmis_r(info_r.targ_tuning,:)];
            end
            
        else
            display(['dref is void for file ' num2str(dlist(dd)) ' and align '  info.align '!'])
        end
        
        %axis tight
        %axis([ -137.5705  205.4758    8.7354   38.7013])
        %dd
        %pause
    end
    
    
    %%%%%%%%%%%%%%%%%%
    %spk
    %plot vmis_r_avg with outliers
    vmis_r_avg=nanmean(vmis_r_alltuning,1);
    vmis_r_var=nanvar(vmis_r_alltuning,1);
    vmis_r_std=nanstd(vmis_r_alltuning,1);
    %moving mean or median.
    %NOTE: See Matlab commands like: movmedian, movmean,  tsmovavg, and medfilt1
    %vmis_r_avg=nanmedian(vmis_r_alltuning,1);
    %vmis_r_avgmov=movmean(vmis_r_alltuning,3,1,'omitnan');
    %vmis_r_avgmovmed=movmedian(vmis_r_alltuning,3,1,'omitnan');
    %plot vmis
    plot_vmis(vmis_r_avg,1,'-',color_vmis,3,info_r,hdlfig,[]);
    
    
%     %%%%%%%%%%%%%%%%%%
%     %LFP
%     %detect outliers
%     vmis_r_nout={};outliers={};
%     for o=1:size(vmis_r_alltuning,2)
%         aux=vmis_r_alltuning(find(vmis_r_alltuning(:,o)<200),o);
%         [vmis_r_nout{o} outliers{o}] = findoutliers(aux);
%         
%         %thresh=6;
%         %aux=vmis_r_alltuning(find(abs(vmis_r_alltuning(:,o))<=thresh),o);
%         %outliers{o}=vmis_r_alltuning(find(abs(vmis_r_alltuning(:,o))>thresh),o);
%         %vmis_r_nout{o}=aux;
%         
%         vmis_r_avg(o)=nanmean(vmis_r_nout{o});
%         
%         %plot outliers
%         out=outliers{o};
%         if ~isempty(out)
%         for oi=1:length(out)
%             plot(out(oi),o,'o','linewidth',2,'MarkerFaceColor','b');
%         end
%         end
%         
%     end
%     %plot vmis
%     plot_vmis(vmis_r_avg,1,'-',color_vmis,3,info_r,hdlfig,[]);
%     %plot_dprimes(vmis_r_avg,1,'-',color_vmis,3,info_r,hdlfig,[]);
%     
    
    %%%%%%%%%%%%%%%%%%
    %line at ch_ref
    hl=line([-1 1] ,[ch_ref ch_ref]);
    
    %hl=line([-5 5] ,[ch_ref ch_ref]); %dprimes
    %set(hl,'Color',colorlist(1,:),'LineStyle','--','Linewidth',1);
    
    %     %standard deviation
    %     vmis_r_std=nanstd(vmis_r_alltuning,1);
    %     chs_r=find(~isnan(vmis_r_avg));
    %     vmis_r_avg_plot=vmis_r_avg(min(chs_r):max(chs_r));
    %     vmis_r_std_plot=vmis_r_std(min(chs_r):max(chs_r));
    %     fill([vmis_r_avg_plot-vmis_r_std_plot fliplr(vmis_r_avg_plot+vmis_r_std_plot)],[min(chs_r):max(chs_r) max(chs_r):-1:min(chs_r)], 1,'facecolor','blue','edgecolor','none','facealpha', 0.3);
    
    
        %%%%%%%%%%%%%%%%%%%
        % 95% confidence interval
        %case of missing value in chs_r (because of not enough data point)
%         chs_r=find(~isnan(vmis_r_avg));
%         [vmiss imiss]=find(chs_r(2:end)-chs_r(1:end-1)>1);
%     
       
        %find channel range
        chs_r=find(~isnan(vmis_r_avg));
        [vmiss imiss]=find(chs_r(2:end)-chs_r(1:end-1)>1);
        
        %consider only the last consecutive channels
        %     %if ~isempty(imiss),min_ch=chs_r(max(imiss)+1);else min_ch=chs_r(1);end
        %min_ch=chs_r(1);
        %     if ~isempty(imiss),max_ch=chs_r(imiss-1);else max_ch=max(chs_r);end
        %max_ch=chs_r(imiss(1)-1);
        
        min_ch=chs_r(1);
        max_ch=max(chs_r);
        if length(imiss)==1 & imiss(1)<10
            min_ch=chs_r(imiss(1)+1);
        elseif length(imiss)==1 & imiss(1)>=10
            max_ch=chs_r(imiss(1)-1);
        elseif length(imiss)==2 & imiss(1)<10 & imiss(2)>10
            min_ch=chs_r(imiss(1)+1);
            max_ch=chs_r(imiss(2)-1);
        elseif length(imiss)==2 & imiss(1)<10 & imiss(2)<10
            min_ch=chs_r(imiss(2)+1);
        elseif length(imiss)==3 & imiss(2)<10 & imiss(3)>10
            min_ch=chs_r(imiss(2)+1);
            max_ch=chs_r(imiss(3)-1);
        end

        
        
        ind=0;vmis_r_ci=[];
        for ch=min_ch:max_ch,
            ind=ind+1;
    
            %alltuning data
            aux=(vmis_r_alltuning(find(~isnan(vmis_r_alltuning(:,ch))),ch));
            %after removing outliers (for LFP)
            %aux=vmis_r_nout{ch};
    
            if numel(aux)<=1,
                vmis_r_ci(ind,:)=[vmis_r_avg(ch) ; vmis_r_avg(ch)];
            else
                vmis_r_ci(ind,:) = bootci(2000,{@mean,aux},'type','per');
            end
        end
        fill([vmis_r_ci(:,1)' fliplr(vmis_r_ci(:,2)')],[min_ch:max_ch max_ch:-1:min_ch], 1,'facecolor',color_conf,'edgecolor','none','facealpha', 0.3);
    
    
    %%%%%%%%%%%%%%%%
    %variance (measure of quality of alignment)
    %vmis_r_var=nanvar(vmis_r_alltuning,1);
    %mvar(realign)=nanmean(vmis_r_var(lims(1):lims(2)));
    title({title_al})% ; ['variance=' num2str(mvar(realign))]})
    
    %%
    %%%%%%%%%%%%%%%
    %axis
    %lims=findlimits(vmis_r_avg);
    switch realign
        case 1
            lims=[6 24];
        case 2
            lims=[16 32];%[16 33];%[12 35];
            
        case 3
            lims=[12 34];
    end
    axis([-1 1 lims(1) lims(2)]);
    yyaxis left
    set(gca,'Xtick',[-1:0.2:1],'Ytick',[lims(1):2:lims(2)+1],'Yticklabel',[-8:2:10])
    axis([-1 1 lims(1) lims(2)]);
    yyaxis right
    set(gca,'Ytick',[lims(1):2:lims(2)+1],'Yticklabel',[1.2:-0.3:-1.5])
    %repeat?? to display second axis
    axis([-1 1 lims(1) lims(2)]);
    yyaxis right
    set(gca,'Ytick',[lims(1):2:lims(2)+1],'Yticklabel',[1.2:-0.3:-1.5])
    ylabel('Depth (mm)')
    
    %%
    %set axis for bursts avg
    lims=[16 32]
    yyaxis left
    axis([0 180 lims(1) lims(2)])
    set(gca,'Xtick',[0:20:180])
    yyaxis right
    axis([0 180 lims(1) lims(2)])
    set(gca,'Ytick',[lims(1):2:lims(2)],'Yticklabel',[1.2:-0.3:-1.2])
    
    %%
    
    %     %%%%%%%%%%%%%%%%%%
    %     %dprimes
    %     axis([-5 5 lims(1) lims(2)]);
    %     yyaxis left
    %     set(gca,'Xtick',[-5:1:5],'Ytick',[lims(1):2:lims(2)+1],'Yticklabel',[-8:2:10])
    %     axis([-5 5 lims(1) lims(2)]);
    %     yyaxis right
    %     set(gca,'Ytick',[lims(1):2:lims(2)+1],'Yticklabel',[1.2:-0.3:-1.5])
    %     %repeat?? to display second axis
    %     axis([-5 5 lims(1) lims(2)]);
    %     yyaxis right
    %     set(gca,'Ytick',[lims(1):2:lims(2)+1],'Yticklabel',[1.2:-0.3:-1.5])
    %     ylabel('Depth (mm)')
    
    
    
    %     %%%%%%%%%%%%%%%%%%
    %     %plot all dref
    %     figure;hold on;
    %     plot(1:size(vmis_list,2),ch_ref_list,'o-');
    %     xlabel('Session');ylabel('CSD reference channel')
    %     axis([0 30 1 16])
    
    
%         %%%%%%%%%%%%%%%%%%
%         %save vmis for stats
%         %spk
%         %file=['Stats/vmis_' info.datafile(20:31)];
%         %vmis=vmis_r_alltuning;
%     
%         %LFP
%         file=['Stats/vmislfp_' info.datafile(20:31)];
%         vmis=vmis_r_nout;
%     
%         %option
%         %file=[file '_matched'];
%     
%         file
%         save(file, 'vmis' );


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%stats (mean , CI ....)
stats=[];
stats(:,1)=vmis_r_avg(min_ch:max_ch);
stats(:,2:3)=vmis_r_ci(:,1:2); 

min_ch

display(['Across depths mean=' num2str(mean(stats(:,1)))]);
display(['Across depths meanCI inf=' num2str(mean(stats(:,2)))]);
display(['Across depths meanCI sup=' num2str(mean(stats(:,3)))]);

[valstats istats]=max(stats(:,1));
display(['Across channel of max=' num2str(istats)]);
display(['Across depths max=' num2str(valstats)]);
display(['Across depths maxCI inf=' num2str(stats(istats,2))]);
display(['Across depths maxCI sup=' num2str(stats(istats,3))]);

[valstats istats]=min(stats(:,1));
display(['Across channel of min=' num2str(istats)]);
display(['Across depths min=' num2str(valstats)]);
display(['Across depths minCI inf=' num2str(stats(istats,2))]);
display(['Across depths minCI sup=' num2str(stats(istats,3))]);

    
end

% %%
% %%%%%%%%%%%%%%%%%%
% %save figs
% if savefigs
%     %saveas(hdlfig,[save_path info.datafile(1:2) '_VMI_align_' num2str(realign) .' figtype],figtype);
%     print([save_path info.datafile(1:2) '_VMIs_align_' num2str(realign) '.' figtype],'-depsc','-painters','-loose',gcf)
% end;







%%
%MISC
% %jitter analysis of sensitivity of dref estimation
% compute mean/var of pair-wise rmse for each alignment: before alignment,after alignment, with jitter

% nbjit=100;mvar=[];
% rdlist=[1:16];
% for jit=1:nbjit
%     for rd=rdlist
%
%         %figvmicsd=figure('Position',[1 100 scrsz(3)-100 scrsz(4)-200]);
%
%         %without re-alignment
%         al=1;
%         info.align=alignlist{al};
%         %hdlfig=subplot(2,3,2);hold on;
%         for dd=1:size(vmis_list,2),
%             vmis=vmis_list(al,dd).vmis;
%             info=vmis_list(al,dd).info;
%          %   plot_vmis(vmis,info.targ_tuning,dd,1,info,hdlfig,[]);
%         end
%         %title({[info.datafile(1:18) ' (' num2str(dd) ' files) rd:' num2str(rd)]});
%
%
%         %with re-alignment
%         vmis_listr=[];
%         for realign=1:3 %assuming alignlist={'targ' 'sacc'}
%             switch realign
%                 case 1
%                     al=1;
%                     title_al='targ (csdfeat 1)';
%                 case 2
%                     al=1;
%                     title_al='targ (csdfeat 2)';
%                 case 3
%                     al=2;
%                     title_al='sacc';
%             end
%             info.align=alignlist{al};
%
%             vmis_r_alltuning=[];
%             %hdlfig=subplot(2,3,3+realign);hold on;
%             for dd=1:size(vmis_list,2),
%                 vmis=vmis_list(al,dd).vmis;
%                 info=vmis_list(al,dd).info;
%                 %dref
%                 switch realign
%                     case 1
%                         dref=info.csdfeat_avg_targ(1)+ rd-1;%( rd-1+round(randn*5) );
%                     case 2
%                         dref=info.csdfeat_avg_targ(2)+rd-1;%( rd-1+round(randn*5) );
%                     case 3
%                         dref=info.csdfeat_avg_sacc+rd-1;%( rd-1+round(randn*5) );
%                 end
%
%                 %realign vmis
%                 if ~isempty(dref)
%                     [vmis_r info_r ch_ref]=get_vmis_aligndepth(vmis,dref,info);
%
%                     %plot realigned vmis
%                     %plot_vmis(vmis_r,info_r.targ_tuning,dd,1,info_r,hdlfig,[]);
%
%                     %vmis_listr
%                     %vmis_r_list(al,dd).vmis=vmis_r;
%                     %vmis_r_list(al,dd).info=info_r;
%
%                     %vmis_r_avg
%                     if isempty(vmis_r_alltuning),
%                         vmis_r_alltuning=vmis_r(info_r.targ_tuning,:);
%                     else
%                         vmis_r_alltuning=[vmis_r_alltuning ; vmis_r(info_r.targ_tuning,:)];
%                     end
%
%                 else
%                     display(['dref is void for file ' num2str(dlist(dd)) ' and align '  info.align '!'])
%                 end
%             end
%
%             %plot vmis_r_avg
%             vmis_r_avg=nanmean(vmis_r_alltuning,1);
%             %lims=findlimits(vmis_r_avg);
%             %plot_vmis(vmis_r_avg,1,dd,3,info_r,hdlfig,[]);
%             %hl=line([-1 1] ,[ch_ref ch_ref]);
%             %set(hl,'Color',colorlist(1,:),'LineStyle','-','Linewidth',1);
%             %axis([-1 1 lims(1) lims(2)]);
%
%             %variance
%             vmis_r_var=nanvar(vmis_r_alltuning,1);
%             mvar(jit,rd-min(rdlist)+1,realign)=mean(vmis_r_var(lims(1):lims(2)));
%
%             %title({title_al ; ['variance=' num2str(mvar(realign))]})
%
%
%         end
%         %pause
%         rd
%         %close(figvmicsd)
%     end
%     jit
% end
%
% %plot results of jitter analysis on variance
% figure;hold on;
% for al=1:3
%     plot([0:size(mvar,2)-1],mean(squeeze(mvar(:,:,al)),1),'linewidth',2,'color',colorlist(al,:))
% end
% xlabel('Jitter (nb channels)');ylabel('Variance of realignment')
% legend({'targ (csdfeat 1)','targ (csdfeat 2)','sacc'},'Location','NorthWest')
% title({[info.datafile(1:18) ' (' num2str(dd) ' files) nb repetitions:' num2str(nbjit)]});
%
%
%
