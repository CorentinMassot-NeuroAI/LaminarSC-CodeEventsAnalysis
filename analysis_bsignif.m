%function analysis_bsignif

%function analysis_bsignif
%   analyse significance of visual and movement bursts
%
% see also compute_vmi
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 06/22/2017 last modified 06/22/2017
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
%alignlist={'targ' 'sacc'};
%alignlist={'targ_pburst_ch' 'sacc'};

%algo for targ_pburst
algo='pburst';%Hanes et al.
%algo='rsburst';%Gourevitch et al.

%thresholds
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

%figure
%fig=figure;hold on;%('Position',[1 100 scrsz(3)-100 scrsz(4)-200]);

data=[];
info=[];
targ_bsignif_list=[];
sacc_bsignif_list=[];
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
    
    %target tuning (after compute_tuning)
    targ_tuning=data(1).offline.targ_tuning;
    
    
    %%%%%%%%%%%%%
    %bsignif based on targ alignment (see compute_bsignif)
    %targ_bsignif=data(1).offline.targ_bsignif;
    %sacc_bsignif=data(1).offline.sacc_bsignif;
    
    %targ_bsignif=data(1).offline.targ_bsignif & data(1).offline.targ_bthresh';
    %sacc_bsignif=data(1).offline.sacc_bsignif & data(1).offline.sacc_bthresh';
    
    %%%%%%%%%%%%%
    %bsignif based on targ_pburstch alignment (see compute_bsignif compute_pburst)
    %targ_pburst
    switch algo
        case 'pburst'
            %targ
            ratios_targ=data(1).offline.targ_pburst_ratio(targ_tuning,:)>thresh_ratios;
            surprises_targ=data(1).offline.targ_pburst_msurprises(targ_tuning,:)>thresh_surprises;
            bsignif_targ=data(1).offline.targ_pburstch_bsignif;
            bthresh_targ=data(1).offline.targ_pburstch_bthresh_trials';
            %bthresh_targ=data(1).offline.targ_pburstch_bthresh';
            
            %sacc
            ratios_sacc=data(1).offline.sacc_pburst_ratio(targ_tuning,:)>thresh_ratios;
            surprises_sacc=data(1).offline.sacc_pburst_msurprises(targ_tuning,:)>thresh_surprises;
            bsignif_sacc=data(1).offline.sacc_bsignif;
            bthresh_sacc=data(1).offline.sacc_bthresh_trials';
            %bthresh_sacc=data(1).offline.sacc_bthresh';
            
        case 'rsburst'
            %targ
            ratios_targ=data(1).offline.targ_rsburst_ratio(targ_tuning,:)>thresh_ratios;
            surprises_targ=data(1).offline.targ_rsburst_msurprises(targ_tuning,:)>thresh_surprises;
            bsignif_targ=data(1).offline.targ_rsburstch_bsignif;
            bthresh_targ=data(1).offline.targ_rsburstch_bthresh_trials';
            %bthresh_targ=data(1).offline.targ_rsburstch_bthresh_';
            
            %sacc
            ratios_sacc=data(1).offline.sacc_rsburst_ratio(targ_tuning,:)>thresh_ratios;
            surprises_sacc=data(1).offline.sacc_rsburst_msurprises(targ_tuning,:)>thresh_surprises;
            bsignif_sacc=data(1).offline.sacc_bsignif;
            bthresh_sacc=data(1).offline.sacc_bthresh_trials';
            %bthresh_sacc=data(1).offline.sacc_bthresh';
    end
    
    %targ_bsignif
    %targ_bsignif=(ratios_targ & surprises_targ & bsignif_targ & bthresh_targ);
    targ_bsignif=(bsignif_targ & bthresh_targ);
    %targ_bsignif=(ratios_targ & surprises_targ);
    %targ_bsignif=(bsignif_targ);
    
    %sacc_bsignif
    %sacc_bsignif=(ratios_sacc & surprises_sacc & bsignif_sacc & bthresh_sacc);
    sacc_bsignif=(bsignif_sacc & bthresh_sacc);
    %sacc_bsignif=(ratios_sacc & surprises_sacc);
    %sacc_bsignif=(bsignif_sacc);
    
    %%%%%%%%%%%%%%%%
    %classification
    vm_list=zeros(1,info.nchannels);vis_list=zeros(1,info.nchannels);mov_list=zeros(1,info.nchannels);
    for ch=1:info.nchannels
        if targ_bsignif(ch)==1 & sacc_bsignif(ch)==1
            vm_list(ch)=1;
        end
        
        if targ_bsignif(ch)==1 & sacc_bsignif(ch)==0
            vis_list(ch)=1;
        end
        
        if targ_bsignif(ch)==0 & sacc_bsignif(ch)==1
            mov_list(ch)=1;
        end
    end
    %update data
    %     data=update_data(0,1,0,data,data_path,info.datafile,['classif_' algo '_vm'],vm_list);
    %     data=update_data(0,1,0,data,data_path,info.datafile,['classif_' algo '_vis'],vis_list);
    %     data=update_data(0,1,0,data,data_path,info.datafile,['classif_' algo '_mov'],mov_list);
    display('NOT SAVED!')
    %     %update_data(1,0,0,data,data_path,info.datafile,[],[]);
    
    
    %%%%%%%%%%%%%%%%
    %alignment using CSD features (after compute_CSDfeature)
    info.csdfeat_avg_targ=data(1).offline.csdfeat_avg_targ;
    info.zs=data(1).offline.csdzs;
    dref=info.csdfeat_avg_targ(2);
    [targ_bsignif_r info_r ch_ref ~]=get_vmis_aligndepth(targ_bsignif,dref,info);
    [sacc_bsignif_r info_r ch_ref ~]=get_vmis_aligndepth(sacc_bsignif,dref,info);
    
    targ_bsignif_list=[targ_bsignif_list ; targ_bsignif_r];
    sacc_bsignif_list=[sacc_bsignif_list ; sacc_bsignif_r];
    
    %    classif.accum_only{dd} = ~isnan(accum_r(:,1)) & isnan(burst_r(:,1));
    %    classif.burst_only{dd} = isnan(accum_r(:,1)) & ~isnan(burst_r(:,1));
    %    classif.accum_burst{dd} = ~isnan(accum_r(:,1)) & ~isnan(burst_r(:,1));
end





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot histos of bsignif
nchannels=48;
nneurons=size(targ_bsignif_list,1);


% %subplot(1,2,1);hold on;
% targ_sum=nansum(targ_bsignif_list);
% plot(targ_sum,1:nchannels,'b')
%
% sacc_sum=nansum(sacc_bsignif_list);
% plot(sacc_sum,1:nchannels,'r')

vm_list=zeros(nneurons,nchannels);
vis_list=zeros(nneurons,nchannels);
mov_list=zeros(nneurons,nchannels);
normvect=zeros(1,nchannels);
for ch=1:nchannels
    for n=1:nneurons
        if targ_bsignif_list(n,ch)==1 & sacc_bsignif_list(n,ch)==1
            vm_list(n,ch)=1;
        end
        
        if targ_bsignif_list(n,ch)==1 & sacc_bsignif_list(n,ch)==0
            vis_list(n,ch)=1;
        end
        
        if targ_bsignif_list(n,ch)==0 & sacc_bsignif_list(n,ch)==1
            mov_list(n,ch)=1;
        end
        
        if targ_bsignif_list(n,ch)>0 | sacc_bsignif_list(n,ch)>0
            normvect(ch)=normvect(ch)+1;
        end
        
    end
end

%%
%%%%%%%%%%%%%%
figure;hold on;
%lims=[10 34];
lims=[16 32];
for t=1:4
    switch t
        case 1
            type='vm';
            vals_list=vm_list;
            color_t='g';
        case 2
            type='vis';
            vals_list=vis_list;
            color_t='b';
            
        case 3
            type='mov';
            vals_list=mov_list;
            color_t='r';

        case 4
            type='total';
            vals_list=vm_list+vis_list+mov_list;
            color_t='k';
    end
    
    %%%%%
    %no normalization
    subplot(1,3,1);hold on;
    vals_sum=nansum(vals_list);
    plot(vals_sum,1:nchannels,color_t)
    axis([0 20 lims(1) lims(2)]);
    title('no normalization')
    ylabel('Channel')
    xlabel('Number of neurons')
    
    %%%%%
    %display number of neurons for each type
    display([type ': '  num2str(sum(vals_sum))]);
    
    
    %%%%%
    %normalization by number of elements
    subplot(1,3,2);hold on;
    normmax=size(vals_list,1);
    vals_sumn=nansum(vals_list,1)/normmax;
    plot(vals_sumn,1:nchannels,color_t)
    axis([0 1 lims(1) lims(2)]);
    title('normalization by max')
    ylabel('Channel')
    xlabel('Number of neurons')
    
    
    %%%%%
    %normalization with normvect
    subplot(1,3,3);hold on;
    vals_sumnv=nansum(vals_list,1)./normvect;
    plot(vals_sumnv,1:nchannels,color_t)
    axis([0 1 lims(1) lims(2)]);
    title('normalization per channel')
    ylabel('Channel')
    xlabel('Number of neurons')
    
    
end







%MISC
% %%
% tried to use plot_histos_depths_v
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %plot histos type of activity
% fig_class_h=figure('Position',[scrsz(3)/4 100 scrsz(3)/2 scrsz(4)-200]);hold on;
% max_val=length(dlist);
% %plot total
% ltotal=max([size(classif.accum_only,2) size(classif.burst_only,2) size(classif.accum_burst,2)]);
% classif.total={};
% for dtot=[1:ltotal],
%     classif.total{dtot}=classif.accum_only{dtot}+classif.burst_only{dtot}+classif.accum_burst{dtot};
% end
%
% total=[];
% total(1,:)=plot_histos_depths_v(classif.total,classif.total,[ ],'k',[0 max_val],fig_class_h);
%
% total(2,:)=plot_histos_depths_v(classif.accum_only,classif.total,[ ],'g',[0 max_val],fig_class_h);
% total(3,:)=plot_histos_depths_v(classif.burst_only,classif.total,[ ],'m',[0 max_val],fig_class_h);
% total(4,:)=plot_histos_depths_v(classif.accum_burst,classif.total,[ ],'b',[0 max_val],fig_class_h);
% %plot_histos_depths_v(classif.nothing,classif.total,[ ],'k',[0 max_val],fig_class_h)
%
% %total(5,:)=plot_histos_depths_v(classif.accum_only_accum_burst_only,classif.total,[ ],'r',[0 max_val],fig_class_h);
%
% % direct 5 of each type
% total(:,4)=total(:,1)/total(1,1)*100;
%
% total
%
%
%
% subplot(1,3,1);
% text(10,8,num2str(round(total,2)))
% %legend('Significant bursts','Build-up only','Burst only','Build-up+Burst','Build-up only and Build-up+Burst','Location','NorthEast')
% legend('Significant bursts','Build-up only','Burst only','Build-up+Burst','Location','NorthEast')
%
% %proportion of each type of activity


