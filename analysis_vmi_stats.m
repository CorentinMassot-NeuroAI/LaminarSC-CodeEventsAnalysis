%function analysis_vmi_stats

%function analysis_vmi_stats
%   Analysis of data averaged over trials recorded with a laminar probe (LMA)
%   Compute stats of VMI index
%
% see also analysis_vmi
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 06/06/2017 last modified 06/06/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% The main statistical analysis was a repeated-measures ANOVA, whose
% assumptions were met for use in this study: (1) The independent variables
% (saliency-level and saccade-goal) were repeated/matched (that is, saliency-level
% and saccade-goal were extracted from within the same experimental sessions).
% (2) The dependent variable (neuronal discharge rate) was continuous. (3) The data
% across all conditions met the normality assumption based on Kolmogorov–
% Smirnov test (Po0.05). (4) The error variance between conditions was equal
% (that is, Mauchly’s sphericity test of equal variances was not violated). Observed
% statistical power was greater than 0.7 for the main interaction between saliencylevel
% and saccade-goal, and greater than 0.93 for all other significant main effects.
% T-tests were one-tailed unless otherwise stated, and based on a priori directional
% hypotheses. Alpha levels for multiple t-tests were Bonferroni-corrected, and
% running statistical tests were corrected using Bonferroni-Holm method for timeseries
% data.


%channels of interest
min_ch=16;
max_ch=32;%33
p=[];h=[];

tnormality=1


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SPK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get vmis or burst spk
%load('Stats/vmis_Delaysacc_vl')
%load('Stats/vmis_Delaysacc_vl_matched')
%for paper l functional organization
%load('Stats/vmis_r_alltuning_VG_vis');vmis=vmis_r_alltuning_VG_vis; %bursts
%load('Stats/vmis_r_alltuning_VG_mov');vmis=vmis_r_alltuning_VG_mov; %bursts
load('Stats/vmis_r_alltuning_VG_vmis');vmis=vmis_r_alltuning_VG_vmis; %vmis

vmis_vg=vmis;

%load('Stats/vmis_DelaySaccMem')
%load('Stats/vmis_DelaySaccMem_matched')
%for paper l functional organization
%load('Stats/vmis_r_alltuning_MG_vis');vmis=vmis_r_alltuning_MG_vis; %bursts
%load('Stats/vmis_r_alltuning_MG_mov');vmis=vmis_r_alltuning_MG_mov; %bursts
load('Stats/vmis_r_alltuning_MG_vmis');vmis=vmis_r_alltuning_MG_vmis; %vmis
vmis_mg=vmis;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normality test (Kolmogorov-Smirnov)
if tnormality
    %vg
    chi=0;
    for ch=min_ch:max_ch %limit test at channels of interest
        chi=chi+1;
        [h(chi),p(chi)] = kstest(vmis_vg(:,ch));
    end
    %output
    h
    p
    
    
    
    %mg
    chi=0;
    for ch=min_ch:max_ch
        chi=chi+1;
        [h(chi),p(chi)] = kstest(vmis_mg(:,ch));
    end
    %output
    h
    p
    
    
    %conclusions
    %most channels do not have a normal distribution
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%non-parametric test

chi=0;
for ch=min_ch:max_ch
    chi=chi+1;
    [p(chi),h(chi),stats] = ranksum(vmis_vg(:,ch),vmis_mg(:,ch));
    %[p(chi),tab,stats] = anova1([vmis_vg(:,ch) vmis_mg(:,ch)]);
end
h
p

display(['Significant difference for channels (0.05):' num2str(find(p<0.05))])
display(['Significant difference for channels (0.05/3):' num2str(find(p<0.05/3))])

%results
%Significant difference for channels:1 5 7

%chi=0;
% for ch=min_ch:max_ch-1
%      chi=chi+1;
%     [p(chi),h(chi),stats] = ranksum(vmis_mg(:,ch),vmis_mg(:,ch+1));
% end
% h
% p





return


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LFP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get vmis lfp
load('Stats/vmislfp_Delaysacc_vl')
vmis_vg=vmis;
load('Stats/vmislfp_DelaySaccMem')
vmis_mg=vmis;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normality test (Kolmogorov-Smirnov)
if tnormality
    %vg
    chi=0;
    for ch=min_ch:max_ch %limit test at channels of interest
        chi=chi+1;
        [h(chi),p(chi)] = kstest(vmis_vg{ch});
    end
    %output
    h
    p
    
    
    
    %mg
    chi=0;
    for ch=min_ch:max_ch
        chi=chi+1;
        [h(chi),p(chi)] = kstest(vmis_mg{ch});
    end
    %output
    h
    p
    
    
    %conclusions
    %most channels do not have a normal distribution
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%non-parametric test

chi=0;
for ch=min_ch:max_ch
    chi=chi+1;
    [p(chi),h(chi),stats] = ranksum(vmis_vg{ch},vmis_mg{ch});
end
h
p

display(['Significant difference for channels:' num2str(find(p<0.05))])

%results:
%Significant difference for channels:8

%chi=0;
% for ch=min_ch:max_ch-1
%      chi=chi+1;
%     [p(chi),h(chi),stats] = ranksum(vmis_mg{ch},vmis_mg{ch+1});
% end
% h
% p


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SPK vs. LFP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%vg
%get vmis spk and lflp
load('Stats/vmis_Delaysacc_vl')
%load('Stats/vmis_Delaysacc_vl_matched')
vmis_spk=vmis;
load('Stats/vmislfp_Delaysacc_vl')
vmis_lfp=vmis;

% %mg
% load('Stats/vmis_DelaySaccMem')
% %load('Stats/vmis_DelaySaccMem_matched')
% vmis_spk=vmis;
% load('Stats/vmislfp_DelaySaccMem')
% vmis_lfp=vmis;
% 

%%%%%%%%%%%%%%%%%%%%%%%%%
%non-parametric test

chi=0;
for ch=min_ch:max_ch
    chi=chi+1;
    [p(chi),h(chi),stats] = ranksum(vmis_spk(:,ch),vmis_lfp{ch});
    %[p(chi),tab,stats] = anova1([vmis_vg(:,ch) vmis_mg(:,ch)]);
end
h
p

display(['Significant difference for channels (0.05):' num2str(find(p<0.05))])
display(['Significant difference for channels (0.05/16):' num2str(find(p<0.05/16))])

%results
%Significant difference for channels (0.05):1   2   3   5   7   8   9  10  11  12  13
%Significant difference for channels (0.05/16):8  10


%zscore data
chi=0;
vmis_spk_z=nan(size(vmis_spk));
vmis_lfp_z={};
for chi=min_ch:max_ch
     ind_notnan=find(~isnan(vmis_spk(:,chi)));
     vmis_spk_z(ind_notnan,chi)=zscore(vmis_spk(ind_notnan,chi));

     vmis_lfp_z{chi}=zscore(vmis_lfp{chi});
end


chi=0;
for ch=min_ch:max_ch
    chi=chi+1;
    [p(chi),h(chi),stats] = ranksum(vmis_spk_z(:,ch),vmis_lfp_z{ch});
    %[p(chi),tab,stats] = anova1([vmis_vg(:,ch) vmis_mg(:,ch)]);
end
h
p

display(['Zscore Significant difference for channels (0.05):' num2str(find(p<0.05))])
display(['Zscore Significant difference for channels (0.05/16):' num2str(find(p<0.05/16))])

%results


%%
%regressions

%spk
channels=1:size(vmis_spk,2);
vmis=[];chs=[];
for s=1:size(vmis_spk,1)
    ind_notnan=find(~isnan(vmis_spk(s,:)));
    ind_notnan=ind_notnan(find(ind_notnan>=min_ch & ind_notnan<=max_ch))
    chs=[chs channels(ind_notnan)];
    vmis=[vmis vmis_spk(s,ind_notnan)];
end

[r,pstats]=corrcoef(chs,vmis);
display(['SPK correlation coefficient:' num2str(r(1,2))])
display(['p=' num2str(pstats(1,2))])


%lfp
channels=1:size(vmis_lfp,2);
vmis=[];chs=[];
i=0;
for chi=min_ch:max_ch
    i=i+1;
    vals=vmis_lfp{chi};
    ind_notnan2=find(~isnan(vals));
    chs=[chs ones(1,length(ind_notnan2))*chi];
    vmis=[vmis vals(ind_notnan2)'];
end

[r,pstats]=corrcoef(chs,vmis);
display(['LFP correlation coefficient:' num2str(r(1,2))])
display(['p=' num2str(pstats(1,2))])


%MISC
% [p,stats] = polyfit(channels,vmis,1);
% p
% yfit=p(2)+channels*p(1);
% 
% figure;hold on;
% plot(vmis,channels,'k');
% plot(yfit,channels,'b');
% line([-1 1] ,[25 25]);
% line([0 0],[min_ch max_ch]);
% axis([-1 1 min_ch max_ch])
% set(gca,'Ytick',[min_ch:max_ch],'Yticklabel',[-8:1:9]);
% ylabel('Channel')
% xlabel('VMI')
% %text(-0.95,max_ch-1,['vmi=' num2str(p(2)) '+ channels*' num2str(p(1))]);
% 
% 
% %     SSresid = sum((squeeze(meanFRlist(2,stimfit)) - yfit).^2);
% %     SStotal = (length(squeeze(meanFRlist(2,stimfit)))-1) * var(squeeze(meanFRlist(2,stimfit)));
% %     display('residual of responses:')
% %     rsq = 1 - SSresid/SStotal
% %


%(inv(stats.R)*inv(stats.R)')*stats.normr^2/stats.df
%r*sqrt(cov(channels,channels)*cov(vmis,vmis))









