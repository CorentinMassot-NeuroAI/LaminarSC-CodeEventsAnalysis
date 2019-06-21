function [figsessions max_plot stats min_ch max_ch fit_vect]=plot_stats_depths_v(latmag,allonset,t_outliers,param,colorline,lims_xaxis,hdlfig,info,datalist,dlist,hdlfigfit,corr)

%function plot_stats_depths
%   Statistical analysis of any measure across depth of data recorded with a
%   laminar probe (LMA)
%
% display measures vertically along the depths
%
% latmag='latency' or 'magnitude' or 'delay_magnitude'
% t_outliers: threshold for outliers
% param: coefficient multiplicator of all values or bin number of delay analysis
% correction: remove list data points
% max_plot: maximum value of
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 01/05/2018 last modified 03/17/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%colorlist
colorlist=get_colorlist;

%channel limits
%lims=[1 48];
%lims=[16 33];
%lims=[16 32];
lims=[14 34];

%nchannels
nchannels=size(allonset{1},1);

ddlist=[1:length(allonset)];

%figsessions=figure;hold on;
figsessions=0;
onsets_all=[];
dd=0;
for d=[1:length(allonset)],
    %d
    dd=dd+1;
    info.datafile=datalist{dlist(d)};
    
    onset_aux=allonset{d};
    switch latmag
        case 'latency'
            onsets=onset_aux(:,1)';
        case 'magnitude'
            onsets=param*onset_aux(:,2)';
        case 'magnitude_nnorm'
            onsets=param*onset_aux(:,3)';
        case 'magnitude_nnorm_targ'
            onsets=param*onset_aux(:,4)';
        case 'magnitude_nnorm_targ_pburst'
            onsets=param*onset_aux(:,5)';
            
        case 'delay_magnitude'
            if param<=size(onset_aux,2)
                onsets=onset_aux(:,param)';
            else
                onsets=[];
            end
        case 'lat_vis'
            onsets=param*squeeze(onset_aux(1,:,1));
            nchannels=size(allonset{1},2);

    end
    
    %lists
    onsets_all=[onsets_all ; onsets];
    
%     %plot fr latency
%     figure(figsessions)
%     subplot(1,1,1);hold on;
%     titlestr={info.datafile};
%     plot(onsets_all(dd,:),1:nchannels,'o-','color',colorlist(ddlist(dd),:),'Linewidth',2);
%     %plot(nchannels:-1:1,onsets,'-','color',colorlist(ddlist(dd),:),'Linewidth',2);
%     %axis([16 33 2*lims_xaxis(1) lims_xaxis(2)]);
%     axis tight;ax=axis;
%     %axis([ax(1) ax(2) lims(1) lims(2) ])
%     axis([lims_xaxis(1) lims_xaxis(2) lims(1) lims(2) ])
%     %set(gca,'Xtick',[lims(1):2:lims(2)+1],'Xticklabel',[-8:2:10]);
%     ylabel('Channel');
%     switch latmag
%         case 'latency'
%             xlabel('Latency (ms)');
%         case 'magnitude'
%             xlabel('FR amplitude (spk/s)')
%     end
    
% %     %ci if latency measures
% %     if strcmp(latmag,'latency')
% %         chs_r=find(~isnan(onset_aux(:,1)));
% %         [vmiss imiss]=find(chs_r(2:end)-chs_r(1:end-1)>1);
% %         chs_r
% %         imiss
% %         
% %         min_ch=chs_r(1);
% %         max_ch=max(chs_r);
% %         if length(imiss)==1 & imiss(1)<10
% %             min_ch=chs_r(imiss(1)+1);
% %         end
% %         if length(imiss)==1 & imiss(1)>=10
% %             max_ch=chs_r(imiss(1)-1);
% %         end
% %         if length(imiss)==2 & imiss(1)<10 & imiss(2)>10
% %             min_ch=chs_r(imiss(1)+1);
% %             max_ch=chs_r(imiss(2)-1);
% %         end
% %         min_ch
% %         max_ch
% %         size(onset_aux)
% %         onset_aux(min_ch:max_ch,3)
% %         onset_aux(min_ch:max_ch,4)
% %         fill([min_ch:max_ch max_ch:-1:min_ch],[onset_aux(min_ch:max_ch,3)' fliplr(onset_aux(min_ch:max_ch,4)')],1,'facecolor',colorlist(ddlist(dd),:),'edgecolor','none','facealpha', 0.3);
   % pause
% %     end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %plot average and ci
% %remove outliers 
% NOTE: better not use this option
if ~isempty(t_outliers)
    display('OUTLIERS HAVE BEEN REMOVED!!')
    onsets_all(find(onsets_all<t_outliers(1) | onsets_all>t_outliers(2)))=nan;
end 


%remove if singleton value (one value, separated by NaNs)
%NOTE: not useful
%ind_1=find(sum(~isnan(onsets_all),1)==1);
% for id=ind_1
%     if isnan(onsets_all(:,id-1)) & isnan(onsets_all(:,id+1)) 
%         onsets_all(:,id)=nan;
%     end
% end

%remove data point with only 1 estimation
%NOTE: mostly should discard first and last data point
ind_1=find(sum(~isnan(onsets_all),1)==1);
for id=ind_1
    onsets_all(:,id)=nan;
end


%average onsets
onsets_avg=nanmean(onsets_all,1);
%onsets_avg=nanmedian(onsets_all,1);


%correction is want to remove one data point
%NOTE: not useful
%onsets_avg(corr)=nan;

figure(hdlfig)
subplot(1,1,1);hold on;
%hdl=plot(onsets_avg,1:nchannels,colorline,'Linewidth',3);
%hdl=plot(onsets_avg,1:nchannels);
hdl=plot(1:nchannels,onsets_avg);
set(hdl,'color',colorline,'Linewidth',3);
xlabel('Channel');
%%max_plot (based on avg)
%max_plot=max(onsets_avg);


switch latmag
    case 'latency'
        ylabel('Latency (ms)');
    case 'magnitude'
        ylabel('FR amplitude (spk/s)')
    case 'delay_magnitude'
        ylabel('Firing rate (spk/s)')
end


%number of data point per channels
nval_ch=[];
for ch=1:size(onsets_all,2)
    nval_ch(ch)=numel((find(~isnan(onsets_all(:,ch)))));
end
%nval_ch
fliplr(nval_ch(lims(1):lims(2)))';
%pause

%find channel range
chs_r=find(~isnan(onsets_avg));
[vmiss imiss]=find(chs_r(2:end)-chs_r(1:end-1)>1);
%chs_r
%imiss
    
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
        
% %compute and display standard error of the mean (SEM)
% onsets_sem=[];
% sem=nanstd(onsets_all(:,min_ch:max_ch),1)./sqrt(nval_ch(min_ch:max_ch));
% onsets_sem(:,1)=onsets_avg(:,min_ch:max_ch)+sem;
% onsets_sem(:,2)=onsets_avg(:,min_ch:max_ch)-sem;
% fill([onsets_sem(:,1)' fliplr(onsets_sem(:,2)')],[min_ch:max_ch max_ch:-1:min_ch],1,'facecolor',colorline,'edgecolor','none','facealpha', 0.3);
    
%compute and display ci
ind=0;onsets_ci=[];
for ch=min_ch:max_ch,
    ind=ind+1;
    
    aux=(onsets_all(find(~isnan(onsets_all(:,ch))),ch));
    
    if numel(aux)<=1,
        onsets_ci(ind,:)=[onsets_avg(ch) ; onsets_avg(ch)];
    else
        onsets_ci(ind,:) = bootci(1000,{@mean,aux},'type','per');
    end
end
fill([min_ch:max_ch max_ch:-1:min_ch],[onsets_ci(:,1)' fliplr(onsets_ci(:,2)')],1,'facecolor',colorline,'edgecolor','none','facealpha', 0.3);

    

%axes
subplot(1,1,1);hold on;
%axis([lims_xaxis(1) lims_xaxis(2) lims(1) lims(2)])
axis([lims(1) lims(2) lims_xaxis(1) lims_xaxis(2)])
set(gca,'Xtick',[lims(1):2:lims(2)],'Xticklabel',[-10:2:10])        

display(['average: ' num2str(mean(onsets_avg(min_ch:max_ch))) ' and variance: ' num2str(var(onsets_avg(min_ch:max_ch)))])
  


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%stats (mean , CI ....)
stats=[];
stats(:,1)=onsets_avg(min_ch:max_ch);
stats(:,2:3)=onsets_ci(:,1:2); 


display(['Across depths mean=' num2str(mean(stats(:,1)))]);
display(['Across depths meanCI inf=' num2str(mean(stats(:,2)))]);
display(['Across depths meanCI sup=' num2str(mean(stats(:,3)))]);

[valstats istats]=min(stats(:,1));
display(['Across channel of min=' num2str(istats)]);
display(['Across depths min=' num2str(valstats)]);
display(['Across depths minCI inf=' num2str(stats(istats,2))]);
display(['Across depths minCI sup=' num2str(stats(istats,3))]);

[valstats istats]=max(stats(:,1));
display(['Across channel of max=' num2str(istats)]);
display(['Across depths max=' num2str(valstats)]);
display(['Across depths maxCI inf=' num2str(stats(istats,2))]);
display(['Across depths maxCI sup=' num2str(stats(istats,3))]);



%return
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fit
%onsets_avg(1,:)
%min(find(~isnan(onsets_avg(1,:))))

%limfit(1)=max(lims(1),min((find(~isnan(onsets_avg(1,min_ch:max_ch)))+min_ch-1) & (find(onsets_ci(:,1)'~=onsets_avg(1,min_ch:max_ch))+min_ch-1)));
limfit(1)=max(lims(1),min(find(~isnan(onsets_avg(1,:)))));
limfit(2)=min(lims(2),max(find(~isnan(onsets_avg(1,:)))));

onsets_avgaux=onsets_avg(1,limfit(1):limfit(2));
xp=find(~isnan(onsets_avgaux));
yp=onsets_avgaux(xp);
x0p=[xp(1) yp(1)];

% %linear fit
% fun = @(x,xdata)x(1)+x(2)*xdata;
% options=optimset('Display','on');
% [fitparams resnorm] = lsqcurvefit(fun,x0p,xp,yp,[],[],options);
% figure(hdlfig)
% subplot(1,1,1);hold on;
% plot(xp+lims(1)-1,fun(fitparams,xp),'k-','linewidth',2)

%regression
npoly=3;
[p resid_fit yfit_avg]=get_regressioncoefs(xp,yp,npoly);
%regressparams_avg(sig,1)=rsq;
%regressparams_avg(sig,2:2+npoly)=p;
%regressparams_avg(sig,2+npoly+1)=numel(xp);
vect=[min(xp):1:max(xp)];
fit_vect=polyval(p,vect);
%plot(vect,fun(fitparams,vect),'k-','linewidth',2)
figure(hdlfigfit);hold on;
hdlfit=plot(vect+limfit(1)-1,fit_vect,'k--','linewidth',2);
axis([lims(1) lims(2) lims_xaxis(1) lims_xaxis(2)])

%max_plot
[max_plot max_ploti]=max(polyval(p,vect));
[min_plot min_ploti]=min(polyval(p,vect));


display(['residuals=' num2str(resid_fit)]); 
%yfit_avg


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%stats (mean , CI ....)
stats=[];
stats(:,1)=onsets_avg(min_ch:max_ch);
min_ch
limfit(1)
%onsets_avg(1,limfit(1):limfit(2))
stats(:,2:3)=onsets_ci(:,1:2); 


display(['Across depths mean=' num2str(mean(stats(:,1)))]);
display(['Across depths meanCI inf=' num2str(mean(stats(:,2)))]);
display(['Across depths meanCI sup=' num2str(mean(stats(:,3)))]);

max_ploti
min_ploti

istats=min_ploti;
valstats=stats(min_ploti,1);
display(['Across channel of min=' num2str(istats)]);
display(['Across depths min=' num2str(valstats)]);
display(['Across depths minCI inf=' num2str(stats(istats,2))]);
display(['Across depths minCI sup=' num2str(stats(istats,3))]);

istats=max_ploti;
valstats=stats(max_ploti,1);
display(['Across channel of max=' num2str(istats)]);
display(['Across depths max=' num2str(valstats)]);
display(['Across depths maxCI inf=' num2str(stats(istats,2))]);
display(['Across depths maxCI sup=' num2str(stats(istats,3))]);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%permutation test to compute p value of cubic fit
yp_fit=polyval(p,xp);
figure;
subplot(1,4,1);hold on;
plot(yp_fit,'r','Linewidth',3)

[rsq_fit rmse_fit] = get_rsquare(yp,yp_fit) %with constant term
[resid_fit ~]=get_residuals(yp,yp_fit)


nshuff=1000 %1000;
rsqshuff_list=[];
rmseshuff_list=[];
residshuff_list=[];
for sh=1:nshuff
    %shuffle x indexes
    yp_shuff=yp(randperm(length(xp)));
    
    %plot new shuffle
    subplot(1,4,1);hold on;
    plot(yp_shuff);
    
    %compute new fit
    npoly=3;
    [p_shuff ~]=get_regressioncoefs(xp,yp_shuff,npoly);
    yp_fit_shuff=polyval(p_shuff,xp);

    %rsq and rmse
    [rsqshuff_list(sh) rmseshuff_list(sh)] = get_rsquare(yp_shuff,yp_fit_shuff);

    %residuals
    [residshuff_list(sh) stats_resid]=get_residuals(yp_shuff,yp_fit_shuff);
end


%histograms of rsq rmse and residuals
subplot(1,4,2);hold on;
edges=[-0.5:0.01:0.5];
hist=histc(rsqshuff_list,edges);
bar(edges,hist,'histc')
plot([rsq_fit rsq_fit],[0 max(hist)],'-b');%plot rsq_fit
xlabel('R squares')

subplot(1,4,3);hold on;
edges=[0:0.1:10];
hist=histc(rmseshuff_list,edges);
bar(edges,hist,'histc')
plot([rmse_fit rmse_fit],[0 max(hist)],'-b');%plot rsq_fit
xlabel('RMSEs')

subplot(1,4,4);hold on;
edges=[-0.5:0.01:0.5];
hist=histc(rsqshuff_list,edges);
bar(edges,hist,'histc')
plot([resid_fit resid_fit],[0 max(hist)],'-b');%plot rsq_fit
xlabel('Residuals')


p=double((1+sum(rsqshuff_list>=rsq_fit))/length(rsqshuff_list))
display(['P value of fit (rsq): ' num2str(p)]);
p=double((1+sum(rmseshuff_list<=rmse_fit))/length(rmseshuff_list))
display(['P value of fit (rmse): ' num2str(p)]);
p=double((1+sum(residshuff_list>=resid_fit))/length(residshuff_list))
display(['P value of fit (resid): ' num2str(p)]);


%pause

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normality test (Kolmogorov-Smirnov)
%H=1 reject null hypothesis : not from a normal distribution

%channels of interest
%min_ch=16;
%max_ch=32;
p=[];h=[];

chi=0;
for ch=min_ch:max_ch %limit test at channels of interest
    chi=chi+1;
    [h(chi),p(chi)] = kstest(onsets_all(:,ch));
end
%output
h
p

pause

%%%%%%%%%%%%%%%%%%%%%%%%%
% non-parametric test (kruskalwallis)
min_ch=16;
max_ch=32;
[p_onsets,table_onsets,stats_onsets] = kruskalwallis(squeeze(onsets_all(:,min_ch:max_ch)));
%[c,m,h,nms] = 
multcompare(stats_onsets,'display','on');
%output
h
p


% %%%%%%%%%%%%%%%%%%%%%%%%%
% % paired ttest
% min_ch=limfit(1);
% max_ch=limfit(2);
% onsets_aux=squeeze(onsets_all(:,min_ch:max_ch));
% 
% nchs=size(onsets_aux,2);
% 
% H=[];P=[];
% for ch=1:nchs-1
%     [H(ch) P(ch)]=ttest(onsets_aux(:,ch), onsets_aux(:,ch+1));
% end
% 
% %output
% H
% P
% 
% size(onsets_aux)
% onsets_aux(:,9)
% onsets_aux(:,end-2)
% onsets_aux(:,1+2)
% [H P]=ttest(onsets_aux(:,9), onsets_aux(:,9+8))
% [H P]=ttest(onsets_aux(:,1), onsets_aux(:,9-8))


pause

% %%%%%%%%%%%%%%%%%%%%%%%%%
% %parametric test
% %multiple anova
% min_ch=16;
% max_ch=32;
% [p_onsets,table_onsets,stats_onsets] = anova1(squeeze(onsets_all(:,min_ch:max_ch)));
% %[c,m,h,nms] = 
% multcompare(stats_onsets,'display','on');
% %output
% h
% p
% 
% pause



%MISC
% %%%%%%%%%%%%%%%%%%%%%%%%%
% % non-parametric test
% %spearman correlation
% min_ch=16;
% max_ch=32;
% [row_s col_s]=size(onsets_all(:,min_ch:max_ch));
% x_corr=repmat([1:col_s],[row_s 1]);
% y_corr=squeeze(onsets_all(:,min_ch:max_ch));
% 
% size(y_corr)
% [R_spear,p_spear] = corr(y_corr,'type','Spearman');
% %output
% R_spear
% p_spear




% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MISC: shuffling also sessions
% %computing p value of fit
% yp_fit=polyval(p,xp);
% figure;
% subplot(1,2,1);hold on;
% plot(yp_fit,'r','Linewidth',3)
% rsq_fit
% 
% 
% %compute shuffled fit
% limfit(1)=max(lims(1),min(find(~isnan(onsets_avg(1,:)))));
% limfit(2)=min(lims(2),max(find(~isnan(onsets_avg(1,:)))));
% 
% nshuff=1000;
% for sh=1:nshuff
%     onsets_allshuff=onsets_all(:,limfit(1):limfit(2));
%     
%     onsets_avgshuff=nanmean(onsets_allshuff,1);
% 
%     xpshuff=find(~isnan(onsets_avgshuff));
%     yp=onsets_avgsuff(xpshuff);
% 
%     yp_shuff=yp(randperm(length(xp)));
%     
%     subplot(1,2,1);hold on;
%     plot(yp_shuff);
%     [rsqshuff_list(sh) stats_resid]=get_residuals(yp_shuff,yp_fit);
%     %     rsqshuff_list(sh)
%     %     stats_resid.yresid
%     %     stats_resid.SSresid
%     %     stats_resid.SStotal
%     %     %pause
% end
% 
% %histo
% subplot(1,2,2);hold on;
% edges=[-3:0.1:3];
% hist=histc(rsqshuff_list,edges);
% bar(edges,hist,'histc')
% plot([rsq_fit rsq_fit],[0 max(hist)],'-b');%plot rsq_fit
% xlabel('Residuals')
% 
% p=double((sum(rsqshuff_list>=rsq_fit)+1)/nshuff);
% display(['P value of fit: ' num2str(p)])


% %%%%%%%%%%%%%%%%%%%%%%%%%
% %spk/lfp
% %non-parametric test
% chi=0;
% for ch=min_ch:max_ch
%     chi=chi+1;
%     [p(chi),h(chi),stats] = ranksum(onsets_spk_all(:,ch),onsets_lfp_all(:,ch));
%     %[h(chi),p(chi),stats] = ttest2(onsets_spk_all(:,ch),onsets_lfp_all(:,ch));
% end
% h
% p
% 
% %display(['Significant difference spk/lfp for channels:' num2str(find(p<0.01) + min_ch-1 -23)])
% display(['Significant difference spk/lfp for channels:' num2str(find(p<0.05))])
% 
% %results:
% %Significant difference for channels:14 (29)


% %spk/csd
% chi=0;
% for ch=min_ch:max_ch
%     chi=chi+1;
%     [p(chi),h(chi),stats] = ranksum(onsets_spk_all(:,ch),onsets_lfp_all(:,ch));
%     %[h(chi),p(chi),stats] = ttest2(onsets_spk_all(:,ch),onsets_csd_all(:,ch));
% end
% h
% p
%
% display(['Significant difference spk/csd for channels:' num2str(find(p<0.01) + min_ch-1 -23)])
%
%
% %lfp/csd
% chi=0;
% for ch=min_ch:max_ch
%     chi=chi+1;
%     [p(chi),h(chi),stats] = ranksum(onsets_lfp_all(:,ch),onsets_csd_all(:,ch));
%     %[h(chi),p(chi),stats] = ttest2(onsets_lfp_all(:,ch),onsets_csd_all(:,ch));
% end
% h
% p
%
% display(['Significant difference lfp/csd for channels:' num2str(find(p<0.01) + min_ch-1 -23)])
%
% %results:
% %Significant difference for channels:14 (29)



% return
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %plot ranked histograms
% color_avgconf=['.b' '.r']
% [nlat nchannels]=size(onsets_spk_all);
% 
% min_ch=16;max_ch=33;
% 
% %%remove outliers
% %onsets_spk_all(find(onsets_spk_all<60 | onsets_spk_all>150))=nan;
% %onsets_lfp_all(find(onsets_lfp_all<60 | onsets_lfp_all>200))=nan;
% 
% %rank order
% onsets_spk_rk=[];onsets_lfp_rk=[];
% for ch=min_ch:max_ch
%     [onsets_spk_rk(:,ch) i_rk]=sort(onsets_spk_all(:,ch));
%     onsets_lfp_rk(:,ch)=onsets_lfp_all(i_rk,ch);
%     
%     %[onsets_lfp_rk(:,ch) i_rk]=sort(onsets_lfp_all(:,ch));
%     %onsets_spk_rk(:,ch)=onsets_spk_all(i_rk,ch);
% end
% 
% figure;hold on;
% for ch=min_ch:max_ch,
%     subplot(5,4,ch-min_ch+1);hold on;
%     if ch==min_ch,title('Lat / ch');end
%     plot(onsets_spk_rk(:,ch),1:nlat,'o','MarkerSize',5,'MarkerFaceColor','b');
%     plot(onsets_lfp_rk(:,ch),1:nlat,'o','MarkerSize',5,'MarkerFaceColor','r');
%     axis([80 150 1 nlat+1]);axis square;
%     xlabel(['ch' num2str(ch-min_ch+1)])
% end
% 
% 
% %rank order by latency differences
% onsets_diff=[];
% for ch=min_ch:max_ch
%     onsets_diff=abs(onsets_spk_all(:,ch)-onsets_lfp_all(:,ch));
%     [val_rk i_rk]=sort(onsets_diff);
%     onsets_spk_rk(:,ch)=onsets_spk_all(i_rk,ch);
%     onsets_lfp_rk(:,ch)=onsets_lfp_all(i_rk,ch);
% end
% 
% figure;hold on;
% for ch=min_ch:max_ch,
%     subplot(5,4,ch-min_ch+1);hold on;
%     if ch==min_ch,title('Lat (sort by diff) / ch');end
%     plot(onsets_spk_rk(:,ch),1:nlat,'o','MarkerSize',5,'MarkerFaceColor','b');
%     plot(onsets_lfp_rk(:,ch),1:nlat,'o','MarkerSize',5,'MarkerFaceColor','r');
%     axis([80 150 1 nlat+1]);axis square;
%     xlabel(['ch' num2str(ch-min_ch+1)])
% end
% 
% 
% return
% %%
% %latency differences vs. snr
% onsets_diff=[];snr_spk=[];snr_lfp=[];
% onsets_diff_rk=[];snr_spk_rk=[];snr_lfp_rk=[];
% varbsl_spk_rk=[];varbsl_lfp_rk=[];
% 
% for ch=min_ch:max_ch
%     
%     %latency difference
%     onsets_diff(:,ch)=(onsets_lfp_all(:,ch)-onsets_spk_all(:,ch));
%     [onsets_diff_rk(:,ch) i_rk]=sort(onsets_diff(:,ch));
%     
%     %latency
%     %[onsets_spk_rk(:,ch) i_rk]=sort(onsets_spk_all(:,ch));
%     %onsets_lfp_rk(:,ch)=onsets_lfp_all(i_rk,ch);
%     
%     %[onsets_lfp_rk(:,ch) i_rk]=sort(onsets_lfp_all(:,ch));
%     %onsets_spk_rk(:,ch)=onsets_spk_all(i_rk,ch);
%     
%     %     %snr
%     %     snr_spk(:,ch)=var_spk_all(:,ch)./varbsl_spk_all(:,ch);
%     %     snr_lfp(:,ch)=var_lfp_all(:,ch)./varbsl_lfp_all(:,ch);
%     %     snr_spk_rk(:,ch)=snr_spk(i_rk,ch);
%     %     snr_lfp_rk(:,ch)=snr_lfp(i_rk,ch);
%     %
%     %     %var
%     %     varbsl_spk_rk(:,ch)=varbsl_spk_all(i_rk,ch);
%     %     varbsl_lfp_rk(:,ch)=varbsl_lfp_all(i_rk,ch);
%     
% end
% 
% 
% figure;hold on;
% for ch=min_ch:max_ch,
%     subplot(5,4,ch-min_ch+1);hold on;
%     if ch==min_ch,title('Lat diff / snr');end
%     plot(onsets_diff_rk(:,ch),snr_spk_rk(:,ch),'o','MarkerSize',5,'MarkerFaceColor','b');
%     plot(onsets_diff_rk(:,ch),snr_lfp_rk(:,ch),'o','MarkerSize',5,'MarkerFaceColor','r');
%     axis([-20 20 0 50]);axis square;
%     xlabel(['ch' num2str(ch-min_ch+1)])
% end
% 
% figure;hold on;
% for ch=min_ch:max_ch,
%     subplot(5,4,ch-min_ch+1);hold on;
%     if ch==min_ch,title('Lat diff / var');end
%     plot(onsets_diff_rk(:,ch),varbsl_spk_rk(:,ch),'o','MarkerSize',5,'MarkerFaceColor','b');
%     plot(onsets_diff_rk(:,ch),varbsl_lfp_rk(:,ch),'o','MarkerSize',5,'MarkerFaceColor','r');
%     axis([-20 20 0 50]);axis square;
%     xlabel(['ch' num2str(ch-min_ch+1)])
% end
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% %MISC
% %for verticla plot
% % set(hdlfit,'Marker','o');
% % direction = [1 1 0];
% % rotate(hdlfit,direction,180)
% 
% % [xplot is]=sort(polyval(p,vect));
% % xplot
% % yplot=vect+lims(1)-1;yplot=yplot(is);
% % yplot
% %hdlfit=plot(yplot,xplot,[colorline '-'],'linewidth',2);
% %set(hdlfit,'Marker','o');
% 
% 
% 
% 
