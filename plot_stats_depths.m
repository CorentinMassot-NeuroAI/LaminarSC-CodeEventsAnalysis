function figsessions=plot_stats_depths(latmag,allonset,t_outliers,coefmult,colorline,lims_xaxis,hdlfig,info,datalist,dlist)

%function plot_stats_depths
%   Statisical analysis of any measure across depth of data recorded with a
%   laminar probe (LMA)
%
% latmag='latency' or 'magnitude';
% t_outliers: threshold for outliers
% coefmult: coefficient multiplicator of all values
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 01/05/2018 last modified 01/05/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%colorlist
colorlist=get_colorlist;

%channel limits
lims=[16 33];%[16 33];

%nchannels
nchannels=size(allonset{1},1);

ddlist=[1:length(allonset)];

%figsessions=figure;hold on;
figsessions=0;
onsets_all=[];
for dd=1:length(allonset),
    info.datafile=datalist{dlist(dd)};
    
    onset_aux=allonset{dd};
    switch latmag
        case 'latency'
            onsets=onset_aux(:,1)';
        case 'magnitude'
            onsets=coefmult*onset_aux(:,2)';
    end
    
    %lists
    onsets_all=[onsets_all ; onsets];
    
%     %plot fr latency
%     figure(figsessions)
%     subplot(1,1,1);hold on;
%     titlestr={info.datafile};
%     plot(1:nchannels,onsets,'-','color',colorlist(ddlist(dd),:),'Linewidth',2);
%     %plot(nchannels:-1:1,onsets,'-','color',colorlist(ddlist(dd),:),'Linewidth',2);
%     %axis([16 33 2*lims_xaxis(1) lims_xaxis(2)]);
%     axis tight;ax=axis;
%     axis([lims(1) lims(2) ax(3) ax(4)])
%     set(gca,'Xtick',[lims(1):2:lims(2)+1],'Xticklabel',[-8:2:10]);
%     xlabel('Channel');
%     switch latmag
%         case 'latency'
%             ylabel('Latency (ms)');
%         case 'magnitude'
%             ylabel('FR amplitude (spk/s)')
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
% %     pause
% %     end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%plot average and ci
%remove outliers
if ~isempty(t_outliers)
    onsets_all(find(onsets_all<t_outliers(1) | onsets_all>t_outliers(2)))=nan;
end 

%remove if only one value
%onsets_all
ind_1=find(sum(~isnan(onsets_all),1)==1);
onsets_all(:,ind_1)=nan;
%onsets_all
%pause
    

%average onsets
onsets_avg=nanmean(onsets_all,1);
%onsets_avg=nanmedian(onsets_all,1);
    
figure(hdlfig)
subplot(1,1,1);hold on;
plot(1:nchannels,onsets_avg,colorline,'Linewidth',3);
%plot(nchannels:-1:1,onsets_avg,colorline,'Linewidth',3);
xlabel('Channel (deepest -> highest)');


switch latmag
    case 'latency'
        ylabel('Latency (ms)');
    case 'magnitude'
        ylabel('FR amplitude (spk/s)')
end



%number of data point per channels
nval_ch=[];
for ch=1:size(onsets_all,2)
    nval_ch(ch)=numel((find(~isnan(onsets_all(:,ch)))));
end
nval_ch

%find channel range
chs_r=find(~isnan(onsets_avg));
[vmiss imiss]=find(chs_r(2:end)-chs_r(1:end-1)>1);
chs_r
imiss
    
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
        
    
%compute ci
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
%fill([max_ch:-1:min_ch min_ch:1:max_ch]+5,[onsets_ci(:,1)' fliplr(onsets_ci(:,2)')],1,'facecolor',colorline,'edgecolor','none','facealpha', 0.3);
%fill(48-[min_ch:max_ch max_ch:-1:min_ch],[onsets_ci(:,1)' fliplr(onsets_ci(:,2)')],1,'facecolor',colorline,'edgecolor','none','facealpha', 0.3);
    
%subplot(2,2,sig);hold on;
subplot(1,1,1);hold on;

axis([lims(1) lims(2) lims_xaxis(1) lims_xaxis(2)])
set(gca,'Xtick',[lims(1):2:lims(2)+1],'Xticklabel',[-8:2:10])        
%set(gca,'Xtick',[lims(1):2:lims(2)],'Xticklabel',[8:-2:-8])

display(['average: ' num2str(mean(onsets_avg(min_ch:max_ch))) ' and variance: ' num2str(var(onsets_avg(min_ch:max_ch)))])
   

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fit
onsets_avgaux=onsets_avg(1,lims(1):lims(2));
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
npoly=4;
[p rsq yfit_avg]=get_regressioncoefs(xp,yp,npoly);
%regressparams_avg(sig,1)=rsq;
%regressparams_avg(sig,2:2+npoly)=p;
%regressparams_avg(sig,2+npoly+1)=numel(xp);
vect=[1:1:length(xp)];
%plot(vect,fun(fitparams,vect),'k-','linewidth',2)
plot(vect+lims(1)-1,polyval(p,vect),[colorline '--'],'linewidth',2)




return
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%stats
%channels of interest
min_ch=16;
max_ch=33;
p=[];h=[];

%normality test (Kolmogorov-Smirnov)
chi=0;
for ch=min_ch:max_ch %limit test at channels of interest
    chi=chi+1;
    [h(chi),p(chi)] = kstest(onsets_all(:,ch));
end
%output
h
p



%%%%%%%%%%%%%%%%%%%%%%%%%
%non-parametric test

%spk/lfp
chi=0;
for ch=min_ch:max_ch
    chi=chi+1;
    [p(chi),h(chi),stats] = ranksum(onsets_spk_all(:,ch),onsets_lfp_all(:,ch));
    %[h(chi),p(chi),stats] = ttest2(onsets_spk_all(:,ch),onsets_lfp_all(:,ch));
end
h
p

%display(['Significant difference spk/lfp for channels:' num2str(find(p<0.01) + min_ch-1 -23)])
display(['Significant difference spk/lfp for channels:' num2str(find(p<0.05))])

%results:
%Significant difference for channels:14 (29)


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



return
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot ranked histograms
color_avgconf=['.b' '.r']
[nlat nchannels]=size(onsets_spk_all);

min_ch=16;max_ch=33;

%%remove outliers
%onsets_spk_all(find(onsets_spk_all<60 | onsets_spk_all>150))=nan;
%onsets_lfp_all(find(onsets_lfp_all<60 | onsets_lfp_all>200))=nan;

%rank order
onsets_spk_rk=[];onsets_lfp_rk=[];
for ch=min_ch:max_ch
    [onsets_spk_rk(:,ch) i_rk]=sort(onsets_spk_all(:,ch));
    onsets_lfp_rk(:,ch)=onsets_lfp_all(i_rk,ch);
    
    %[onsets_lfp_rk(:,ch) i_rk]=sort(onsets_lfp_all(:,ch));
    %onsets_spk_rk(:,ch)=onsets_spk_all(i_rk,ch);
end

figure;hold on;
for ch=min_ch:max_ch,
    subplot(5,4,ch-min_ch+1);hold on;
    if ch==min_ch,title('Lat / ch');end
    plot(onsets_spk_rk(:,ch),1:nlat,'o','MarkerSize',5,'MarkerFaceColor','b');
    plot(onsets_lfp_rk(:,ch),1:nlat,'o','MarkerSize',5,'MarkerFaceColor','r');
    axis([80 150 1 nlat+1]);axis square;
    xlabel(['ch' num2str(ch-min_ch+1)])
end


%rank order by latency differences
onsets_diff=[];
for ch=min_ch:max_ch
    onsets_diff=abs(onsets_spk_all(:,ch)-onsets_lfp_all(:,ch));
    [val_rk i_rk]=sort(onsets_diff);
    onsets_spk_rk(:,ch)=onsets_spk_all(i_rk,ch);
    onsets_lfp_rk(:,ch)=onsets_lfp_all(i_rk,ch);
end

figure;hold on;
for ch=min_ch:max_ch,
    subplot(5,4,ch-min_ch+1);hold on;
    if ch==min_ch,title('Lat (sort by diff) / ch');end
    plot(onsets_spk_rk(:,ch),1:nlat,'o','MarkerSize',5,'MarkerFaceColor','b');
    plot(onsets_lfp_rk(:,ch),1:nlat,'o','MarkerSize',5,'MarkerFaceColor','r');
    axis([80 150 1 nlat+1]);axis square;
    xlabel(['ch' num2str(ch-min_ch+1)])
end


return
%%
%latency differences vs. snr
onsets_diff=[];snr_spk=[];snr_lfp=[];
onsets_diff_rk=[];snr_spk_rk=[];snr_lfp_rk=[];
varbsl_spk_rk=[];varbsl_lfp_rk=[];

for ch=min_ch:max_ch
    
    %latency difference
    onsets_diff(:,ch)=(onsets_lfp_all(:,ch)-onsets_spk_all(:,ch));
    [onsets_diff_rk(:,ch) i_rk]=sort(onsets_diff(:,ch));
    
    %latency
    %[onsets_spk_rk(:,ch) i_rk]=sort(onsets_spk_all(:,ch));
    %onsets_lfp_rk(:,ch)=onsets_lfp_all(i_rk,ch);
    
    %[onsets_lfp_rk(:,ch) i_rk]=sort(onsets_lfp_all(:,ch));
    %onsets_spk_rk(:,ch)=onsets_spk_all(i_rk,ch);
    
    %     %snr
    %     snr_spk(:,ch)=var_spk_all(:,ch)./varbsl_spk_all(:,ch);
    %     snr_lfp(:,ch)=var_lfp_all(:,ch)./varbsl_lfp_all(:,ch);
    %     snr_spk_rk(:,ch)=snr_spk(i_rk,ch);
    %     snr_lfp_rk(:,ch)=snr_lfp(i_rk,ch);
    %
    %     %var
    %     varbsl_spk_rk(:,ch)=varbsl_spk_all(i_rk,ch);
    %     varbsl_lfp_rk(:,ch)=varbsl_lfp_all(i_rk,ch);
    
end


figure;hold on;
for ch=min_ch:max_ch,
    subplot(5,4,ch-min_ch+1);hold on;
    if ch==min_ch,title('Lat diff / snr');end
    plot(onsets_diff_rk(:,ch),snr_spk_rk(:,ch),'o','MarkerSize',5,'MarkerFaceColor','b');
    plot(onsets_diff_rk(:,ch),snr_lfp_rk(:,ch),'o','MarkerSize',5,'MarkerFaceColor','r');
    axis([-20 20 0 50]);axis square;
    xlabel(['ch' num2str(ch-min_ch+1)])
end

figure;hold on;
for ch=min_ch:max_ch,
    subplot(5,4,ch-min_ch+1);hold on;
    if ch==min_ch,title('Lat diff / var');end
    plot(onsets_diff_rk(:,ch),varbsl_spk_rk(:,ch),'o','MarkerSize',5,'MarkerFaceColor','b');
    plot(onsets_diff_rk(:,ch),varbsl_lfp_rk(:,ch),'o','MarkerSize',5,'MarkerFaceColor','r');
    axis([-20 20 0 50]);axis square;
    xlabel(['ch' num2str(ch-min_ch+1)])
end



