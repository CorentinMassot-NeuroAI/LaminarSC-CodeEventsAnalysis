function [total]=plot_stats_depths_v(allvals,allnormvect,t_outliers,colorline,lims_xaxis,hdlfig)

%function plot_stats_depths
%   histograms of any measure across depth of data recorded with a
%   laminar probe (LMA)
%
% display measures vertically along the depths

% t_outliers: threshold for outliers (should not use that)
% allnormvect: structure containing vectors of normalization values if not by default use the number
% of elements in allvals
% total=table containing the total number of values computng in each plotted histo
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 01/05/2018 last modified 03/17/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%colorlist
colorlist=get_colorlist;

%channel limits
lims=[14 32];%[16 32] ;%[16 33];

%nchannels
nchannels=size(allvals{1},1);

vals_all=[];normvect_all=[];
for d=[1:length(allvals)],
    %lists
    vals_all=[vals_all ; allvals{d}'];
 
    %normalization
    normvect_all=[normvect_all ; allnormvect{d}'];
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%histo

% %remove outliers (better not use this option)
% if ~isempty(t_outliers)
%     onsets_all(find(onsets_all<t_outliers(1) | onsets_all>t_outliers(2)))=nan;
% end 

% %remove if singleton value (one value, separated by NaNs)
% ind_1=find(sum(~isnan(onsets_all),1)==1);
% for id=ind_1
%     if isnan(onsets_all(:,id-1)) & isnan(onsets_all(:,id+1)) 
%         onsets_all(:,id)=nan;
%     end
% end

%sum of vals
normmax=size(vals_all,1);
vals_sum=nansum(vals_all,1);

%normalization by number of elements
vals_sumn=nansum(vals_all,1)/normmax;

%normalization with normvect
vals_sumnv=nansum(vals_all,1)./sum(normvect_all,1);

% %correction
% onsets_avg(corr)=nan;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot
figure(hdlfig)
subplot(1,3,1);hold on;
hdl=plot(1:nchannels,vals_sum);
set(hdl,'color',colorline,'Linewidth',3);
xlabel('Channel');
ylabel('Count');

subplot(1,3,2);hold on;
hdl=plot(1:nchannels,vals_sumn);
set(hdl,'color',colorline,'Linewidth',3);
xlabel('Channel');
ylabel(['Ratio (norm:' num2str(normmax) ')']);

subplot(1,3,3);hold on;
hdl=plot(1:nchannels,vals_sumnv);
set(hdl,'color',colorline,'Linewidth',3);
xlabel('Channel');
ylabel(['Ratio (by number of significant bursts)']);


%axes
subplot(1,3,1);hold on;
grid
axis([lims(1) lims(2) lims_xaxis(1) lims_xaxis(2) ])
set(gca,'Xtick',[lims(1):2:lims(2)],'Xticklabel',[-10:2:8])        

subplot(1,3,2);hold on;
grid
axis([lims(1) lims(2) 0 1 ])
set(gca,'Xtick',[lims(1):2:lims(2)],'Xticklabel',[-10:2:8])        

subplot(1,3,3);hold on;
grid
axis([lims(1) lims(2) 0 1 ])
set(gca,'Xtick',[lims(1):2:lims(2)],'Xticklabel',[-10:2:8])        
% subplot(1,3,3);hold on;
% axis([lims(1) lims(2) 0 1 ])

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Total
total=[];

%number of values
aux=vals_sum(lims(1):lims(2));
total(1,1)=nansum(aux);


%normalized by number of elements
aux=vals_sumn(lims(1):lims(2));
size(aux);
naux=numel(find(aux~=0));
total(1,2)=nansum(aux)/naux;

%normalized by normvect
aux=vals_sumnv(lims(1):lims(2));
size(aux);
naux=numel(find(~isnan(aux)));
total(1,3)=nansum(aux)/naux;



