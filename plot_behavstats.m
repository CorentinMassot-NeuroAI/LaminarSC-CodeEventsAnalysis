function plot_behavstats(stats,info,hdlfig,titlestr)

%function plot_behavstats(data,info,hdlfig,titlestr)
%  plot behavioral data statistics
%
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 01/23/2017 last modified 01/23/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%plot
if ~isempty(hdlfig)
    subplot(hdlfig);
else
    figure;
end

%histo
histo=stats';
edges=[0:10:600];
hist=histc(histo,edges);
bar(edges,hist,'histc')

%set(gca,'xtick',xtick_vec,'xticklabel',xticklabel_vec);xlabel('Time (ms)');



%titlestr
if ~isempty(titlestr),
    title(titlestr);
else
    title(['CSD'])
end
