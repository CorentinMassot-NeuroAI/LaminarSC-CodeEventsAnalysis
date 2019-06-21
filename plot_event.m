function plot_event(events,aligntime,range,colorind,hdlfig)

%function plot_event(events,aligntime,range,colorind,hdlfig)
%  plot any event times as vertical lines
%
% event_align: timing of event used to align data
%
%Corentin University of Pittsburgh Montreal 08/22/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plotting gaze position with XY coordinates
%figure1
if ~isempty(hdlfig)
    subplot(hdlfig);
else
    figure;
end

%colorlist
colorlist=get_colorlist;

%range
mintime=range(1);
maxtime=range(2);
minval=range(3);
maxval=range(4);


%events
for ev=1:numel(events)    
    eventtime=events(ev)+aligntime;
    hl=line([eventtime eventtime] ,[minval maxval]);
    set(hl,'Color',colorlist(colorind,:),'LineStyle','--','Linewidth',2);
end
