function plot_events(events,event_align,aligntime,range,hdlfig,disp_legend)

%function plot_events(events,event_align,range,hdlfig,disp_legend)
%  plot events times as vertical lines
%
% event_align: timing of event used to align data
%
%Corentin University of Pittsburgh Montreal 08/22/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TO DO everywhere check plot_ch or plot event times +1ms


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
%fixation onset
eventtime=events.fixon-event_align;
if eventtime>=mintime & eventtime<=maxtime,
    %display('fixation beginning')
    eventtime=eventtime+aligntime;
    hl=line([eventtime eventtime] ,[minval maxval]);
    set(hl,'Color',colorlist(1,:),'LineStyle','--');
end

%target onset
eventtime=events.targ-event_align;
if eventtime>=mintime & eventtime<=maxtime,
    %display('targ onset')
    eventtime=eventtime+aligntime;
    hl=line([eventtime eventtime] ,[minval maxval]);
    set(hl,'Color',colorlist(2,:),'LineStyle','--');
end

%go onset
eventtime=events.go-event_align;
if eventtime>=mintime & eventtime<=maxtime,
    %display('go')
    eventtime=eventtime+aligntime;
    hl=line([eventtime eventtime] ,[minval maxval]);
    set(hl,'Color',colorlist(3,:),'LineStyle','--');
end

%sacc onset
eventtime=events.sacc-event_align;
if eventtime>=mintime & eventtime<=maxtime,
    %display('sacc onset')
    eventtime=eventtime+aligntime;
    hl=line([eventtime eventtime] ,[minval maxval]);
    set(hl,'Color',colorlist(4,:),'LineStyle','--','Linewidth',2);
end

%fixation target
eventtime=events.fixtarget-event_align;
if eventtime>=mintime & eventtime<=maxtime,
    %display('fixation target')
    eventtime=eventtime+aligntime;
    hl=line([eventtime eventtime] ,[minval maxval]);
    set(hl,'Color',colorlist(5,:),'LineStyle','--','Linewidth',2);
end

% %reward
%eventtime=events.reward-event_align;
%if eventtime>=mintime & eventtime<=maxtime,
%eventtime=eventtime+aligntime;
%hl=line([eventtime eventtime] ,[minval maxval]);
% set(hl,'Color',colorlist(6,:),'LineStyle','--');

% %peak (from plot_gazedata)
% eventtime=events.peak-event_align;
% if eventtime>=mintime & eventtime<=maxtime,
%     eventtime=eventtime+aligntime;
%     hl=line([eventtime eventtime] ,[minval maxval]);
%     set(hl,'Color',colorlist(7,:),'LineStyle','--');
% end


if disp_legend,
    legend('fixon','targ','go','sacc','fixtarget','peak','location','North');
end

