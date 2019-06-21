function plot_gazedata(gazepos,events,event_align,wind,info,hdlfig1,titlestr1,hdlfig2,titlestr2)

%function plot_gazedata(gazepos,events,event_align,wind,info,hdlfig1,titlestr1,hdlfig2,titlestr2)
%  plot list of trials
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 07/07/2016 last modified 01/17/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plotting gaze position with XY coordinates
%figure1
if ~isempty(hdlfig1)
    subplot(hdlfig1);
else
    scrsz = get(groot,'ScreenSize');
    figure('Position',[scrsz(1)+400 scrsz(2)+300 scrsz(3)/2 scrsz(4)/2]);
    hdlfig1=subplot(1,2,1);hold on;
end

%colorlist
colorlist=get_colorlist;

%plot gaze position
plot(gazepos(1,events.fixon:events.go),gazepos(2,events.fixon:events.go),'Linewidth',1,'color',colorlist(1,:));%blue
plot(gazepos(1,events.go:events.fixtarget),gazepos(2,events.go:events.fixtarget),'Linewidth',2,'color',colorlist(4,:));%red

%plot markers
%beginning
%plot(gazepos(1,1),gazepos(2,1),'s','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',colorlist(6,:));
%end
%plot(gazepos(1,end),gazepos(2,end),'s','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',colorlist(8,:));
%fixation onset
plot(gazepos(1,events.fixon),gazepos(2,events.fixon),'s','MarkerSize',5,'MarkerFaceColor',colorlist(1,:));%blue
%target onset
plot(gazepos(1,events.targ),gazepos(2,events.targ),'s','MarkerSize',5,'MarkerFaceColor',colorlist(2,:));%cyan
%go onset
plot(gazepos(1,events.go),gazepos(2,events.go),'s','MarkerSize',5,'MarkerFaceColor',colorlist(3,:));%magenta
%sacc onset
plot(gazepos(1,events.sacc),gazepos(2,events.sacc),'s','MarkerSize',5,'MarkerFaceColor',colorlist(4,:));%red
%fixation target
plot(gazepos(1,events.fixtarget),gazepos(2,events.fixtarget),'s','MarkerSize',5,'MarkerFaceColor',colorlist(5,:));%green

axis tight;

if ~isempty(titlestr1),
    title(titlestr1);
    % else
    %     title({info.datafile ; [info.align ' t' num2str(info.targ) ' #trials:' num2str(info.ntrials)]});
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting X and Y traces across time
if ~isempty(hdlfig2)
    subplot(hdlfig2);hold on;
    % else
    %     hdlfig2=subplot(1,2,2);hold on;
    % end
    
    %colorlist
    colorlist=get_colorlist;
    
    %plot
    %gazepos = resample(gazepos,3,1);
    plot([1:size(gazepos,2)],gazepos(1,:)','Linewidth',2,'color',colorlist(1,:));
    plot([1:size(gazepos,2)],gazepos(2,:)','Linewidth',2,'color',colorlist(4,:));
    
    
    %axes
    axis tight;ax=axis;
    if ~isempty(wind),
        minval=min(min(gazepos(:,events.fixon:events.fixtarget+50)))-3;
        maxval=max(max(gazepos(:,events.fixon:events.fixtarget+50)))+3;
        axis([event_align+wind(1) event_align+wind(2) minval maxval]);
    else
        minval=min(gazepos(:))-0.5;maxval=max(gazepos(:))+0.5;
        axis([ax(1) ax(2) minval-1 maxval+1]);
    end
    
    %tick marks zeroed on event_align
    ax=axis;
    mintime=wind(1);maxtime=wind(2);step=50;%(maxtime-mintime)/5
    vec=[ax(1):step:ax(2)];vectime=[mintime:step:maxtime];
    if ~isempty(find(vec==event_align))
        %if isempty(event_align) | ~isempty(find(vec==event_align))
        xtick_vec=vec;
        xticklabel_vec=vectime;
    else
        al_ind=min(find(vec>event_align));
        xtick_vec=[vec(1:al_ind-1) event_align vec(al_ind:end)]
        xticklabel_vec=[vectime(1:al_ind-1) 0 vectime(al_ind:end)]
    end
    set(gca,'xtick',xtick_vec,'xticklabel',xticklabel_vec);xlabel('Time (ms)');
    ylabel('Position (degree)');
    
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot events
    %always first plot event at aligntime
    %if ~isempty(event_align)
    hl=line([event_align event_align] ,[minval maxval]);
    set(hl,'Color',colorlist(1,:),'LineStyle','-','Linewidth',1);
    if ~isempty(events),
        plot_events(events,event_align,event_align,[mintime maxtime minval maxval],hdlfig2,0);
    end
    %end
    
    %legend
    %legend('X','Y','fixon','targ','go','sacc','fixtarget','location','North');
    
  
    %title
    if ~isempty(titlestr2),
        title(titlestr2);
    else
        title({info.datafile ; [info.align ' t' num2str(info.targ) ' #trials:' num2str(info.ntrials)]});
    end;
    
    
end


