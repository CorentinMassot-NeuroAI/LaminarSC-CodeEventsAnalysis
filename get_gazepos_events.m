function events_new=get_gazepos_events(gazepos,events,event_align,wind,info,display_gazepos,hdlfig1,hdlfig2)

%function events_new=get_gazepos_events(gazepos,events,event_align,wind,info,disp)
%  get new events from gaze position signals
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 08/24/2017 last modified 08/24/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%minval maxval
if ~isempty(wind),
    minval=min(min(gazepos(:,events.fixon:events.fixtarget+50)))-1;
    maxval=max(max(gazepos(:,events.fixon:events.fixtarget+50)))+1;
else
    minval=min(gazepos(:))-0.5;
    maxval=max(gazepos(:))+0.5;
end

%colorlist
colorlist=get_colorlist;


if display_gazepos
    
    %Plotting gaze position with XY coordinates
    %scrsz = get(groot,'ScreenSize');
    %figure('Position',[1 100 scrsz(3)-100 scrsz(4)-200]);
                
    subplot(hdlfig1);hold on;
    
 
    %plot gaze position
    plot(gazepos(1,events.fixon:events.go),gazepos(2,events.fixon:events.go),'Linewidth',1,'color',colorlist(1,:));
    plot(gazepos(1,events.go:events.fixtarget),gazepos(2,events.go:events.fixtarget),'Linewidth',2,'color',colorlist(4,:));
    
    %plot markers
    %beginning
    %plot(gazepos(1,1),gazepos(2,1),'s','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',colorlist(6,:));
    %end
    %plot(gazepos(1,end),gazepos(2,end),'s','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',colorlist(8,:));
    %fixation onset
    plot(gazepos(1,events.fixon),gazepos(2,events.fixon),'s','MarkerSize',5,'MarkerFaceColor',colorlist(1,:));
    %target onset
    plot(gazepos(1,events.targ),gazepos(2,events.targ),'s','MarkerSize',5,'MarkerFaceColor',colorlist(2,:));
    %go onset
    plot(gazepos(1,events.go),gazepos(2,events.go),'s','MarkerSize',5,'MarkerFaceColor',colorlist(3,:));
    %sacc onset
    plot(gazepos(1,events.sacc),gazepos(2,events.sacc),'s','MarkerSize',5,'MarkerFaceColor',colorlist(4,:));
    %fixation target
    plot(gazepos(1,events.fixtarget),gazepos(2,events.fixtarget),'s','MarkerSize',5,'MarkerFaceColor',colorlist(5,:));
    
    axis([-35 35 -35 35]);grid;
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Plotting X and Y traces across time
    subplot(hdlfig2);hold on;
    
    %colorlist
    colorlist=get_colorlist;
    
    %plot
    %gazepos = resample(gazepos,3,1);
    plot([1:size(gazepos,2)],gazepos(1,:)','Linewidth',2,'color',colorlist(1,:));
    plot([1:size(gazepos,2)],gazepos(2,:)','Linewidth',2,'color',colorlist(4,:));
    
    
    %axes
    axis tight;ax=axis;
    if ~isempty(wind),
        axis([event_align+wind(1) event_align+wind(2) minval maxval]);
    else
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
    %title({info.datafile ; [info.align ' #trials:' num2str(info.ntrials)]});
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extracting new events for sacc onset and fixation target

%position of final fixation
finalfix(1) = gazepos(1,events.fixtarget);
finalfix(2) = gazepos(2,events.fixtarget);
events_new.finalfix=finalfix;


%velocity peaks
windp=[0 50];
gazepos_dx = smooth(diff(gazepos(1,:)),5,'moving')*10;
gazepos_dy = smooth(diff(gazepos(2,:)),5,'moving')*10;

end_dx=min([events.sacc+windp(2) length(gazepos_dx)]);
[peakv_gx peakt_gx]=max(abs(gazepos_dx(events.sacc+windp(1):end_dx)));
[peakv_gy peakt_gy]=max(abs(gazepos_dy(events.sacc+windp(1):end_dx)));
peakt_gx=peakt_gx+events.sacc+windp(1)-1;
peakt_gy=peakt_gy+events.sacc+windp(1)-1;

if peakv_gx>=peakv_gy
    peak(1,:)=[peakt_gx peakv_gx];
    peak(2,:)=[peakt_gy peakv_gy];
else
    peak(1,:)=[peakt_gy peakv_gy];
    peak(2,:)=[peakt_gx peakv_gx];
end

%events_new
events_new.peak=peak;


%plot
if display_gazepos
    plot([1:size(gazepos,2)-1],gazepos_dx,'Linewidth',1,'color',colorlist(2,:));
    plot([1:size(gazepos,2)-1],gazepos_dy,'Linewidth',1,'color',colorlist(5,:));
    
    hl=line([peak(1,1) peak(1,1)] ,[minval maxval]);
    set(hl,'Color',colorlist(2,:),'LineStyle','--','Linewidth',2);
    
    hl=line([peak(2,1) peak(2,1)] ,[minval maxval]);
    set(hl,'Color',colorlist(5,:),'LineStyle','--','Linewidth',1);
    
    title('Saccade velocity peaks')
end



