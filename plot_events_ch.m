function plot_events_ch(events_ch,index,vshift,range,info,hdlfig,titlestr,linestyle,linewidth,color)

%function [range]=plot_events_ch(events_ch,index,vshift,info,hdlfig,titlestr,linestyle,linewidth,color)
%  plot events for each individual channel
%
% index: index of channels to display if missing channels
% vshift: verticla shift between channels (control scaling)
%
% see also plot_trials
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 11/07/2016 last modified 06/21/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure
if ~isempty(hdlfig)
    subplot(hdlfig);
else
    figure;
end


%range
maxtime=range(1);maxtime=range(2);
minval=range(3);maxval=range(4);


%plot
nchannels=size(events_ch,1);

%color
if isempty(color)
    colorlist=get_colorlist;
else
    colorlist(1:nchannels,:)=color;
end


%vshift
if isempty(vshift)
    vshift=maxval/4;
end

%index
if isempty(index),index=[1:nchannels];end;



%adjust timing
events_ch(:,1)=events_ch(:,1)+info.aligntime;


if nchannels~=1
for ch=1:nchannels
    chi=index(ch);
    if events_ch(chi,1)~=0
        hl=line([events_ch(chi,1) events_ch(chi,1)] ,[vshift*chi+events_ch(chi,2)-vshift/2 vshift*chi+events_ch(chi,2)+vshift/2]);
        set(hl,'Color',colorlist(chi,:),'LineStyle',linestyle,'Linewidth',linewidth);
    end
end
else
    if events_ch(1,1)~=0
        hl=line([events_ch(1,1) events_ch(1,1)] ,[minval maxval]);
        set(hl,'Color',colorlist(1,:),'LineStyle',linestyle,'Linewidth',linewidth);
    end
end

if ~strcmp(titlestr,'n')
    if ~isempty(titlestr)
        title(titlestr);
    else
        title({info.datafile ; [info.align ' t' num2str(info.targ) ' #trials:' num2str(info.ntrials)]});
    end
end


