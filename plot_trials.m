function [range vshift]=plot_trials(trials,var,index,vshift,events,event_align,info,hdlfig,titlestr,linetype,linewidth)

%function [range]=plot_trials(trials,var,index,vshift,events,event_align,info,hdlfig,titlestr,linetype,linewidth)
%  plot list of trials
%
% var: variance if trials is a list of averaged trials
% index: index of channels to display if missing channels
% vshift: vertical shift between channels (control scaling)
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 07/07/2016 last modified 01/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure
if ~isempty(hdlfig)
    subplot(hdlfig);
else
    figure;
end

%colorlist
colorlist=get_colorlist;

%linetype
if isempty(linetype)
    linetype='-';
end

%linewidth
if isempty(linewidth)
    linewidth=2;
end


%plot
if info.nchannels~=1,
    nchannels=size(trials,1);
    if isempty(vshift)
        vshift=max(abs(trials(:)))/4;
    end
    
    %index
    if isempty(index),index=[1:nchannels];end;
    
    for ch=1:nchannels
        chi=index(ch);
        
        if vshift~=0,
            plot(trials(ch,:)+vshift*chi,linetype,'Linewidth',linewidth,'color',colorlist(chi,:));
            if ~isempty(var)
                plot(trials(ch,:)+vshift*chi+var(ch,:),'Linewidth',1,'color',colorlist(chi,:));
                plot(trials(ch,:)+vshift*chi-var(ch,:),'Linewidth',1,'color',colorlist(chi,:));
            end
            %axes
            minval=min(trials(index(1),:));maxval=vshift*(chi)+max(trials(ch,:));
            %%minval=min(trials(index(1),:));maxval=vshift*(index(1))+max(trials(ch,:));
            
        elseif vshift==0
            plot(trials(ch,:),linetype,'Linewidth',linewidth,'color',colorlist(chi,:));
            if ~isempty(var)
                plot(trials(ch,:)+var(ch,:),'Linewidth',1,'color',colorlist(chi,:));
                plot(trials(ch,:)-var(ch,:),'Linewidth',1,'color',colorlist(chi,:));
            end
            %axes
            minval=min(trials(:));maxval=max(trials(:));
        end
        
        

    end
    
else
    %only 1 channel
    plot(trials,'Linewidth',2,'color',info.color);
    
    if ~isempty(var)
        plot(trials+var,'Linewidth',1,'color',info.color);
        plot(trials-var,'Linewidth',1,'color',info.color);
    end
    
    %axes
    minval=min(trials(:));maxval=max(trials(:));
    
end



axis tight;ax=axis;mintime=0-info.aligntime;maxtime=ax(2)-info.aligntime;step=50;%(maxtime-mintime)/5
vec=[ax(1):step:ax(2)];
vectime=[mintime:step:maxtime];
al_ind=min(find(vec>info.aligntime+1));

if ~isempty(find(vec==info.aligntime+1)) %+1 because vec starts at 1
    xtick_vec=vec;
    xticklabel_vec=vectime;
elseif ~isempty(find(vec==info.aligntime)) %one case once where vec started by 0 (?)
    xtick_vec=vec;
    xticklabel_vec=vectime;
elseif al_ind==2
    xtick_vec=[info.aligntime vec(al_ind:end)];
    xticklabel_vec=[0 vectime(al_ind:end)];
else
    %vec
    %info.aligntime
    xtick_vec=[vec(1:al_ind-1) info.aligntime vec(al_ind:end)];
    %xtick_vec
    xticklabel_vec=[vectime(1:al_ind-1) 0 vectime(al_ind:end)];
end

% figure;plot(xtick_vec)
% size(xtick_vec)
% size(xticklabel_vec)
% pause


set(gca,'xtick',xtick_vec,'xticklabel',xticklabel_vec);xlabel('Time (ms)');
set(gca,'ytick',[vshift:vshift:length(info.chmap)*vshift],'yticklabel',info.chmap);ylabel('Channel')
%set(gca,'ytick',[vshift:vshift:nchannels*vshift],'yticklabel',depths(chmap));ylabel('Depth (mm)');

%range
range=[mintime maxtime minval maxval];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot events
%always first plot event at aligntime
hl=line([info.aligntime info.aligntime] ,[minval maxval]);
set(hl,'Color',colorlist(1,:),'LineStyle','-','Linewidth',1);
if ~isempty(events),
    plot_events(events,event_align,info.aligntime,range,hdlfig,0);
end


if ~strcmp(titlestr,'n')
    if ~isempty(titlestr)
        title(titlestr);
    else
        title({info.datafile ; [info.align ' t' num2str(info.targ) ' #trials:' num2str(info.ntrials)]});
    end
end


