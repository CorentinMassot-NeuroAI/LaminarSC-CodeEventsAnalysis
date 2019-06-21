function plot_vmis(vmis,targ,type,colorind,width,info,hdlfig,titlestr);

%function plot_vmis(vmis,targ,type,colorind,width,info,hdlfig,titlestr);
%  plot vmis
%
% see also compute_vmi
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 10/16/2016 last modified 01/19/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%figure1
if ~isempty(hdlfig)
    subplot(hdlfig);
else
    figure;
end

%colorlist
colorlist=get_colorlist;

%vmis
if strcmp(targ,'all'),
    for t=1:size(vmis,1)
        plot(vmis(t,:),[1:info.nchannels],'linewidth',1,'color',colorlist(t,:));
    end
else
    switch type
        case '-'
            plot(vmis(targ,:),[1:info.nchannels],'linewidth',2,'color',colorlist(colorind,:),'linewidth',width);
        case 'o'
            plot(vmis(targ,:),[1:info.nchannels],'o','linewidth',2,'MarkerFaceColor',colorlist(colorind,:),'linewidth',width);
    end
    
end

%vertical line at vmi=0
hl=line([0 0] ,[1 info.nchannels]);
set(hl,'Color',colorlist(1,:),'LineStyle','--','Linewidth',1);

%axis
axis([-1 1 1 info.nchannels]);
xlabel('VMI');ylabel('Channel');

%ticks
set(gca,'Xtick',[-1:0.2:1])

%titlestr
if ~isempty(titlestr),
    title(titlestr);
end;


