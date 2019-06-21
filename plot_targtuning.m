function plot_targtuning(alltrials,targs_ind,info,hdlfig,titlestr)

%function plot_targtuning(alltrials,targs_ind,hdlfig,titlestr)
%   Compute the tuning of each channel recorded with a laminar probe (LMA)
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 06/03/2016 last modified 01/10/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ntargs=size(alltrials,1);
nchannels=size(alltrials{1},1);
tuning_ch=zeros(nchannels,ntargs);

%to change the order of display of the targets and center the tuning
targs_index_display=[1:ntargs];%[1 8 7 6 5 4 3 2];

for ts=1:ntargs
   %get trials of target
   trials=alltrials{ts};
   
   %computing trial averages
   trials_avg=[];
   for ch=1:nchannels
       %average across trials (takes into account NaN values)
       trials_avg(ch,:)=nanmean(squeeze(trials(ch,:,:)));
       
       %%with normalization
       %trials_avg(ch,:)=mean(squeeze(trials(ch,:,:)));-mean(squeeze(trials(ch,:,:)));
   
       %peak of average firing rate
       tuning_ch(ch,targs_index_display(ts))=max(trials_avg(ch,:));
       %tuning_ch(ch,s)=min(trials_avg(ch,:));
   end
end


%figure
if ~isempty(hdlfig)
    subplot(hdlfig);
else
    figure;hold on;
end

%colorlist
colorlist=get_colorlist;


% %plot
% for ch=1:nchannels
%     plot(1:ntargs,tuning_ch(ch,targs_ind),'-o','Linewidth',2,'color',colorlist(ch,:));
% end
% set(gca,'xtick',[1:length(targs_ind)],'xticklabel',targs_ind);xlabel('Target index');                
% ylabel('Response');
% if ~isempty(titlestr),
%     title(titlestr)
% else
%     title({info.datafile ; info.align});
% end

%plot with shift
vshift=max(tuning_ch(:))/16;
for ch=1:nchannels
    plot(1:ntargs,tuning_ch(ch,targs_ind)+vshift*ch,'-o','Linewidth',2,'color',colorlist(ch,:));
end
axis tight;
set(gca,'xtick',[1:length(targs_ind)],'xticklabel',targs_ind);xlabel('Target index');                
set(gca,'ytick',[vshift:vshift:length(info.chmap)*vshift],'yticklabel',info.chmap);ylabel('Channel number')
if ~isempty(titlestr),
    title(titlestr)
else
    title({info.datafile ; info.align});
end


%MISC
% figtuning_targ=figure;hold on;
% for it=targs_ind;
%    % plot(1:nchannels,tuning_ch(1:nchannels,it),'-o');
%     plot(tuning_ch(1:nchannels,it),1:nchannels,'-o');
% pause
% end
% %set(gca,'xticklabel',1:nchannels);xlabel('Channels');ylabel('Response');
% set(gca,'yticklabel',1:nchannels);xlabel('Response');ylabel('Channels');
% title([datalist(d) [align ' t' num2str(it)]]);
% 
