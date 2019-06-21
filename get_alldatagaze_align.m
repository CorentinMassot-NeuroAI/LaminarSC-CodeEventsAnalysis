function [allgazepos,allevents]=get_alldatagaze_align(data,seltrials,info,targslist)

%[allgazepos,events]=get_alldatagaze_align(data,seltrials,info,targslist)
%   get all eye gaze data for each trial and each target
%
% Warning: the conditions on the selected trials have to be the same as in
% get_alltrials_align
%
% seltrials:see also get_seltrials
% wind: extracted windows over which signals will be analyzed
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 07/30/2016 last modified 01/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calibration parameters
%NOTE retrieve and use numbers for each session instead
gh=1;%0.8;
bh=0;
gv=1;%0.8;
bv=0;

%ratio of data sampling over eye tracker sampling
samp_ratio=1;


ntargs=size(targslist,1);
allgazepos=cell(ntargs,1);
allevents=cell(ntargs,1);
ntrials=zeros(ntargs,1);
ntrials_d=zeros(ntargs,1);

for t=seltrials
    
    %get target index for the trial
    pos=data(t).offline.targpos;
    it=find(pos(1,1)==targslist(:,1) & pos(1,2)==targslist(:,2));
    
    %get gaze data and calibrate
    gazepos=data(t).gazePosition;
    %resampling
    gazeposaux=[];
    gazeposaux(1,:)=resample(gazepos(1,:),samp_ratio,1);
    gazeposaux(2,:)=resample(gazepos(2,:),samp_ratio,1);
    gazepos=gazeposaux;
    
    %%%aligntime=0;
    gazepos(1,:)=gh*(gazepos(1,:)+bh);
    gazepos(2,:)=gv*(gazepos(2,:)+bv);
    
    %gazedelayhvel
    %gazedelayvvel
    
    %creating list of trials (taking into account different sizes)
    if ~isempty(allgazepos{it}),allgazeposaux=allgazepos{it};else allgazeposaux={};end;
    ntrials(it)=ntrials(it)+1;
    %if ntrials(it)==1, tlen=size(gazepos,2);else tlen=min(tlen,size(gazepos,2));end
    allgazeposaux{ntrials(it)}=gazepos(:,:);
    allgazepos{it}=allgazeposaux;
    
    % colorlist=get_colorlist;
    % figure;hold on;
    % plot(gazepos(1,:),gazepos(2,:),'Linewidth',2,'color',colorlist(1,:));
    % plot(gazepos(1,1:100),gazepos(2,1:100),'Linewidth',2,'color',colorlist(3,:));
    % pause
    
    %gazepos(1,:)=fliplr(gazepos(1,:));
    %plot(gazepos(1,1:200),gazepos(2,1:200),'Linewidth',2,'color',colorlist(4,:));
    %plot(1:size(gazepos,2),gazepos(1,:),'Linewidth',2,'color',colorlist(4,:));
    %pause
    
    
    %events
    events=[];
    
    %sampling correction
    %data(t).stateTransitions(2,:)=data(t).stateTransitions(2,:);%*1/samp_ratio;
    
    %fixation onset
    events.fixon=get_event(data(t),'fixonCode');
    
    %target onset
    events.targ=get_event(data(t),'targCode');
    
    %go onset
    events.go=get_event(data(t),'goCode');
    
    %sacc onset
    events.sacc=get_event(data(t),'sacc','rpt');
    
    %fixation target
    events.fixtarget=get_event(data(t),'fixtarget');
    
    %%reward
    %events.reward=get_event(data(t),'reward');
    
    %creating list of events
    if ~isempty(allevents{it}),alleventsaux=allevents{it};else alleventsaux={};end;
    alleventsaux{ntrials(it)}=events;
    allevents{it}=alleventsaux;
    
    
end

%extracting window of analysis %and ensuring all trials have same length
%tlen
%tlen=min(tlen, 5000);
% for it=1:ntargs,
%     allgazeposaux=allgazepos{it};
%     if ~isempty(wind)
%WARNING change the following
%         allgazeposaux=allgazeposaux{:,:,wind(1)+aligntime:wind(2)+aligntime};
%     %else
%     %    allgazeposaux=allgazeposaux(:,:,1:tlen);
%     end
%     allgazepos{it}=allgazeposaux;
% end



%info
display(['Gaze ' info.align])
display(['# of trials: ' num2str(ntrials')])
display(['# of trials with NaN values: ' num2str(ntrials_d')])
display(['Duration of each trial: ' num2str(size(allgazepos{1},3))])

