function [alltrials aligntime lut_trials]=get_alltrials_align(data,seltrials,wind,signal,info,targslist,sigma_FR,dispnan)

%[alltrials aligntime]=get_alltrials_align(data,seltrials,wind,signal,info,targslist,sigma_FR,dispnan)
%   get all trials for each channels from data recorded with a laminar probe (LMA)
%       spikes
%       LFP
%
% seltrials: see select_trials
% wind: extracted windows over which signals will be analyzed
% dispnan: 1 if want to display trial where there are NaN values
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 07/06/2016 last modified 01/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%timevec
if ~isempty(wind)
    timevec=[wind(1):wind(2)];
else
    timevec=[];
end

%seltrials
if isempty(seltrials)
    seltrials=[1:numel(data)];
end


%initializations
nchs=numel(data(1).spikeTimestamps);
ntargs=size(targslist,1);
alltrials=cell(ntargs,1);
ntrials=zeros(ntargs,1);
ntrials_d=zeros(ntargs,1);
ntrials_n=zeros(ntargs,1);
%if dispnan,fignan=figure;hold off;end;

%trials loo-up table
lut_trials=cell(ntargs,1);


%get data only for selected_trials
for t=seltrials
    
    %get target index for the trial
    pos=data(t).offline.targpos;
    it=find(pos(1,1)==targslist(:,1) & pos(1,2)==targslist(:,2));
    
    %get data according to signal
    option='';
    switch info.align
        %             case 'no'
        %                 switch signal
        %                     case 'fr'
        %                         trialsize=size(data(t).lfp,2)+20;
        %                         trial=zeros(nchs,trialsize);
        %                         for ch=1:nchs
        %                             spiketimes=round(1000*data(t).spikeTimestamps{ch});
        %                             spiketimes=spiketimes(find(spiketimes~=0));
        %                             trial(ch,spiketimes)=1;
        %                         end
        %                     case 'lfp'
        %                         trial=data(t).lfp;
        %                 end
        
        case 'targ'
            code='targCode';
            
        case 'go'
            code='goCode';
            
        case 'sacc'
            code='sacc';
            option='rpt'; %use .rpt file (EyeMove)
%             if t==seltrials(1),data(t).offline.sacc_pburst_ch_align
%            end
            
        case 'peak'
            code='peak';
            %option='offline'; %use offline data
            
        case 'targ_pburst'
            code='targ_pburst';
            
        case 'targ_pburst_ch'
            code='targ_pburst_ch';  
            if t==seltrials(1),data(t).offline.targ_pburst_ch_align
            end
            
        case 'targ_rsburst'
            code='targ_rsburst';
            
        case 'targ_rsburst_ch'
            code='targ_rsburst_ch';            
            if t==seltrials(1),data(t).offline.targ_rsburst_ch_align
            end

        case 'sacc_pburst'
            code='sacc_pburst';
            option='rpt'; %use .rpt file (EyeMove) for sacc
            if t==seltrials(1),data(t).offline.sacc_pburst_ch_align
            end
            
        case 'sacc_rsburst'
            code='sacc_rsburst';
            option='rpt'; %use .rpt file (EyeMove) for sacc
            if t==seltrials(1),data(t).offline.sacc_rsburst_ch_align
            end
    end
    
    %get event and snippet of signal
    %NOTE: there should not be missing event because select only successful trials
    %t
    %data(t).offline.targ_pburst_trial.b_begin
    event_align = get_event(data(t),code,option);
    %pause
% %     data(t).offline
%      lfpfeat=data(t).offline.lfpfeat.lfpfeat_sacc_1;
%      event_align=event_align+round(lfpfeat(1))
% %     %pause
    

    if size(event_align,2)==1 
        trial = get_snippet(signal,data(t),event_align,timevec);
    elseif size(event_align,2)>1 
        trial = get_snippet_ch(signal,data(t),event_align,timevec);
    else
        trial=NaN(nchs,1);
        ntrials_n(it)=ntrials_n(it)+1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %signal filtering
    %     switch signal
    %         case 'fr'
    %             %compute firing rate of spiking activity
    %             for ch=1:nchs
    %                 %gaussian filtering
    %                 %trial(ch,:) = filter_FR(trial(ch,:),'gauss',1000,sigma_FR);
    %
    %                 %epsp filtering
    %                 trial(ch,:) = filter_FR(trial(ch,:),'epsp',1000,1,6);
    %
    %                 %Chronux toolbox
    %                 %[V,t,Err] = evoked(data,Fs,win,width,plt,err)
    %
    %             end
    %
    %         case 'lfp'
    %             %                 %remove 60Hz and harmonics
    %             %                 %bandpass digital filter design
    %             %                 [b60,a60]=butter(2,[(60-5)/500,(60+5)/500],'stop');
    %             %                 %[b60,a60]=butter(8,[(60)/500],'low');
    %             %                 %h = fvtool(b60,a60); % to visualize filter
    %             %                 %[b120,a120]=butter(2,[(120-2.5)/500,(120+2.5)/500],'stop');
    %             %                 %h = fvtool(b120,a120); %to visualize filter
    %             %                 trial60=zeros(size(trial));trial120=zeros(size(trial));
    %             %                 for ch=1:nchs
    %             %                     trial60(ch,:) = filtfilt(b60,a60,double(trial(ch,:)));
    %             %                     %trial120(ch,:) = filtfilt(b120,a120,double(trial(ch,:)));
    %             %
    %             %                     %figure;plot([1:length(trial)],trial(ch,:),'b');hold on
    %             %                     %plot([1:length(trial)],trial60(ch,:),'r');hold on
    %             %                     %%plot([1:length(trial)],trial120(ch,:),'g');hold on
    %             %                     %pause
    %             %
    %             %                     trial(ch,:)=trial60(ch,:);
    %             %
    %             %                 end
    %     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
    %check if NaN values in trial
    s=isnan(trial);
    if sum(s(:))>=1
        ntrials_d(it)=ntrials_d(it)+1;
        %             if dispnan
        %                 figure(fignan);imagesc(~isnan(trial));axis ij; grid;
        %                 colorbar;%('Ticks',[0,1]);%,'TickLabels',{'NaN','Value'});
        %                 xlabel('Time (ms)');set(gca,'ytick', [1:16]);ylabel('Channels');
        %                 title(['trial ' num2str(t) ' contains NaN values.']);
        %                 display(['trial ' num2str(t) ' contains NaN values.']);
        %                 pause(0.1)
        %             end
    end
    
    %creating list of trials (taking into account different sizes)
    if ~isempty(alltrials{it}),alltrialsaux=alltrials{it};else alltrialsaux=[];end;
    ntrials(it)=ntrials(it)+1;
    if ntrials(it)==1, tlen=size(trial,2);else tlen=min(tlen,size(trial,2));end
    
    %remapping of channels using info.chmap
    alltrialsaux(:,ntrials(it),1:tlen)=trial(info.chmap(:),1:tlen);
    %alltrialsaux(:,ntrials(it),1:tlen)=trial(:,1:tlen);
    
   
    %update alltrials
    alltrials{it}=alltrialsaux;
    
    %update look-up table
    lut_trials{it}=[lut_trials{it} t];
    
end


%aligntime (time of beginning of snippet in original signal)
if ~isempty(wind)
    aligntime=abs(min(timevec));
else
    aligntime=1;
end

%info
if dispnan
    display([signal ' ' info.align])
    display(['# of trials per channels: ' num2str(ntrials')])
    display(['# of trials with NaN values per channels: ' num2str(ntrials_d')])
    display(['# of trials with NaN values due to missing event: ' num2str(ntrials_n')])
    display(['Duration of each trial: ' num2str(size(alltrials{1},3))])
end



