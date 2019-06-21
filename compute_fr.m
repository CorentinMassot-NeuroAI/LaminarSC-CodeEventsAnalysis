%function compute_fr

%function compute_fr
%   Compute firing rate from laminar data
%
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 03/08/2017 last modified 03/08/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set paths
[root_path data_path save_path]=set_paths;

%screen size
scrsz = get(groot,'ScreenSize');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
%print figures and save data
savedata=0;
savefigs=0;
figtype='epsc2';%'png';%'epsc2';

%alignement
%alignlist={'no' 'targ' 'go' 'sacc'};
alignlist={'no'};

%gaussian filtering
sigma_FR=6;
filt=['gauss_' num2str(sigma_FR)];


%epsp filtering
param2=20;%6; %20ms drecrease rate is the right value (See Thompson et al. Schall 1996)
filt=['epsp_' num2str(param2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get data
datalist=load_data_gandhilab(data_path);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analyzing data
dlist=get_dlist

for d=dlist%1:numel(datalist)
    data=[];
    %get data and info
    info.datafile=datalist{d};
    load ([data_path info.datafile]);
    display(info.datafile)
    
    %getting channel mapping and discard selected bad channels
    [info.chmap info.nchannels info.depths]=get_chmap(data(1).info.electrode{2},[]);
    
    info.align=alignlist{1};
    
    %get all trials
    for t=[1:numel(data)]
        
        %get spiketimes for the whole trial
        size_trial=size(data(t).lfp,2);
        trialsize=size_trial+20;

        trial=zeros(info.nchannels,trialsize);
        
        for ch=1:info.nchannels
            trial_spkt=zeros(1,trialsize);
            spiketimes=round(1000*data(t).spikeTimestamps{ch});
            spiketimes_sel=spiketimes(find(spiketimes~=0 & spiketimes<=size_trial));
            trial_spkt(spiketimes_sel)=1;
       
            %compute firing rate of spiking activity
            %gaussian filtering
            %trial(ch,:) = filter_FR(trial_spkt(:),'gauss',1000,sigma_FR);
            
            %epsp filtering
            trial(ch,:) = filter_FR(trial_spkt(:),'epsp',1000,1,param2);
            
            %Chronux toolbox
            %[V,t,Err] = evoked(data,Fs,win,width,plt,err)
            
        end
        
        %Update data
        eval(['data(t).offline.fr_' filt '=trial;']); 
        
    end
    
    %save updated data
    update_data(1,0,0,data,data_path,info.datafile,[],[]);
end

