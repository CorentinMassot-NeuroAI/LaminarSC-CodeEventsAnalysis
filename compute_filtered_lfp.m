%function compute_filtered_lfp

%function compute_filtered_lfp
%   Compute a filtered version of the lfp from laminar data
%
% see also compute_lfp
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 05/12/2018 last modified 05/12/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set paths
[root_path data_path save_path]=set_paths;

%screen size
scrsz = get(groot,'ScreenSize');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters

%filters of LFP signal
filt=['low50Hz'];%low100Hz
 
srate=1000/2
[b_low,a_low]=butter(4,(50)/srate,'low');        % 4 BW (Ripple)
% %h = fvtool(b_lfp,a_lfp);                 % Visualize filter
% [b_60,a_60]=butter(2,[57.5/srate,62.5/srate],'stop');        % Bandpass digital filter design
% %h = fvtool(b_60,a_60);                 % Visualize filter
% [b_120,a_120]=butter(2,[(120-2.5)/srate,(120+2.5)/srate],'stop');        % Bandpass digital filter design
% %h = fvtool(b120,a120);                 % Visualize filter


%display
disp=0;

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
    
    %get all trials
    trial_lfp=[];
    for t=[1:numel(data)]
        t
        trial_lfp=data(t).lfp;
        
        trialsize=size(trial_lfp,2);
        trial_lfpf=zeros(info.nchannels,trialsize);
        for ch=1:info.nchannels
            
            %filter lfp
            trial_chf=trial_lfp(ch,:);
            %cutoff
            trial_chf = filter(b_low,a_low,trial_chf);
            %             %60Hz
            %             trial_chf = filter(b_60,a_60,trial_chf);
            %             %120Hz
            %             trial_chf = filter(b_120,a_120,trial_chf);
     
            trial_lfpf(ch,:)=trial_chf;
            
            if disp
                fig=figure;hold on
                plot(trial_lfp(ch,:),'b');
                plot(trial_lfpf(ch,:),'r');
                pause
                close(fig)
            end
        end
        
        %Update data
        %WARNING: change filt name!!!
        eval(['data(t).offline.lfp_' filt '=trial_lfpf;']); 
        
    end
    
    %save updated data
    %display('NOT SAVED!')
    update_data(1,0,0,data,data_path,info.datafile,[],[]);
    
    %pause
end

