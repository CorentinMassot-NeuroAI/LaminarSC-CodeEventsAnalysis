%function compute_tuning

%function compute_tuning
%   Compute tuning from laminar data
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 10/16/2016 last modified 01/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%set paths
[root_path data_path save_path]=set_paths;

%screen size
scrsz = get(groot,'ScreenSize');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
%alignement
%alignlist={'no' 'targ' 'go' 'sacc'};
alignlist={'sacc'};

%window of analysis
%wind=[100 170];%targ
wind=[-25 50];%sacc

%sigma FR
sigma_FR=6;

%display
disp=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get datalist
datalist=load_data_gandhilab(data_path);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analyzing data
dlist=get_dlist;

data=[];
info=[];
for d=dlist%1:numel(datalist)
    
    %get data and info
    info.datafile=datalist{d};
    load ([data_path info.datafile]);
    display(info.datafile)
    
    if disp
        %getting channel mapping and discard selected bad channels
        [info.chmap info.nchannels info.depths]=get_chmap(data(1).info.electrode{2},[]);
    end
    
    %getting trial type
    info.trialtype=data(1).sequence(1);
    %getting list of targets
    targslist=data(1).offline.targslist;
    %targets index
    targs_ind=get_targsindex(targslist,info);
        
    %select trials
    seltrials=get_seltrials(data,'rpt');
    
    
    if disp
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %get spk data for tuning
        info.align=alignlist{1};       
        [alltrials_spk_tuning info.aligntime]=get_alltrials_align(data,seltrials,[],'fr',info,targslist,sigma_FR,1);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %display all targets
    figtrials=figure('Position',[scrsz(1)+300 300 scrsz(3)-500 scrsz(4)-500]);
    hdlfig=subplot(1,2,1);hold on;
    display_alltargets(targslist,info,hdlfig);
    
    if disp
        %compute target tuning
        hdlfig=subplot(1,2,2);hold on;
        titlestr={info.datafile ; ['Target tuning ' info.align ]};
        plot_targtuning(alltrials_spk_tuning,targs_ind,info,hdlfig,titlestr);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %selection of target tuning
    %yn=[];
    %while ~strcmp(yn,'y')
    targ_tuning=input(['select tuning target (' num2str(length(targslist)) ' total):']);
    %targ_tuning=6
    
    %yn=input('sure (y/n)?','s');
    %end
    display(['Selected target tuning: ' num2str(targ_tuning)])
    
    
    %%
    %%%%%%%%%%%%%%%%%%
    %Update data: tuning
    update_data(1,1,0,data,data_path,info.datafile,'targ_tuning',targ_tuning);
    
    %close(figtrials)
end

