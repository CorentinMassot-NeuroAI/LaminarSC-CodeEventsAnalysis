%function targslist=compute_targsinfo

%targslist=compute_targsinfo
% get targets information: list of positions and sort targets.
%
% see also compute_tuning
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 01/03/2017 last modified 01/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set paths
[root_path data_path save_path]=set_paths;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get data
datalist=load_data_gandhilab(data_path);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analyzing data
dlist=get_dlist;


info=[];
for d=dlist%1:numel(datalist)
    
    %get data and info
    info.datafile=datalist{d};
    load ([data_path info.datafile]);
    display(info.datafile)
    
    %Extracting coordinates of target windows
    %Tind = 1;% only one target displayed at a time params.Tind;
    targ_ind = 2;%second target displayed is target (after fixation) params.Tind;
    
    %coordinates of all targets
    targs=[];targpos=[];targslist=[];
    for i=1:length(data)
        %trialNum = str2num(data(i).trialNum);
        targpos(1) = data(i).targets.window(targ_ind,1)';
        targpos(2) = data(i).targets.window(targ_ind,2)';
        
        %in case of unsuccessful trial without target1 appearing
        if (targpos(1)==0 & targpos(2)==0)
            targpos=NaN(2,1);
        end
            
        %Update data: targpos
        data(i).offline.targpos=targpos;
        
        %targs
        targs(i,:)=targpos;
  
    end
    
    %list of targets without repetitions
    indextargs=find(~isnan(targs(:,1)) & ~isnan(targs(:,2)));
    targslist = unique(targs(indextargs,:),'rows');
    
    %sort targets index according to polar coordinates
    [th, r]=cart2pol(targslist(:,1),targslist(:,2));
    th=th*180/pi;th(th<0)=th(th<0)+360;
    
    %sort according to radius first
    [rs irs]=sort(r);
    targslist=targslist(irs,:);
    [ths iths]=sort(th(irs));
    targslist=targslist(iths,:);
    
    %sort according to angular position first
    %[ths iths]=sort(th);
    %targslist=targslist(iths,:);
    %[rs irs]=sort(r(iths));
    %targslist=targslist(irs,:);
    
    %Update data: targslist
    update_data(1,1,0,data,data_path,info.datafile,'targslist',targslist);
    
    %check
    %clear('data'); 
    %load ([data_path info.datafile]);
    data(1).offline
    targslist
    clear('data'); 
    %pause
end

%MISC
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %target index for each trial
% targindcell = cell(size(targslist,1),1);
% for ut=1:size(targslist,1)
%     targindcell{ut,1} = (xtargs(:,Tind) == targslist(ut,1)) &...
%         (ytargs(:,Tind) == targslist(ut,2));
% end


