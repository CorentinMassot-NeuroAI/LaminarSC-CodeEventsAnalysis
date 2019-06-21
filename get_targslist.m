function targslist=get_targslist(data)
%function targslist=get_targslist(data)
% get targets information: list of positions and sort targets.
%
% see also compute_tuning
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 01/03/2017 last modified 01/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Extracting coordinates of target windows
%Tind = 1;% only one target displayed at a time params.Tind;
targ_ind = 2;%second target displayed is target (after fixation) params.Tind; 

%coordinates of all targets
xtargs = [];ytargs = [];
for i=1:length(data)
    %trialNum = str2num(data(i).trialNum);
    posTarg(1) = data(i).targets.window(targ_ind,1)';
    posTarg(2) = data(i).targets.window(targ_ind,2)';
    %Update data: postTarg
    %update_data(0,1,0,data,data_path,info.datafile,'posTarg',posTarg);
    
    targs(i,:)=posTarg;
end

%list of targets without repetitions
targslist = unique(targs,'rows');


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


%MISC
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %target index for each trial
% targindcell = cell(size(targslist,1),1);
% for ut=1:size(targslist,1)
%     targindcell{ut,1} = (xtargs(:,Tind) == targslist(ut,1)) &...
%         (ytargs(:,Tind) == targslist(ut,2));
% end


