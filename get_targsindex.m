function targs_ind=get_targsindex(targs_list,info)

%function targs_ind=get_targsindex(targs_list,info)
%   get targets index sorted according to x position
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created  07/06/2016 last modified 01/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NOTE: targs_ind is the index of the argets as stored in targslist
%then targs_ind is organized according to the order chosen in this program
% for example targs_ind can be [4 5 6 1 2 3] to follow the order chosen


%target index
if strcmp(info.trialtype,'DelaySacc-8T-R');targs_ind=[1:8];
elseif strcmp(info.trialtype,'DelaySacc-8T-L');targs_ind=[1:8];
else
    %sort according to x position
    [val targs_ind]=sort(targs_list(:,1));
    targs_ind=targs_ind';
end

%display(['Index of targets: ' num2str(targs_ind)])

