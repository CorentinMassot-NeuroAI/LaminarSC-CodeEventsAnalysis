function [H]=get_testsignifci(dataci1,dataci2,overlap)

%function [H]=get_testsignifci(data1,data2,p_level)
%   compute significance test between two datasets using their confidence
%   intervals
%  originally for comparing two sets of bootstrapped estimations
%
% size(dataci1)=1 3 (mean / mean+ci above / mean-ci below) averaged from
% bootstrapped data
%
% overlap: amount of time bin that ci can overlap and still being
% considered significantrly different.
%
% H=0 if no singificance
% H=1 if singificant difference and mean of dataci1 < mean of dataci2 
% H=-1 if singificant difference and mean of dataci1 > mean of dataci2 
%
% see also get_testsignif analysis_onset_buildup
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 02/14/2018 last modified 02/14/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dataci1
%dataci2

dataci1=sort(dataci1);
dataci2=sort(dataci2);

%sign of H
alpha=0;

if dataci1(2)<=dataci2(2)
    ci1=dataci1(3);
    ci2=dataci2(1);
    alpha=1;
else
    ci1=dataci2(3);
    ci2=dataci1(1);
    alpha=-1;
end    

%ci1
%ci2

if ~isempty(ci1) & ~isempty(ci2) & ~isnan(ci1) & ~isnan(ci2) 

    %if ci1<ci2+overlap*ci1/100
    if ci1-ci2<overlap*(ci1+ci2)/100
        
        H=alpha*1;
    else
        H=0;
    end
else
    H=nan;
end

