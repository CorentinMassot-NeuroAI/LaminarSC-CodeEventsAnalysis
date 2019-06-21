function [H p]=get_testsignif(data1,data2,p_level)

%function meanci=get_testsignif(data1,data2,p_level)
%   compute significance test between two datasets 
%  originally for comparing two sets of bootstrapped estimations
%
% see also analysis_onset_buildup
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 01/18/2018 last modified 01/18/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind=find(~isnan(data1));
data1_aux=data1(ind);
[ind data2_aux]=find(~isnan(data2));
data2_aux=data2(ind);


if ~isempty(data1_aux) & ~isempty(data2_aux)
    %normality test (Kolmogorov-Smirnov)
    [h_ks1 p_ks1] = kstest(data1_aux);
    [h_ks2 p_ks2] = kstest(data2_aux);
    
    
    %Mann-Whitney test (ranksum)
    [p h stats] = ranksum(data1_aux,data2_aux);
    
    %H
%     h_ks1
%     h_ks2
    H=p<p_level;
    
else
    H=nan;
    data1_aux'
    data2_aux'
%    pause
end

