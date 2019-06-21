function [H]=compute_dprime(data1,data2,level)

%function [H]=compute_dprime(data1,data2,level)
%   compute dprime between 2 sets of data
%
%
%Corentin Pitt Pittsburgh 11/09/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m1=mean(data1);
m2=mean(data2);
std1=std(data1);
std2=std(data2);

dprime=(m1-m2)/(sqrt(0.5*(std1^2+std2^2)));

if abs(dprime)>=level, 
    H=1;
else
    H=0;
end