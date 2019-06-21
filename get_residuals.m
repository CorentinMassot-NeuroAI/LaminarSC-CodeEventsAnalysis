function [resid stats]=get_regressioncoefs(y,yfit)


%function [p rsq]=get_regressioncoefs(x,y)
%   get regression coefficients
%
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 11/08/2017 last modified 11/08/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Residuals
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
resid = 1 - SSresid/SStotal;

stats.yresid=yresid;
stats.SSresid=SSresid;
stats.SStotal=SStotal;
