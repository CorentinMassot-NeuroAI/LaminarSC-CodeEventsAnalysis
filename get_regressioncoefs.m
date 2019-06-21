function [p resid yfit]=get_regressioncoefs(x,y,n)


%function [p rsq]=get_regressioncoefs(x,y)
%   get regression coefficients
%
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 11/08/2017 last modified 11/08/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%regression spk
p = polyfit(x,y,n);%p(1) is the slope and p(2) is the intercept 

%Residuals
yfit = polyval(p,x);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
resid = 1 - SSresid/SStotal;


