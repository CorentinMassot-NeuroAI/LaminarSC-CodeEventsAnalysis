function fitangles=get_fitangles(fitparams)

%function fitangles=get_fitangles(fitparams)
%   compute angles from ftiparams obtained from a linear fit
%
%
% see also get_inflection_2pwlr
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 12/29/2017 last modified 12/29/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%loop on list of fitparams
nboot=size(fitparams,1);
fitangles=nan(nboot,1);
for b=1:nboot
    fitangles(b)=atan(fitparams(b,2))*180/pi;
end
