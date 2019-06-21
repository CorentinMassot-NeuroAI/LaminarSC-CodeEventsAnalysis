function [align]=get_features_alignlist(feature)

%[align]=get_features_alignlist(feature)
%   get all alignment dpeending on feature to be analyzed
%
% see also analysis_spklfp_munoz

% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 09/18/2017 last modified 09/18/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch feature
    case 'buildup'
        align='sacc';
    case ''