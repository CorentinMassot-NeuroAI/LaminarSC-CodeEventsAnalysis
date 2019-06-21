function [chmap nchannels depths]=get_chmap(channelOrder,discard)

%function get_chmap(channelOrder,discard)
%  get channel mapping for LMA probe recording and discard bad channels for
%  the session
%
% NOTE: first channel in chmap is the DEEPEST contact on the probe 
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 11/03/2015 last modified 01/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%WARNING: hack for keep analyzing bad channels
if length(channelOrder)<16
channelOrder=[1:16]
end

%estimated depths of each channel in mm
depths_all=[2.35:-0.15:0.1];
depths=depths_all(channelOrder);



%chmap is channelOrder if all channels are good and discard is []
chmap=channelOrder;


%discarding bad channels
indexes=[1:length(chmap)];
for b=1:length(discard)
    indexes=indexes(find(indexes~=discard(b)));
end
if ~isempty(discard),display(['bad channel: ' num2str(discard)]);end;

%chmap 
chmap=chmap(indexes);

%number of channels
nchannels=length(chmap);

%depths
depths=depths_all(indexes);

