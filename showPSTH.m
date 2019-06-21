function [psth bounds] = showPSTH(spikeTimes,bounds,sigma,color,display)
%function psth = showPSTH(spikeTimes,bounds,sigma)
%
%  spikeTimes is a cell of spikeTimes (in ms)
%  bounds in ms
%  spikeTimes is a cell array and all values inside are pooled
%
%  sigma is the standard deviation of the gaussian smoothing window
%
% modifed from showPSTH Corentin CNBC 03/17/2013

sampling = 1000;% i.e1/0.001

spikeTimes = reshape(spikeTimes,numel(spikeTimes),1);
ntrials = length(spikeTimes);

%cell2mat
for i = 1:ntrials
    spikeTimes{i} = spikeTimes{i}(:);
end
spikeTimes = cell2mat(spikeTimes);

%bounds
if isempty(bounds)
    if max(spikeTimes)>0,
        bounds=[0 max(spikeTimes)];
    else
        bounds=[0 1];
    end        
end

%psth histogram
psth = hist(spikeTimes,bounds(1):bounds(2));

%smoothing
psth = gaussSmooth(psth,sigma);

%normalization
psth = psth*sampling/ntrials;

%plot
if display,
    if sigma==1
        plot(psth)
    else
        plot(psth,color,'Linewidth',1)
    end
end