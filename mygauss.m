function smoothedvector=mygauss(vector,sigma)

%mygauss(vector)
%   fitler with a gaussian filter
%
%Corentin Pitt Pittsburgh 09/29/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gaussian window
gaussfilter = gausswin(sigma/2);
gaussfilter = double(gaussfilter/sum(gaussfilter)); % Normalize.

%figure;plot(gaussfilter)

% filter vector
smoothedvector = conv(double(vector),gaussfilter,'same');

% % plot
% figure
% hold on;
% plot(double(vector), 'b-', 'linewidth', 1);
% plot(smoothedvector, 'r-', 'linewidth', 1);