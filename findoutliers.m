function  [noutliers outliers i_out ]=findoutliers(data);

% [noutliers i_out outliers]=findoutliers(data)
%
%   find outliers in vector of values
%   noutliers: vector without outliers
%   outliers: vector of outliers
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 06/01/2017 last modified 06/06/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%remove nan values
nanvals=isnan(data);
data=data(find(~nanvals));

%median
median_data=median(data);


%median absolute deviation MAD
mad= median(abs(data-median_data))


%find outliers (>3*mad)
coef=10;
%outliers and abs(outliers) have to be >1
%i_out=find(abs(data)>1);
i_out=find((data>median_data+coef*mad | data<median_data-coef*mad)& abs(data)>1);
outliers=data(i_out);

%noutliers
%i_nout=find(abs(data)<=1);
i_nout=find((data<=median_data+coef*mad & data>=median_data-coef*mad) & abs(data)<=1);
noutliers=data(i_nout);