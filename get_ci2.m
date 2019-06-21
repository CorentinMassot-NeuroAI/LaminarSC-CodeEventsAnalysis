function meanci=get_ci2(data,ci_level)

%function meanci=get_ci2(data,ci_level)
%   compute ci from data with an operation other than mean
%   and take mean of bootstrapped value as the mean estimation
%
% data: all values obtained through resampling (bootstrapping)
% ci_level: % of confidence interval (e.g. 0.95 for 95%)
%
% Excerpt from "Bootstrap confidence intervals" (Jeremy Orloff and Jonathan Bloom)
% Let’s walk through a summary of the steps needed to answer the question.
% 1.
% Data: x1,...,x272
% 2.
% Data median: xmedian = 240
% ? ??
% 3. Find the median x of a bootstrap sample x1,...,x Repeat 1000 times.
% median 272.
% 4. Compute the bootstrap differences
% ?
% ??
% = x
% median ? xmedian
% Put these 1000 values in order and pick out the .95 and .05 critical values, i.e. the 50th and 950th biggest values. Call these ?.?95 and ??
% .05.
% 5. The bootstrap principle says that we can use ?? and ?? as estimates of ?.95 and ?.05.
% .95 .05
% So our estimated 90% bootstrap confidence interval for the median is
% [xmedian ? ?.?05,xmedian ? ?.?95]

%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 12/28/2017 last modified 12/28/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%remove nan values

data_nnan=data(find(~isnan(data)));

meanci=nan(1,3);
if ~isempty(data_nnan) & length(data_nnan)>0.5*length(data)
    %mean data_nnan
    value=mean(data_nnan);
    %value=median(data_nnan);
    %value_std=std(data_nnan);
    
    %difference bootstrapped values - value original sample
    diffval=data_nnan-value;
    
    %sort difference
    diffval_s=sort(diffval);
    
    %take only ci_level
    alpha=(1-ci_level)/2;
    ind_low=floor(length(diffval_s)*alpha);
    ind_high=ceil(length(diffval_s)*(1-alpha));
    
    %get CI from sorted differences
    %NOTE: add instead of - because negative value
    meanci(1)=value;
    if meanci(1)<0
        meanci(2)=value+diffval_s(ind_low);
        meanci(3)=value+diffval_s(ind_high);
    else
        meanci(2)=value-diffval_s(ind_low);
        meanci(3)=value-diffval_s(ind_high);
    end
    %meanci(2)=value+value_std
    %meanci(3)=value-value_std
end

%MISC
% %remove nan values
% data_nnan=data(find(~isnan(data)));
%
% %alpha
% alpha=1-ci_level;
%
% %T-score for level
% ts = tinv([1-alpha/2],length(data_nnan)-1);
%
% %mean data_nnan
% mean_data_nnan=mean(data_nnan);
%
% %Standard Error
% sem = std(data_nnan)/sqrt(length(data_nnan));
%
%
% %mean and confidence intervals
% meanci = [mean_data_nnan ts*sem];
