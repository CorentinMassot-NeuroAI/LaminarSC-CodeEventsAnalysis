function [csd_prof_r info_r csdsample_ref]=get_csd_prof_aligndepth(csd_prof,dref,info)

%function [csd_prof_aldref ref]=get_csd_prof_aligndepth(csd_prof,dref,info)
%   get csd_prof realigned according to dref from data recorded with a laminar probe (LMA)
%
% see also get_vmis_aligndepth 
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 05/16/2017 last modified 05/19/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ntarg ncsdsamples]=size(csd_prof);

%realign csd_prof
csd_prof_r=NaN(ntarg,ncsdsamples*3);
shift=ncsdsamples/2-dref;
csd_prof_r(:,ncsdsamples+shift+1:2*ncsdsamples+shift)=csd_prof(:,:);

%pause
%csdsample of ref
csdsample_ref=ncsdsamples+ncsdsamples/2;


%update info ncsdsamples
info_r=info;
info_r.ncsdsamples=length(csd_prof_r);

