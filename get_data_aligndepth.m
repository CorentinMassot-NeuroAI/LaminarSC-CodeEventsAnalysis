function [data_r info_r ch_ref dref_conv]=get_data_aligndepth(data,dref,info,coef)

%function [data_r info_r ch_ref dref_conv]=get_data_aligndepth(data,dref,info,coef)
%   get data matrix realigned according to dref from data recorded with a laminar probe (LMA)
%
% use zs_iCSD saved from plot_CSD and compute_iCSD
%
% coef= [] fort data vector -1 for image
%
%see also get_vmis_aligndepth
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 01/09/2018 last modified 01/09/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nchannels ncols]=size(data);

%number of csd samples (see compute_iCSD)
load zs_iCSD
ncsdsamples=length(zs_iCSD);

%coef vector or data
if isempty(coef)
    coef=1;
end


%convert dref into channel number (round to closest channel)
%WARNING this is an approximation, best to convert data into depths(mm)
dref_conv=round(dref*(nchannels+2)/ncsdsamples);

%%%%
%NOTE: add +1 after finding discrepencies with figure CSD alignment from paper
% however does not change alignment but just 1 channel shift relative to 0
%dref*(nchannels+2)/ncsdsamples
%dref_conv=round(dref*(nchannels+2)/ncsdsamples)+1
%%%&

%realign data
data_r=NaN(nchannels*3,ncols);
shift=coef*(nchannels/2-dref_conv);
data_r(nchannels+shift+1:2*nchannels+shift,:)=data;

%pause
%channel of ref
ch_ref=nchannels+nchannels/2;


%update info nchannels
info_r=info;
info_r.nchannels=size(data_r,1);

%MISC
%el_pos = flip(info.depths(info.chmap))*1e-3;
%out_zs = el_pos_with_ends(1):(el_pos_with_ends(N+2)-el_pos_with_ends(1))/(num_out_zs-1):el_pos_with_ends(N+2);
