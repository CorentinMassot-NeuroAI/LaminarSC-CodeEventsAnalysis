function [ch_ref]=get_chref(nchannels,dref)

%function [ch_ref]=get_chref(nchannels,dref)
%   get ch_rf according to dref from data recorded with a laminar probe (LMA)
%
% use zs_iCSD saved from plot_CSD and compute_iCSD
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 11/05/2016 last modified 11/05/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%number of csd samples (see compute_iCSD)
load zs_iCSD
ncsdsamples=length(zs_iCSD);

%convert dref into channel number (round to closest channel)
%WARNING this is an approximation, best to convert vmis into depths(mm)
ch_ref=round(dref*(nchannels+2)/ncsdsamples);

%%%%
%NOTE: add +1 after finding discrepencies with figure CSD alignment from paper
% however does not change alignment but just 1 channel shift relative to 0
%dref*(nchannels+2)/ncsdsamples
%dref_conv=round(dref*(nchannels+2)/ncsdsamples)+1
%%%&


%MISC
%el_pos = flip(info.depths(info.chmap))*1e-3;
%out_zs = el_pos_with_ends(1):(el_pos_with_ends(N+2)-el_pos_with_ends(1))/(num_out_zs-1):el_pos_with_ends(N+2);
