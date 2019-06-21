function [trials_csd depths_ch]=get_csdchannels(csd,info)

%function [trials_csd depths_ch]=get_csdchannels(csd,info)
%   get channles equivalent from CSD plot from data recorded with a laminar probe (LMA)
%
% use zs_iCSD saved from plot_CSD and compute_iCSD
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 09/06/2016 last modified 09/06/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of csd samples (see compute_iCSD)
load zs_iCSD
ncsdsamples=length(zs_iCSD);

%convert channel into depth (round to closest channel)
%%%%%%%%%%
%WARNING this is an approximation, best to convert vmis into depths(mm)
%NOTE: add +1 or 0.5 after finding discrepencies with figure CSD alignment from paper
% however does not change alignment but just 1 channel shift relative to 0
%dref*(nchannels+2)/ncsdsamples
%dref_conv=round(dref*(nchannels+2)/ncsdsamples)+1
%%%%%%%%%%
depths_ch=[];
for ch=1:info.nchannels+2
    %depths_ch=round(dref*(nchannels+2)/ncsdsamples);
    %depths_ch(info.nchannels+2+1-ch)=round((ch)*ncsdsamples/(info.nchannels+2));
    %center better
    depths_ch(info.nchannels+2+1-ch)=round((ch-0.5)*ncsdsamples/(info.nchannels+2));
end


%extract csd at each depth
for ch=1:info.nchannels
    trials_csd(ch,:)=csd(depths_ch(ch+1),:);%compensate for additional channel from csd
end


%MISC
%el_pos = flip(info.depths(info.chmap))*1e-3;
%out_zs = el_pos_with_ends(1):(el_pos_with_ends(N+2)-el_pos_with_ends(1))/(num_out_zs-1):el_pos_with_ends(N+2);
