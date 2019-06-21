function datalist=load_data_gandhilab(data_path)

%analysis_trials_LMA
%   Analysis of data tirals-by-trials recorded with a laminar probe (LMA)
%   in gandhilab
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 11/03/2015 last modified 11/03/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%read directory
files=dir(data_path);

datalist={};
dateslist=[];
ff=0;
for f=1:size(files),
    file=files(f).name;
    if length(file)>3,
        ff=ff+1;
        dateslist(ff)=str2num([file(17:18) file(13:16)]);
        datalist{ff}=file;
    end
end

%sort files according to date
[~,ind]=sort(dateslist);
datalist=datalist(ind);


% %get data SFN2015
% datalist={
%     'bl_sc_031115_mcell_spikelfp_cSC'; %
%     'bb_sc_031315_mcell_spikelfp_cSC'; %lowfreq lfps
%     'bb_sc_070915_mcell_spikelfp_cSC'; %low ok
%     %'bb_sc_071215_mcell_spikelfp_cSC'; %lowfreq lfps
%     'bl_sc_071415_mcell_spikelfp_cSC';
%     'bl_sc_072315_1_mcell_spikelfp_cSC';
%     'bl_sc_072315_2_mcell_spikelfp_cSC';
%     'bb_sc_080415_mcell_spikelfp_cSC'; %low good
%     }


%get data Data_SC
% datalist={
%     'bl_lSCTrack_033116_DelaySacc-A';
%     'bl_lSCTrack_033116_DelaySaccMem-A';
%     'bl_lSCTrack_033116_DelaySacc-A_vl';
%     'bl_lSCTrack_033116_DelaySaccMem-A_vl';
%     'bl_lSCTrack_033116_DelaySacc-8T-R';
% 
%     'bl_lSCTrack_033016_DelaySacc-A';
%     'bl_lSCTrack_033016_DelaySaccMem-A';
%     'bl_lSCTrack_033016_DelaySacc-A_vl';
%     'bl_lSCTrack_033016_DelaySaccMem-A_vl';
%     'bl_lSCTrack_033016_DelaySacc-8T-R';
%    }