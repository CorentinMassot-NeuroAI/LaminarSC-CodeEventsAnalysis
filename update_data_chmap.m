%function update_data_chmap

%update_data_chmap
%   update chmap data
%
%Corentin Pitt Pittsburgh 10/14/2015 10/16/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parameters

%paths
data_path=['C:\Users\massotc@upmc.edu\Work\NeuroPITT\Data\SFN2015Data_aligned\'];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get LFP data 
datalist={
    'bl_sc_031115_mcell_spikelfp_cSC';
    'bb_sc_031315_mcell_spikelfp_cSC'; %lowfreq lfps
    'bb_sc_070915_mcell_spikelfp_cSC';
    'bb_sc_071215_mcell_spikelfp_cSC'; %lowfreq lfps
    'bl_sc_071415_mcell_spikelfp_cSC';
    'bl_sc_072315_1_mcell_spikelfp_cSC';
    'bl_sc_072315_2_mcell_spikelfp_cSC';
    'bb_sc_080415_mcell_spikelfp_cSC'; 
    }


chmaplist{1}=flip([9:13,15:16,1:6,8]);
chmaplist{2}=flip([9:16,1:6,8]);
chmaplist{3}=flip([9:13,15:16,1:6,8]);
chmaplist{4}=flip([9:16,1:6,8]);
chmaplist{5}=flip([9:13,15:16,1:6,8]);
chmaplist{6}=flip([9:16,1:8]);
chmaplist{7}=flip([9:16,1:8]);
chmaplist{8}=flip([9:16,1:8]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analyse LFPs

for d=1:numel(datalist)
    load ([data_path datalist{d}]);
    display(datalist{d})
    for i=1:numel(data)
    data(i).chmap = chmaplist{d};
    end
  save([data_path datalist{d}], 'data')
    
end
    
