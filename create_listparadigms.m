%function create_listparadigms(dlist)

%create_listparadigms(dlist)
%  cretae list of indices for all the different paradigms in dlist
%
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 01/22/2017 last modified 01/22/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set paths
[root_path data_path save_path]=set_paths;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get data
datalist=load_data_gandhilab(data_path);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analyzing data
dlist=get_dlist;

data=[];info=[];

vg=[];
mg=[];
vg_A=[];
mg_A=[];
step=[];
gap=[];
for d=dlist
    %get data and info
    info.datafile=datalist{d};
    load ([data_path info.datafile]);
    display(info.datafile)
    
    %getting trial type
    trialtype=char(data(1).sequence{1});
    trialtype_supp=char(data(1).sequence{4});
    
    switch trialtype
        case 'DelaySacc'
            vg=[vg d];
        case 'NoPuff_Overlap'
            vg=[vg d];
        case 'DelaySaccMem'
            mg=[mg d];
        case 'DelaySacc-A'
            vg_A=[vg d];
        case 'DelaySaccMem-A'
            mg_A=[mg d];
        case 'Step_basic'
            step=[step d];
        case 'Gap'
            if isempty(trialtype_supp)
                gap=[gap d];
            end
    end

end

vg
mg
step
gap



