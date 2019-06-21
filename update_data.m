function data=update_data(savedata,updatedata,loaddata,data,data_path,datafile,field,update)

%data=update_data(savedata,updatedata,loaddata,data,data_path,datafile,field,update)
%
%   update data
%   if field and update are indicated then each trial can be updated with
%   the same information
%   if field==[], then each trial can updated individually depending on
%   the value of data
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 10/16/2016 last modified 10/16/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if updatedata,
    display(['updating data: offline.' field])
    
    %load data
    if loaddata
        display('loading data...')
        load ([data_path datafile]);
    else
        if isempty(data)
            display('Warning: data is empty')
            pause
        end
    end
        
    %update field
    if ~isempty(field),
        for d=1:numel(data)
            eval(['data(d).offline.' field '=update;']);
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save
if savedata
    display(['saving updated data: offline.' field])
    save([data_path datafile], 'data')
end