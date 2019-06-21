function event_align=get_eventalign(events,align)

%function event_alignign=get_eventalign(events,align)
%   get event timing corresponding to the chosen alignement
%
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 08/22/2016 last modified 08/22/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch align
    case 'no'
        event_align=0;       
    case 'targ'
        event_align=events.targ;
    case 'go'
        event_align=events.go;
    case 'sacc'
        event_align=events.sacc;
    case 'fixtarget'
        event_align=events.fixtarget; 
end


