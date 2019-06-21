function [lims]=findlimits(vect)

%function [lims]=findlimits(vect)
%  find the first non-zero extremities defining the limits of a vector 
%
%
%Corentin Pitt Pittsburgh 10/17/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


l=length(vect);
lims=zeros(1,2);
for i=1:length(vect)
i
    if ~isnan(vect(i)) & vect(i)~=0 & lims(1)==0
        lims(1)=i;
    elseif (vect(i)==0 | isnan(vect(i))) & lims(1)~=0 & lims(2)==0,
        lims(2)=i-1;
    elseif (vect(i)~=0 & ~isnan(vect(i))) & lims(1)~=0 & lims(2)==0 & i==length(vect)
        lims(2)=i;
    end
end