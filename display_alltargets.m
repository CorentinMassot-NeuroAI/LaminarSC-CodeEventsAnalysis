function display_alltargets(targs_list,info,hdlfig)

%function targs_ind=display_alltargets(targs_list,info,hdlfig)
%   display all targets in targs_list
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 07/06/2016 last modified 07/06/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%figure

if ~isempty(hdlfig)
    subplot(hdlfig);
else
    figure;hold on;
end

ntargs=size(targs_list,1);
plot(targs_list(:,1),targs_list(:,2),'o','MarkerFaceColor','k');
for it=1:ntargs
    text(targs_list(it,1),targs_list(it,2),['  ' num2str(it)]);
end
xmax=max(targs_list(:,1));ymax=max(targs_list(:,2));
%xmax=xmax+0.1*xmax;ymax=ymax+0.1*ymax;
xmax=2*xmax;ymax=2*ymax;
axis([-xmax xmax -ymax ymax]);
line([0 0] ,[-ymax ymax]);
line([-xmax xmax] ,[0 0]);
xlabel('Horizontal positions');ylabel('Vertical positions');
%title({info.datafile ; 'Targets'})
title('Targets')
grid

