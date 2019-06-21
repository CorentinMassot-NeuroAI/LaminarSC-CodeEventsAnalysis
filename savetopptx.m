function savetopptx(hdlfig,file,figtype,titlename)

%function savetopptx(hdlfig,file,figtype,titlename)
%   save figure to pptx file
%
%
% see also exportToPPTX

% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 12/17/2017 last modified 12/07/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%save tmp figure
saveas(hdlfig,file,figtype);

%create slide
slideNum = exportToPPTX('addslide',hdlfig);
fprintf('Added slide %d\n',slideNum);

%add figure
exportToPPTX('addpicture',file,'scale','maxfixed');
%add title
exportToPPTX('addtext',titlename,'Position',[0 0 10 10],'FontSize',20);
