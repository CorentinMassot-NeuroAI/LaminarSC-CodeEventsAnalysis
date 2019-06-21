function plot_histo(data,type,info,hdlfig,titlestr)

%function plot_histo(data,type,info,hdlfig,titlestr)
%  plot histo of data
%
%
%Corentin Pitt Pittsburgh 11/03/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure
if ~isempty(hdlfig)
    subplot(hdlfig);
else
    figure;
end

%colorlist
colorlist=get_colorlist;

meanval=mean(data);
maxval=max(abs(data));

%edges=[-maxval:0.05:maxval];
edges=[-maxval:1:maxval];
[nhist binhist]=histc(data,edges);
hdl_b=bar(edges,nhist,'histc');
%set(hdl_b(1),'facecolor','k');

line([meanval meanval],[0 max(nhist)]);
%axis([-maxval maxval 0 max(nhist)])
axis([-maxval maxval 0 max(nhist)])
xlabel(type)
grid;grid minor;

if ~isempty(titlestr),
    title(titlestr);
else
    title(info.datafile);
end;


