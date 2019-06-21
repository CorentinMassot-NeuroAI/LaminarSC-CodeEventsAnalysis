function [r_output p_output]=plot_regress(data1,data2,type,ind,info,hdlfig,titlestr)

%function r=plot_regress(data1,data2,type,texth,info,hdlfig,titlestr)
%  plot regession between 2 sets of data
%
% ind: index of datafile
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 11/03/2016 last modified 06/21/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure
if ~isempty(hdlfig)
    subplot(hdlfig);
    hold on;
else
    figure;
end


%colorlist
colorlist=get_colorlist;

%NaN values
vals=find(~isnan(data1) & ~isnan(data2));
data1=data1(vals);
data2=data2(vals);

%plot
%plot(data1,data2,'-','color',colorlist(ind,:));
for ch=1:size(data1,2)
    %plot(data1(ch),data2(ch),'o','MarkerSize',5,'MarkerFaceColor',colorlist(ch,:));
    plot(data1(ch),data2(ch),'o','MarkerSize',5,'MarkerFaceColor',colorlist(ind,:));
end
%plot(data1(1),data2(1),'o','MarkerSize',10,'MarkerFaceColor',colorlist(ind,:));
%plot(data1(end),data2(end),'o','MarkerSize',10,'MarkerEdgeColor',colorlist(ind,:));



% %regression
r_output=[];
p_output=[];

[p,polystats] = polyfit(data1,data2,1);
yfit=p(2)+data1*p(1);
plot(data1,yfit,'color',colorlist(ind,:));

%                 SSresid = sum((data2 - yfit).^2);
%                 SStotal = (length(data2)-1) * var(data2);
%                 display('residual of responses:')
%                 rsq = 1 - SSresid/SStotal

%correlation coefficients also used for stats analysis
[r,pstats]=corrcoef(data1,data2);
%output
r_output=r(1,2);
p_output=pstats(1,2);

%display
maxplot=max(max(data1(:)),max(data2(:)))+10;
minplot=min(min(data1(:)),min(data2(:)))-10;

text(minplot+1,maxplot-ind*2,['r=' num2str(round(r_output,3)) ' p=' num2str(p_output)]);

hl=line([minplot maxplot] ,[minplot maxplot]);
set(hl,'Color','k','LineStyle','--','Linewidth',1);
xlabel([type ' spk']);ylabel([type ' lfp']);
axis([minplot maxplot minplot maxplot])
%grid;grid minor;


if ~isempty(titlestr),
    title(titlestr);
else
    title(info.datafile);
end;


