function [inflect,resnorm,resnorm_fit,fitparams]=get_inflection_2pwlr(trial,anchor_b,anchor_e,min_bins)

%function [inflect,resnorm,resnorm_fit,fitparams]=get_inflection_2pwlr(trial,anchor_b,anchor_e)
%   get inflection point of activity in a trial between anchor_b and
%   anchor_e
%
% use parfor and matlabpool
% to get pool for 4 hours type in terminal: myPool = parpool('IdleTimeout', 4*60); 
% to close pool: delete(gcp('noCreate'));
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 11/30/2017 last modified 11/30/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%debug
debug=0;

%two piece-wise linear fit

trial_sel=trial(anchor_b:anchor_e);

%linear function
fun = @(x,xdata)x(1)+x(2)*xdata;
options=optimset('Display','off');

if debug,
    figtrial=figure;plot(trial_sel);
    hdlreg=figure;hold on;
end;

reg_min=min_bins;%min nb time bins for each regression
reg_max=length(trial_sel)-reg_min;
ln=length(reg_min:reg_max);
resnorm_list1=nan(ln,1);
resnorm_list2=nan(ln,1);
fitparams_list1=nan(ln,2);
fitparams_list2=nan(ln,2);
parfor k=reg_min:1:reg_max
%for k=reg_min:reg_max
    
    %linear piece 1
    trial_p1=trial_sel(1:k);
    xp1=[1:length(trial_p1)];
    %x0p1=[xp1(1) trial_p1(1)];
    x0p1=[xp1(end) trial_p1(end)];%NOTE: changing x0p1 does not change fit
    %free regression
    [fitparams1 resnorm1] = lsqcurvefit(fun,x0p1,xp1,trial_p1,[],[],options);
    %     %fix last point
    %     display('2PWLR using fix point!')
    %     p = polyfix(xp1,trial_p1,1,x0p1(1),x0p1(2));%use polyfix function
    %     yp1=polyval(p,xp1);
    %     %compute resnorm and fitparams
    %     resnorm1=sum((yp1-trial_p1).^2);
    %     fitparams11=yp1(1);
    %     fitparams12=(yp1(end)-yp1(1))/(xp1(end)-xp1(1));
    
    
    %linear piece 2
    trial_p2=trial_sel(k:end);
    xp2=[1:length(trial_p2)];
    %x0p2=[xp2(1) trial_p2(1)];
    x0p2=[xp2(1) trial_p2(1)];
    %free regression
    [fitparams2 resnorm2] = lsqcurvefit(fun,x0p2,xp2,trial_p2,[],[],options);
    %     %fix first point
    %     display('2PWLR using fix point!')
    %     p = polyfix(xp2,trial_p2,1,x0p2(1),x0p2(2));%use polyfix function
    %     yp2=polyval(p,xp2);
    %     %compute resnorm and fitparams
    %     resnorm2=sum((yp2-trial_p2).^2);
    %     fitparams21=yp2(1);
    %     fitparams22=(yp2(end)-yp2(1))/(xp2(end)-xp2(1));
    
    
%     if debug
%         figure(hdlreg)
%         subplot(1,2,1);hold on;
%         plot(xp1,trial_p1,'k-',xp1,fun(fitparams1,xp1),'r-',x0p1(1),x0p1(2),'ro')
%         %plot(xp1,fun(fitparams1,xp1),'g-',x0p12(1),x0p12(2),'go')
%         %plot(xp1,yp1,'g-',x0p1(1),x0p1(2),'go')
%         title(['angle: ' num2str(atan(fitparams1(2))*180/pi) 'deg'])
%         axis([0 50 -20 200])
%         hold off
%         subplot(1,2,2);hold on;
%         plot(xp2,trial_p2,'k-',xp2,fun(fitparams2,xp2),'b-',x0p2(1),x0p2(2),'bo')
%         %plot(xp2,yp2,'g-',x0p2(1),x0p2(2),'go')
%         title(['angle: ' num2str(atan(fitparams2(2))*180/pi) 'deg'])
%         axis([0 100 -20 200])
%         hold off
%         pause
%     end
%     
    resnorm_list1(k)=resnorm1;
    resnorm_list2(k)=resnorm2;
    fitparams_list1(k,:)=fitparams1;
    fitparams_list2(k,:)=fitparams2;
    %     %using fix point
    %     fitparams_list1(k,:)=[fitparams11 fitparams12];
    %     fitparams_list2(k,:)=[fitparams21 fitparams22];

    
end
resnorm_list=[resnorm_list1  resnorm_list2];
fitparams_list=[fitparams_list1  fitparams_list2];

%replace zeros by nans
resnorm_list(find(resnorm_list(:,:)==0))=nan;


%find min residuals
inflect=nan(1,3);
resnorm=nan(1,1);
resnorm_fit=nan(1,2);
fitparams=nan(1,4);
if sum(~isnan(resnorm_list(:,1)))>0 & sum(~isnan(resnorm_list(:,2)))>0 %trial_sel was long enough
    
    resnormaux=resnorm_list(:,1)+resnorm_list(:,2); 
    resnorm=nan(1,anchor_b-1+length(resnormaux));
    resnorm(1,anchor_b:end)=resnormaux';
    
    [k_v k_ind]=nanmin(resnorm);

    inflect(1)=k_ind;
    inflect(2)=trial_sel(k_ind-anchor_b+1);
    inflect(3)=k_v;

    resnorm_fit=resnorm_list(k_ind-anchor_b+1,:);
    fitparams=fitparams_list(k_ind-anchor_b+1,:);
    %size(fitparams)
    
    if debug
        figure(figtrial)
        hold on;
        %subplot(1,2,1);hold on;
        trial_p1=trial_sel(1:k_ind-anchor_b+1);
        xp1=[1:length(trial_p1)];
        x0p1=[xp1(end) trial_p1(end)];
        %plot(xp1,trial_p1,'k-',xp1,fun(fitparams(1:2),xp1),'r-',x0p1(1),x0p1(2),'ro')
        plot(xp1,fun(fitparams(1:2),xp1),'b-',x0p1(1),x0p1(2),'ko')
        %title(['angle: ' num2str(atan(fitparams1(2))*180/pi) 'deg'])
        %axis([0 50 -20 200])
        %hold off
        
        %subplot(1,2,2);hold on;
        trial_p2=trial_sel(k_ind-anchor_b+1:end);
        xp2=[1:length(trial_p2)];
        x0p2=[xp2(1) trial_p2(1)];
        %plot(xp2,trial_p2,'k-',xp2,fun(fitparams(3:4),xp2),'b-',x0p2(1),x0p2(2),'bo')
        plot(xp2+k_ind-anchor_b,fun(fitparams(3:4),xp2),'r-',x0p2(1)+k_ind-anchor_b,x0p2(2),'ko')
        %title(['angle: ' num2str(atan(fitparams2(2))*180/pi) 'deg'])
        %axis([0 100 -20 200])
        %hold off
        pause
    end
    
    
end



%%
%MISC
% %%
% %spline fit
% % Breaks interpolated from data
% y=trial(latency(1):sig_minmax(1));
% x=[1:length(y)]+latency(1)-1;
% pp1 = splinefit(x,y,1,10);  % 11 breaks, 10 pieces
%
%
% % Plot
% figure(1)
% y1 = ppval(pp1,x);
% plot(x,y,'-',x,y1);
% plot(x(1:end-1),diff(y)*10,'b--');
% plot(x(1:end-1),diff(y1)*10,'r--');
% %axis([0,2*pi,-2.5,2.5]), grid on
% %legend('data','41 breaks, 40 pieces','11 breaks, 10 pieces')
% %title('EXAMPLE 1: Breaks and pieces')

