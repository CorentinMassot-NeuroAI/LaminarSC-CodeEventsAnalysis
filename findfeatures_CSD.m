function  [csdfeat csd_prof]=findfeatures_CSD(find_feat,display,auto,csd,zs,info,hdlfig);

% [csdfeat csd_prof]=findfeatures_CSD(find_feat,display,auto,csd,zs,info.hdlfig);
%
%   find features in CSD plot
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 10/15/2016 last modified 10/15/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%region depending on align
switch info.align
    case 'targ'
        region=[info.aligntime+140,1,150,length(zs)-1];
    case 'sacc'
        region=[info.aligntime+5,1,20,length(zs)-1];
end

%get patch of CSD
%csd_patch=flip(csd(region(2):region(2)+region(4),region(1):region(1)+region(3)));%using plot_CSD
csd_patch=(csd(region(2):region(2)+region(4),region(1):region(1)+region(3)));%using plot_csdtrials

%get CSD profile
%csd_prof=mean(csd_patch,2);%using plot_CSD
csd_prof=mean(flip(csd_patch),2); %using plot_csdtrials


if display
    %plot csd and region
    figcsdfeat=figure;
    % subplot(1,2,1);hold on;
    % plot_CSD(flip(csd),[1:d],1,1,0)
    % hdlrect=rectangle('Position', region);
    % set(hdlrect,'linewidth',3)
    
    hdlfig=subplot(1,2,1);hold on;
    %plot_CSD(csd_patch,flip(zs),1,1,[])
    plot_csdtrials('csd',csd_patch,zs,[],[],[],[],info,hdlfig,[]);
    colormap(colormap_redblackblue)
    title({info.datafile ; info.align})
    
    subplot(1,2,2);hold on;
    plot(csd_prof,[1:length(csd_prof)])
    line([0 0] ,[0 length(zs)]);
    axis tight;ax=axis;
    
end


csdfeat=[];
if find_feat,
    
    %heuristically
    switch info.align
        case 'targ'
            [peaksourceval peaksource]=min(csd_prof);
            upcross=min(find(csd_prof(peaksource:end)>0)-1);
            if ~isempty(upcross)
                csdfeat(1)=upcross+peaksource;
            else
                csdfeat(1)=length(csd_prof);
            end
            downcross=max(find(csd_prof(1:peaksource)>0));
            if ~isempty(downcross)
                csdfeat(2)=downcross+1;
            else
                csdfeat(2)=1;
            end
            hl=line([ax(1) ax(2)],[csdfeat(1) csdfeat(1)]);
            set(hl,'linewidth',2);
            hl=line([ax(1) ax(2)],[csdfeat(2) csdfeat(2)]);
            set(hl,'linewidth',2);
        case 'sacc'
            f=0;
            delta=10;
            csdfeataux=[];
            for z=1+delta:length(csd_prof)-delta
                if sum(csd_prof(z-delta:z)>0)>0.99*delta & sum(csd_prof(z:z+delta)<0)>0.99*delta
                    f=f+1;
                    csdfeataux(f)=z;
                    line([ax(1) ax(2)],[csdfeataux(f) csdfeataux(f)]);
                end
            end
            if ~isempty(csdfeataux)
                [peaksourceval peaksource]=min(csd_prof);
                distsource=abs(csdfeataux-peaksource);
                [valmin indmin]=min(distsource);
                csdfeat(1)=csdfeataux(indmin);
                hl=line([ax(1) ax(2)],[csdfeat(1) csdfeat(1)]);
                set(hl,'linewidth',2);
            else
                csdfeat=[];
            end
    end
    
    %by hand
    csdfeat
    if ~auto
        load zs_iCSD
        ncsdsamples=length(zs_iCSD);

        ok=0;
        ok = str2num(input('Are you ok with the auto features? 1/0:','s'));
        
        if ~ok
            switch info.align
                case 'targ'
                    while ~ok
                        [x,y] = ginput(2);
                        if x(:)<0
                            csdfeat(1:2)=[];
                        else
                            csdfeat(1)=round(y(1));
                            csdfeat(2)=round(y(2));
                            
                            hl=line([ax(1) ax(2)],[csdfeat(1) csdfeat(1)]);
                            set(hl,'linewidth',2);
                            hl=line([ax(1) ax(2)],[csdfeat(2) csdfeat(2)]);
                            set(hl,'linewidth',2);
                        end
                        csdfeat
                        
                        %%%%%%%%%%%%%%%%%%%%%%
                        %display chosen depth (bottom sink) as channel number
                        %alignment of onset using CSD features (after compute_CSDfeature)
                        %info.csdfeat_avg_targ=data(1).offline.csdfeat_avg_targ;
                        %info.zs=data(1).offline.csdzs;
                        dref=csdfeat(2);
                        
                        %convert dref into channel number (round to closest channel)
                        %WARNING this is an approximation, best to convert data into depths(mm)
                        dref_conv=round(dref*(info.nchannels+2)/ncsdsamples)
                        pause
                        
                        %%%%%%%%%%%%%%%%%%%%%%
            
                        ok = str2num(input('Are you ok with the features? 1/0:','s'));
                    end
                case 'sacc'
                    while ~ok
                        [x,y] = ginput(1);
                        if x(:)<0
                            csdfeat(1)=[];
                        else
                            csdfeat(1)=round(y);
                            
                            hl=line([ax(1) ax(2)],[csdfeat(1) csdfeat(1)]);
                            set(hl,'linewidth',2);
                        end
                        csdfeat
                        ok = str2num(input('Are you ok with the features? 1/0:','s'));
                    end
            end
        end
    end
end

if display
    close(figcsdfeat)
end
