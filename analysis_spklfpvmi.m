 %function analysis_spklfpvmi

%function analysis_spklfpvmi
%   PLot spiking activity,lfp and VMi from laminar data
%
% see also compute_vmi
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 10/18/2016 last modified 01/22/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set paths
[root_path data_path save_path]=set_paths;

%screen size
scrsz = get(groot,'ScreenSize');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
%print figures and save data
savedata=0;
savefigs=0;
figtype='png';%'epsc2';

%alignement
%alignlist={'no' 'targ' 'go' 'sacc'};
alignlist={'targ' 'sacc'};
%alignlist={'targ' 'peak'};
    
%windows of analysis
wind_targ=[0 250];%[-50 350];
wind_sacc=[-100 150];%[-200 200];
wind_peak=[-100 150];%[-200 200];

wind_bsl=[-50 50];

% %adaptive window
% [p,polystats] = polyfit([5 20],[55 100],1);


% %vshift
% vshift_spk=100;
% vshift_lfp=20;%29;

%sigma FR
sigma_FR=6;


%start .pptx file
savepptx=0;
if savepptx
    isOpen  = exportToPPTX();
    if ~isempty(isOpen),
        % If PowerPoint already started, then close first and then open a new one
        exportToPPTX('close');
    end
    
    exportToPPTX('new','Dimensions',[12 6], ...
        'Title','SC dataset', ...
        'Author','Corentin', ...
        'Subject','Automatically generated PPTX file form output of analysis_spklfpvmi.m', ...
        'Comments',' ');

    %tmp filename
    file=[save_path 'spklfp_tmp' '.' figtype];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get data
datalist=load_data_gandhilab(data_path);

%colorlist
colorlist=get_colorlist;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create list of vmis
dlist=get_dlist

hdlfigallvmis=figure;hold on;
data=[];info=[];
dd=0;
for d=64%dlist
    
    %counter
    dd=dd+1;
    
    %get data and info
    info.datafile=datalist{d};
    load ([data_path info.datafile]);
    display(info.datafile)
    
    %getting channel mapping and discard selected bad channels
    [info.chmap info.nchannels info.depths]=get_chmap(data(1).info.electrode{2},[]);
    %getting trial type
    info.trialtype=data(1).sequence(1);
    %getting list of targets
    targslist=data(1).offline.targslist;
    %targets index
    targs_ind=get_targsindex(targslist,info);
    
    %target tuning (after compute_tuning)
    info.targ_tuning=data(1).offline.targ_tuning;
    %info.targ_tuning=targ_tuning;

    %select trials
    seltrials=get_seltrials(data,'rpt');
    
    %VMI index (after compute_vmi)
    vmis=data(1).offline.vmi.vmis_peak_bsl_bsignif;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Data aligned on target and saccade onset
    for al=1:numel(alignlist)
        info.align=alignlist{al};
        switch info.align
            case 'targ'
                [alltrials_spk_targ aligntime_targ]=get_alltrials_align(data,seltrials,wind_targ,'fr',info,targslist,sigma_FR,1);
                [alltrials_lfp_targ aligntime_targ]=get_alltrials_align(data,seltrials,wind_targ,'lfp',info,targslist,sigma_FR,1);

                %baseline
                [alltrials_spk_bsl aligntime_bsl]=get_alltrials_align(data,seltrials,wind_bsl,'fr',info,targslist,sigma_FR,1);
                [alltrials_lfp_bsl aligntime_bsl]=get_alltrials_align(data,seltrials,wind_bsl,'lfp',info,targslist,sigma_FR,1);
                
            case 'sacc'
                [alltrials_spk_sacc aligntime_sacc]=get_alltrials_align(data,seltrials,wind_sacc,'fr',info,targslist,sigma_FR,1);
                [alltrials_lfp_sacc aligntime_sacc]=get_alltrials_align(data,seltrials,wind_sacc,'lfp',info,targslist,sigma_FR,1);
            
            case 'peak'
                [alltrials_spk_peak aligntime_peak]=get_alltrials_align(data,seltrials,wind_peak,'fr',info,targslist,sigma_FR,1);
                [alltrials_lfp_peak aligntime_peak]=get_alltrials_align(data,seltrials,wind_peak,'lfp',info,targslist,sigma_FR,1);
        end
    end
 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %analysis of trials for each target
    allspk_vmi_avg=zeros(2,length(targslist),info.nchannels);
    for tg=info.targ_tuning;%[2:6]%%targs_ind,
        
        %target index
        info.targ=tg;
        
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         %Compute size of adaptive window for analysis
        %         pos=targslist(tg,:);
        %         amp=sqrt(pos(1)^2+pos(2)^2);
        %         if amp>5,
        %             wadapt=floor(p(2)+amp*p(1));
        %         else
        %             wadapt=ceil(p(2)+5*p(1));
        %         end
        %
        %         wind_targ_vmi=[110 110+wadapt];
        %         wind_targ_bsl=[50-wadapt 50];
        %         wind_sacc_vmi=[-25 -25+wadapt];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figtrials=figure('Position',[1 100 scrsz(3)-100 scrsz(4)-200]);
        
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         %%
        %         %plot vmi
        %         hdlfig=subplot(1,3,1);hold on;
        %         titlestr={info.datafile ; [' t(tuning)' num2str(info.targ)]};
        %         plot_vmis(vmis,info.targ_tuning,'-',1,2,info,hdlfig,titlestr);

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %plot spk of targ and sacc
        vshift=80;trials_al=[];
        for al=1:2
            info.align=alignlist{al};
            switch info.align
                case 'targ'
                    info.aligntime=aligntime_targ;
                    %wind_vmi=wind_targ_vmi; 
                    trials_spk=alltrials_spk_targ{tg};
                case 'sacc'
                    info.aligntime=aligntime_sacc;
                    %wind_vmi=wind_sacc_vmi;
                    trials_spk=alltrials_spk_sacc{tg};
                case 'peak'
                    info.aligntime=aligntime_peak;
                    %wind_vmi=wind_peak_vmi;
                    trials_spk=alltrials_spk_peak{tg};
            end
            %baseline
            %wind_bsl=wind_targ_bsl;
            trials_spk_bsl=alltrials_spk_bsl{tg};
            
            [info.nchannels info.ntrials info.triallen]=size(trials_spk);
            %compute average trials
            [trials_spk_avg trials_spk_var]=get_trials_avg(trials_spk);
            %remove trials with amplitude that is too small
            [trials_spk_avgc index_spk_c]=clean_trials(trials_spk_avg,'fr');
            %compute average trials of baseline
            [trials_spk_bsl_avg trials_spk_bsl_var]=get_trials_avg(trials_spk_bsl);
            %normalize average trials by baseline
            trials_spk_avgcn=get_trials_normalized(trials_spk_avgc,trials_spk_bsl_avg,'fr',info);
            trials_al(al,:,:)=trials_spk_avgcn;
            
            %vshift
            vshift=min([vshift max(max(abs(squeeze(trials_al(al,:,:))),[],2))/4]);
        end

        for al=1:2
            info.align=alignlist{al};
            
            %hdlfig=subplot(1,3,al+1);hold on;
            hdlfig=subplot(2,2,al);hold on;
            if strcmp(info.align,'targ')
                info.aligntime=aligntime_targ;
                titlestr={info.datafile ; ['FR ' info.align ' (unit=' num2str(round(vshift)) 'spk/s) t' num2str(info.targ) ' #trials:' num2str(info.ntrials) ]};
            elseif strcmp(info.align,'sacc')
                info.aligntime=aligntime_sacc;
                titlestr={'' ; ['FR ' info.align ]};
            elseif strcmp(info.align,'peak')
                info.aligntime=aligntime_peak;
                titlestr={'' ; ['FR ' info.align ]};
            end
            plot_trials(squeeze(trials_al(al,:,:)),[],index_spk_c,vshift,[],[],info,hdlfig,titlestr);
           
            %             %plot wind_vmi limits
            %             plot_event(wind_vmi,info.aligntime,range,1,hdlfig);
            %             %plot wind_bsl limits
            %             if strcmp(info.align,'targ')
            %                 plot_event(wind_bsl,info.aligntime,range,3,hdlfig);
            %             end
            
            %same axis
            if al==1,
                axis tight;
                ax1=axis;
            else
                axis tight;
                ax2=axis;
                axis([ax2(1) ax2(2) ax1(3) ax1(4)]);
            end
            grid;
            
            %pause
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %plot lfp of targ and sacc
        vshift=30;trials_al=[];
        for al=1:2
            info.align=alignlist{al};
            switch info.align
                case 'targ'
                    trials_lfp=alltrials_lfp_targ{tg};
                    info.aligntime=aligntime_targ;
                    %wind_vmi=wind_targ_vmi;
                case 'sacc'
                    trials_lfp=alltrials_lfp_sacc{tg};
                    info.aligntime=aligntime_sacc;
                    %wind_vmi=wind_sacc_vmi;
                case 'peak'
                    trials_lfp=alltrials_lfp_peak{tg};
                    info.aligntime=aligntime_peak;
                    %wind_vmi=wind_sacc_vmi;
            end
            %baseline
            trials_lfp_bsl=alltrials_lfp_bsl{tg};
            
            [info.nchannels info.ntrials info.triallen]=size(trials_lfp);
            %compute average trials
            [trials_lfp_avg trials_lfp_var]=get_trials_avg(trials_lfp);
            %remove trials with amplitude that is too small
            [trials_lfp_avgc index_lfp_c]=clean_trials(trials_lfp_avg,'lfp');
            %compute average trials of baseline
            [trials_lfp_bsl_avg trials_lfp_bsl_var]=get_trials_avg(trials_lfp_bsl);
            %normalize average trials by baseline
            trials_lfp_avgcn=get_trials_normalized(trials_lfp_avgc,trials_lfp_bsl_avg,'fr',info);
            trials_al(al,:,:)=trials_lfp_avgcn;
            
            %vshift
            vshift=min([vshift max(max(abs(squeeze(trials_al(al,:,:))),[],2))/4]);
            
        end

        for al=1:2
            info.align=alignlist{al};
            
            %hdlfig=subplot(2,3,al+4);hold on;
            hdlfig=subplot(2,2,al+2);hold on;
            if strcmp(info.align,'targ')
                info.aligntime=aligntime_targ;
                titlestr=['LFP ' info.align ' (unit=' num2str(round(vshift)) 'mV)'];
            elseif strcmp(info.align,'sacc')
                info.aligntime=aligntime_sacc;
                titlestr=['LFP ' info.align];
            elseif strcmp(info.align,'peak')
                info.aligntime=aligntime_peak;
                titlestr=['LFP ' info.align];
            end
            plot_trials(squeeze(trials_al(al,:,:)),[],index_lfp_c,vshift,[],[],info,hdlfig,titlestr);
            %plot wind
            %plot_event(wind_vmi,info.aligntime,range,1,hdlfig);
            
            %same axis
            if al==1,
                axis tight;
                ax1=axis;
            else
                axis tight;
                ax2=axis;
                axis([ax2(1) ax2(2) ax1(3) ax1(4)]);
            end
            grid;
        end
        
        if savepptx,
            %save tmp figure
            saveas(hdlfig,file,figtype);
            %create slide
            slideNum = exportToPPTX('addslide',hdlfig);
            fprintf('Added slide %d\n',slideNum);
            
            %add figure
            exportToPPTX('addpicture',file,'scale','maxfixed');
        end
        
        pause
    
    end

    %pause
    close(figtrials)
       
end


if savepptx
    %close .pptx
    newFile = exportToPPTX('saveandclose','SCdataset')
end

