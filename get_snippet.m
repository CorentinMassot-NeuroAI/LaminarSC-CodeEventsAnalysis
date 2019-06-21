function snippet=get_snippet(signal,trial,event_align,timevec)

%function snippet=get_snippet(signal,trial,event_align,timevec)
%   get snippet of spiketimes in timevec based on event
%
% see also compute_fr
%
% Corentin Massot based on Uday Jagadisan's code
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 01/09/2017 last modified 01/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nchs
nchs=numel(trial.spikeTimestamps);

%timevec
if isempty(timevec)
    timevec=[1:size(trial.lfp,2)];
end

snippet=nan(nchs,length(timevec));

if event_align~=0
    switch signal
        case 'spk'
            for ch=1:nchs
                spiketimes=round(1000*trial.spikeTimestamps{ch});
                spiketimes_align = spiketimes - event_align;
                ind = ismember(spiketimes_align,timevec);
                
                %NOTE:hist (spike may be = 2 or 3 because of round to nearest integer)
                snippet(ch,:)= hist(spiketimes_align(ind),timevec);
            end
            
        case 'fr'
            %get fr (see also compute_fr)
            %trial_fr=trial.offline.fr_epsp_6;
            trial_fr=trial.offline.fr_epsp_20;%20ms is the correct value (Thompson et al. Schall, 1996)
            
            for ch=1:nchs
                pre = event_align+min(timevec);
                post = event_align+max(timevec);
                if pre<0;
                    %NOTE: if pre<0, adjust wind lim instead of risking alignment error
                    display('Error: inferior lim of wind is outside lfp signal!')
                    pause
                elseif pre>0 & post<=size(trial_fr,2)
                    snippet(ch,:) = trial_fr(ch,event_align+timevec);
                elseif post>size(trial_fr,2)
                    %display('Warning: lfp signal shorter than wind limit!')
                    lfpaux=[trial_fr(ch,:),NaN(1,post-size(trial_fr,2))];
                    snippet(ch,:) = lfpaux(event_align+timevec);
                end
            end
            
        case 'lfp'
            %LFP from Ripple filtered lfp
            %trial_lfp=trial.lfp;
            %LFP offline 60Hz filtered (see also compute_lfp) 
            trial_lfp=trial.offline.lfp_60120Hz;
            %LFP offline filtered (see also compute_filtered_lfp) 
            %trial_lfp=trial.offline.lfp_low50Hz;%100Hz
            
            for ch=1:nchs
                pre = event_align+min(timevec);
                post = event_align+max(timevec);
                if pre<0;
                    %NOTE: if pre<0, adjust wind lim instead of risking alignment error
                    display('Error: inferior lim of wind is outside lfp signal!')
                    pause
                elseif pre>0 & post<=size(trial_lfp,2)
                    snippet(ch,:) = trial_lfp(ch,event_align+timevec);
                elseif post>size(trial_lfp,2)
                    %display('Warning: lfp signal shorter than wind limit!')
                    lfpaux=[trial_lfp(ch,:),NaN(1,post-size(trial_lfp,2))];
                    snippet(ch,:) = lfpaux(event_align+timevec);
                end
            end
            
        case 'raw'
            samp=30;%30KHz
            for ch=1:nchs
                pre = samp*(event_align+min(timevec));
                post = samp*(event_align+max(timevec));
                timings=samp*event_align+ [samp*min(timevec):1:samp*max(timevec)];
                
                if pre<0;
                    %NOTE: if pre<0, adjust wind lim instead of risking alignment error
                    display('Error: inferior lim of wind is outside raw signal!')
                    pause
                elseif pre>0 & post<=size(trial.raw,2)
                    snippet(ch,:) = trial.raw(ch,timings);
                elseif post>size(trial.raw,2)
                    %display('Warning: raw signal shorter than wind limit!')
                    rawpaux=[trial.raw(ch,:),NaN(1,post-size(trial.raw,2))];
                    snippet(ch,:) = rawpaux(timings);
                end
            end
    end
end
