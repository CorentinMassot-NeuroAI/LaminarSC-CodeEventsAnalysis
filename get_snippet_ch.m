function snippet=get_snippet_ch(signal,trial,event_align,timevec)

%function snippet=get_snippet_ch(signal,trial,event_align,timevec)
%   get snippet of spiketimes in timevec based on event that is different
%   for each trial
%
% see also compute_fr
%
% Corentin Massot based on Uday Jagadisan's code
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh  
% created 10/17/2017 last modified 10/17/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nchs
nchs=numel(trial.spikeTimestamps);

%timevec
if isempty(timevec)
    timevec=[1:size(trial.lfp,2)];
end

snippet=nan(nchs,length(timevec));

switch signal
    case 'spk'
        for ch=1:nchs
            if event_align(ch)~=0
                spiketimes=round(1000*trial.spikeTimestamps{ch});
                spiketimes_align = spiketimes - event_align(ch);
                ind = ismember(spiketimes_align,timevec);
                
                %NOTE:hist (spike may be = 2 or 3 because of round to nearest integer)
                snippet(ch,:)= hist(spiketimes_align(ind),timevec);
            end
        end
        
    case 'fr'
        %get fr (see also compute_fr)
        trial_fr=trial.offline.fr_epsp_6;
        
        for ch=1:nchs
            if event_align(ch)~=0
                pre = event_align(ch)+min(timevec);
                post = event_align(ch)+max(timevec);
                if pre<0;
                    %NOTE: if pre<0, adjust wind lim instead of risking alignment error
                    display('Error: inferior lim of wind is outside lfp signal!')
                    pause
                elseif pre>0 & post<=size(trial_fr,2)
                    snippet(ch,:) = trial_fr(ch,event_align(ch)+timevec);
                elseif post>size(trial_fr,2)
                    %display('Warning: lfp signal shorter than wind limit!')
                    lfpaux=[trial_fr(ch,:),NaN(1,post-size(trial_fr,2))];
                    snippet(ch,:) = lfpaux(event_align(ch)+timevec);
                end
            end
        end
        
    case 'lfp'
        for ch=1:nchs
            if event_align(ch)~=0
                pre = event_align(ch)+min(timevec);
                post = event_align(ch)+max(timevec);
                if pre<0;
                    %NOTE: if pre<0, adjust wind lim instead of risking alignment error
                    display('Error: inferior lim of wind is outside lfp signal!')
                    pause
                elseif pre>0 & post<=size(trial.lfp,2)
                    snippet(ch,:) = trial.lfp(ch,event_align(ch)+timevec);
                elseif post>size(trial.lfp,2)
                    %display('Warning: lfp signal shorter than wind limit!')
                    lfpaux=[trial.lfp(ch,:),NaN(1,post-size(trial.lfp,2))];
                    snippet(ch,:) = lfpaux(event_align(ch)+timevec);
                end
            end
        end
        
    case 'raw'
        samp=30;%30KHz
        for ch=1:nchs
            if event_align(ch)~=0
                pre = samp*(event_align(ch)+min(timevec));
                post = samp*(event_align(ch)+max(timevec));
                timings=samp*event_align(ch)+ [samp*min(timevec):1:samp*max(timevec)];
                
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
            else
                snippet(ch,:)=nan;
            end
        end
end
