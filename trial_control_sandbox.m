%% setup

% timing
total_time = 1 ; % total time to run converted from minutes to seconds.
total_time = total_time * 60; % converted from minutes to seconds.
dt = 0.05; % sampling rate. 0.05 is much faster than you probably need.

% states
states = repmat([1;2],100,1); % lazy way to repeating numbers.

prior_states = [];
trials = [];
iState = 0; % initial state;

% make a figure to track everything
figure(101)
xlim([.5 .75])
ylim([.45 .55])

for iT = 0:dt:total_time
    eTime = sprintf('%2.1f', iT); % time in a nice string format.
    
    if iT > 30
        if iState == 0 % move to a state;
            iState = datasample(states, 1);
            while length(prior_states) >= 3 &&  sum(prior_states(end-2:end) == iState) ==3 % check if the same state has appeared three times in a row.
                iState = datasample(states, 1);
            end
        end
        
        %update plot
        clf
        text(.5, .5, ['Time: ' eTime 's'], 'FontSize', 24)
        text(.5, .25, ['State: ' num2str(iState)], 'FontSize', 24)
        drawnow
        
        if iState == 1
            licks = 0; % initialize licks.
            
            %            play a tone;
            
            
            %check for licks (during tone)% when writing this make sure it
            %logs the number and the time.
            
            
            % do a thing like fire valve;
            
            
            % log it and return to 0 state
            prior_states(end+1) = iState;
            % save the output as 'trials' with a bunch of informaiton.
            trials{end+1} = [];
            trials{end}.state = iState;
            trials{end}.licks = licks;
            % add in other details.
            
            
            iState = 0; % return to 0 state.
            
        elseif iState == 2
            licks = 0; % initialize licks.
            
            %play a tone;
            
            
            %check for licks%
            
            
            % do a thing like fire valve;
            
            
            % log it and return to 0 state
            prior_states(end+1) = iState;
            % save the output as 'trials' with a bunch of informaiton.
            trials{end+1} = [];
            trials{end}.state = iState;
            trials{end}.licks = licks;
            % add in other details.
            
            
            iState = 0; % return to 0 state.
        end
        
    else
        
        clf
        text(.5, .5, ['Time: ' eTime 's'], 'FontSize', 24)
        text(.5, .25, 'State: waiting...', 'FontSize', 24)
        drawnow
    end
    
end % main loop


%% save things.

save(['Trials_' datestr(date, 'yyyy_mm_dd') '_' datestr(now, 'HH_MM')], 'trials')

