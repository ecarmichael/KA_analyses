%% setup

% timing
total_time = .5 ; % total time to run converted from minutes to seconds.
total_time = total_time * 60; % converted from minutes to seconds.

% states
states = repmat([1;2],100,1); % lazy way to repeating numbers.

prior_states = [];
trials = [];
iState = 0; % initial state;

%tones
Fs=44100;                       %CD quality - also conveniently divisible by 30 and 25
stim_dur=5;               %duration of puretone in seconds
T=0:1/Fs:stim_dur-1/Fs;

%tone 1
t1_freq =8000; 
Y1 = sin(2*pi*t1_freq*T);
%tone2
t2_freq = 1000; 
Y2 = sin(2*pi*t2_freq*T);

% make a figure to track everything
figure(101)
xlim([.5 .75])
ylim([.45 .55])

t0 = clock; 

while etime(clock, t0) < total_time
    this_time = sprintf('%2.1f', etime(clock, t0)); % time in a nice string format.
    
    if etime(clock, t0) > 10
        if iState == 0 % move to a state;
            iState = datasample(states, 1);
            while length(prior_states) >= 3 &&  sum(prior_states(end-2:end) == iState) ==3 % check if the same state has appeared three times in a row.
                iState = datasample(states, 1);
            end
        end
        
        %update plot
        clf
        text(.5, .5, ['Time: ' this_time 's'], 'FontSize', 24)
        text(.5, .25, ['State: ' num2str(iState)], 'FontSize', 24)
        drawnow
        
        if iState == 1
            tone_t = clock; % get elapsed time. 
            licks = []; % initialize licks.
            
            % play a tone; log the time. 
            sound(Y1, Fs)
            
           while etime(clock, tone_t) < 5
                
                % lick check function here. check for licks (during tone)% when writing this make sure it logs the time. time points can be counted later and can also be added to the plot.
                clf
                text(.5, .5, ['Time: ' this_time 's'], 'FontSize', 24)
                text(.5, .25, ['State: Tone 1'], 'FontSize', 24, 'color', 'b')
                drawnow
            end
     
            
            
            %check for licks
            
            
            % do a thing like fire valve; log the time. 
            valve_t = clock; 
            
            %check for post-valve licks.  Could also be done later but this
            %way you can plot them in a different color or track
            %correct/incorrect responses in real time. 
            post_licks = []; 
            
            % log it and return to 0 state
            prior_states(end+1) = iState;
            % save the output as 'trials' with a bunch of informaiton.
            trials{end+1} = [];
            trials{end}.state = iState;
            trials{end}.licks = licks;
            trials{end}.tone_time = tone_t; 
            trials{end}.valve_t = valve_t; 
            % add in other details.
            
            iState = 0; % return to 0 state.
            
        elseif iState == 2
            licks = []; % initialize licks.
            
            % play a tone; log the time. 
            tone_t = clock; 
            tic; sound(Y2, Fs); toc; % seems to take .33 seconds to play the tone. 

            while etime(clock, tone_t) < 5
                % lick check function here. check for licks (during tone)% when writing this make sure it logs the time. time points can be counted later and can also be added to the plot.
                clf
                text(.5, .5, ['Time: ' this_time 's'], 'FontSize', 24)
                text(.5, .25, ['State: Tone 2'], 'FontSize', 24, 'color', 'r')
                drawnow
            end
            
            
            % do a thing like fire valve; log the time. 
            valve_t = clock; 
            
            %check for post-valve licks
            post_licks = []; 
            
            % log it and return to 0 state
            prior_states(end+1) = iState;
            % save the output as 'trials' with a bunch of informaiton.
            trials{end+1} = [];
            trials{end}.state = iState;
            trials{end}.licks = licks;
            trials{end}.tone_time = tone_t; 
            trials{end}.valve_t = valve_t; 
            % add in other details.
            
            
            iState = 0; % return to 0 state.
        end
        
    else
        
        clf
        text(.5, .5, ['Time: ' this_time 's'], 'FontSize', 24)
        text(.5, .25, 'State: waiting...', 'FontSize', 24)
        drawnow
    end
    
end % main loop


%% save things.

% save useful variables with dates and times.  Probably good to add a
% subect ID or save it in a specific directory. 
save(['Trials_' datestr(date, 'yyyy_mm_dd') '_' datestr(now, 'HH_MM')], 'trials')

