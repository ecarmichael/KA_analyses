%% setup

% set mouse ID
%mouseID = '_'; % which mouse is being run
%date = datestr(now,'yyyy-mm-dd');
%disp(['Running ' mouseID ' on ' date])
%%
% initialize arduino (edited by KAA), need to include capacitive sensor
% library path inside the ardunio() function, included library for capsense
a_board = arduino('COM4', 'Due');
%%
% timing
total_time = 1.5; % total time to run converted from minutes to seconds.
total_time = total_time * 60; % converted from minutes to seconds.

% states
states = repmat([1;2],100,1); % lazy way to repeating numbers.

prior_states = [];
trials = [];
iState = 0; % initial state;

%tones
Fs=44100/10;      %CD quality - also conveniently divisible by 30 and 25
stim_dur=5;               %duration of puretone in seconds
T=0:1/Fs:stim_dur-1/Fs;

%tone 1
t1_freq =1200;
Y1 = sin(2*pi*t1_freq*T);
%tone2
t2_freq = 800;
Y2 = sin(2*pi*t2_freq*T);

% make a figure to track everything
figure(101)
xlim([.5 .75])
ylim([.45 .55])

t0 = clock;
t_state = t0;  % get the time at the end of the state, used to know when to start the next state.
%%
while etime(clock, t0) < total_time
    this_time = sprintf('%2.1f', etime(clock, t0)); % time in a nice string format.
    
    % I want the capacitance sensor to run here (continuous
    % sampling of the capacitance, and if cap_value > x, log a lick event
    % sample code below (line 49 and 50; 51 and 52)
    while readDigitalPin(a_board, 'D10') ~=0
    sprintf('Cap1: %d', cap_1_value)
    end
    while readDigitalPin(a_board, 'D8') ~=0
    sprintf('%d', cap_2_value)
    end
    
    
    if etime(clock, t_state) > 10 % change this to the initial time.
        if iState == 0 % move to a state;
            iState = datasample(states, 1);
            while length(prior_states) >= 3 &&  sum(prior_states(end-2:end) == iState) ==3 % check if the same state has appeared three times in a row.
                iState = datasample(states, 1);
            end
        end
        
        if iState == 1
            tone_t = clock; % get elapsed time.
            licks_1 = []; %initialize licks
            licks_2 = [];
            % play a tone; log the time. Fire laser
            sound(Y1, Fs)
            
            freq = 4; % in Hz
            t_laser_start = clock;
            t_laser_end = 5;
            
            p_start = clock;
            p1 = 1; p2 = 1;
            while etime(clock, t_laser_start) < t_laser_end
                % check the caps
                if cap_1_value > 0.5
                    licks_1(end+1) = etime(clock, t0);
                    disp(['lick1:' num2str(licks_1(end))])
                end
                
                if cap_2_value > 0.5
                    licks_2(end+1) = etime(clock, t0);
                    disp(['lick2:' num2str(licks_2(end))])
                end
                
                if etime(clock, p_start) < (1/freq)/2 %pulse on
                    if p1 == 1
                        writeDigitalPin(a_board, 'D6', 1) % turn the TTL on
                        fprintf('on \n') % print on first instance
                        p1 = 0;
                    end
                elseif etime(clock, p_start) >= (1/freq)/2 && etime(clock, p_start) < (1/freq) % pulse off
                    if p2 == 1
                        fprintf('off \n')
                        writeDigitalPin(a_board, 'D6', 0); % turn TTl off
                        p2 = 0;
                    end
                elseif etime(clock, p_start) >= (1/freq) % reset the pulse cycle.
                    p_start = clock;
                    p1 = 1; p2 = 1;
                    disp('cycle')
                end
            end
            writeDigitalPin(a_board, 'D6', 0); % turn TTl off
            
            
            %while etime(clock, tone_t) < 5
            
            %   clf
            %   text(.5, .5, ['Time: ' this_time 's'], 'FontSize', 24);
            %   text(.5, .4, ['Trial: ' sprintf('%2.0f',length(trials)+1)  ], 'FontSize', 24);
            %   text(.5, .3, ['State: Tone 1 ' sprintf('%2.1f', etime(clock, tone_t))], 'FontSize', 24, 'color', 'b');
            %   drawnow
            %end
            
            % do a thing like fire valve; log the time.
            % edited by KAA to send TTL via arduino, D2 is the 2nd digital
            % pin, will be connected to Outcome 1
            valve_t = clock;
            if etime(clock, valve_t) < 0.1
                writeDigitalPin(a_board, 'D2', 0)
                %             pause(0.1)
            elseif etime(clock, valve_t) > 0.1 && etime(clock, valve_t) <= 0.5
                writeDigitalPin(a_board, 'D2', 1)
                %             pause(0.5)
            elseif etime(clock, valve_t) < 0.5
                writeDigitalPin(a_board, 'D2', 0)
                %             pause(5.0)
            end
            post_licks_1 = []; %trying to sample post-tone licks during 5 second reward bin (see pause in line above)
            post_licks_2 = [];
            if cap_1_value > 0.5
                post_licks_1(end+1) = etime(clock, t0);
            end
            if cap_2_value > 0.5
                post_licks_2(end+1) = etime(clock, t0);
            end
            
            % log it and return to 0 state
            prior_states(end+1) = iState;
            % save the output as 'trials' with a bunch of informaiton.
            trials{end+1} = [];
            trials{end}.state = iState;
            trials{end}.licks_1 = licks_1;
            trials{end}.licks_2 = licks_2;
            trials{end}.post_licks_1 = post_licks_1;
            trials{end}.post_licks_2 = post_licks_2;
            trials{end}.tone_time = tone_t;
            trials{end}.valve_t = valve_t;
            % add in other details.
            
            iState = 0; % return to 0 state.
            t_state = clock;  % get the time at the end of the state, used to know when to start the next state.
        elseif iState == 2
            
            % play a tone; log the time.
            tone_t = clock;
            licks_1 = []; %initialize licks
            licks_2 = [];
            tic; sound(Y2, Fs); toc; % seems to take between 0.01 and .33 seconds to play the tone.
            
            freq = 4; % in Hz
            t_laser_start = clock;
            t_laser_end = 5;
            
            p_start = clock;
            p1 = 1; p2 = 1;
            while etime(clock, t_laser_start) < t_laser_end
                % check the caps
                if cap_1_value > 0.5
                    licks_1(end+1) = etime(clock, t0);
                    disp(['lick1:' num2str(licks_1(end))])
                end
                
                if cap_2_value > 0.5
                    licks_2(end+1) = etime(clock, t0);
                    disp(['lick2:' num2str(licks_2(end))])
                end
                
                if etime(clock, p_start) < (1/freq)/2 %pulse on
                    if p1 == 1
                        writeDigitalPin(a_board, 'D6', 1) % turn the TTL on
                        fprintf('on \n') % print on first instance
                        p1 = 0;
                    end
                elseif etime(clock, p_start) >= (1/freq)/2 && etime(clock, p_start) < (1/freq) % pulse off
                    if p2 == 1
                        fprintf('off \n')
                        writeDigitalPin(a_board, 'D6', 0); % turn TTl off
                        p2 = 0;
                    end
                elseif etime(clock, p_start) >= (1/freq) % reset the pulse cycle.
                    p_start = clock;
                    p1 = 1; p2 = 1;
                    disp('cycle')
                end
            end
            writeDigitalPin(a_board, 'D6', 0); % turn TTl off
            
            %while etime(clock, tone_t) < 5
            
            %    clf
            %    text(.5, .5, ['Time: ' this_time 's'], 'FontSize', 24)
            %    text(.5, .4, ['Trial: ' sprintf('%2.0f',length(trials)+1)  ], 'FontSize', 24)
            %    text(.5, .3, ['State: Tone 2 ' sprintf('%2.1f', etime(clock, tone_t))], 'FontSize', 24, 'color', 'r')
            %    drawnow
            %end
            
            
            % do a thing like fire valve; log the time.
            % edited by KA to send TTL via arduino, pin D4 sends TTL to
            % Outcome 2
            %             valve_t = clock;
            %             writeDigitalPin(a_board, 'D4', 0)
            %             pause(0.1)
            %             writeDigitalPin(a_board, 'D4', 1)
            %             pause(0.5)
            %             writeDigitalPin(a_board, 'D4', 0)
            %             pause(5.0)
            valve_t = clock;
            if etime(clock, valve_t) < 0.1
                writeDigitalPin(a_board, 'D4', 0)
                %             pause(0.1)
            elseif etime(clock, valve_t) > 0.1 && etime(clock, valve_t) <= 0.5
                writeDigitalPin(a_board, 'D4', 1)
                %             pause(0.5)
            elseif etime(clock, valve_t) < 0.5
                writeDigitalPin(a_board, 'D4', 0)
                %             pause(5.0)
            end
            post_licks_1 = []; % same issue about post-tone reward bin lick logging
            post_licks_2 = [];
            if cap_1_value > 0.5
                post_licks_1(end+1) = etime(clock, t0);
            end
            if cap_2_value > 0.5
                post_licks_2(end+1) = etime(clock, t0);
            end
            
            % log it and return to 0 state
            prior_states(end+1) = iState;
            % save the output as 'trials' with a bunch of informaiton.
            trials{end+1} = [];
            trials{end}.state = iState;
            trials{end}.licks_1 = licks_1;
            trials{end}.licks_2 = licks_2;
            trials{end}.post_licks_1 = post_licks_1;
            trials{end}.post_licks_2 = post_licks_2;
            trials{end}.tone_time = tone_t;
            trials{end}.valve_t = valve_t;
            % add in other details.
            
            iState = 0; % return to 0 state.
            t_state = clock;  % get the time at the end of the state, used to know when to start the next state.
        end
        
    else
        clf
        text(.5, .5, ['Time: ' this_time 's'], 'FontSize', 24)
        text(.5, .4, ['Trial: ' sprintf('%2.0f',length(trials))  ], 'FontSize', 24)
        text(.5, .3, 'State: waiting...', 'FontSize', 24)
        drawnow
    end
    
end % main loop
disp('all done!')

%% save things.

% save useful variables with dates and times.  Probably good to add a
% subect ID or save it in a specific directory.
save(['Trials_' date '_' datestr(now, 'HH_MM')], 'trials')
disp([mouseID ' data saved!'])

