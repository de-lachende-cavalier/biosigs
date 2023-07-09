function [r_peaks_amps, r_peaks_idxs, delay] = pan_tompkins(acc_i, f, graph)
%PAN_TOMPKINS This function applies the Pan-Tompkins algorithm, given the
%ith component of the acceleration and the frequency f.
%   
%   The third argument (graph) is a boolean indicating whether we want to
%   graph the procedure or not. 
%   The third output argument returned (delay) returns
%   the number of samples the signal is delayed by due to the filtering.
%
%   The initial algorithm was the one provided Hooman Sedghamiz (hooman.sedghamiz@gmail.com), 
%   which I cleaned up and edited following the tips in the Pan-Tompkins++ (Imtiaz et al.) paper 
%   and some heuristic observations of my own.

    % Initialisation
    if ~isvector(acc_i)
        error("acceleration must be a row or column vector");
    end
    
    if nargin < 3
        graph = 0; % don't plot by default
    end
    
    acc_i = acc_i(:); % vectorize
    
    delay = 0;
    skip = 0; % becomes 1 when a T wave is detected
    m_selected_RR = 0;
    mean_RR_int = 0;
    ser_back = 0; 
    
    ax = zeros(1,6);
    
    % Bandpass filter
    [b, a] = butter(3, [5 18]/(f/2), 'bandpass');
    acc_bf = filtfilt(b, a, acc_i);
    acc_bf = acc_bf / max(abs(acc_bf));
    
    if graph
        figure
        ax(1) = subplot(3, 2, [1 2]);
        plot(acc_i)
        axis tight
        title("Raw signal")
        
        ax(3) = subplot(3, 2, 3);
        plot(acc_bf)
        axis tight
        title("Bandpass filtered")
    end
    
    % Derivative filter
    % H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2)) 
    int_c = (5 - 1)/(f * 1/40);
    b = interp1(1:5, [1 2 0 -2 -1] .* (1/8) * f, 1:int_c:5);
    
    acc_df = filtfilt(b, 1, acc_bf);
    acc_df = acc_df/max(acc_df);
    
    if graph
        ax(4) = subplot(3, 2, 4);
        plot(acc_df)
        axis tight
        title("Derivative filtered")
    end
    
    % Squaring 
    acc_s = acc_df .^ 2;
    
    if graph
        ax(5) = subplot(325);
        plot(acc_s)
        axis tight
        title("Squared")
    end

    % Smoothing with flattop window
    N = round(0.06 * f);
    flattop = flattopwin(N);
    acc_ft = conv(flattop, acc_s);
    
    % Moving average calculation
    acc_ma = conv(acc_ft, ones(1, round(0.150 * f))/round(0.150 * f));
    delay = delay + round(0.150 * f)/2;
    
    if graph
        ax(6) = subplot(3, 2, 6);
        plot(acc_ma)
        axis tight
        title(["MWI", "Noise (black), Adaptive thresh (green), Signal (red), R-peaks (circles)"])
        axis tight
    end
    
    % Fiducial marks (select peaks that are at least 231ms apart) => max
    % 260 bpm
    [pks, locs] = findpeaks(acc_ma, "MinPeakDistance", round(0.231 * f));
    
    % Second initialization
    Npks = length(pks);
    
    qrs_count = zeros(1, Npks);
    qrs_i = zeros(1, Npks);           
    r_peaks_idxs = zeros(1, Npks);       
    r_peaks_amps = zeros(1, Npks);
    
    % Noise and signal buffers
    noise_count_buf = zeros(1, Npks);
    noise_i = zeros(1, Npks);
    SIGL_buf = zeros(1, Npks);
    NOISL_buf = zeros(1, Npks);
    SIGL_buf1 = zeros(1, Npks);
    NOISL_buf1 = zeros(1, Npks);
    THRS_buf1 = zeros(1, Npks);
    THRS_buf = zeros(1, Npks);
    
    % TRAINING PHASE 

    % Initialize using 2 seconds of original signal
    THR1 = max(acc_ma(1:2*f)) * 1/3; % 0.25 of the max amplitude 
    THR2 = mean(acc_ma(1:2*f)) * 1/2; % 0.5 of the mean signal is considered to be noise
    SPK = THR1;
    NPK = THR2;
    
    % Initialise bandpass filter thresholds (2 seconds of bandpass signal)
    THR1_1 = max(acc_bf(1:2*f)) * 1/3; % 0.25 of the max amplitude 
    THR2_1 = mean(acc_bf(1:2*f)) * 1/2; 
    SPK1 = THR1_1; % signal level in bandpass filtered
    NPK1 = THR2_1; % noise level in bandpass filtered
    
    % Threshold and decision rules determination
    beat_count = 0; % raw beats
    beat_count_filt = 0; % filtered beats
    noise_count = 0; % noise counter

    for i = 1:Npks  
        % locate the possible peaks in the bandpass filtered signal 
        if locs(i) - round(0.150 * f) >= 1 && locs(i) <= length(acc_bf)
            [y_i, x_i] = max(acc_bf(locs(i) - round(0.150 * f):locs(i)));
        else
            if i == 1
                [y_i, x_i] = max(acc_bf(1:locs(i)));
            ser_back = 1;
            elseif locs(i) >= length(acc_bf)
                [y_i, x_i] = max(acc_bf(locs(i) - round(0.150 * f):end));
            end       
        end   
    
        % number of detected beats > 8 => adjust based on mean RR interval 
        if beat_count > 8        
            mean_RR_int = mean(diff(qrs_i(beat_count - 8:beat_count))); % calculate RR interval
            cur_RR_int = qrs_i(beat_count) - qrs_i(beat_count - 1); % latest RR
        
            if cur_RR_int <= 0.92 * mean_RR_int || cur_RR_int >= 1.16 * mean_RR_int
                % lower threshold
                THR1 = 0.5 * (THR1);
                THR1_1 = 0.5 * (THR1_1);               
            else
                m_selected_RR = mean_RR_int;
            end        
        end
    
        if m_selected_RR
            test_m = m_selected_RR;
        elseif mean_RR_int && m_selected_RR == 0
            test_m = mean_RR_int;   
        else
            test_m = 0;
        end
    
        if test_m
            if (locs(i) - qrs_i(beat_count)) >= round(1.66 * test_m) % it shows a QRS is missed 
                [pks_temp, locs_temp] = max(acc_ma(qrs_i(beat_count) + round(0.231*f):locs(i) - round(0.231*f))); % search back and locate the max in this interval
                locs_temp = qrs_i(beat_count) + round(0.231*f) + locs_temp - 1; % location 
                
                if pks_temp > THR2
                    beat_count = beat_count + 1;
                    qrs_count(beat_count) = pks_temp;
                    qrs_i(beat_count) = locs_temp;
                
                    % locate R peak in filtered signal
                    if locs_temp <= length(acc_bf)
                        [y_i_t, x_i_t] = max(acc_bf(locs_temp - round(0.150 * f):locs_temp));
                    else
                        [y_i_t, x_i_t] = max(acc_bf(locs_temp - round(0.150 * f):end));
                    end
                    
                    % set band pass threshold
                    if y_i_t > THR2_1 
                        beat_count_filt = beat_count_filt + 1;
                        r_peaks_idxs(beat_count_filt) = locs_temp - round(0.150 * f) + (x_i_t - 1); % save index of bandpass 
                        r_peaks_amps(beat_count_filt) = y_i_t; % save amplitude of bandpass 
                        SPK1 = 0.25*y_i_t + 0.75*SPK1; % when found with the second thres 
                    end
                    
                    SPK = 0.25*pks_temp + 0.75*SPK; % when found with the second threshold             
                end             
            end
        end
            
        % find noise and QRS peaks
        if pks(i) >= THR1      
            if beat_count >= 3
                % if no QRS in 360ms of the previous QRS see if T wave
                if (locs(i) - qrs_i(beat_count)) <= round(0.360 * f)
                  cur_slope = mean(diff(acc_ma(locs(i) - round(0.075 * f):locs(i)))); % mean slope of current R wave
                  prev_slope = mean(diff(acc_ma(qrs_i(beat_count) - round(0.075 * f):qrs_i(beat_count)))); % mean slope of previous R wave

                  if abs(cur_slope) <= abs(0.6 * (prev_slope)) % slope less then 0.6 of previous R
                     noise_count = noise_count + 1;
                     noise_count_buf(noise_count) = pks(i);
                     noise_i(noise_count) = locs(i);
                     
                     skip = 1; % T wave identification
                     
                     % adjust noise levels
                     SPK = 0.125*y_i + 0.875*SPK;
                     SPK1 = 0.125*y_i + 0.875*SPK1;

                     NPK1 = 0.125*y_i + 0.875*NPK1;
                     NPK = 0.125*pks(i) + 0.875*NPK; 
                  else
                     skip = 0;
                  end
                end
            end
        
        % skip is 1 when a T wave is detected
        if skip == 0    
            beat_count = beat_count + 1;
            qrs_count(beat_count) = pks(i);
            qrs_i(beat_count) = locs(i);
        
            % check bandpass filter threshold
            if y_i >= THR1_1  
              beat_count_filt = beat_count_filt + 1;
              
              if ser_back 
                 r_peaks_idxs(beat_count_filt) = x_i; % save index of bandpass 
              else
                 r_peaks_idxs(beat_count_filt) = locs(i) - round(0.150*f) + (x_i - 1); % save index of bandpass 
              end
              r_peaks_amps(beat_count_filt) = y_i; % save amplitude of bandpass 
              SPK1 = 0.125*y_i + 0.875*SPK1; % adjust threshold for bandpass filtered sig
            end
           
            SPK = 0.125*pks(i) + 0.875*SPK; % adjust signal level
        end
          
        elseif (THR2 <= pks(i)) && (pks(i) < THR1)
            NPK1 = 0.125*y_i + 0.875*NPK1; % adjust Noise level in filtered sig
            NPK = 0.125*pks(i) + 0.875*NPK; % adjust Noise level in MVI       
        
        elseif pks(i) < THR2
            noise_count = noise_count + 1;
            noise_count_buf(noise_count) = pks(i);
            noise_i(noise_count) = locs(i);    
            NPK1 = 0.125*y_i + 0.875*NPK1; % noise level in filtered signal    
            NPK = 0.125*pks(i) + 0.875*NPK; % adjust Noise level in MVI     
        end
           
        % adjust the threshold with SNR for original signal
        if NPK ~= 0 || SPK ~= 0
            THR1 = NPK + 0.25*(abs(SPK - NPK));
            THR2 = 0.4 * (THR1);
        end
        
        % ad hoc adjustements (else we have dimensionality issues)
        if isempty(NPK1); NPK1 = 0; end
        if isempty(SPK1); SPK1 = 0; end

        % adjust the threshold with SNR for bandpassed signal
        if NPK1 ~= 0 || SPK1 ~= 0
            THR1_1 = NPK1 + 0.25*(abs(SPK1 - NPK1));
            THR2_1 = 0.4 * (THR1_1);
        end
        
        % track thresholds of smoothed signal
        SIGL_buf(i) = SPK;
        NOISL_buf(i) = NPK;
        THRS_buf(i) = THR1;
        
        % track thresholds of filtered signal
        SIGL_buf1(i) = SPK1;
        NOISL_buf1(i) = NPK1;
        THRS_buf1(i) = THR1_1;
        
        % reset parameters
        skip = 0;                                                   
        ser_back = 0;    
    end
    
    % length adjustement
    r_peaks_idxs = r_peaks_idxs(1:beat_count_filt);
    r_peaks_amps = r_peaks_amps(1:beat_count_filt);
    qrs_count = qrs_count(1:beat_count);
    qrs_i = qrs_i(1:beat_count);

    % Final adjustement
    % delete peaks whose amplitude is too low
    mean_amp = mean(r_peaks_amps);
    min_amp = mean_amp / 2; % the denominator was choosen heuristically

    to_delete = find(r_peaks_amps < min_amp);

    r_peaks_idxs(to_delete) = [];
    r_peaks_amps(to_delete) = [];

    % delete peaks that are too close
    distances = diff(r_peaks_idxs); % they're all > 0 and increasing
    mean_dist = mean(distances);
    min_dist = mean_dist / 1.75; % the denominator was choosen heuristically

    [~, m] = max(r_peaks_amps); % position of max peak => almost certainly an R-peak

    to_delete = zeros(size(r_peaks_idxs)); % to store indexes to delete
    % rightward traversal
    for r = m:length(r_peaks_idxs) - 1
        cur_dist = abs(r_peaks_idxs(r) - r_peaks_idxs(r + 1));
        if cur_dist <= min_dist
            if r_peaks_amps(r) > r_peaks_amps(r + 1)
                to_delete(r) = r + 1;
            else
                to_delete(r) = r;
            end
        end
    end

    % leftward traversal
    for l = m:-1:2
        cur_dist = abs(r_peaks_idxs(l) - r_peaks_idxs(l - 1));
        if cur_dist <= min_dist
            if r_peaks_amps(l) > r_peaks_amps(l - 1)
                to_delete(l) = l - 1;
            else
                to_delete(l) = l;
            end
        end
    end

    to_delete = nonzeros(to_delete);

    r_peaks_idxs(to_delete) = [];
    r_peaks_amps(to_delete) = [];

    % Final plot
    if graph
        hold on
        scatter(qrs_i, qrs_count, "m")
        
        hold on
        plot(locs, NOISL_buf, "--k", "LineWidth", 2)
        
        hold on
        plot(locs, SIGL_buf, "--r", "LineWidth", 2)
        
        hold on
        plot(locs, THRS_buf, "--g", "LineWidth", 2)
        
        if any(ax)
            ax(~ax) = []; 
            linkaxes(ax, "x")
            zoom on
        end

        figure
        az(1) = subplot(3, 1, 1);
        plot(acc_bf)
        title("R-peaks on bandpass filtered signal")
        axis tight

        hold on
        scatter(r_peaks_idxs, r_peaks_amps, "m")
        hold on
        plot(locs, NOISL_buf1, "LineWidth", 2, "Linestyle", "--", "color", "k")
        hold on
        plot(locs, SIGL_buf1, "LineWidth", 2, "Linestyle", "-.", "color", "r")
        hold on
        plot(locs, THRS_buf1, "LineWidth", 2, "Linestyle", "-.", "color", "g")

        az(2) = subplot(3, 1, 2);
        plot(acc_ma)
        title(["R-peaks on MWI signal", "Noise (black), Signal (red) and Adaptive thresh (green)"])
        axis tight
        
        hold on
        scatter(qrs_i, qrs_count, "m")
        hold on
        plot(locs, NOISL_buf, "LineWidth", 2, "Linestyle", "--", "color", "k")
        hold on
        plot(locs, SIGL_buf, "LineWidth", 2, "Linestyle", "-.", "color", "r")
        hold on
        plot(locs, THRS_buf, "LineWidth", 2, "Linestyle", "-.", "color", "g")
        
        az(3) = subplot(3, 1, 3);
        plot(acc_i - mean(acc_i))
        title("Detected pulse train");
        axis tight
        line(repmat(r_peaks_idxs, [2 1]), ...
            repmat([min(acc_i - mean(acc_i))/2; max(acc_i - mean(acc_i))/2], size(r_peaks_idxs)), ...
            "LineWidth", 2.5, "LineStyle", "-.", "Color", "r");
        linkaxes(az, "x")
        zoom on
    end
end
 








