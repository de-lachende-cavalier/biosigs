function [rr, imfs] = estimate_RR(acc_i, t_axis, f, max_n_imfs, method, graph)
%ESTIMATE_RR Estimates RR given the ith component of acceleration, the time
%axis, the sampling frequency, the maximum number of IMFs desired and the
%method to utilize.
%
%   The final boolean argument (graph) is used to tell the function whether it should 
%   plot the intermediate steps or not.

    if nargin < 6
        graph = 0; % don't plot by default
    end

    Nt = length(t_axis);
    imfs = eemd(acc_i, Nt, max_n_imfs);    

    if graph
        plot_imfs(imfs)
    end

    % The last IMF clearly isolates breathing (for "vert_xiphoid", accel_y
    % and max_n_imfs = 6)
    resp_sig = imfs(:, end);

    if graph
        figure
        subplot(2, 1, 1);
        plot(t_axis, resp_sig)
        axis tight
        title("Respiratory signal")
    end
    
    if method == "fft"
        f_axis = (0:Nt - 1)/Nt * f;
        resp_fft = fft(resp_sig);
        spect = abs(resp_fft);

        target_freqs = f_axis(spect >= floor(max(spect)));
        rr = round(target_freqs(1) * 60);

        if graph
            subplot(2, 1, 2);
            plot(f_axis - f/2, fftshift(spect))
            axis tight
            title("Respiratory signal spectrum")
        end
        return

    elseif method == "estrada"
        [~, locs] = findpeaks(resp_sig); 
        min_peak_dist = mean(diff(locs)) / 2; % to mitigate mode mixing

        [~, locs] = findpeaks(resp_sig, 'MinPeakDistance', min_peak_dist);
        peak_times = t_axis(locs);
        
        est_breath_time = mean(diff(peak_times));
        rr = round(60 / est_breath_time);

        if graph
            subplot(2, 1, 2);
            findpeaks(resp_sig, 'MinPeakDistance', min_peak_dist)
            grid off
            axis tight
            title("Respiratory signal peaks")
        end
        return
        
    else
        error("The method was not recognized (available ones are 'fft' and 'estrada').")
    end   
end

