function hr = estimate_HR(acc_i, T, f, graph)
%estimate_HR Estimates the heart rate given the ith acceleration component, 
%the final time point T and the frequency f.
%
%   The graph boolean argument indicates whether one wants to graph the
%   intermediate steps or not (during the Pan-Tompkins procedure).

    if nargin < 4
        graph = 0; % no graphing by default
    end

    [~, r_peaks_idxs, ~] = pan_tompkins(acc_i, f, graph);
    hr = round(length(r_peaks_idxs)/T * 60);
end
