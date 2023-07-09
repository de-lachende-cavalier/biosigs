clearvars; close all; clc; warning("off", "all") % to ignore the 'VariableNamingRule' warning

%% Setup
base = "../data"; % base directory for the data
specific = "vert_xiphoid/hrv/"; % specific subdirectory
file = ""; % file of interest (if we desire a specific file)

%% Preprocess data
accel_data = get_data(fullfile(base, specific, file));

f_resample = 100; % Hz
clean_data = preprocess_data(accel_data, size(accel_data, 2), f_resample);

%% Plot
% plot_data(size(clean_data, 2), clean_data)

%% HR and RR estimation
% i = 5;
% hr_est = estimate_HR(clean_data{i}.accel_z, clean_data{i}.time(end), f_resample, 1);
% 
% max_n_imfs = 6;
% [rr_est, resp_imfs] = estimate_RR(clean_data{i}.accel_y, clean_data{i}.time, f_resample, max_n_imfs, "fft", 0);

%% HRV 
% Useful references: https://en.wikipedia.org/wiki/Heart_rate_variability
% and the Landreani et al. paper

% ms5 is an outlier (and mostly ruins the boxplot)!
% we can see it by uncommenting the following two lines (and commenting the
% two below)
% ms_data = clean_data(:, 1:5);
% rest_data = clean_data(:, 6:end);

ms_data = clean_data(:, 1:4);
rest_data = clean_data(:, 6:end - 1);

rmssds_ms = zeros(size(ms_data, 2), 1);
sdnns_ms = zeros(size(ms_data, 2), 1);
for i = 1:size(ms_data, 2)
    [~, r_peaks_idxs_ms, ~] = pan_tompkins(ms_data{i}.accel_z, f_resample, 0);

    rr_series_ms = diff(r_peaks_idxs_ms);
    
    rmssds_ms(i) = sqrt(sum(rr_series_ms .^ 2) / length(rr_series_ms));
    % technically, NN intervals should only contain normal R-peaks but practically there's 
    % no real difference between RR and NN intervals:
    % https://psychology.stackexchange.com/questions/16076/what-is-the-difference-between-rr-intervals-and-nn-intervals-in-hrv-data
    % also Muaremi et al. (especially section 3.3.2)
    sdnns_ms(i) = std(rr_series_ms);
    
%     x = linspace(1, length(r_peaks_idxs_ms) - 1, length(rr_series_ms));
%     figure
%     plot(x, rr_series_ms)
end

rmssds_rest = zeros(size(rest_data, 2), 1);
sdnns_rest = zeros(size(rest_data, 2), 1);
for i = 1:size(rest_data, 2)
    [~, r_peaks_idxs_rest, ~] = pan_tompkins(rest_data{i}.accel_z, f_resample, 0);

    rr_series_rest = diff(r_peaks_idxs_rest);
    
    rmssds_rest(i) = sqrt(sum(rr_series_rest .^ 2) / length(rr_series_rest));
    sdnns_rest(i) = std(rr_series_rest);
    
%     x = linspace(1, length(r_peaks_idxs_rest) - 1, length(rr_series_rest));
%     figure
%     plot(x, rr_series_rest)
end

rmssds = [rmssds_rest, rmssds_ms];
figure
boxplot(rmssds, 'Labels', {'REST', 'MS'})
title("RMSSDs")

sdnss = [sdnns_rest, sdnns_ms];
figure
boxplot(sdnss, 'Labels', {'REST', 'MS'})
title("SDNNs")