clearvars; close all; clc; warning("off", "all") % to ignore the 'VariableNamingRule' warning

%% Setup
fs_phone = 198.9; % Hz (from https://phyphox.org/sensordb/)

base = "../data"; % base directory for the data
specific = "vert_xiphoid/"; % specific subdirectory
file = ""; % file of interest (if we desire a specific file)

%% Preprocess data
accel_data = get_data(fullfile(base, specific, file));

% % (optional) check sampling rate reported by phyphox
% Ts = mean(diff(accel_data{1}.time));
% fs_est = 1 / Ts;
% assert(abs(fs_est - fs_phone) < 0.05)

f_resample = 100; % Hz
clean_data = preprocess_data(accel_data, size(accel_data, 2), f_resample);

%% Plot
% plot_data(1, clean_data{2})

%% HR estimation
i = 2;
hr_est = estimate_HR(clean_data{i}.accel_z, clean_data{i}.time(end), f_resample, 1);

%% RR estimation
max_n_imfs = 6; % this works well for "vert_xiphoid/" and accel_y
[rr_est, resp_imfs] = estimate_RR(clean_data{i}.accel_y, clean_data{i}.time, f_resample, max_n_imfs, "fft", 1);