clearvars; close all; clc; warning("off", "all") % to ignore the 'VariableNamingRule' warning

%% Setup
base = "../data"; % base directory for the data
specific = "vert_xiphoid/misc/"; % specific subdirectory
file = "standing.csv"; % file of interest (if we desire a specific file)

% N.B.: The following code assumes we're operating on a single table (i.e. that
%       we have specified a file)!

%% Preprocess data
accel_data = get_data(fullfile(base, specific, file));

f_resample = 100; % Hz
clean_data = preprocess_data(accel_data, 1, f_resample);

%% Plot
plot_data(1, clean_data)

%% Additional preprocessing => remove initial and final parts (they're noise deriving from setting up and removing the phone)
first_noise = 10 * f_resample; % 10s
last_noise = 5 * f_resample; % 5s

clean_data(1:first_noise, :) = [];
clean_data(end:-1:end - last_noise, :) = [];
% adapt the time axis to make it start from 0 and end when it's supposed to
clean_data.time = clean_data.time - clean_data.time(1);

%% Plot
plot_data(1, clean_data)

%% HR estimation
hr_est = estimate_HR(clean_data.accel_z, clean_data.time(end), f_resample, 1);

%% RR estimation
max_n_imfs = 6;
[rr_est, resp_imfs] = estimate_RR(clean_data.accel_y, clean_data.time, f_resample, max_n_imfs, "fft", 1);