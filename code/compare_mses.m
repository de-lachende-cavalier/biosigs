clearvars; close all; clc; warning("off", "all") % to ignore the 'VariableNamingRule' warning

%% Setup
f_resample = 100; % Hz

base = "../data"; 
specific = ["vert_xiphoid/", "vert_chest/", "vert_stomach/", ...
    "hor_xiphoid/", "hor_chest/", "hor_stomach/"];
file = "";

ref_rrs = repmat(8, length(specific), 3); % counted manually, attempted to mantain it constant throughout the measurements
ref_hrs = [[57, 58, 61]; ... % read from the pulse oximeter
           [54, 53, 58]; ...
           [56, 59, 57]; ...
           [57, 58, 58]; ...
           [54, 51, 59]; ...
           [59, 57, 60]];


rr_methods = ["fft", "estrada"];
max_n_imfs = 6;

%% MSE estimation
for rm = 1:length(rr_methods)
    fprintf("[+] rr_method = '%s'\n", rr_methods(rm))
    for i = 1:length(specific)
        accel_data = get_data(fullfile(base, specific(i), file));
        clean_data = preprocess_data(accel_data, size(accel_data, 2), f_resample);
        
        [mse_hr, mse_rr] = compute_mse(clean_data, f_resample, ref_hrs(i, :), ref_rrs(i, :), max_n_imfs, rr_methods(rm));
        fprintf("\t%s:", specific(i))
        fprintf("\n\t\tHR MSE: %.4f\n", mse_hr);
        fprintf("\t\tRR MSE: %.4f\n", mse_rr);
    end
    fprintf("\n")
end
% => best when "vert_xiphoid/" and rr_method = "fft" ("estrada" is still pretty close, but "fft" seems more robust)


