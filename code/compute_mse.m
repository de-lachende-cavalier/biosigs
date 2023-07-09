function [mse_hr, mse_rr] = compute_mse(data, f, ref_hr, ref_rr, max_n_imfs, rr_method)
%COMPUTE_MSE Computes the mse obtained from the data, given the signal
%data, the frequency f, reference HR and RR, the maximum number of IMFs desired 
%and the method to utilize for RR estimation.
%
%   This function only operates on single cells (i.e. on single data sources). 
%   This decision has been made to mantain generality and keep the code DRY.
%   To check the mse for multiple sources, then, it is only sufficient
%   to run it on each cell and then compare the accuracies (e.g. compare_accuracies.m).

    if iscell(data) && size(data, 1) == 1
        if size(ref_hr, 2) ~= size(data, 2)
            error("Not enough reference values for HR: expected %d, got %d.", ...
                size(data, 2), size(ref_hr, 2))
        end
    
        if size(ref_rr, 2) ~= size(data, 2)
            error("Not enough reference values for RR: expected %d, got %d.", ...
                size(data, 2), size(ref_rr, 2))
        end
    
        hr_ests = zeros(1, size(data, 2));
        rr_ests = zeros(1, size(data, 2));
    
        for i = 1:size(data, 2)
            hr_ests(i) = estimate_HR(data{i}.accel_z, data{i}.time(end), f);
            rr_ests(i) = estimate_RR(data{i}.accel_y, data{i}.time, f, max_n_imfs, rr_method);
        end
        
        N = length(ref_hr);
        mse_hr = sum((hr_ests - ref_hr).^2) / N;
        mse_rr = sum((rr_ests - ref_rr).^2) / N;
    end
end

