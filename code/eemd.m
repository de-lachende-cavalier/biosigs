function eemd_imfs = eemd(clean_data_ii, Nt, max_n_imfs)
%EEMD Implements the EEMD algorithm, given the ith acceleration component of a data point,
%the number of time steps and the maximum number of IMFs desired.
%
%   The hyperparameters used have mainly been adapted from Xinmiao et al.
    
    std = 0.2;
    N = 100;
    wgn = std * randn(Nt, 100);
    
    noised_acc_ii = clean_data_ii + wgn;
    
    imfs = zeros(Nt, max_n_imfs, N);
    for i = 1:N
        % we don't care about residuals
        [imfs(:, :, i), ~] = emd(noised_acc_ii(:, i), 'MaxNumIMF', max_n_imfs);
    end
    
    eemd_imfs = mean(imfs, 3);
end

