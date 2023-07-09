function plot_imfs(imf_array)
%PLOT_IMFS A utility function to quickly visualize the IMFs
%extracted by the EMD procedure.

    figure
    for i = 1:size(imf_array, 2)
        subplot(size(imf_array, 2), 1, i)
        plot(imf_array(:, i))
        title(sprintf("IMF %d", i))
    end

end

