function clean_data = preprocess_data(accel_data, n_data_sources, f_resample)
%PREPROCESS_DATA Preprocesses the raw acceleration data and prepares it for
%further computation, given the data itself (as read from the files), the number of data sources, 
%and the resampling frequency.

    g = 9.784; % m/s^2

    if istable(accel_data)
        % we only have one data source
        
        % phyphox doesn't remove g
        accel_data.accel_z = accel_data.accel_z - g;

        % Resample the signal (Estrada et al. suggest 100 Hz as resampling
        % frequency)
        [res_acc_x, res_t_axis] = resample(accel_data.accel_x, accel_data.time, f_resample);
        res_acc_y = resample(accel_data.accel_y, accel_data.time, f_resample);
        res_acc_z = resample(accel_data.accel_z, accel_data.time, f_resample);
    
        clean_data = table('size', [length(res_acc_x), size(accel_data, 2)], ...
            'VariableTypes', repmat("double", size(accel_data, 2), 1), ...
            'VariableNames', accel_data.Properties.VariableNames);
    
        clean_data.time = res_t_axis;
        clean_data.accel_x = res_acc_x;
        clean_data.accel_y = res_acc_y;
        clean_data.accel_z = res_acc_z;
        return
    end
    
    % we have multiple data sources => same procedure as above for each
    % table
    clean_data = cell(size(accel_data));
    for d = 1:n_data_sources 
        accel_data{d}.accel_z = accel_data{d}.accel_z - g;

        [res_acc_x, res_t_axis] = resample(accel_data{d}.accel_x, accel_data{d}.time, f_resample);
        res_acc_y = resample(accel_data{d}.accel_y, accel_data{d}.time, f_resample);
        res_acc_z = resample(accel_data{d}.accel_z, accel_data{d}.time, f_resample);
    
        clean_data{d} = table('size', [length(res_acc_x), size(accel_data{d}, 2)], ...
            'VariableTypes', repmat("double", size(accel_data{d}, 2), 1), ...
            'VariableNames', accel_data{d}.Properties.VariableNames);
    
        clean_data{d}.time = res_t_axis;
        clean_data{d}.accel_x = res_acc_x;
        clean_data{d}.accel_y = res_acc_y;
        clean_data{d}.accel_z = res_acc_z;
    end
end

