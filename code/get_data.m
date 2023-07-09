function data = get_data(given_path)
%GET_DATA Returns the data available at a certain path.
%
%   The path can indicate a directory, in which case all the data found in
%   it (i.e. all the .csv files) will be collected and returned as a
%   cell of tables, or can contain a filename, in which case only the data
%   in that file will be read.
%
%   N.B.: When only a directory is specified the search doesn't proceed
%   recursively (i.e. the subdirectories of the directory won't be checked,
%   even if they may contain other .csv files)

    [filepath, ~, ext] = fileparts(given_path);
    
    if (ext ~= "")
        % have been given a full path (i.e. a specific .csv file)
        data = readtable(given_path);
        
        % for easier access later
        data = renamevars(data, ...
            ["Time_s_", "AccelerationX_m_s_2_", "AccelerationY_m_s_2_", "AccelerationZ_m_s_2_"], ...
            ["time", "accel_x", "accel_y", "accel_z"]);
        return
    end
    
    cur_dir = pwd; % so that we return to proper directory once the function is done
    
    cd(filepath)
    csv_files = dir("*.csv");
    
    data = cell(1, length(csv_files));
    for i = 1:length(csv_files)
        data{i} = readtable(csv_files(i).name);

        data{i} = renamevars(data{i}, ...
            ["Time_s_", "AccelerationX_m_s_2_", "AccelerationY_m_s_2_", "AccelerationZ_m_s_2_"], ...
            ["time", "accel_x", "accel_y", "accel_z"]);
    end   

    cd(cur_dir)
end