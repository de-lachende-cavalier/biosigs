function plot_data(n_data_sources, data)
%PLOT_DATA Plots the data given the number of data sources and the data itself.

    if n_data_sources == 1
        % not very elegant...
        figure
        
        subplot(3, 1, 1)
        plot(data.time, data.accel_x)
        title("Acceleration in the x direction")
        xlabel("$t\ [s]$", 'Interpreter','latex')
        ylabel("$a\ [\frac{m}{{s}^2}]$", 'Interpreter','latex')
        
        subplot(3, 1, 2)
        plot(data.time, data.accel_y)
        title("Acceleration in the y direction")
        xlabel("$t\ [s]$", 'Interpreter','latex')
        ylabel("$a\ [\frac{m}{{s}^2}]$", 'Interpreter','latex')
    
        subplot(3, 1, 3)
        plot(data.time, data.accel_z)
        title("Acceleration in the z direction")
        xlabel("$t\ [s]$", 'Interpreter','latex')
        ylabel("$a\ [\frac{m}{{s}^2}]$", 'Interpreter','latex')

        return
    end

    for d = 1:n_data_sources
        figure
        
        subplot(3, 1, 1)
        plot(data{d}.time, data{d}.accel_x)
        title("Acceleration in the x direction")
        xlabel("$t\ [s]$", 'Interpreter','latex')
        ylabel("$a\ [\frac{m}{{s}^2}]$", 'Interpreter','latex')
        
        subplot(3, 1, 2)
        plot(data{d}.time, data{d}.accel_y)
        title("Acceleration in the y direction")
        xlabel("$t\ [s]$", 'Interpreter','latex')
        ylabel("$a\ [\frac{m}{{s}^2}]$", 'Interpreter','latex')
    
        subplot(3, 1, 3)
        plot(data{d}.time, data{d}.accel_z)
        title("Acceleration in the z direction")
        xlabel("$t\ [s]$", 'Interpreter','latex')
        ylabel("$a\ [\frac{m}{{s}^2}]$", 'Interpreter','latex')
    end
end

