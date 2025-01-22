using FFTW
using DataFrames
using CSV
using Plots
plotly()

function create_interactive_plot(conv_fac_atm, conv_fac_oce, t, z_atm, z_oce)
    # Create grids for plotting
    z = vcat(z_oce, z_atm)

    # Concatenate the 2D temperature arrays for both atmosphere and ocean
    combined_data = vcat(conv_fac_oce, conv_fac_atm)  # Concatenate vertically (combine atmosphere and ocean data)

    println(size(combined_data))
    println(size(t))
    println(size(z))
    # Create the surface plot using PlotlyJS for interactivity
    plot_surface = Plots.surface(
        t,  # Time data for x-axis
        z,  # Depth data for y-axis
        combined_data,  # Temperature data for z-axis
        label="Conv. Factor",  # Label for the plot
        xlabel="Time",        # x-axis label
        ylabel="Depth",       # y-axis label
        zlabel="Temperature", # z-axis label
        title="Atmosphere and Ocean - Time vs Depth",  # Plot title
    )

    # Display the plot
    display(plot_surface)
end

function main()
    # Get the list of files in the directories
    atm_files = readdir("atm_temps/")
    oce_files = readdir("oce_temps/")

    # Initialize lists to store the data frames
    atm_data_frames = []
    oce_data_frames = []
    times = []
    atm_z_range = []
    oce_z_range = []

    # Read the ATM temperature files and store the data frames
    for (i, atm_file) in enumerate(atm_files)
        file_path = joinpath("atm_temps", atm_file)
        df = CSV.read(file_path, DataFrame)
        if i == 1
            times = parse.(Float64, names(df)[2:end])  # Convert column names to Float64 if they are numeric strings
            atm_z_range = df[:, 1]
        end
        temps = Matrix(df[:, 2:end])
        push!(atm_data_frames, temps)
    end

    # Read the OCE temperature files and store the data frames
    for oce_file in oce_files
        file_path = joinpath("oce_temps", oce_file)
        df = CSV.read(file_path, DataFrame)
        if oce_z_range == []
            oce_z_range = df[:, 1]
            # print(oce_z_range)
        end
        temps = Matrix(df[:, 2:end])
        push!(oce_data_frames, temps)
    end

    # Convert the data frames to 3D arrays by stacking them along the third axis
    data_atm = cat(atm_data_frames..., dims=3)
    data_oce = cat(oce_data_frames..., dims=3)

    # Print the shape of the resulting 3D arrays
    # println(size(data_atm))
    # println(size(data_oce))
    opt_sol_atm = data_atm[:, :, end]
    opt_sol_oce = data_oce[:, :, end]

    error_atm = 0
    error_oce = 0
    fourier_error_atm = 0
    fourier_error_oce = 0
    conv_fac_atm_real = []
    conv_fac_oce_real = []
    conv_fac_atm_fourier = []
    conv_fac_oce_fourier = []
    for i in 1:length(data_atm[1, 1, :])
        pre_error_atm = error_atm
        pre_error_oce = error_oce
        error_atm = data_atm[:, :, i] - opt_sol_atm
        error_oce = data_oce[:, :, i] - opt_sol_oce

        pre_fourier_error_atm = fourier_error_atm
        pre_fourier_error_oce = fourier_error_oce
        fourier_error_atm = hcat(map(fft, eachrow(error_atm))...)'
        fourier_error_oce = hcat(map(fft, eachrow(error_oce))...)'
        if i > 1
            # Adding real convergence factors
            # print(size(error_atm ./ pre_error_atm))
            push!(conv_fac_atm_real, abs.(error_atm ./ (pre_error_atm .+ 1)))
            push!(conv_fac_oce_real, abs.(error_oce ./ (pre_error_oce .+ 1)))

            # Adding fourier convergence factors, to compare to analysis.
            push!(conv_fac_atm_fourier, abs.(fourier_error_atm ./ (pre_fourier_error_atm .+ 1)))
            push!(conv_fac_oce_fourier, abs.(fourier_error_oce ./ (pre_fourier_error_oce .+ 1)))
        end
    end
    create_interactive_plot(conv_fac_atm_real[1], conv_fac_oce_real[1], times, atm_z_range, oce_z_range)

end

# Run the script
main()