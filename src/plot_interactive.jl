using CSV
using DataFrames
using Plots
plotly()

function create_interactive_plot(data_atm, data_oce, t, z_atm, z_oce)
    # Create grids for plotting
    z = vcat(z_oce, z_atm)

    for i in 1:size(data_atm, 3)
        # Extract the data for the current time step
        atm_data = data_atm[:, :, i]
        oce_data = data_oce[:, :, i]

        # Concatenate the 2D temperature arrays for both atmosphere and ocean
        combined_data = vcat(oce_data, atm_data)  # Concatenate vertically (combine atmosphere and ocean data)

        # Create the surface plot using PlotlyJS for interactivity
        plot_surface = Plots.surface(
            t,  # Time data for x-axis
            z,  # Depth data for y-axis
            combined_data,  # Temperature data for z-axis
            label="Temperature",  # Label for the plot
            xlabel="Time",        # x-axis label
            ylabel="Depth",       # y-axis label
            zlabel="Temperature", # z-axis label
            title="Atmosphere and Ocean - Time vs Depth",  # Plot title
        )

        # Display the plot
        display(plot_surface)
    end
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
            print(oce_z_range)
        end
        temps = Matrix(df[:, 2:end])
        push!(oce_data_frames, temps)
    end

    # Convert the data frames to 3D arrays by stacking them along the third axis
    data_atm = cat(atm_data_frames..., dims=3)
    data_oce = cat(oce_data_frames..., dims=3)

    # Print the shape of the resulting 3D arrays
    println(size(data_atm))
    println(size(data_oce))

    create_interactive_plot(data_atm, data_oce, times, atm_z_range, oce_z_range)
end

# Run the script
main()
