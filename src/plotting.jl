using PlotlyJS
using Plotly
using CSV
using DataFrames
using FFMPEG

# data = CSV.read("atm_temps/iter_1.csv", DataFrame, header=false)
# print(data)

# file = open("atm_temps/iter_1.csv", "r")
# println(read(file, Char, 100))  # Read first 100 characters
# close(file)
atm_files = readdir("atm_temps/")
oce_files = readdir("oce_temps/")

num_iters = length(atm_files)
atm_data_frames = []
oce_data_frames = []
for atm_file in atm_files
    file_path = joinpath(@__DIR__, "atm_temps/", atm_file)
    df = CSV.read(file_path, DataFrame)
    push!(atm_data_frames, df)
end

for oce_file in oce_files
    file_path = joinpath(@__DIR__, "oce_temps/", oce_file)
    df = CSV.read(file_path, DataFrame)
    push!(oce_data_frames, df)
end

data_atm = cat([Matrix(df) for df in atm_data_frames]..., dims=3)
data_oce = cat([Matrix(df) for df in oce_data_frames]..., dims=3)

# Assuming 'data_atm' is your 3D data matrix (e.g., size (20003, 200, num_frames))
# We'll create a 3D surface plot for a specific frame (e.g., the first frame)

z = data_atm[:, :, 1]  # Extract the first 2D slice/frame from the 3D matrix

# Create the x and y grid
x = 1:size(z, 1)
y = 1:size(z, 2)
surface = Surface(x, y, z)

# # Create meshgrid (this can be done using a broadcasting operation)
# X, Y = meshgrid(x, y)

# Create the plot
surface = plot(
    [surface],
    type="surface",  # Make it a surface plot
    colorscale="Viridis",  # Color scale for the surface
    layout=Layout(
        title="3D Surface Plot of Temperature",  # Plot title
        scene=attr(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Temperature"
        )
    )
)

# Show the plot
display(surface)

# # using HDF5
# using Plots;
# using CSV
# using DataFrames
# using FFMPEG

# # data = CSV.read("atm_temps/iter_1.csv", DataFrame, header=false)
# # print(data)

# # file = open("atm_temps/iter_1.csv", "r")
# # println(read(file, Char, 100))  # Read first 100 characters
# # close(file)
# atm_files = readdir("atm_temps/")
# oce_files = readdir("oce_temps/")

# num_iters = length(atm_files)
# atm_data_frames = []
# oce_data_frames = []
# for atm_file in atm_files
#     file_path = joinpath(@__DIR__, "atm_temps/", atm_file)
#     df = CSV.read(file_path, DataFrame)
#     push!(atm_data_frames, df)
# end

# for oce_file in oce_files
#     file_path = joinpath(@__DIR__, "oce_temps/", oce_file)
#     df = CSV.read(file_path, DataFrame)
#     push!(oce_data_frames, df)
# end

# data_atm = cat([Matrix(df) for df in atm_data_frames]..., dims=3)
# data_oce = cat([Matrix(df) for df in oce_data_frames]..., dims=3)

# println(size(data_atm))
# println(size(data_atm))
# println(size(data_atm))
# # println(data_atm)
# num_frames = size(data_atm, 3)

# # Create a grid for plotting
# x = 1:size(data_atm, 1)  # Rows
# y = 1:size(data_atm, 2)  # Columns
# z = data_atm[:, :, 1]
# print(size(x))
# print(size(y))
# print(size(z))
# surface(x, y, z, color=:viridis, xlabel="X", ylabel="Y", zlabel="Temperature")

# Generate 2D grids for X and Y
# Initialize the video
# anim = @animate for i in 1:num_frames
#     z = data_atm[:, :, i]  # Extract the i-th frame (2D matrix)
#     print(size(z))
#     surf(x, y, z, color=:viridis, xlabel="X", ylabel="Y", zlabel="Temperature")
# end

# Save the animation to a video file
# gif(anim, "output_animation.gif", fps=10)
# println("Saved video to $mp4_path")



# # List all the files in the directory
# file_paths = readdir("./checkpoint/")

# # Separate the files into "Atmos" and "Ocean"
# atm_files = filter(x -> occursin("Atmos", x), file_paths)
# oce_files = filter(x -> occursin("Ocean", x), file_paths)

# # Sort each list by the number after the last underscore
# # Sort by both extracted numbers (first by the first number, then by the second number)

# # Define the method to extract and parse numbers from the filename
# function extract_numbers(filename::String)
#     m = match(r"_(\d+)_(\d+)\.hdf5$", filename)
#     return (parse(Int, split(m.match, "_")[2]), parse(Int, split(m.match, "_")[3][1:end-5]))    # Parse both numbers from the match and return them as a tuple
# end

# # Now use the new method in the sorting
# sorted_atm_files = sort(atm_files, by=extract_numbers)
# sorted_oce_files = sort(oce_files, by=extract_numbers)

# # Print the sorted file paths
# println(sorted_atm_files)
# println(sorted_oce_files)

# oce_temps_list = []
# atm_temps_list = []
# for atm_file in sorted_atm_files
#     file_path = joinpath(@__DIR__, "oce_temps", atm_file)
#     file = h5open(file_path)
#     data = read(file)
#     atm_temps = data["fields"]["model_state"]["atm"]
#     push!(atm_temps_list, atm_temps)
# end
# for oce_file in sorted_oce_files
#     file_path = joinpath(@__DIR__, "checkpoint", oce_file)
#     file = h5open(file_path)
#     data = read(file)
#     oce_temps = data["fields"]["model_state"]["oce"]
#     push!(oce_temps_list, oce_temps)
# end

# print(atm_temps_list)
