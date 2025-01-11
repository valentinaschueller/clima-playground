using HDF5
using Plots

#Define file paths
file_path_atm_0 = joinpath(@__DIR__, "checkpoint", "checkpoint_HeatEquationAtmos_0.hdf5")
file_path_atm_3500 = joinpath(@__DIR__, "checkpoint", "checkpoint_HeatEquationAtmos_3500.hdf5")
file_path_oce_0 = joinpath(@__DIR__, "checkpoint", "checkpoint_HeatEquationOcean_0.hdf5")
file_path_oce_3500 = joinpath(@__DIR__, "checkpoint", "checkpoint_HeatEquationOcean_3500.hdf5")

# Open the HDF5 files
file_atm_0 = h5open(file_path_atm_0, "r")
file_atm_3500 = h5open(file_path_atm_3600, "r")
file_oce_0 = h5open(file_path_oce_0, "r")
file_oce_3500 = h5open(file_path_oce_3600, "r")

# Read data
data_atm_0 = read(file_atm_0)
data_atm_3500 = read(file_atm_3600)
data_oce_0 = read(file_oce_0)
data_oce_3500 = read(file_oce_3600)

# Get temperature profile
temps_atm_0 = data_atm_0["fields"]["model_state"]["atm"]
temps_atm_3500 = data_atm_3500["fields"]["model_state"]["atm"]
temps_oce_0 = data_oce_0["fields"]["model_state"]["oce"]
temps_oce_3500 = data_oce_3500["fields"]["model_state"]["oce"]
println(temps_atm_3500)
println(temps_oce_3500)

# Plot temperature profile
y_atm = range(0, 1, length(temps_atm_3500))
y_oce = range(-1, 0, length(temps_oce_3500))
p = plot(temps_atm_3500, y_atm, xlabel="temperature", ylabel="index", label="atmosphere")
plot!(temps_oce_0, y_oce, label="ocean")
display(p)

# Access a specific dataset
# data = read(file["dataset_name"])  # Replace "dataset_name" with the actual name

# Close the file when done
close(file_atm_0)
close(file_atm_3500)
close(file_oce_0)
close(file_oce_3500)
