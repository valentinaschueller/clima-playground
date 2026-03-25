module SeaIce

import ClimaCore.MatrixFields: @name
import clima_playground: SimulationParameters, get_diagnostic, get_prognostic_data!

include("ice/shared.jl")
include("ice/constant.jl")
include("ice/temperature_feedback.jl")
include("ice/thickness_feedback.jl")

end