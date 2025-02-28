using CSV
using DataFrames
using Plots
using LaTeXStrings
using GLM

# Read the CSV file
file_path = "/home/namiko/Projects/LysisInhibitionJulia/agentbasedmodel/Culture_Figures_Paper/LIC_thresholdDependence_initialB_(5.0e7)_initialP_(1.0e7).csv"
data = CSV.read(file_path, DataFrame)

# Extract the threshold and collapse_time columns
threshold = data.threshold
collapse_time = data.collapse_time

# Define the linear fitting line
slope = 2.95
y_intercept = -4.22
extended_threshold = vcat(0, threshold)
linear_fit = slope .* extended_threshold .+ y_intercept

df = DataFrame(threshold = threshold, collapse_time = collapse_time)

# Fit a linear model
linear_model = lm(@formula(collapse_time ~ threshold), df)

# Get the fitted values
fitted_values = predict(linear_model)


# Create the plot
scatter(threshold, collapse_time, label = "", xlabel = L"Threshold $M_{LOR}$", ylabel = "Collapse Time (min)", title = "", xlims = (0, maximum(threshold)), ylims = (0, maximum(collapse_time)))
plot!(extended_threshold, linear_fit, label = "", line = :dash, color = :red)
#plot!(threshold, fitted_values, label = "", line = :dash, color = :red)

# Save the plot
figures_dir = "Culture_Figures_Paper"
mkpath(figures_dir)
savefig(joinpath(figures_dir,"threshold_vs_collapse_time.pdf"))
println(coeftable(linear_model))
