using CSV
using DataFrames
using Plots

# Load the CSV files
P_diff_matrix = CSV.read("Bode_Figures_Paper/P_diff_matrix_eta(2.0e-9)_burst(150)_lysis_timer(250).csv", DataFrame)
time_df = CSV.read("Bode_Figures_Paper/time_eta(2.0e-9)_burst(150)_lysis_timer(250).csv", DataFrame)
api_df = CSV.read("Bode_Figures_Paper/api_burst(150)_lysis_timer(250).csv", DataFrame)

# Define the API values to plot
api_values = [2, 5, 10, 15,20, 25]

# Initialize the plot
p = plot(xlabel="Time (minutes)", ylabel="Phage flux")

# Loop through each API value
for chosen_api in api_values
    # Find the index of the chosen API
    api_index = findfirst(x -> abs(x - chosen_api) < 1e-6, api_df[:, 1])

    if isnothing(api_index)
        println("API value $chosen_api not found in the DataFrame.")
        continue
    end

    # Extract the corresponding P_diff and time values
    P_diff = collect(P_diff_matrix[api_index, :])
    time_values = time_df[:, 1]  # Assuming time values are in the first column

    # Add the series to the plot
    plot!(p, time_values, P_diff, label="API=$chosen_api")
end

# Add annotation "(b)" at the top left of the plot
#annotate!(p, 0.01, 0.98, text("(b)", :left, 12, :black), :relative)
annotate!(p, time_df[1, 1], 6500000, text("(b)", :left, 12, :black))


# Save the plot to the specified file path
savefig("Bode_Figures_Paper/FixedAPI_eta(2.0e-9)_burst(150)_lysis_timer(250).pdf")