using CSV
using DataFrames
using Plots
using LaTeXStrings
using GLM

# Read the CSV file
file_path = "figure_files_space_samples/prob_survive_lysis_inhibition(false)_LO(false)_LIC(false).csv"
dataFFF = CSV.read(file_path, DataFrame)
file_path = "figure_files_space_samples/prob_survive_lysis_inhibition(true)_LO(false)_LIC(true).csv"
dataTFT = CSV.read(file_path, DataFrame)
file_path = "figure_files_space_samples/prob_survive_lysis_inhibition(true)_LO(true)_LIC(true).csv"
dataTTT = CSV.read(file_path, DataFrame)
file_path = "figure_files_space_samples/prob_survive_lysis_inhibition(true)_LO(true)_LIC(false).csv"
dataTTF = CSV.read(file_path, DataFrame)
file_path = "figure_files_space_samples/prob_survive_lysis_inhibition(true)_LO(false)_LIC(false).csv"
dataTFF = CSV.read(file_path, DataFrame)


# Extract the threshold and collapse_time columns
BFFF = dataFFF.InitialBacteria
PFFF = dataFFF.Prob_survive
BTFT = dataTFT.InitialBacteria
PTFT = dataTFT.Prob_survive
BTTT = dataTTT.InitialBacteria
PTTT = dataTTT.Prob_survive
BTTF = dataTTF.InitialBacteria
PTTF = dataTTF.Prob_survive
BTFF = dataTFF.InitialBacteria
PTFF = dataTFF.Prob_survive



# Create the plot
colors = ["#0072B2", "#E69F00", "#56B4E9", "#D55E00", "#009E73"]


plot(BTTT, PTTT, seriestype = :scatter, label = "LO, LIN, LINC", xlabel = "Initial uninfected bacteria (cells)", ylabel = "Probability to survive", title = "(f)", legend = :bottomright, color = colors[1], xlims = (0, 24), ylims = (0, 1), size = (600, 600), legendfontsize = 12, guidefontsize = 16, tickfontsize = 16)
plot!(BTTT, PTTT, seriestype = :line, label = "", color=colors[1])

plot!(BTFT, PTFT, seriestype = :scatter, label = "No LO, LIN, LINC", color=colors[2])
plot!(BTFT, PTFT, seriestype = :line, label = "", color=colors[2])


plot!(BTTF, PTTF, seriestype = :scatter, label = "LO, LIN, No LINC", color=colors[3])
plot!(BTTF, PTTF, seriestype = :line, label = "", color=colors[3])


plot!(BFFF, PFFF, seriestype = :scatter, label = "No LO, No LIN, No LINC", color=colors[4])
plot!(BFFF, PFFF, seriestype = :line, label = "", color=colors[4])

plot!(BTFF, PTFF, seriestype = :scatter, label = "No LO, LIN, No LINC", color=colors[5])
plot!(BTFF, PTFF, seriestype = :line, label = "", color=colors[5])






# Save the plot
figures_dir = "figure_files_space_samples"
mkpath(figures_dir)
figure_file_path = joinpath(figures_dir, "prob_survive_comparison.pdf")
savefig(figure_file_path)