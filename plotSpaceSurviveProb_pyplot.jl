using CSV
using DataFrames
using PyPlot

# Read data
figures_dir = "figure_files_space_phageBC"
smallburst_dir = "figure_files_space_phageBC_smallburst"

function read_data(path)
    CSV.read(path, DataFrame)
end

dataFFF = read_data(joinpath(figures_dir, "prob_survive_lysis_inhibition(false)_LO(false)_LIC(false).csv"))
dataTFT = read_data(joinpath(figures_dir, "prob_survive_lysis_inhibition(true)_LO(false)_LIC(true).csv"))
dataTTT = read_data(joinpath(figures_dir, "prob_survive_lysis_inhibition(true)_LO(true)_LIC(true).csv"))
dataTTF = read_data(joinpath(figures_dir, "prob_survive_lysis_inhibition(true)_LO(true)_LIC(false).csv"))
dataTFF = read_data(joinpath(figures_dir, "prob_survive_lysis_inhibition(true)_LO(false)_LIC(false).csv"))
dataFTT = read_data(joinpath(figures_dir, "prob_survive_lysis_inhibition(false)_LO(true)_LIC(true).csv"))
dataTTTH = read_data(joinpath(smallburst_dir, "prob_survive_lysis_inhibition(true)_LO(true)_LIC(true).csv"))

# Extract columns
BFFF, PFFF = dataFFF.InitialBacteria, dataFFF.Prob_survive
BTFT, PTFT = dataTFT.InitialBacteria, dataTFT.Prob_survive
BTTT, PTTT = dataTTT.InitialBacteria, dataTTT.Prob_survive
BTTF, PTTF = dataTTF.InitialBacteria, dataTTF.Prob_survive
BTFF, PTFF = dataTFF.InitialBacteria, dataTFF.Prob_survive
BFTT, PFTT = dataFTT.InitialBacteria, dataFTT.Prob_survive
BTTTH, PTTTH = dataTTTH.InitialBacteria, dataTTTH.Prob_survive

# Colors (Okabe–Ito)
colors = ["#0072B2", "#E69F00", "#56B4E9", "#D55E00", "#009E73", "#F0E442"]

# Figure and axes
fig, ax = subplots(figsize=(6, 6))
ax.set_xlabel("Initial uninfected bacteria (cells)")
ax.set_ylabel("Probability to survive")
ax.set_xlim(0, 15)
ax.set_ylim(0, 1)
ax.tick_params(labelsize=16)

# Ensure axes/frame render behind data (use Python dict values())
for sp in ax[:spines][:values]()
    sp[:set_zorder](-1)
end
ax.set_axisbelow(true)

# Helper to plot both scatter and line with explicit zorder
function add_series(ax, x, y; color, label, fill=true)
    if fill
        ax.scatter(x, y; s=36, c=color, edgecolors=color, linewidths=0.8, zorder=3, label=label)
    else
        ax.scatter(x, y; s=36, facecolors="white", edgecolors=color, linewidths=1.5, zorder=3, label=label)
    end
    ax.plot(x, y; color=color, zorder=2, label="_nolegend_")
end

add_series(ax, BTTT, PTTT; color=colors[1], label="LO, LIN, LORO")
add_series(ax, BFTT, PFTT; color=colors[2], label="LO, No LIN, LORO")
add_series(ax, BTFT, PTFT; color=colors[3], label="No LO, LIN, LORO")
add_series(ax, BTTF, PTTF; color=colors[4], label="LO, LIN, No LORO")
add_series(ax, BFFF, PFFF; color=colors[5], label="No LO, No LIN, No LORO")
add_series(ax, BTFF, PTFF; color=colors[6], label="No LO, LIN, No LORO")
add_series(ax, BTTTH, PTTTH; color=colors[1], label="LO, LIN, LORO, half burst", fill=false)

ax.legend(loc="upper left", fontsize=8)

# Save
mkpath(figures_dir)
outfile = joinpath(figures_dir, "prob_survive_comparison_withhalfburst_pyplot.pdf")
savefig(outfile, bbox_inches="tight")
println("Saved: " * outfile)
