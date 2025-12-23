using CSV, DataFrames, PlotlyJS

figures_dir = "Culture_Figures_Paper"
#figures_dir = "Culture_Figures_Paper_smallLORtime"

# Load each timeseries from its CSV file
df_S1 = CSV.read(joinpath(figures_dir, "timeseries_S1.csv"), DataFrame)
df_I1 = CSV.read(joinpath(figures_dir, "timeseries_I1.csv"), DataFrame)
df_P1 = CSV.read(joinpath(figures_dir, "timeseries_P1.csv"), DataFrame)
df_S2 = CSV.read(joinpath(figures_dir, "timeseries_S2.csv"), DataFrame)
df_I2 = CSV.read(joinpath(figures_dir, "timeseries_I2.csv"), DataFrame)
df_P2 = CSV.read(joinpath(figures_dir, "timeseries_P2.csv"), DataFrame)
df_S3 = CSV.read(joinpath(figures_dir, "timeseries_S3.csv"), DataFrame)
df_I3 = CSV.read(joinpath(figures_dir, "timeseries_I3.csv"), DataFrame)
df_P3 = CSV.read(joinpath(figures_dir, "timeseries_P3.csv"), DataFrame)

colors = ["#0072B2", "#E69F00", "#009E73"]

# Recreate traces
tracesB1 = scatter(x = df_S1.time, y = df_S1.Stimeseries, mode = "lines", line = attr(width = 2, dash = "dash", color = colors[1]), name="Sensitive, Mₗₒᵣ=∞")
tracesI1 = scatter(x = df_I1.time, y = df_I1.Itimeseries, mode = "lines", line = attr(width = 2, dash = "dash", color = colors[2]), name="Infected, Mₗₒᵣ=∞")
tracesP1 = scatter(x = df_P1.time, y = df_P1.Ptimeseries, mode = "lines", line = attr(width = 2, dash = "dash", color = colors[3]), name="Phage, Mₗₒᵣ=∞")
tracesB2 = scatter(x = df_S2.time, y = df_S2.Stimeseries, mode = "lines", line = attr(width = 2, dash = "solid", color = colors[1]), name="Sensitive, Mₗₒᵣ=150")
tracesI2 = scatter(x = df_I2.time, y = df_I2.Itimeseries, mode = "lines", line = attr(width = 2, dash = "solid", color = colors[2]), name="Infected, Mₗₒᵣ=150")
tracesP2 = scatter(x = df_P2.time, y = df_P2.Ptimeseries, mode = "lines", line = attr(width = 2, dash = "solid", color = colors[3]), name="Phage, Mₗₒᵣ=150")
tracesB3 = scatter(x = df_S3.time, y = df_S3.Stimeseries, mode = "lines", line = attr(width = 2, dash = "dot", color = colors[1]), name="Sensitive, Phage removal")
tracesI3 = scatter(x = df_I3.time, y = df_I3.Itimeseries, mode = "lines", line = attr(width = 2, dash = "dot", color = colors[2]), name="Infected, Phage removal")
tracesP3 = scatter(x = df_P3.time, y = df_P3.Ptimeseries, mode = "lines", line = attr(width = 2, dash = "dot", color = colors[3]), name="Phage, Phage removal")

# ...existing code...

# --- Plot (a): log y-axis ---
layout_log = Layout(
    title = attr(
        text = "(a)",
        font = attr(size = 24),
        x = 0.01,
        y = 0.9,
        xanchor = "left",
        yanchor = "top"
    ),
    xaxis = attr(
        title = attr(text = "Time (minutes)", font = attr(size = 24)),
        linecolor = "black",
        linewidth = 2,
        ticks = "inside",
        showline = true,
        mirror = true,
        range = [0, 500],
        tickfont = attr(size = 20)
    ),
    yaxis = attr(
        title = attr(text = "Bacteria /ml or Phage /ml", font = attr(size = 24)),
        linecolor = "black",
        linewidth = 2,
        ticks = "inside",
        type = "log",
        showline = true,
        mirror = true,
        tickfont = attr(size = 20),
        range = [log10(1e4), log10(1e11)],
        tickvals = [1e4, 1e6, 1e8, 1e10],
        ticktext =  ["10⁴", "10⁶", "10⁸", "10¹⁰"] #["10$(Char(0x2070 + 4))", "10$(Char(0x2070 + 6))", "10$(Char(0x2070 + 8))","10$(Char(0x00B9)Char(0x2070))"]
    ),
        legend = attr(
        x = 0.5,           
        y = 0.25,         
        xanchor = "left",
        yanchor = "middle",
        bgcolor = "rgba(255,255,255,0.7)",  # optional: semi-transparent background
        bordercolor = "black",
        borderwidth = 1,
        font = attr(size = 12)
    ),
    showlegend = true,
    plot_bgcolor = "rgba(0,0,0,0)",
    paper_bgcolor = "rgba(0,0,0,0)"
)

plot_log = plot([tracesB1, tracesI1, tracesP1, tracesB2, tracesI2, tracesP2, tracesB3, tracesI3, tracesP3], layout_log)

# --- Plot (b): linear y-axis ---
layout_linear = Layout(
    title = attr(
        text = "(b)",
        font = attr(size = 24),
        x = 0.01,
        y = 0.9,
        xanchor = "left",
        yanchor = "top"
    ),
    xaxis = attr(
        title = attr(text = "Time (minutes)", font = attr(size = 24)),
        linecolor = "black",
        linewidth = 2,
        ticks = "inside",
        showline = true,
        mirror = true,
        range = [0, 500],
        tickfont = attr(size = 20)
    ),
    yaxis = attr(
        title = attr(text = "Bacteria /ml or Phage /ml", font = attr(size = 24)),
        linecolor = "black",
        linewidth = 2,
        ticks = "inside",
        type = "linear",
        showline = true,
        mirror = true,
        tickfont = attr(size = 20),
        range = [3e7, 1.7e8],
        tickvals = [5e7, 1e8, 1.5e8, 2e8],
        ticktext = ["0.5×10⁸", "1.0×10⁸", "1.5×10⁸", "2.0×10⁸"]  
    ),
    showlegend = false,
    plot_bgcolor = "rgba(0,0,0,0)",
    paper_bgcolor = "rgba(0,0,0,0)"
)

plot_linear = plot([tracesB1, tracesI1, tracesP1, tracesB2, tracesI2, tracesP2, tracesB3, tracesI3, tracesP3], layout_linear)


# Save log-scale plot
savefig(plot_log, joinpath(figures_dir, "Population_timeseries_log.pdf"))

# Save linear-scale plot
savefig(plot_linear, joinpath(figures_dir, "Population_timeseries_linear.pdf"))