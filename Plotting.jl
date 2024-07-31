using JLD2
using Plots

#I will now try to create a struct
mutable struct State
    Bstate::UInt16 # 1 to growth_timer is uninfected, growth_timer+1 to growth_timer+lysis_timer is infected
    # assign a random integer number between 1 and growth_timer to each bacteria
    Istate::UInt16  # Count the number of infection in total
    Pstate::Float64 # recording time spent in infected state, to compute number of produced phages for an infected bacteria, proportional to time spent in infected state (minus eclipse)
    LORstate::Bool # boolean, True if it is in Lysis from without resistant state
end

# Create directories if they do not exist
data_dir = "data_files_struct_culture"
figures_dir = "figure_files_struct_culture"
mkpath(data_dir)
mkpath(figures_dir)

# Load the data from the JLD2 file
data_file_path = joinpath(data_dir, "population_data_lysis_timer(100).jld2")
@load data_file_path time Btimeseries Itimeseries Ptimeseries lysis_time_record states time_step record_time_step final_time volume growth_rate nutrient lysis_rate burst_rate eclipse growth_timer lysis_timer eta lysis_inhibition lysis_inhibition_timer lysis_from_without lysis_from_without_phage lo_resistance lo_resistance_timer li_collapse li_collapse_phage


# Create the plot
# Filter out zero or negative values
function filter_positive(time, series)
    positive_indices = findall(x -> x > 0, series)
    return time[positive_indices], series[positive_indices]
end

time_B, Btimeseries_filtered = filter_positive(time, Btimeseries/volume)
time_I, Itimeseries_filtered = filter_positive(time, Itimeseries/volume)
time_P, Ptimeseries_filtered = filter_positive(time, Ptimeseries/volume)

# Create the plot with log scale on y-axis
plot(time_B, Btimeseries_filtered, label="Bacteria", linewidth=2) #, yscale=:log10)
plot!(time_I, Itimeseries_filtered, label="Infected Bacteria", linewidth=2)
#plot!(time_P, Ptimeseries_filtered, label="Phage", linewidth=2)


# Set labels and title
xlabel!("Time (minutes)")
ylabel!("Population")
title!("Population Dynamics")

# Show legend
#plot!(legend=:topright)

figure_file_path = joinpath(figures_dir, "population_dynamics_plot_lysis_timer(100).pdf")
savefig(figure_file_path)
