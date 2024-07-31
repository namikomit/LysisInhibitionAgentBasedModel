using JLD2
using Plots

# Load the data from the JLD2 file
@load "population_data.jld2" time2 Btimeseries2 Itimeseries2 Ptimeseries2 irecord2
# Create the plot
plot(time2[1:irecord2], Btimeseries2[1:irecord2], label="Bacteria", linewidth=2)
plot!(time2[1:irecord2], Itimeseries2[1:irecord2], label="Infected Bacteria", linewidth=2)
plot!(time2[1:irecord2], Ptimeseries2[1:irecord2], label="Phage", linewidth=2)

# Set labels and title
xlabel!("Time (minutes)")
ylabel!("Population")
title!("Population Dynamics")

# Show legend
#plot!(legend=:topright)


savefig("population_dynamics_plot.png")



# Calculate the difference in Phage levels
P_diff = diff(Ptimeseries2[1:irecord2])

# Create a new figure for the difference plot
plot(size=(800, 480))

# Plot the difference in Phage levels
plot!(time2[2:irecord2] .+ 15, P_diff, label="Difference in Phage Level")

# Label the axes
xlabel!("Time (minutes)")
ylabel!("Difference in Phage Level")

savefig("phage_difference_plot.png")
