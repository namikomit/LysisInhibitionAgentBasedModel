# This is a Julia file

# Add your code here
using Random
using Plots

# Main function
function simulate_population_agents(Bstate, Pstate, Istate, LORstate, time_step, record_time_step, final_time, bacteria, phage, infected, bacteriaever, carrying_capacity, 
    growth_rate, lysis_rate, burst_rate, eclipse, growth_timer, lysis_timer, eta; lysis_inhibition=false, lysis_inhibition_timer=5, 
    lysis_from_without=false, lysis_from_without_phage=10, lo_resistance=false, lo_resistance_time=5, li_collapse=false, li_collapse_phage=100)
    # Add your code here
    
    println("started")
    """
    Simulates a population dynamics with individual bacteria having growth timer and lysis timer. 
    Simulate V ml. The number of phages (integer)
    ...
    """
    volume = carrying_capacity / 1e9  # volume in ml

    timenow = 0.0
    nrecord = Int(final_time / record_time_step)
    time = zeros(nrecord)
    Btimeseries = zeros(nrecord)
    Itimeseries = zeros(nrecord)
    Ptimeseries = zeros(nrecord)
    irecord = 0
    time[irecord + 1] = timenow  # Julia arrays are 1-indexed
    Btimeseries[irecord + 1] = bacteria
    Itimeseries[irecord + 1] = infected
    Ptimeseries[irecord + 1] = phage
    grate = growth_rate * growth_timer
    lrate = lysis_rate * lysis_timer
    println(volume)
    while timenow < final_time
        timenow += time_step
        # Simulation logic here
        if timenow >= irecord * record_time_step
            irecord += 1
            time[irecord + 1] = timenow
            Btimeseries[irecord + 1] = bacteria
            Itimeseries[irecord + 1] = infected
            Ptimeseries[irecord + 1] = phage
        end
    end

    for j in 1:length(Bstate)
        if !LORstate[j] && lysis_from_without
            if Istate[j] > lysis_from_without_phage
                Bstate[j] = 0
                infected -= 1
                bacteria -= 1
                phagenew += max(Int(burst_rate * (Pstate[j] - eclipse)), 0)
            end
        end
        if LORstate[j] || !lysis_from_without
            if lysis_inhibition
                Bstate[j] = max(growth_timer + 1, Bstate[j] - lysis_inhibition_timer * phageinfect)
            end
            if li_collapse && Istate[j] > li_collapse_phage
                Bstate[j] = 0
                infected -= 1
                bacteria -= 1
                phagenew += max(Int(burst_rate * (Pstate[j] - eclipse)), 0)
            end
        end
        if Bstate[j] < growth_timer
            if rand() < growth_rate * time_step
                Bstate[j] += 1
            end
        elseif Bstate[j] == growth_timer
            if rand() < growth_rate * time_step
                Bstate[j] = 1
                bacteria += 1
                bacteriaever += 1
                if bacteriaever < carrying_capacity
                    Bstate[bacteriaever] = 1
                else
                    println("Carrying capacity reached")
                    timenow = final_time
                    break
                end
            end
        end
    end
    return time, Btimeseries, Itimeseries, Ptimeseries
end


# Define initial states and parameters
Bstate = [1, 2, 3]  # Example initial states for bacteria
Pstate = [0, 0, 0]  # Example initial states for phages
Istate = [0, 0, 0]  # Example initial states for infected bacteria
LORstate = [false, false, false]  # Example initial states for lysis or resistance

time_step = 0.1
record_time_step = 1.0
final_time = 10.0
bacteria = 3
phage = 0
infected = 0
bacteriaever = 3
carrying_capacity = 1000

growth_rate = 0.1
lysis_rate = 0.05
burst_rate = 50
eclipse = 1
growth_timer = 3
lysis_timer = 5
eta = 0.1

# Optional parameters
lysis_inhibition = false
lysis_inhibition_timer = 5
lysis_from_without = false
lysis_from_without_phage = 10
lo_resistance = false
lo_resistance_time = 5
li_collapse = false
li_collapse_phage = 100

# Call the main function
time, Btimeseries, Itimeseries, Ptimeseries = simulate_population_agents(
    Bstate, Pstate, Istate, LORstate, time_step, record_time_step, final_time, 
    bacteria, phage, infected, bacteriaever, carrying_capacity, growth_rate, 
    lysis_rate, burst_rate, eclipse, growth_timer, lysis_timer, eta; 
    lysis_inhibition=lysis_inhibition, lysis_inhibition_timer=lysis_inhibition_timer, 
    lysis_from_without=lysis_from_without, lysis_from_without_phage=lysis_from_without_phage, 
    lo_resistance=lo_resistance, lo_resistance_time=lo_resistance_time, 
    li_collapse=li_collapse, li_collapse_phage=li_collapse_phage
)

# Plot the results with labels and title
plot(time, Btimeseries, label="Bacteria", xlabel="Time", ylabel="Population", title="Population Dynamics")
plot!(time, Itimeseries, label="Infected")
plot!(time, Ptimeseries, label="Phage")

savefig("population_dynamics_plot.png")


# Display the plot
#display(plot(time, Btimeseries, label="Bacteria", xlabel="Time", ylabel="Population", title="Population Dynamics"))
#display(plot!(time, Itimeseries, label="Infected"))
#display(plot!(time, Ptimeseries, label="Phage"))
