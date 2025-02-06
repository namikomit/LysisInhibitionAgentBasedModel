# This is a Julia file

# Add your code here
using Random
#using Plots
using Distributions
using JLD2
using CSV
using DataFrames
using PlotlyJS



#Now try to create a struct
mutable struct State
    Bstate::UInt16 # 1 to growth_timer is uninfected, growth_timer+1 to growth_timer+lysis_timer is infected
    # assign a random integer number between 1 and growth_timer to each bacteria
    Istate::UInt16  # Count the number of infection in total
    Pstate::Float64 # recording time spent in infected state, to compute number of produced phages for an infected bacteria, proportional to time spent in infected state (minus eclipse)
    LORstate::Bool # boolean, True if it is in Lysis from without resistant state
end

# Main function
function simulate_population_agents(states::Vector{State}, time_step, record_time_step, final_time, bacteria, phage, infected, volume, 
    growth_rate, nutrient, lysis_rate, burst_rate, eclipse, growth_timer, lysis_timer, eta; lysis_inhibition=false, lysis_inhibition_timer=5, 
    lysis_from_without=false, lysis_from_without_phage=10, lo_resistance=false, lo_resistance_timer=5, li_collapse=false, li_collapse_phage=100)
    # Add your code here
    
    println("started")
    """
    Simulates a population dynamics with individulal bacteria has growth timer and lysis timer. 
    Simulate V ml. The number of phages (integer)

    Parameters:
    - time_step: The time step for evaluating the dynamics.
    -record_time_step: The time step for recording the dynamics.
    - final_time: The final time to simulate.
    - bacteria: The initial number of bacteria.
    - infected: The initial number of infected bacteria.
    - phage: The initial number of phage.
    - carrying_capacity: The carrying capacity of the environment. 
    system volume is calulcated from this assuming 1ml has carrying capacity of 10^9 cells. 
    - growth_rate: The growth rate of the bacteria.
    - lysis_rate: 1/The latency time of the phage production.
    - burst_rate The burst rate of the phage.
    - ecliplse: The ecliplse time of the phage.
    - growth_timer: Number of the growth timer of the bacteria.
    - lysis_timer: Number of the lysis timer of upon infection.
    - eta: The adsorption rate of the phage.
    - lysis_inhibition: Whether the lysis inhibition is present. Default is False.
    
    Returns:
    - time: A list of times.
    - Btimeseries: A list of the number of bacteria at each record time step.
    - Itimeseries: A list of the number of infected bacteria at each record time step.
    - Ptimeseries: A list of the number of phage at each record time step.
    ...
    """

    timenow = 0.
    time = Float64[]
    Btimeseries = Int[]
    Itimeseries = Int[]
    Ptimeseries = Int[]
    push!(time,timenow)  
    push!(Btimeseries,bacteria)
    push!(Itimeseries,infected)
    push!(Ptimeseries,phage)
    grate = growth_rate * growth_timer
    lrate = lysis_rate * lysis_timer
    lysis_time_record = Float64[]
    println(volume)
    itime=Int64(0)
    while timenow < final_time
        lo_new_phage=0 
        lysis_new_phage=0
        new_bacteria=0
        phageinfect_array= rand(Poisson(eta*phage/volume*time_step), bacteria)
        random_numbers = rand(bacteria)
        for i in 1:bacteria
            #println("phage infection", phageinfect_array[i])
            if phageinfect_array[i] > 0
                #Only execute the infection actions when it actually hallens!
                states[i].Istate += phageinfect_array[i]
                if states[i].Bstate <= growth_timer
                #first time infection
                    states[i].Bstate = growth_timer+1   
                end 
                if lysis_from_without
                    
                    #The lysis from without can happen
                    check_LO = !states[i].LORstate * states[i].Istate  #If true for LOR state, Istate will be multiplied by 0. 
                    if check_LO >= lysis_from_without_phage
                        #println("lysis from without")
                        states[i].Bstate = 0
                        lo_new_phage =0 #+= max(Int(round(burst_rate * (states[i].Pstate - eclipse))), 0)
                        # Append the values of Pstate that match mask_LO to lysis_time_record
                        append!(lysis_time_record, states[i].Pstate)                     
                    end
                end
                if li_collapse
                    if states[i].Istate > li_collapse_phage
                        #println("lysis inhibition collapse: ", states[i].Istate, " threshold: ",li_collapse_phage)
                        states[i].Bstate = 0
                        lo_new_phage += max(Int(round(burst_rate * (states[i].Pstate - eclipse))), 0)
                        # Append the values of Pstate that match mask_LI to lysis_time_record
                        append!(lysis_time_record, states[i].Pstate)
                    end
                end
                if lysis_inhibition
                    #println("lysis inhibition")
                    check_LI = true
                    if  lo_resistance & !states[i].LORstate
                        #if LOR can happen, then the lysis inhibition should start only after LOR is established
                        check_LI = false
                    end
                    #The lysis inhibition can happen, but only the ones that did not do lysis from without 
                    if states[i].Bstate == 0
                        check_LI = false
                    end
                    if check_LI
                        states[i].Bstate = max(growth_timer+1, states[i].Bstate - lysis_inhibition_timer*phageinfect_array[i])
                    end    
                end
            end
            # Update Bstate elements less than growth_timer with probability grate * time_step
            #println("growth")
            if states[i].Bstate < growth_timer
                states[i].Bstate += random_numbers[i] < (grate * time_step)
            elseif states[i].Bstate == growth_timer
                if random_numbers[i] < (grate * time_step)
                    #println("new bacteria")    
                    new_bacteria += 1
                    states[i].Bstate = 1
                end
            # Update Bstate elements greater than growth_timer and less than growth_timer + lysis_timer with probability lrate * time_step
            #println("lysis")
            elseif states[i].Bstate < growth_timer + lysis_timer
                #lysis actions
                # Add time_step to all masked elements of Pstate
                states[i].Pstate += time_step
                if lo_resistance & !states[i].LORstate
                    states[i].LORstate = states[i].Bstate - growth_timer-1 > lo_resistance_timer
                end
                states[i].Bstate += random_numbers[i] < (lrate * time_step)
            elseif states[i].Bstate == growth_timer + lysis_timer
                if random_numbers[i] < (lrate * time_step)
                    states[i].Bstate = 0
                    lysis_new_phage += max(Int(round(burst_rate * (states[i].Pstate - eclipse))), 0)
                    # Append the values of Pstate that match mask_lysis to lysis_time_record
                    append!(lysis_time_record, states[i].Pstate)
                end 
            end 
        end


        #update of phage, infected, and bacteria should be done here at the end
        #println("update")
        phage += lo_new_phage + lysis_new_phage
        phage=max(0,phage-sum(phageinfect_array))
        # Filter out elements for no cell
        filter!(states -> states.Bstate != 0, states)
        if new_bacteria>0
        #Add the new bacteria, number new_bacteria
            for _ in 1:new_bacteria
                states = push!(states, State(1, 0, 0.0, false))
            end
        end
        bacteria=length(states)

        if bacteria > nutrient
            println("Nutrient exhausted")
            return (time, Btimeseries, Itimeseries, Ptimeseries, lysis_time_record)
        end
        #println(bacteria, " bacteria ", sum(Bstate .>=0), " sum ", length(Bstate))
        #println("Pstate ", length(Pstate), " Istate ", length(Istate), " LORstate ", length(LORstate))
        itime+=1
        timenow = time_step*itime
        if timenow - time[end] >= record_time_step
            push!(time,timenow)  
            push!(Btimeseries,bacteria)
            infected = sum([state.Bstate > growth_timer for state in states])
            push!(Itimeseries,infected)
            push!(Ptimeseries,phage)
            println(time[end], " ", Btimeseries[end], " ", Itimeseries[end], " ", Ptimeseries[end])
        end
        #println("bacteria", bacteria, timenow, final_time)
    end

    return time, Btimeseries, Itimeseries, Ptimeseries, lysis_time_record, states, bacteria, phage
end


#Here we define the system parameters.
#We start with the simulation done in the Julia's thesis of different MSOI
growth_rate = 2.0/60. #per minute
lysis_rate = 1.0/27.0  #per minute
growth_timer = 10 #max growth timer
eclipse = 15    #eclipse time in minutes
burst_size = 150 #burst size
burst_rate=burst_size/((1/lysis_rate)-eclipse)

lysis_timer = 250 #max lysis timer
lysis_inhibition=true
lysis_inhibition_timer=Int(round(5*(lysis_timer*lysis_rate)))
lysis_from_without=true
lysis_from_without_phage=50
lo_resistance=true
lo_resistance_timer=Int(round(10*(lysis_timer*lysis_rate)))
li_collapse=true
li_collapse_phage=100
time_step=0.01
eta=2.e-9

#Now set the initial condition and run the simulation. 
record_time_step = 1 #minutes
P_diff_matrix = []
time_save_list = []
P_diff_error   = []

lysis_timer_flag=true
# Create directories if they do not exist
figures_dir = "Culture_Figures_Paper"
mkpath(figures_dir)

volume = 0.001 #ml
nutrient = Int(round(1.e9*volume)) #cells/ml, growth rate does not depends on it but growth stops if bacteria number reach nutrient
bacteria = Int(round(1e8*volume)) #cells
infected= 0
final_time = 6*60 # minutes
phage = Int(round(2e7*volume)) # pfu
# Generate random values
initial_values = rand(1:growth_timer, bacteria-infected)
append!(initial_values, [growth_timer + 1 for _ in 1:infected])
#make it into an array with default values
states = [State(initial_values[i], 0, 0.0, false) for i in 1:bacteria]


# Create directories if they do not exist
figures_dir = "Culture_Figures_Paper"
mkpath(figures_dir)
#Simulation without no interference

time, Btimeseries, Itimeseries, Ptimeseries, lysis_time_record, states, bacteria, phage = simulate_population_agents(
    states, time_step, record_time_step, final_time, 
    bacteria, phage, infected, volume, growth_rate, nutrient,
    lysis_rate, burst_rate, eclipse, growth_timer, lysis_timer, eta; 
    lysis_inhibition=lysis_inhibition, lysis_inhibition_timer=lysis_inhibition_timer, 
    lysis_from_without=lysis_from_without, lysis_from_without_phage=lysis_from_without_phage, 
    lo_resistance=lo_resistance, lo_resistance_timer=lo_resistance_timer, 
    li_collapse=li_collapse, li_collapse_phage=li_collapse_phage)


# Create the plot
# Filter out zero or negative values
function filter_positive(time, series)
    positive_indices = findall(x -> x > 0, series)
    return time[positive_indices], series[positive_indices]
end

time_B, Btimeseries_filtered = filter_positive(time, Btimeseries/volume)
time_I, Itimeseries_filtered = filter_positive(time, Itimeseries/volume)
time_P, Ptimeseries_filtered = filter_positive(time, Ptimeseries/volume)

# Create the timeseries plot using PlotlyJS
trace_B = scatter(x = time_B, y = Btimeseries_filtered, mode = "lines", name = "Bacteria", line = attr(width = 2))
trace_I = scatter(x = time_I, y = Itimeseries_filtered, mode = "lines", name = "Infected Bacteria", line = attr(width = 2))
trace_P = scatter(x = time_P, y = Ptimeseries_filtered, mode = "lines", name = "Phage", line = attr(width = 2))

layout_timeseries = Layout(
    xaxis = attr(title = "Time (minutes)"),
    yaxis = attr(title = "Bacteria or Phage (/ml)", type = "log",tickformat = ".0e" )  # Set y-axis to logarithmic scale
)

plot_timeseries = Plot([trace_B, trace_I,trace_P], layout_timeseries)

# Save the timeseries plot
savefig(plot_timeseries, joinpath(figures_dir, "Population_LIC($li_collapse)_LICT($li_collapse_phage).pdf"))

#layout_Ptimeseries = Layout(
#    xaxis = attr(title = "Time (minutes)"),  
#    yaxis = attr(title = "Phage (/ml)", type = "log")  # Set y-axis to logarithmic scale
#)

#plot_timeseries = Plot([trace_P], layout_Ptimeseries)

# Save the timeseries plot
#savefig(plot_timeseries, joinpath(figures_dir, "population_dynamics_LIC($li_collapse)_LICT($li_collapse_phage).pdf"))



# Create the histogram using PlotlyJS
trace_histogram = histogram(x = lysis_time_record, nbinsx = 200, name = "Lysis Time Histogram")

layout_histogram = Layout(
    xaxis = attr(title = "Lysis Time"),
    yaxis = attr(title = "Frequency"),
    title = "Histogram of Lysis Time"
)

plot_histogram = Plot([trace_histogram], layout_histogram)

# Save the histogram plot
savefig(plot_histogram, joinpath(figures_dir, "lysis_time_histogram_LIC($li_collapse)_LICT($li_collapse_phage).pdf"))


data_file_path = joinpath(figures_dir, "Btimeseries_filtered_LIC($li_collapse)_LICT($li_collapse_phage).csv")
CSV.write(data_file_path, DataFrame(time = time_B, Btimeseries = Btimeseries_filtered))

data_file_path = joinpath(figures_dir, "Itimeseries_filtered_LIC($li_collapse)_LICT($li_collapse_phage).csv")
CSV.write(data_file_path, DataFrame(time = time_I, Itimeseries = Itimeseries_filtered))

data_file_path = joinpath(figures_dir, "Ptimeseries_filtered_LIC($li_collapse)_LICT($li_collapse_phage).csv")
CSV.write(data_file_path, DataFrame(time = time_P, Ptimeseries = Ptimeseries_filtered))