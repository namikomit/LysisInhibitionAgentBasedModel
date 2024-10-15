# This is a Julia file

# Add your code here
using Random
using Plots
using Distributions
using JLD2


#I will now try to create a struct
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
lysis_rate = 1.0/25.0  #per minute
growth_timer = 10 #max growth timer
lysis_timer = 200 #max lysis timer, 4 timer is 1 minute
eclipse = 15    #eclipse time in minutes
burst_size = 100 #burst size
burst_rate=burst_size/((1/lysis_rate)-eclipse)
eta = 1e-9  #adsorption rate per ml/min
lysis_inhibition=true
lysis_inhibition_timer=4*10*2
lysis_from_without=true
lysis_from_without_phage=50
lo_resistance=true
lo_resistance_timer=4*5
li_collapse=true
li_collapse_phage=100
time_step=0.01

#Now set the initial condition and run the simulation. 
record_time_step = 1 #minutes

culture_growth=false
#I will now make 2 versions of the simulation, one with MSOI and another is culture growth
if culture_growth
    volume = 0.001 #ml
    nutrient = Int(round(1.e9*volume)) #cells/ml, growth rate does not depends on it but growth stops if bacteria number reach nutrient
    bacteria = Int(round(1e8*volume)) #cells
    infected= 0
    final_time = 10*60 # minutes
    phage = Int(round(2e7*volume)) # pfu
    # Generate random values
    initial_values = rand(1:growth_timer, bacteria-infected)
    append!(initial_values, [growth_timer + 1 for _ in 1:infected])
    #make it into an array with default values
    states = [State(initial_values[i], 0, 0.0, false) for i in 1:bacteria]


    # Create directories if they do not exist
    data_dir = "data_files_struct_culture_Default"
    figures_dir = "figure_files_struct_culture_Default"
    mkpath(data_dir)
    mkpath(figures_dir)


    time, Btimeseries, Itimeseries, Ptimeseries, lysis_time_record, states, bacteria, phage = simulate_population_agents(
        states, time_step, record_time_step, final_time, 
        bacteria, phage, infected, volume, growth_rate, nutrient,
        lysis_rate, burst_rate, eclipse, growth_timer, lysis_timer, eta; 
        lysis_inhibition=lysis_inhibition, lysis_inhibition_timer=lysis_inhibition_timer, 
        lysis_from_without=lysis_from_without, lysis_from_without_phage=lysis_from_without_phage, 
        lo_resistance=lo_resistance, lo_resistance_timer=lo_resistance_timer, 
        li_collapse=li_collapse, li_collapse_phage=li_collapse_phage)
    


    data_file_path = joinpath(data_dir, "population_data_lysis_timer($lysis_timer).jld2")
    @save  data_file_path time Btimeseries Itimeseries Ptimeseries lysis_time_record states time_step record_time_step final_time volume growth_rate nutrient lysis_rate burst_rate eclipse growth_timer lysis_timer eta lysis_inhibition lysis_inhibition_timer lysis_from_without lysis_from_without_phage lo_resistance lo_resistance_timer li_collapse li_collapse_phage 


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

    figure_file_path = joinpath(figures_dir, "population_dynamics_plot_lysis_timer($lysis_timer).pdf")
    savefig(figure_file_path)

    plot(size=(800, 480))

    # Create the histogram
    histogram(lysis_time_record, bins=200, label="lysis time histogram", xlabel="lysis time", ylabel="Frequency", title="Histogram of lysis time")
    figure_file_path = joinpath(figures_dir, "lysis_time_histogram_lysis_timer($lysis_timer).pdf")
    savefig(figure_file_path)

else
    volume = 0.01 #ml
    nutrient = Int(round(1.e9*volume)) #cells/ml, growth rate does not depends on it but growth stops if bacteria number reach nutrient
    bacteria = Int(round(2e7*volume)) #cells
    infected=Int(round(1e7*volume))
    si_duration=3. #minutes
    msoi=0.4 
    #deltaP=P_0*(1-exp(-eta*(B/volume)*s)
    #The total number of phages adsorbed, if P_0 is per ml, then dP is also per ml
    #Then msoi=deltaP/(B/volume)=   P_0(1-exp(-eta*(B/volume)*si_duration))/(B/volume)
    eta0=5e-9
    P0 = msoi * ((Float64(bacteria) / volume)) / (1 - exp(-eta0 * Float64(bacteria) / volume * si_duration))
    si_time = 15. # minutes
    final_time = si_time # minutes
    #eta=eta*0.3


    # Generate random values
    initial_values = rand(1:growth_timer, bacteria-infected)
    append!(initial_values, [growth_timer + 1 for _ in 1:infected])
    #make it into an array with default values
    states = [State(initial_values[i], 0, 0.0, false) for i in 1:bacteria]


    # Create directories if they do not exist
    data_dir = "data_files_struct_growth_paper_eta1e-9"
    figures_dir = "figure_files_struct_growth_paper_eta1e-9"
    mkpath(data_dir)
    mkpath(figures_dir)

    phage =0

    time, Btimeseries, Itimeseries, Ptimeseries, lysis_time_record, states, bacteria, phage = simulate_population_agents(
        states, time_step, record_time_step, final_time, 
        bacteria, phage, infected, volume, growth_rate, nutrient,
        lysis_rate, burst_rate, eclipse, growth_timer, lysis_timer, eta; 
        lysis_inhibition=lysis_inhibition, lysis_inhibition_timer=lysis_inhibition_timer, 
        lysis_from_without=lysis_from_without, lysis_from_without_phage=lysis_from_without_phage, 
        lo_resistance=lo_resistance, lo_resistance_timer=lo_resistance_timer, 
        li_collapse=li_collapse, li_collapse_phage=li_collapse_phage)
    
    #   27.553249 seconds (668.70 k allocations: 30.180 GiB, 15.08% gc time, 0.71% compilation time)
    phage = Int(round(P0 * volume)) # pfu
    print(phage)
    final_time = si_duration     #minutes
    time, Btimeseries, Itimeseries, Ptimeseries, lysis_time_record, states, bacteria, phage = simulate_population_agents(
        states, time_step, record_time_step, final_time, 
        bacteria, phage, infected, volume, growth_rate, nutrient,
        lysis_rate, burst_rate, eclipse, growth_timer, lysis_timer, eta; 
        lysis_inhibition=lysis_inhibition, lysis_inhibition_timer=lysis_inhibition_timer, 
        lysis_from_without=lysis_from_without, lysis_from_without_phage=lysis_from_without_phage, 
        lo_resistance=lo_resistance, lo_resistance_timer=lo_resistance_timer, 
        li_collapse=li_collapse, li_collapse_phage=li_collapse_phage)
    


    phage = 0
    final_time = 40 # minutes
    #eta=0.0
    #println(length(Bstate), length(Pstate), length(Istate), length(LORstate), bacteria)

    time2, Btimeseries2, Itimeseries2, Ptimeseries2, lysis_time_record2, states, bacteria, phage = simulate_population_agents(
        states, time_step, record_time_step, final_time, 
        bacteria, phage, infected, volume, growth_rate, nutrient,
        lysis_rate, burst_rate, eclipse, growth_timer, lysis_timer, eta; 
        lysis_inhibition=lysis_inhibition, lysis_inhibition_timer=lysis_inhibition_timer, 
        lysis_from_without=lysis_from_without, lysis_from_without_phage=lysis_from_without_phage, 
        lo_resistance=lo_resistance, lo_resistance_timer=lo_resistance_timer, 
        li_collapse=li_collapse, li_collapse_phage=li_collapse_phage)
    

    append!(lysis_time_record2, lysis_time_record)
    # Save the data

    data_file_path = joinpath(data_dir, "population_data_lysis_timer($lysis_timer)_lysis_inhibition_timer($lysis_inhibition_timer)_MSOI$(msoi).jld2")
    @save  data_file_path time2 Btimeseries2 Itimeseries2 Ptimeseries2 lysis_time_record2 states time_step record_time_step final_time volume growth_rate nutrient lysis_rate burst_rate eclipse growth_timer lysis_timer eta lysis_inhibition lysis_inhibition_timer lysis_from_without lysis_from_without_phage lo_resistance lo_resistance_timer li_collapse li_collapse_phage
    #@save "population_data_lysis_timer($lysis_timer)_MSOI$(msoi).jld2" begin
    #    time2, Btimeseries2, Itimeseries2, Ptimeseries2, irecord2, 
    #    Bstate, Pstate, Istate, LORstate, time_step, record_time_step, 
    #    final_time, volume, growth_rate, lysis_rate, burst_rate, 
    #    eclipse, growth_timer, lysis_timer, eta, lysis_inhibition, 
    #    lysis_inhibition_timer, lysis_from_without, lysis_from_without_phage, 
    #    lo_resistance, lo_resistance_time, li_collapse, li_collapse_phage
    #end

    # Create the plot
    plot(time2 .+ 18, Btimeseries2, label="Bacteria", linewidth=2)
    plot!(time2 .+ 18, Itimeseries2, label="Infected Bacteria", linewidth=2)
    plot!(time2 .+ 18, Ptimeseries2, label="Phage", linewidth=2)

    # Set labels and title
    xlabel!("Time (minutes)")
    ylabel!("Population")
    title!("Population Dynamics")

    # Show legend
    #plot!(legend=:topright)

    figure_file_path = joinpath(figures_dir, "population_dynamics_plot_lysis_timer($lysis_timer)_lysis_inhibition_timer($lysis_inhibition_timer)_MSOI$(msoi).pdf")
    savefig(figure_file_path)



    # Calculate the difference in Phage levels
    P_diff = diff(Ptimeseries2)

    # Create a new figure for the difference plot
    plot(size=(800, 480))

    # Plot the difference in Phage levels
    plot!(time2[2:end] .+ 18, P_diff, label="Difference in Phage Level")

    # Label the axes
    xlabel!("Time (minutes)")
    ylabel!("Difference in Phage Level")

    figure_file_path = joinpath(figures_dir, "phage_difference_plot_lysis_timer($lysis_timer)_lysis_inhibition_timer($lysis_inhibition_timer)_MSOI$(msoi).pdf")
    savefig(figure_file_path)

    plot(size=(800, 480))

    # Create the histogram
    histogram(lysis_time_record2, bins=200, label="lysis time histogram", xlabel="lysis time", ylabel="Frequency", title="Histogram of lysis time")
    figure_file_path = joinpath(figures_dir, "lysis_time_histogram_lysis_timer($lysis_timer)_lysis_inhibition_timer($lysis_inhibition_timer)_MSOI$(msoi).pdf")
    savefig(figure_file_path)
end
