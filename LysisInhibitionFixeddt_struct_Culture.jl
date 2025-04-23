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
    growth_rate, nutrient, lysis_rate, burst_rate, beta_max, eclipse, growth_timer, lysis_timer, eta; lysis_inhibition=false, lysis_inhibition_timer=5, 
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
                        lo_new_phage += max(min(Int(round(burst_rate * (states[i].Pstate - eclipse))), beta_max), 0)
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
            if states[i].Bstate > 0
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
                        lysis_new_phage += max(min(Int(round(burst_rate * (states[i].Pstate - eclipse))), beta_max), 0)
                    # Append the values of Pstate that match mask_lysis to lysis_time_record
                        append!(lysis_time_record, states[i].Pstate)
                    end 
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

        if(timenow>400)
            for j in 1:length(states)
                if states[j].Bstate <= growth_timer
                    println("Index: $j, Bstate: $(states[j].Bstate)")
                end
            end
        end

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
            push!(Btimeseries, bacteria)
            infected = sum([state.Bstate > growth_timer for state in states])
            push!(Itimeseries,infected)
            push!(Ptimeseries,phage)
            println(time[end], " ", Btimeseries[end], " ", Itimeseries[end], " ", Ptimeseries[end],"BS",Btimeseries[end]-Itimeseries[end], " BI", bacteria-infected )
        end
    end
    
    return time, Btimeseries, Itimeseries, Ptimeseries, lysis_time_record, states, bacteria, phage
end


#Here we define the system parameters.
#We start with the simulation done in the Julia's thesis of different MSOI
global growth_rate = 2.0/60. #per minute
lysis_rate = 1.0/27.0  #per minute
growth_timer = 10 #max growth timer
eclipse = 15    #eclipse time in minutes
global burst_size = 150 #burst size
global burst_rate=burst_size/((1/lysis_rate)-eclipse)
global beta_max=500

lysis_timer = 250 #max lysis timer
lysis_inhibition=true
lysis_inhibition_timer=Int(round(5*(lysis_timer*lysis_rate)))
lysis_from_without=true
lysis_from_without_phage=100
lo_resistance=true
lo_resistance_timer=Int(round(10*(lysis_timer*lysis_rate)))
li_collapse=true
global li_collapse_phage = 150
time_step=0.1
eta=2.e-9

#Now set the initial condition and run the simulation. 
record_time_step = 0.01 #minutes
lysis_timer_flag=true
# Create directories if they do not exist
figures_dir = "Culture_Figures_Paper"
mkpath(figures_dir)

volume = 0.001 #ml



#Simulation of a culture
collapse_threshold_test=false
if(collapse_threshold_test)
    collapse_threshold = [100, 125, 150, 175, 200]
    collapse_time = Float64[0, 0, 0, 0, 0]
    collapse_definition_population = 1.
    initialB=5e7
    initialP=1e7
    for i in range(1, stop=length(collapse_threshold))
        global li_collapse_phage = collapse_threshold[i]
        nutrient = Int(round(1.e9*volume)) #cells/ml, growth rate does not depends on it but growth stops if bacteria number reach nutrient
        bacteria = Int(round(initialB*volume)) #cells
        infected= 0
        phage = Int(round(initialP*volume)) # pfu
# Generate random values
        initial_values = rand(1:growth_timer, bacteria-infected)
        append!(initial_values, [growth_timer + 1 for _ in 1:infected])
#make it into an array with default values
        states = [State(initial_values[i], 0, 0.0, false) for i in 1:bacteria]

        final_time = 10*60 # minutes

time, Btimeseries, Itimeseries, Ptimeseries, lysis_time_record, states, bacteria, phage = simulate_population_agents(
    states, time_step, record_time_step, final_time, 
    bacteria, phage, infected, volume, growth_rate, nutrient,
    lysis_rate, burst_rate, beta_max, eclipse, growth_timer, lysis_timer, eta; 
    lysis_inhibition=lysis_inhibition, lysis_inhibition_timer=lysis_inhibition_timer, 
    lysis_from_without=lysis_from_without, lysis_from_without_phage=lysis_from_without_phage, 
    lo_resistance=lo_resistance, lo_resistance_timer=lo_resistance_timer, 
    li_collapse=li_collapse, li_collapse_phage=li_collapse_phage)

    index_below_threshold = findfirst(x -> x < collapse_definition_population, Btimeseries)  
    println("index_below_threshold: ", index_below_threshold)
    println(length(Btimeseries),length(time))
    if index_below_threshold != nothing
        collapse_time[i] = time[index_below_threshold]
    else
        collapse_time[i] = final_time
    end
end
# Create the timeseries plot using PlotlyJS
trace = scatter(x = collapse_threshold, y = collapse_time, mode = "markers", line = attr(size=10, color="blue"))

layout_collapsetime = Layout(
    title = attr(
        text = "(b)",
        font = attr(size = 20),  # Font size for the title
        x = 0.01,  # X position of the title (0 to 1)
        y = 0.9,  # Y position of the title (0 to 1)
        xanchor = "left",  # Anchor the title at the center
        yanchor = "top"  # Anchor the title at the top
    ),
    xaxis = attr(
        title = attr(
            text = "LI collapse threshold",
            font = attr(size = 20)  # Font size for the x-axis label
        ),
        range = [0, maximum(collapse_threshold)]*1.1,
        linecolor = "black",
        linewidth = 2,
        ticks = "inside",  # Add ticks inside the plot
        showline = true,  # Show line on the bottom x-axis
        mirror = true,  # Mirror the axis lines on the top and right
        tickfont = attr(size = 20)  # Font size for the x-axis ticks
    ),
    yaxis = attr(
        title = attr(
            text = "Collapse time (min)",
            font = attr(size = 20)  # Font size for the x-axis label
        ), 
        range = [0, maximum(collapse_time)]*1.1,
        linecolor = "black",
        linewidth = 2,
        ticks = "inside",  # Add ticks inside the plot
        showline = true,  # Show line on the left y-axis
        mirror = true,  # Mirror the axis lines on the top and right
        tickfont = attr(size = 20)  # Font size for the x-axis ticks
    ),
    plot_bgcolor = "rgba(0,0,0,0)",  # Transparent plot background
    paper_bgcolor = "rgba(0,0,0,0)",  # Transparent paper background
    legend = "false"
)

plot_timeseries = Plot([trace], layout_collapsetime)

# Save the timeseries plot
savefig(plot_timeseries, joinpath(figures_dir, "LIC_thresholdDependence_initialB_($initialB)_initialP_($initialP).pdf"))

data_file_path = joinpath(figures_dir, "LIC_thresholdDependence_initialB_($initialB)_initialP_($initialP).csv")
CSV.write(data_file_path, DataFrame(threshold = collapse_threshold, collapse_time = collapse_time))
else
    #andiphage timing 
    all_timesB = []
    all_timesI = []
    all_timesP = []
    all_timesS = []
    all_Btimeseries = []
    all_Itimeseries = []
    all_Ptimeseries = []
    all_Stimeseries = []

    
    global li_collapse_phage = 150
    antiphage_timing_list = [180.,  360.]
   
    for i in range(2, stop=3)
        global burst_size = 150 #burst size
        global burst_rate=burst_size/((1/lysis_rate)-eclipse)
        initialB=5e7
        initialP=1e7
        final_time = 10*60 # minutes
        nutrient = Int(round(1.e9*volume)) #cells/ml, growth rate does not depends on it but growth stops if bacteria number reach nutrient
        bacteria = Int(round(initialB*volume)) #cells
        infected= 0
        phage = Int(round(initialP*volume)) # pfu
# Generate random values
        initial_values = rand(1:growth_timer, bacteria)
#make it into an array with default values
        states = [State(initial_values[i], 0, 0.0, false) for i in 1:bacteria]
        if(i==1)
            time_antiphage=final_time
            global li_collapse_phage = 100000
        elseif(i==2)
            time_antiphage=final_time
            global li_collapse_phage = 150
        else    
            time_antiphage=180
            global li_collapse_phage = 150
        end
   

    time, Btimeseries, Itimeseries, Ptimeseries, lysis_time_record, states, bacteria, phage = simulate_population_agents(
        states, time_step, record_time_step, time_antiphage, 
        bacteria, phage, infected, volume, growth_rate, nutrient,
        lysis_rate, burst_rate, beta_max, eclipse, growth_timer, lysis_timer, eta; 
        lysis_inhibition=lysis_inhibition, lysis_inhibition_timer=lysis_inhibition_timer, 
        lysis_from_without=lysis_from_without, lysis_from_without_phage=lysis_from_without_phage, 
        lo_resistance=lo_resistance, lo_resistance_timer=lo_resistance_timer, 
        li_collapse=li_collapse, li_collapse_phage=li_collapse_phage)

    Stimeseries=Btimeseries-Itimeseries
    #println("Bacteria: ", Btimeseries[end], " Infected: ", Itimeseries[end], " Sensitive: ", Stimeseries[end])
    #sleep(10)
final_time = final_time-time_antiphage 
if(final_time>0)
    phage = 0
burst_rate=0 #no phage production

    time2, Btimeseries2, Itimeseries2, Ptimeseries2, lysis_time_record2, states, bacteria, phage = simulate_population_agents(
        states, time_step, record_time_step, final_time, 
        bacteria, phage, infected, volume, growth_rate, nutrient,
        lysis_rate, burst_rate, beta_max, eclipse, growth_timer, lysis_timer, eta; 
        lysis_inhibition=lysis_inhibition, lysis_inhibition_timer=lysis_inhibition_timer, 
        lysis_from_without=lysis_from_without, lysis_from_without_phage=lysis_from_without_phage, 
        lo_resistance=lo_resistance, lo_resistance_timer=lo_resistance_timer, 
        li_collapse=li_collapse, li_collapse_phage=li_collapse_phage)
        Stimeseries2=Btimeseries2-Itimeseries2

    new_time = vcat(time, time2 .+ time_antiphage)
    new_Btimeseries = vcat(Btimeseries, Btimeseries2)
    new_Itimeseries = vcat(Itimeseries, Itimeseries2)
    new_Ptimeseries = vcat(Ptimeseries, Ptimeseries2)
    new_Stimeseries = vcat(Stimeseries, Stimeseries2)
else
    new_time = time
    new_Btimeseries = Btimeseries
    new_Itimeseries = Itimeseries
    new_Ptimeseries = Ptimeseries
    new_Stimeseries = Stimeseries
end

# Create the plot
#println("Bacteria: ", new_Btimeseries[end], " Infected: ", new_Itimeseries[end], " Sensitive: ", new_Stimeseries[end])
#sleep(10)

time_B, Btimeseries_filtered = new_time, new_Btimeseries/volume
time_I, Itimeseries_filtered = new_time, new_Itimeseries/volume
time_P, Ptimeseries_filtered = new_time, new_Ptimeseries/volume
time_S, Stimeseries_filtered = new_time, new_Stimeseries/volume

if(i==2)
#nonzero_indices = findall(x -> x != 0, Stimeseries)
#println("Non-zero indices: ", nonzero_indices)
#nonzero_indices = findall(x -> x != 0, new_Stimeseries)
#println("Non-zero indices: ", nonzero_indices)
#nonzero_indices = findall(x -> x != 0, Stimeseries_filtered)
#println("Non-zero indices: ", nonzero_indices)
exit()
end
#if(i==1)
#trace=scatter(x = time_B, y = Btimeseries_filtered, mode = "lines", name = "Bacteria $(antiphage_timing_list[i])", line = attr(width = 2))
#layout_timeseries = Layout(
#    xaxis = attr(title = "Time (minutes)"),
#    yaxis = attr(title = "Bacteria (/ml)", type = "log", tickformat = ".0e")  # Set y-axis to logarithmic scale
#)

#plot_timeseries = plot(trace, layout_timeseries)  # Use `plot` instead of `Plot`

# Save the timeseries plot
#savefig(plot_timeseries, joinpath(figures_dir, "Population_LIC($li_collapse)_LICT($li_collapse_phage)_antiphage_test.pdf"))
#end

push!(all_timesB, time_B)
push!(all_Btimeseries, Btimeseries_filtered)
push!(all_timesI, time_I)
push!(all_Itimeseries, Itimeseries_filtered)
push!(all_timesP, time_P)
push!(all_Ptimeseries, Ptimeseries_filtered)
push!(all_timesS, time_S)
push!(all_Stimeseries, Stimeseries_filtered)
end


# Create the timeseries plot using PlotlyJS
colors = ["#0072B2", "#E69F00", "#009E73"]
i=1
tracesB1=scatter(x = all_timesS[1], y = all_Stimeseries[1], mode = "lines", line = attr(width = 1,  dash = "dash", color =  colors[1]),showlegend=false)
tracesI1=scatter(x = all_timesI[i], y = all_Itimeseries[i], mode = "lines", line = attr(width = 1,  dash = "dash", color =  colors[2]),showlegend=false)
tracesP1=scatter(x = all_timesP[i], y = all_Ptimeseries[i], mode = "lines", line = attr(width = 1,  dash = "dash", color =  colors[3]),showlegend=false)
i=2
tracesB2=scatter(x = all_timesS[i], y = all_Stimeseries[i], mode = "lines", line = attr(width = 2,  dash = "solid", color =  colors[1]),showlegend=false)
tracesI2=scatter(x = all_timesI[i], y = all_Itimeseries[i], mode = "lines", line = attr(width = 2,  dash = "solid", color =  colors[2]),showlegend=false)
tracesP2=scatter(x = all_timesP[i], y = all_Ptimeseries[i], mode = "lines", line = attr(width = 2,  dash = "solid", color =  colors[3]),showlegend=false)
i=3
tracesB3=scatter(x = all_timesS[i], y = all_Stimeseries[i], mode = "lines", line = attr(width = 1,  dash = "dot", color =  colors[1]),showlegend=false)
tracesI3=scatter(x = all_timesI[i], y = all_Itimeseries[i], mode = "lines", line = attr(width = 1,  dash = "dot", color =  colors[2]),showlegend=false)
tracesP3=scatter(x = all_timesP[i], y = all_Ptimeseries[i], mode = "lines", line = attr(width = 1,  dash = "dot", color =  colors[3]),showlegend=false)
    #tracesB=scatter(x = all_timesB[1], y = all_Btimeseries[1], mode = "lines", name = "Bacteria $(antiphage_timing_list[1])", line = attr(width = 2))
    #tracesB2=scatter(x = all_timesB[2], y = all_Btimeseries[2], mode = "lines", name = "Bacteria $(antiphage_timing_list[2])", line = attr(width = 2))
    layout_timeseries = Layout(
        title = attr(
            text = "(a)",
            font = attr(size = 24),  # Font size for the title
            x = 0.01,  # X position of the title (0 to 1)
            y = 0.9,  # Y position of the title (0 to 1)
            xanchor = "left",  # Anchor the title at the center
            yanchor = "top"  # Anchor the title at the top
        ),
        xaxis = attr(
            title = attr(
                text = "Time (minutes)",
                font = attr(size = 24)  # Font size for the x-axis label
            ),
            linecolor = "black",
            linewidth = 2,
            ticks = "inside",  # Add ticks inside the plot
            showline = true,  # Show line on the bottom x-axis
            mirror = true,   # Mirror the axis lines on the top and right
            tickfont = attr(size = 20)  # Font size for the x-axis ticks
        ),
        yaxis = attr(
            title = attr(
                text = "Bacteria /ml or Phage /ml",
                font = attr(size = 24)  # Font size for the x-axis label
            ),
            linecolor = "black",
            linewidth = 2,
            ticks = "inside",  # Add ticks inside the plot
            type = "log", 
            tickformat = ".0e",
            showline = true,  # Show line on the left y-axis
            mirror = true,  # Mirror the axis lines on the top and right
            tickfont = attr(size = 20),  # Font size for the x-axis ticks
            range = [log10(1e-1), log10(1e11)],  # Set y-axis range in log scale
            #tickvals = [1e4, 1e6, 1e8,1e10],  # Set specific tick values
            ticktext = ["10$(Char(0x2070 + 4))", "10$(Char(0x2070 + 6))", "10$(Char(0x2070 + 8))","10$(Char(0x00B9)Char(0x2070))"]  # Custom tick labels
        ),
        plot_bgcolor = "rgba(0,0,0,0)",  # Transparent plot background
        paper_bgcolor = "rgba(0,0,0,0)",  # Transparent paper background
        legend = "false"
        )
    #test=[]
    #push!(test, tracesB)
    #push!(test, tracesB2)
    #test=[tracesB, tracesB2]

    #plot_timeseries = plot([tracesB1, tracesI1, tracesP1, tracesB2, tracesI2,tracesP2, tracesB3, tracesI3,tracesP3], layout_timeseries)  # Use `plot` instead of `Plot`
    plot_timeseries = Plot([tracesB1, tracesB2,  tracesB3], layout_timeseries)
    # Save the timeseries plot
    savefig(plot_timeseries, joinpath(figures_dir, "Population_LIC($li_collapse)_LICT($li_collapse_phage)_antiphage.pdf"))

# Create the scatter plot
trace = scatter(x = all_timesS[2], y = all_Stimeseries[2], mode = "lines", line = attr(color = "blue", width = 2))

# Define the layout
layout = Layout(
    title = "Sensitive Bacteria Over Time",
    xaxis = attr(title = "Time (minutes)"),
    yaxis = attr(title = "Sensitive Bacteria (/ml)",type = "log",range = [log10(1e-1), log10(1e11)], tickformat = ".0e")
)

# Combine the trace and layout into a plot
plot_s = Plot([trace], layout)
savefig(plot_s, joinpath(figures_dir, "test.pdf"))


#layout_Ptimeseries = Layout(
#    xaxis = attr(title = "Time (minutes)"),  
#    yaxis = attr(title = "Phage (/ml)", type = "log")  # Set y-axis to logarithmic scale
#)

#plot_timeseries = Plot([trace_P], layout_Ptimeseries)

# Save the timeseries plot
#savefig(plot_timeseries, joinpath(figures_dir, "population_dynamics_LIC($li_collapse)_LICT($li_collapse_phage).pdf"))



# Create the histogram using PlotlyJS
#trace_histogram = histogram(x = lysis_time_record, nbinsx = 200, name = "Lysis Time Histogram")

#layout_histogram = Layout(
#    xaxis = attr(title = "Lysis Time"),
#    yaxis = attr(title = "Frequency"),
#    title = "Histogram of Lysis Time"
#)

#plot_histogram = Plot([trace_histogram], layout_histogram)

# Save the histogram plot
#savefig(plot_histogram, joinpath(figures_dir, "lysis_time_histogram_LIC($li_collapse)_LICT($li_collapse_phage)_antiphage($time_antiphage).pdf"))


#data_file_path = joinpath(figures_dir, "Btimeseries_filtered_LIC($li_collapse)_LICT($li_collapse_phage)_antiphage($time_antiphage).csv")
#CSV.write(data_file_path, DataFrame(time = time_B, Btimeseries = Btimeseries_filtered))

#data_file_path = joinpath(figures_dir, "Itimeseries_filtered_LIC($li_collapse)_LICT($li_collapse_phage)_antiphage($time_antiphage).csv")
#CSV.write(data_file_path, DataFrame(time = time_I, Itimeseries = Itimeseries_filtered))

#data_file_path = joinpath(figures_dir, "Ptimeseries_filtered_LIC($li_collapse)_LICT($li_collapse_phage)_antiphage($time_antiphage).csv")
#CSV.write(data_file_path, DataFrame(time = time_P, Ptimeseries = Ptimeseries_filtered))
end