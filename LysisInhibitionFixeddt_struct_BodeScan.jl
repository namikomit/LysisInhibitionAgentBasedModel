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


# Read the CSV file
file_path = "BodeData/BodeFig2D.csv"
data = CSV.read(file_path, DataFrame; header=false)

# Assign custom names to the columns
rename!(data, [:time_from_infection, :phage_flux])

data.phage_flux = data.phage_flux ./maximum(data.phage_flux)

# Extract the columns
#x = data.time_from_infection
#y = data.phage_flux

# Create the plot using PlotlyJS
#trace = scatter(x=x, y=y, mode="lines+markers")
#layout = Layout(title="Phage Flux Over Time", xaxis_title="Time from Infection", yaxis_title="Phage Flux")
#plot = Plot([trace], layout)

# Save the plot as an HTML file to preserve interactivity
#savefig(plot, "BodeFig2D_plot.html")

# Optionally, save the plot as a static image (e.g., PNG)
#savefig(plot, "BodeFig2D_plot.png")


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
eta0=2.e-9

#Now set the initial condition and run the simulation. 
record_time_step = 1 #minutes
P_diff_matrix = []
time_save_list = []
P_diff_error   = []

lysis_timer_flag=true
# Create directories if they do not exist
figures_dir = "Bode_Figures_Paper"
mkpath(figures_dir)
if(lysis_timer_flag)
    lysistimer_try = [50,100,150,200,250,300]
    for i in 1:length(lysistimer_try)
    #for i in 1:1
        # Set the initial conditions
        lysis_timer = lysistimer_try[i] #max lysis timer
        lysis_inhibition_timer=Int(round(5*(lysis_timer*lysis_rate)))
        lo_resistance_timer=Int(round(10*(lysis_timer*lysis_rate)))
        volume = 0.01 #ml
        nutrient = Int(round(1.e9*volume)) #cells/ml, growth rate does not depends on it but growth stops if bacteria number reach nutrient
        bacteria = Int(round(2e7*volume)) #cells
        infected=bacteria
        #bacteria=bacteria+infected
        si_duration=3. #minutes
        eta=eta0
    
        P0 = 0
        si_time = 15. # minutes
        final_time = si_time # minutes
     
            # Generate random values
    initial_values = rand(1:growth_timer, bacteria-infected)
    append!(initial_values, [growth_timer + 1 for _ in 1:infected])
    #make it into an array with default values
    states = [State(initial_values[i], 1, 0.0, false) for i in 1:bacteria]

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
    eta=0.0  # because phages were continuously removed
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

    P_diff = diff(Ptimeseries2)
    P_diff = P_diff./maximum(P_diff)
    time_save = time2[2:end] .+ 18
        # Save P_diff and time_save for each iteration
        push!(P_diff_matrix, P_diff)
        push!(time_save_list, time_save)

        #Calculate the deviation from the data
        Pdifferror=0
        for j in 1:length(data.time_from_infection)
            for k in 1:length(time_save)-1
                if data.time_from_infection[j]> time_save[k] && data.time_from_infection[j]< time_save[k+1]
                    tdiff=(time_save[k+1]-time_save[k])
                    P_diff_estimate = P_diff[k]*(data.time_from_infection[j]-time_save[k]) + (P_diff[k+1]*(time_save[k+1]-data.time_from_infection[j]))/tdiff
                    Pdifferror+=(P_diff_estimate-data.phage_flux[j])^2
                    break
                end
            end
        end
        push!(P_diff_error, Pdifferror)
    end
        # Extract the columns
        x1 = data.time_from_infection
        y1 = data.phage_flux
    
        # Create the plot using PlotlyJS
        trace1 = scatter(x=x1, y=y1, mode="lines+markers", name="Phage Flux Data")
        traces = [trace1]
    for i in 1:length(P_diff_matrix)
    trace = scatter(x=time_save_list[i], y=P_diff_matrix[i], mode="lines+markers", name="LI timer $(lysistimer_try[i])") #, error $(P_diff_error[i])")
    push!(traces, trace)
    end
    
# Create layout with transparent background, custom legend position, ticks, and lines on top and right
layout = Layout(
    xaxis = attr(
        title = "Time",
        linecolor = "black",
        linewidth = 2,
        ticks = "inside",  # Add ticks inside the plot
        showline = true,  # Show line on the bottom x-axis
        mirror = true  # Mirror the axis lines on the top and right
    ),
    yaxis = attr(
        title = "Phage Flux",
        linecolor = "black",
        linewidth = 2,
        ticks = "inside",  # Add ticks inside the plot
        showline = true,  # Show line on the left y-axis
        mirror = true  # Mirror the axis lines on the top and right
    ),
    plot_bgcolor = "rgba(0,0,0,0)",  # Transparent plot background
    paper_bgcolor = "rgba(0,0,0,0)",  # Transparent paper background
    legend = attr(
        x = 0.7,  # X position of the legend (0 to 1)
        y = 0.9,  # Y position of the legend (0 to 1)
        bgcolor = "rgba(255, 255, 255, 0.5)",  # Background color of the legend
        bordercolor = "black",  # Border color of the legend
        borderwidth = 1  # Border width of the legend
    )
)
        plot = Plot(traces, layout)
    
    
    # Optionally, save the plot as a static image (e.g., PNG)
    figure_file_path = joinpath(figures_dir, "BodeFig2D_comparison_plot_lt$(1/lysis_rate).pdf")
        savefig(plot, figure_file_path)
else
api = collect(0:0.25:25)
for i in 1:length(api)
    println("api", api[i])
    volume = 0.01 #ml
    nutrient = Int(round(1.e9*volume)) #cells/ml, growth rate does not depends on it but growth stops if bacteria number reach nutrient
    bacteria = Int(round(2e7*volume)) #cells
    infected=bacteria
    #bacteria=bacteria+infected
    si_duration=3. #minutes
    eta=eta0

    P0 = api[i] * ((Float64(bacteria) / volume)) 
    si_time = 15. # minutes
    final_time = si_time # minutes


    # Generate random values
    initial_values = []
    append!(initial_values, [growth_timer + 1 for _ in 1:infected])
    #make it into an array with default values
    states = [State(initial_values[i], 1, 0.0, false) for i in 1:bacteria]

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

#I want to get the histogram of the superinfection after 3 min for given MSOI
#Istate_values = [state.Istate for state in states]

# Create a histogram of Istate values using PlotlyJS
#trace_hist = histogram(x=Istate_values, nbinsx=20)
#layout_hist = Layout(title="Histogram of Istate", xaxis_title="Istate", yaxis_title="Frequency")
#plot_hist = Plot([trace_hist], layout_hist)

# Save the histogram as an HTML file to preserve interactivity
#savefig(plot_hist, "Istate_histogram.html")

# Optionally, save the histogram as a static image (e.g., PNG)
#savefig(plot_hist, "Istate_histogram.png")

    phage = 0
    final_time = 40 # minutes
    eta=0.0  # because phages were continuously added
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


    # Create the plot
    #plot(time2 .+ 18, Btimeseries2, label="Bacteria", linewidth=2)
    #plot!(time2 .+ 18, Itimeseries2, label="Infected Bacteria", linewidth=2)
    #plot!(time2 .+ 18, Ptimeseries2, label="Phage", linewidth=2)

    # Set labels and title
    #xlabel!("Time (minutes)")
    #ylabel!("Population")
    #title!("Population Dynamics")

    # Show legend
    #plot!(legend=:topright)

    #figure_file_path = joinpath(figures_dir, "population_dynamics_plot_lysis_timer($lysis_timer)_lysis_inhibition_timer($lysis_inhibition_timer)_MSOI$(msoi).pdf")
    #savefig(figure_file_path)



    # Calculate the difference in Phage levels
    P_diff = diff(Ptimeseries2)

    push!(P_diff_matrix, P_diff)

    # Create a new figure for the difference plot
    #plot(size=(800, 480))

    # Plot the difference in Phage levels
    #plot!(time2[2:end] .+ 18, P_diff, label="Difference in Phage Level")

    # Label the axes
    #xlabel!("Time (minutes)")
    #ylabel!("Difference in Phage Level")

    if i==1 
        time_save = time2[2:end] .+ 18
        push!(time_save_list,time_save)
    end

# Create the histogram
#histogram(lysis_time_record2, bins=200, label="lysis time histogram", xlabel="lysis time", ylabel="Frequency", title="Histogram of lysis time")
#figure_file_path = joinpath(figures_dir, "lysis_time_histogram_lysis_timer($lysis_timer)_lysis_inhibition_timer($lysis_inhibition_timer)_MSOI$(msoi).pdf")
#savefig(figure_file_path)
println("done")
end


# Convert P_diff_matrix to a 2-dimensional array
P_diff_matrix_2d =  hcat(P_diff_matrix...)
P_diff_matrix_2d = P_diff_matrix_2d'

# Prepare data for 3D plot
time = time_save_list[1]  # Assuming all time_save arrays are the same

println(P_diff_matrix_2d)
println(api)
array_size = size(P_diff_matrix_2d)
println("The size of P_diff_matrix_2d is: ", array_size)
# Ensure dimensions match for heat map
println("length(time): ", length(time))
println("size(P_diff_matrix_2d, 1): ", size(P_diff_matrix_2d, 1))
println("size(P_diff_matrix_2d, 2): ", size(P_diff_matrix_2d, 2))
println("length(api): ", length(api))
# Create the 3D plot
# Ensure dimensions match for 3D plot
#println("length", length(time),size(P_diff_matrix, 1), size(P_diff_matrix, 2), length(msoi_values))

#gr()
#trace = surface(x = time, y = api, z = P_diff_matrix)
#layout = Layout(title = "3D Surface Plot", scene = attr(xaxis_title = "Time (minutes)", yaxis_title = "MSOI", zaxis_title = "Phage production"))
#plot = Plot([trace], layout)
#contour(x = time, y = api, z = P_diff_matrix, contours = attr(show=true, start=0, end=0, size=0.5, coloring="heatmap"))
# Create the 3D plot using gr backend
#plot(time, msoi_values, P_diff_matrix, st = :surface, xlabel = "Time (minutes)", ylabel = "MSOI", zlabel = "Difference in Phage Level", title = "3D Plot of Phage Difference")
# Save the 3D plot as PDF
#savefig("contour_plot.png")
heatmap_trace = heatmap(x = time, y = api, z = P_diff_matrix_2d)
contour_trace = contour(
    x = time,
    y = api,
    z = P_diff_matrix_2d,
    line = attr(color = "black"),
    contours = attr(showlabels = true, labelfont = attr(size = 12, color = "black")),
    showscale = false  # Hide the contour legend
)
layout = Layout(xaxis_title = "Time (minutes)", yaxis_title = "API")
plot = Plot([heatmap_trace, contour_trace], layout)

#trace = heatmap(x = time, y = api, z = P_diff_matrix_2d) #, contours = attr(show=true, size=0.5, coloring="heatmap"))
#layout = Layout(title = "Phage flux", xaxis_title = "Time (minutes)", yaxis_title = "API")
#plot = Plot([trace], layout)

# Save the plot as an HTML file to preserve interactivity
#savefig(plot, "heat_and_contour_plot_eta($eta0)_burst($burst_size)_lysis_timer($lysis_timer).html")

# Optionally, save the plot as a static image (e.g., PNG)
figure_file_path = joinpath(figures_dir, "heat_and_contour_plot_eta($eta0)_burst($burst_size)_lysis_timer($lysis_timer).pdf")
savefig(plot, figure_file_path)
#savefig(plot, "heat_and_contour_plot_eta($eta).eps")
#savefig(plot, "heat_and_contour_plot_eta($eta).svg")

# Save P_diff_matrix_2d, time, and api to CSV files
P_diff_matrix_df = DataFrame(Matrix(P_diff_matrix_2d), :auto)
data_file_path = joinpath(figures_dir, "P_diff_matrix_eta($eta0)_burst($burst_size)_lysis_timer($lysis_timer).csv")
CSV.write(data_file_path, P_diff_matrix_df)

time_df = DataFrame(time = time)
data_file_path = joinpath(figures_dir, "time_eta($eta0)_burst($burst_size)_lysis_timer($lysis_timer).csv")
CSV.write(data_file_path, time_df)

api_df = DataFrame(api = api)
data_file_path = joinpath(figures_dir, "api_burst($burst_size)_lysis_timer($lysis_timer).csv")
CSV.write(data_file_path, api_df)
end
