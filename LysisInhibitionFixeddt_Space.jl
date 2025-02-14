# This is a Julia file

# Add your code here
using Random
using Plots
using Distributions
using JLD2


#I will now try to create a struct
mutable struct SState
    #This is state of the site. It can be empty or with bacteria, and it also carries the information of phages
    Bstate::UInt16 # 0: Empty, 1 to growth_timer is uninfected, growth_timer+1 to growth_timer+lysis_timer is infected
    # assign a random integer number between 1 and growth_timer to each bacteria
    Istate::UInt32  # Count the number of infection in total
    Pstate::Float64 # recording time spent in infected state, to compute number of produced phages for an infected bacteria, proportional to time spent in infected state (minus eclipse)
    LORstate::Bool # boolean, True if it is in Lysis from without resistant state
    Phage::Int64 # number of phages at the site
end

function custom_mod(i, lattice_size)
    result = mod(i , lattice_size)
    return result == 0 ? lattice_size : result
end

# Main function
function simulate_space_agents(states::Vector{SState}, time_step, record_time_step, final_time, 
    growth_rate, lysis_rate, burst_rate, eclipse, growth_timer, lysis_timer, eta, hop_rate, push_distance; lysis_inhibition=false, lysis_inhibition_timer=5, 
    lysis_from_without=false, lysis_from_without_phage=10, lo_resistance=false, lo_resistance_timer=5, li_collapse=false, li_collapse_phage=100,
    li_collapse_recovery=false, licR_rate=1.0/10.0)
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

    time_series = Float64[]
    Btimeseries = []
    Itimeseries = []
    Ptimeseries = []
    LORtimeseries = []
    Phagetimeseries = []
    totalLO=0

    timenow = 0.
    lattice_size=length(states)
    println("lattice size: ", lattice_size)
    grate = growth_rate * growth_timer
    lrate = lysis_rate * lysis_timer
    probhop = hop_rate * time_step
    lysis_time_record = Float64[]
    itime=Int64(0)
    new_Phage = zeros(Int32, lattice_size)

    push!(time_series, timenow)
    push!(Btimeseries, [state.Bstate for state in states])
    push!(Itimeseries, [state.Istate for state in states])
    push!(Ptimeseries, [state.Pstate for state in states])
    push!(LORtimeseries, [state.LORstate for state in states])
    push!(Phagetimeseries, [state.Phage for state in states])
    while timenow < final_time
        #phage diffusion first 
        new_Phage .= 0
        for i in 1:lattice_size
            jp=custom_mod(i+1, lattice_size)
            jm=custom_mod(i-1+lattice_size, lattice_size)
            if states[i].Phage > 0
                r = rand(Multinomial(states[i].Phage, [probhop, probhop, 1-2*probhop]))
                new_Phage[jp] += r[1]  # Hop to the right
                new_Phage[jm] += r[2]  # Hop to the left
                new_Phage[i] += r[3]  # Stay in the same position
            end 
        end
        for i in 1:lattice_size
            states[i].Phage = new_Phage[i]
        end
        for _ in 1:lattice_size
            i = rand(1:lattice_size)
            if states[i].Bstate > 0
                #println("i: ", i)
                phageinfect= rand(Poisson(eta*states[i].Phage*time_step))   
                #println("phage infection", phageinfect)
                if phageinfect> 0
                    #Only execute the infection actions when it actually hallens!
                    states[i].Istate += phageinfect
                    states[i].Phage = max(0, states[i].Phage - phageinfect)
                    if states[i].Bstate <= growth_timer
                    #first time infection
                        states[i].Bstate = growth_timer+1   
                    end 
                    if lysis_from_without
                    
                        #The lysis from without can happen
                        check_LO = !states[i].LORstate * states[i].Istate  #If true for LOR state, Istate will be multiplied by 0. 
                        if check_LO > lysis_from_without_phage
                            #println("lysis from without")
                            states[i].Bstate = 0
                            states[i].Phage += max(Int(round(burst_rate * (states[i].Pstate - eclipse))), 0)
                            #Append the values of Pstate that match mask_LO to lysis_time_record
                            append!(lysis_time_record, states[i].Pstate)  
                            states[i].LORstate = false
                            states[i].Istate = 0
                            states[i].Pstate = 0.0
                            totalLO+=1               
                        end
                    end
                    if li_collapse
                        if states[i].Istate > li_collapse_phage
                            #println("lysis inhibition collapse")
                            states[i].Bstate = 0
                            states[i].Phage += max(Int(round(burst_rate * (states[i].Pstate - eclipse))), 0)
                            # Append the values of Pstate that match mask_LI to lysis_time_record
                            append!(lysis_time_record, states[i].Pstate)
                            states[i].LORstate = false
                            states[i].Istate = 0
                            states[i].Pstate = 0.0
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
                            states[i].Bstate = max(growth_timer+1, states[i].Bstate - lysis_inhibition_timer*phageinfect)
                        end    
                    end
                end

                # Update Bstate elements less than growth_timer with probability grate * time_step
                #println("growth")
                if states[i].Bstate < growth_timer
                    states[i].Bstate += rand()< (grate * time_step)
                elseif states[i].Bstate == growth_timer
                    if rand() < (grate * time_step)
                        #println("new bacteria")    
                        #Now the new bacteria needs to push the other phage. Also I should introduce the pushing distance. 
                        ip = false
                        im = false
                        jj=0
                        k_now=0
                        jsign=0
                        for k in 1:push_distance
                            jp=custom_mod(i+k, lattice_size)
                            jm=custom_mod(i-k+lattice_size, lattice_size)
                            #println("jp: ", jp, "jm: ", jm, "state ",states[jp].Bstate, " ",states[jm].Bstate)
                            if states[jp].Bstate == 0
                                ip=true
                            end
                            if states[jm].Bstate ==0
                                im=true
                            end
                            if ip || im
                                if ip&&im
                                    if rand(0:1)==0
                                        jj=jp
                                        jsign=1
                                    else    
                                        jj=jm
                                        jsign=-1
                                    end
                                elseif ip
                                    jj=jp
                                    jsign=1
                                else
                                    jj=jm
                                    jsign=-1
                                end
                                k_now=k
                                break
                            end
                        end
                        if jj>0
                            #there is a place to push. 
                            #for the furthest site, I need to make sure phage is not overridden. 
                            states[i].Bstate = 1
                            jk1=custom_mod(jj-(jsign)+lattice_size, lattice_size)
                            states[jj].Bstate = states[jk1].Bstate
                            states[jj].Istate = states[jk1].Istate
                            states[jj].Pstate = states[jk1].Pstate
                            states[jj].LORstate = states[jk1].LORstate
                            states[jj].Phage += states[jk1].Phage
                            if k_now >1
                                for kk in 2:k_now
                                    jk=jk1
                                    jk1=custom_mod(jj-(jsign*kk)+lattice_size, lattice_size)
                                    states[jk] = states[jk1]
                                end
                            end
                            if jk1!=i
                                println("error")
                            end
                            states[i]=SState(1, 0, 0.0, false, 0)
                        end 
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
                    states[i].Bstate += rand() < (lrate * time_step)
                elseif states[i].Bstate == growth_timer + lysis_timer
                    if rand() < (lrate * time_step)
                        states[i].Bstate = 0
                        states[i].Phage += max(Int(round(burst_rate * (states[i].Pstate - eclipse))), 0)
                        # Append the values of Pstate that match mask_lysis to lysis_time_record
                        append!(lysis_time_record, states[i].Pstate)
                        states[i].LORstate = false
                        states[i].Istate = 0
                        states[i].Pstate = 0.0
                    end 
                end 
            end
        end


        #update time and record status when needed

        itime+=1
        timenow = time_step*itime
        if timenow - time_series[end] >= record_time_step
            push!(time_series, timenow)
            push!(Btimeseries, [state.Bstate for state in states])
            push!(Itimeseries, [state.Istate for state in states])
            push!(Ptimeseries, [state.Pstate for state in states])
            push!(LORtimeseries, [state.LORstate for state in states])
            push!(Phagetimeseries, [state.Phage for state in states])
            non_zero_Bstate_count = count(state -> state.Bstate != 0, states)
            infected_Bstate_count = count(state -> state.Bstate > growth_timer, states)
            totalphage=sum([state.Phage for state in states])
            println(timenow, "bacteria count: ", non_zero_Bstate_count, "infected count: ", infected_Bstate_count, "phage count: ", totalphage, "total LO: ", totalLO) 
            if non_zero_Bstate_count==0
                println("All bacteria are dead")
                timenow=final_time
                break
            end
        end
        non_zero_Bstate_count = count(state -> state.Bstate != 0, states)
        if non_zero_Bstate_count==0
            println("All bacteria are dead")
            timenow=final_time
            break
        end
        
    end

    return time, Btimeseries, Itimeseries, Ptimeseries, LORtimeseries, Phagetimeseries, lysis_time_record, states
end

function custom_color_profile(value, growth_timer, lysis_timer)
    if value == 0
        return RGB(0, 0, 0)  # Black
    elseif value <= growth_timer
        blue_intensity = value / growth_timer
        return RGB(0, 0, blue_intensity)  # Blue, getting brighter
    elseif value <= growth_timer + lysis_timer
        yellow_intensity = (value - growth_timer) / lysis_timer
        #return RGB(yellow_intensity, yellow_intensity, yellow_intensity)  # Yellow, getting brighter
        return RGB(1, 1, 0)  # Yellow
    else
        return RGB(1, 1, 1)  # White for values beyond the specified range
    end
end

function apply_custom_colormap(Btimeseries_2d, growth_timer, lysis_timer)
    heatmap_colors = [custom_color_profile(value, growth_timer, lysis_timer) for value in Btimeseries_2d]
    heatmap(heatmap_colors, color=:auto)
end

#Here we define the system parameters.
#We start with the simulation done in the Julia's thesis of different MSOI
growth_rate = 2.0/60. #per minute
lysis_rate = 1.0/27.0  #per minute
growth_timer = 10 #max growth timer
lysis_timer = 250 #max lysis timer, 4 timer is 1 minute
eclipse = 15    #eclipse time in minutes
burst_size = 150 #burst size
burst_rate=burst_size/((1/lysis_rate)-eclipse)
eta = 200  #adsorption rate per box
#If it is 10^-9 ml/min and the volume is 1 micron^3, then this will be 10^3 micron^3 / minute. 
# 200 could mean I am assuming the volume is 5 micron^3. 

lysis_from_without=true
lysis_from_without_phage=50
lo_resistance=true
lo_resistance_timer=Int(round(10*(lysis_timer*lysis_rate)))
lysis_inhibition=true
lysis_inhibition_timer=Int(round(5*(lysis_timer*lysis_rate)))
li_collapse=true
li_collapse_phage=100
hop_rate=50
#This is equal to the phage diffusion constant in the unit of lattice constant^2/min. 
#Phage diffusion is 4 micron^2/sec = 4*60 micron^2/min in water. 
#Lattice constant is 1.7 micron (gives 5 micron^3 volume) then hop rae is about 80. 
#I make is smaller since the available space is small in a densely packed environment.
#Though I am not changing it between free sites and filled site. 
#This matters if I start thinking about the dead body. 
lattice_size=500
push_distance=lattice_size

time_step=min(0.2/hop_rate, 1/eta)
if(time_step*hop_rate>0.5)
    println("Hop rate is too high")
end 


#Now set the initial condition and run the simulation. 
record_time_step = 1500 #minutes

# Initialize the states

iB_series = []
Alive_series = []
#Infected_series = []
#Phage_series = []
#nsamples=100
#for ibacteria in 1:20
#    nalive=0
#    for isamples in 1:nsamples
#        println("ibacteria: ", ibacteria, "isamples: ", isamples)
    states = [SState(0, 0, 0.0, false, 0) for i in 1:lattice_size]
#bacteria = ibacteria #cells
bacteria=10
infected= 0
final_time = 30*25 # minutes
# Generate random values
initial_values = rand(1:growth_timer, bacteria)

for i in 1:bacteria
     states[Int(floor(lattice_size/2 - bacteria/2)) + i].Bstate = growth_timer # rand(1:growth_timer)
end
 states[Int(floor(lattice_size/2 - bacteria/2))].Phage = Int(floor(lysis_from_without_phage / 2))
 states[Int(floor(lattice_size/2 + bacteria/2)) + 1].Phage = Int(floor(lysis_from_without_phage / 2))
# Create directories if they do not exist


time_series, Btimeseries, Itimeseries, Ptimeseries, LORtimeseries, Phagetimeseries, lysis_time_record, states = simulate_space_agents(
        states, time_step, record_time_step, final_time, 
        growth_rate, lysis_rate, burst_rate, eclipse, growth_timer, lysis_timer, eta, hop_rate, push_distance; 
        lysis_inhibition=lysis_inhibition, lysis_inhibition_timer=lysis_inhibition_timer, 
        lysis_from_without=lysis_from_without, lysis_from_without_phage=lysis_from_without_phage, 
        lo_resistance=lo_resistance, lo_resistance_timer=lo_resistance_timer, 
        li_collapse=li_collapse, li_collapse_phage=li_collapse_phage)
    
        non_zero_Bstate_count = count(state -> state.Bstate != 0, states)
        if non_zero_Bstate_count>0
            nalive+=1
            if non_zero_Bstate_count<0.9*lattice_size
                println("unclear ", ibacteria, " ", non_zero_Bstate_count)
            end
        end
    #end
        #infected_Bstate_count = count(state -> state.Bstate > growth_timer, states)
        #totalphage=sum([state.Phage for state in states])
        #push!(iB_series, ibacteria)
        #push!(Alive_series, nalive)
        #push!(Infected_series, infected_Bstate_count)
        #push!(Phage_series, totalphage)
#end
#data_dir = "data_files_space_ctirical"
figures_dir = "figure_files_space"
#mkpath(data_dir)
mkpath(figures_dir)


#plot(iB_series, Alive_series, label="Alive", linewidth=2) #, yscale=:log10)
#plot!(iB_series, Infected_series, label="Infected Bacteria", linewidth=2)
#plot!(time_P, Ptimeseries_filtered, label="Phage", linewidth=2)

# Set labels and title
#xlabel!("Initial Bacteria")
#ylabel!("Probability to be alive")
#title!("Critical size")

# Show legend
#plot!(legend=:topright)

#figure_file_path = joinpath(figures_dir, "Ctirical_size_LO($lysis_from_without)_LI($lysis_inhibition)_LIC($li_collapse).pdf")
#savefig(figure_file_path)
        

#data_file_path = joinpath(data_dir, "Ctirical_size_LO($lysis_from_without)_LI($lysis_inhibition)_LIC($li_collapse).jld2")
#@save  iB_series Alive_series Infected_series Phage_series 

# Create the plot
Btimeseries_2d = hcat(Btimeseries...)'

# Print the type and size of the new 2D array
println("Type of Btimeseries_2d: ", typeof(Btimeseries_2d))
println("Size of Btimeseries_2d: ", size(Btimeseries_2d))
# Find the maximum value in Btimeseries_2d
max_value = maximum(Btimeseries_2d)

# Print the maximum value
println("Maximum value in Btimeseries_2d: ", max_value)

# Assuming Btimeseries_2d is your 2D array of 0 and positive integers
binary_Btimeseries_2d = Btimeseries_2d .> 0
binaryP_Btimeseries_2d = Btimeseries_2d .> growth_timer

# Convert the boolean array to an integer array
binary_Btimeseries_2d = Int.(binary_Btimeseries_2d+binaryP_Btimeseries_2d)
#custom_colors = cgrad([:black, :yellow], [0, 1])

heatmap(binary_Btimeseries_2d, xlabel="position", ylabel="time", title="Heatmap of Btimeseries")

#heatmap(Btimeseries_2d,  xlabel="position", ylabel="time", title="Heatmap of Btimeseries")
#apply_custom_colormap(Btimeseries_2d, growth_timer, lysis_timer)

# Save the heatmap to a file
figure_file_path = joinpath(figures_dir, "Btimeseries_heatmap_lysis_timer($lysis_timer).pdf")
savefig(figure_file_path)

LORtimeseries_2d = hcat(LORtimeseries...)'
heatmap(LORtimeseries_2d, xlabel="position", ylabel="time", title="Heatmap of LORtimeseries")
figure_file_path = joinpath(figures_dir, "LORtimeseries_heatmap_lysis_timer($lysis_timer).pdf")
savefig(figure_file_path)

Phagetimeseries_2d = hcat(Phagetimeseries...)'
heatmap(LORtimeseries_2d, xlabel="position", ylabel="time", title="Heatmap of Phagetimeseries")
figure_file_path = joinpath(figures_dir, "Phagetimeseries_heatmap_lysis_timer($lysis_timer).pdf")
savefig(figure_file_path)

Itimeseries_2d = hcat(Itimeseries...)'
heatmap(Itimeseries_2d, xlabel="position", ylabel="time", title="Heatmap of Itimeseries")
figure_file_path = joinpath(figures_dir, "Itimeseries_heatmap_lysis_timer($lysis_timer).pdf")
savefig(figure_file_path)

Ptimeseries_2d = hcat(Ptimeseries...)'
heatmap(Ptimeseries_2d, xlabel="position", ylabel="time", title="Heatmap of Ptimeseries")
figure_file_path = joinpath(figures_dir, "Ptimeseries_heatmap_lysis_timer($lysis_timer).pdf")
savefig(figure_file_path)

