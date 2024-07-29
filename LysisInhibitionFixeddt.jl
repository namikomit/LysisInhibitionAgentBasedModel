# This is a Julia file

# Add your code here
using Random
using Plots
using Distributions

# Main function
function simulate_population_agents(Bstate, Pstate, Istate, LORstate, time_step, record_time_step, final_time, bacteria, phage, infected, bacteriaever, carrying_capacity, 
    growth_rate, lysis_rate, burst_rate, eclipse, growth_timer, lysis_timer, eta; lysis_inhibition=false, lysis_inhibition_timer=5, 
    lysis_from_without=false, lysis_from_without_phage=10, lo_resistance=false, lo_resistance_time=5, li_collapse=false, li_collapse_phage=100)
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
    volume = carrying_capacity / 1e9  # volume in ml

    timenow = 0.0
    nrecord = Int(round(final_time / record_time_step))
    time = zeros(nrecord)
    Btimeseries = zeros(nrecord)
    Itimeseries = zeros(nrecord)
    Ptimeseries = zeros(nrecord)
    irecord = 1
    time[irecord] = timenow  # Julia arrays are 1-indexed
    Btimeseries[irecord] = bacteria
    Itimeseries[irecord] = infected
    Ptimeseries[irecord] = phage
    grate = growth_rate * growth_timer
    lrate = lysis_rate * lysis_timer
    println(volume)
    while timenow < final_time
        phagegone=0
        phagenew=0
        for j in 1:bacteriaever
            if Bstate[j]>0
                phageinfect= rand(Poisson(eta*phage/volume*time_step))
                if phageinfect > 0
                    #phage infection actions
                    #println("phage infection")
                    phagegone += phageinfect
                    Istate[j] += phageinfect
                    if Bstate[j]<= growth_timer
                        #This was the first infection, so it moves to the infected state
                        #Note that multiple infection does not affect lysis timer since it is anyway zero
                        infected+=1
                        Bstate[j] = growth_timer + 1
                    end
                    if !LORstate[j] && lysis_from_without
                        #The lysis from without can happen
                        if Istate[j] > lysis_from_without_phage
                            Bstate[j] = 0
                            infected -= 1
                            bacteria -= 1
                            phagenew += max(Int(round(burst_rate * (Pstate[j] - eclipse))), 0)
                        end
                    elseif LORstate[j] || !lysis_from_without
                        #lysis from without cannot happen
                        if lysis_inhibition
                            #lysis inbition can happen
                            Bstate[j] = max(growth_timer + 1, Bstate[j] - lysis_inhibition_timer * phageinfect)
                        end
                        if li_collapse && Istate[j] > li_collapse_phage
                            #lysis inbibition collapse happens
                                Bstate[j] = 0
                                infected -= 1
                                bacteria -= 1
                                phagenew += max(Int(round(burst_rate * (Pstate[j] - eclipse))), 0)
                        end
                    end
                end
                if Bstate[j] < growth_timer
                    #growth actions
                    #println("growth")
                    if rand() < grate * time_step
                        println(j)
                        Bstate[j] += 1
                        println(Bstate)
                    end
                elseif Bstate[j] == growth_timer
                    #can divide
                    #println("division")
                    if rand() < grate * time_step
                        Bstate[j] = 1
                        bacteria += 1
                        bacteriaever += 1
                        println("bacteriaever ", bacteriaever)
                        if bacteriaever < carrying_capacity
                            Bstate[bacteriaever] = 1
                        else
                            println("Carrying capacity reached")
                            timenow = final_time
                            break
                        end
                    end
                elseif Bstate[j] < growth_timer+lysis_timer
                    #lysis actions
                    #println("lysis actions")
                    Pstate[j]+=time_step
                    if lo_resistance
                    #LO resistance can happen
                        if !LORstate[j] && Pstate[j] > lo_resistance_time
                            LORstate[j] = true
                        end
                    end
                    if rand() < lrate * time_step
                        Bstate[j]+=1
                    end
                elseif Bstate[j] == growth_timer+lysis_timer
                    #can lyse
                    #println("lysis")
                    if rand() < lrate * time_step
                        #lysis actions
                        Bstate[j] = 0
                        infected -= 1
                        bacteria -= 1
                        phagenew += max(Int(round(burst_rate * (Pstate[j] - eclipse))), 0)
                    end
                end
            end
        end
        phage=max(0,phage+phagenew-phagegone)
        timenow += time_step
        if timenow - time[irecord] >= record_time_step
            irecord += 1
            time[irecord] = timenow
            Btimeseries[irecord] = bacteria
            Itimeseries[irecord] = infected
            Ptimeseries[irecord] = phage
            println(irecord, " ", time[irecord], " ", Btimeseries[irecord], " ", Itimeseries[irecord], " ", Ptimeseries[irecord])
        end
    end

    return time, Btimeseries, Itimeseries, Ptimeseries, irecord
end


#Here we define the system parameters.
#We start with the simulation done in the Julia's thesis of different MSOI
growth_rate = 0 #per minute
lysis_rate = 1/23  #per minute
growth_timer = 10 #max growth timer
lysis_timer = 23 #max lysis timer
eclipse = 15    #eclipse time in minutes
burst_size = 100 #burst size
burst_rate=burst_size/(1/lysis_rate-eclipse)
eta = 5e-9  #adsorption rate per ml/min
lysis_inhibition=true
lysis_inhibition_timer=5
lysis_from_without=false
lysis_from_without_phage=10
lo_resistance=false
lo_resistance_time=5
li_collapse=false
li_collapse_phage=80
time_step=0.01

#Now set the initial condition and run the simulation. 
record_time_step = 1 #minutes

factor=0.01
carrying_capacity = Int(round(1e9*factor)) #max cells
bacteria = Int(round(2e5)) #cells
infected=Int(round(1e5))
si_duration=3. #minutes
msoi=0. #"=P_0(1-exp(-eta*B*si_duration))/(eta*B^2)"
P0 = msoi * (eta * ((bacteria + infected) / factor)^2) / (1 - exp(-eta * Float64(bacteria + infected) / factor * si_duration))
phage = Int(round(P0 * factor)) # pfu
si_time = 15. # minutes
final_time = si_time # minutes

Bstate = zeros(Int, carrying_capacity) # 1 to growth_timer is uninfected, growth_timer+1 to growth_timer+lysis_timer is infected
# assign a random integer number between 1 and growth_timer to each bacteria
Pstate = zeros(Float64, carrying_capacity) # recording time spent in infected state, to compute number of produced phages for an infected bacteria, proportional to time spent in infected state (minus eclipse)
Istate = zeros(Int, carrying_capacity) # Count the number of infection in total
LORstate = falses(carrying_capacity) # boolean, True if it is in Lysis from without resistant state

Bstate[1:(bacteria-infected)] .= rand(1:growth_timer, bacteria-infected)
Bstate[(bacteria-infected+1):bacteria] .= growth_timer + 1
bacteriaever = bacteria

println(phage)

# Call the main function
time, Btimeseries, Itimeseries, Ptimeseries, irecord = simulate_population_agents(
    Bstate, Pstate, Istate, LORstate, time_step, record_time_step, final_time, 
    bacteria, phage, infected, bacteriaever, carrying_capacity, growth_rate, 
    lysis_rate, burst_rate, eclipse, growth_timer, lysis_timer, eta; 
    lysis_inhibition=lysis_inhibition, lysis_inhibition_timer=lysis_inhibition_timer, 
    lysis_from_without=lysis_from_without, lysis_from_without_phage=lysis_from_without_phage, 
    lo_resistance=lo_resistance, lo_resistance_time=lo_resistance_time, 
    li_collapse=li_collapse, li_collapse_phage=li_collapse_phage
)


phage = 0
mask = (Bstate .>= 1) .& (Bstate .<= growth_timer)
bacteria = sum(mask)
mask = (Bstate .>= growth_timer + 1)
infected = sum(mask)
final_time = 60 # minutes

time2, Btimeseries2, Itimeseries2, Ptimeseries2, irecord2 = simulate_population_agents(
    Bstate, Pstate, Istate, LORstate, time_step, record_time_step, final_time, 
    bacteria, phage, infected, bacteriaever, carrying_capacity, growth_rate, 
    lysis_rate, burst_rate, eclipse, growth_timer, lysis_timer, eta; 
    lysis_inhibition=lysis_inhibition, lysis_inhibition_timer=lysis_inhibition_timer, 
    lysis_from_without=lysis_from_without, lysis_from_without_phage=lysis_from_without_phage, 
    lo_resistance=lo_resistance, lo_resistance_time=lo_resistance_time, 
    li_collapse=li_collapse, li_collapse_phage=li_collapse_phage
)

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



