#Poisson distribution of average 82
using Distributions
       # for formatted printing    

λ = 82  # Mean of the Poisson distribution


poisson_dist = Poisson(λ)


# Example usage of the analytical function
global sumofprob = 0
for k in 0:100
    global sumofprob+= pdf(poisson_dist, k)
    println(k, " ", sumofprob)
    if(sumofprob>0.5)
        println("The 50th percentile is ", k)
        break
    end
end
