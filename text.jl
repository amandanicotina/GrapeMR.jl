t = Template(;
           user="amandanicotina",
           authors=["Amanda Nicotina", "David Dodel"],
           plugins=[
               License(name="MIT"),
               Git(),
               GitHubActions(),
           ],
       )

using GrapeMR
using Flux
using Zygote

# Perform forward pass to calculate cost
cost_function = cost_functions["Euclidean Norm"](isoInit)
cost_value, back = Zygote.forward(cost_function, iso)

# Print the cost value
println("Cost Value: ", cost_value)

# Perform backward pass to calculate gradients
grads = Zygote.gradient(iso -> cost_function(iso), isoInit)

# Print the gradients
println("Gradients: ", grads)


using DifferentialEquations
using Flux
using Zygote


# Define the Bloch equations
function bloch_equations!(du, u, p, t)
    γ, T1, T2, ω0, B1_func = p
    
    Mx, My, Mz = u
    
    # Evaluate B1 at the current time
    B1 = B1_func(t)
    
    # Bloch equations in matrix form
    A = [-1/T2 ω0 -B1;
         -ω0 -1/T2 -B1;
         0 0 -1/T1]
    
    dM = A * [Mx; My; Mz]
    
    du[1] = dM[1]
    du[2] = dM[2]
    du[3] = dM[3]
end

# Define a time-dependent B1 function
function time_dependent_B1(t)
    return 2π * sin(2π * 1e6 * t)
end

# Cost function
function euclidean_norm(iso::Magnetization)
    mag = iso.magnetization[1]
    Mᴬx = mag[2, end]
    Mᴬy = mag[3, end]
    Mᴬz = mag[4, end]

    J = sqrt(Mᴬx^2 + Mᴬy^2 + Mᴬz^2)

    return J
end

# Function to perform Bloch simulation and calculate cost
function simulate_and_calculate_cost(params)
    # Set up parameters
    γ, T1, T2, ω0, B1_func = params
    
    # Set initial conditions
    u0 = [0.0, 0.0, 1.0]
    
    # Set time span
    tspan = (0.0, 5.0e-3)
    
    # Define the problem
    prob = ODEProblem(bloch_equations!, u0, tspan, [γ, T1, T2, ω0, B1_func])
    
    # Solve the problem
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
    
    # Extract final magnetization
    mag = Magnetization(sol.u)
    
    # Calculate cost using the provided cost function
    cost = euclidean_norm(mag)
    
    return cost
end

# Set initial parameter values
initial_params = [42.58e6, 1.0, 0.1, 2π * 10e6, time_dependent_B1]

# Perform forward pass to calculate cost
cost_value = simulate_and_calculate_cost(initial_params)

# Print the cost value
println("Cost Value: ", cost_value)

# Perform backward pass to calculate gradients
grads = Zygote.gradient(params -> simulate_and_calculate_cost(params)[1], initial_params)

# Print the gradients
println("Gradients: ", grads)


using DifferentialEquations
using Flux
using Zygote

struct Magnetization
    magnetization::Vector{Vector{Float64}}
end
# Define the Bloch equations
function bloch_equations!(du, u, p, t)
    γ, T1, T2, ω0, B1_func = p
    
    Mx, My, Mz = u
    
    # Evaluate B1 at the current time
    B1 = B1_func(t)
    
    # Bloch equations in matrix form
    A = [-1/T2 ω0 -B1;
         -ω0 -1/T2 -B1;
         0 0 -1/T1]
    
    dM = A * [Mx; My; Mz]
    
    du[1] = dM[1]
    du[2] = dM[2]
    du[3] = dM[3]
end

# Define a time-dependent B1 function
function time_dependent_B1(t)
    return 2π * sin(2π * 1e6 * t)
end

# Cost function
function euclidean_norm(iso::Magnetization)
    mag = iso.magnetization[1]
    Mᴬx = mag[2, end]
    Mᴬy = mag[3, end]
    Mᴬz = mag[4, end]

    J = sqrt(Mᴬx^2 + Mᴬy^2 + Mᴬz^2)

    return J
end

# Function to perform Bloch simulation and calculate cost
function simulate_and_calculate_cost(params)
    # Set up parameters
    γ, T1, T2, ω0, B1_func = params
    
    # Set initial conditions
    u0 = [0.0, 0.0, 1.0]
    
    # Set time span
    tspan = (0.0, 5.0e-3)
    
    # Define the problem
    prob = ODEProblem(bloch_equations!, u0, tspan, [γ, T1, T2, ω0, B1_func])
    
    # Solve the problem
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
    
    # Extract final magnetization
    mag = Magnetization(sol.u)
    
    # Calculate cost using the provided cost function
    cost = euclidean_norm(mag)
    
    return cost
end

# Set initial parameter values
initial_params = [42.58e6, 1.0, 0.1, 2π * 10e6, time_dependent_B1]

# Perform forward pass to calculate cost
cost_value = simulate_and_calculate_cost(initial_params)

# Print the cost value
println("Cost Value: ", cost_value)

# Perform backward pass to calculate gradients
grads = Zygote.gradient(params -> simulate_and_calculate_cost(params)[1], initial_params)

# Print the gradients
println("Gradients: ", grads)
