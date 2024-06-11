t = Template(;
           user="amandanicotina",
           authors=["Amanda Nicotina", "David Dodel"],
           plugins=[
               License(name="MIT"),
               Git(),
               GitHubActions(),
           ],
       )



##########################################################################################


no_threads = [295.8, 295.9, 289.6, 297.7, 294.3]
multi_threads = [191.5 212.4 190.7 194.0 191.7]
speed_up = [1.54 1.39 1.52 1.53 1.54]

using DataFrames

# Create the DataFrame
data = DataFrame(
    TimeSingleThread = [295.8, 295.9, 289.6, 297.7, 294.3],
    TimeMultiThread = [191.5, 212.4, 190.7, 194.0, 191.7],
    SpeedUp = [1.54, 1.39, 1.52, 1.53, 1.54])

# Display the table
println(data)
using Plots

# Data
input_sizes = [1000, 5000, 10000, 25000, 50000]
time_single_thread = [295.8, 295.9, 289.6, 297.7, 294.3]
time_multi_thread = [191.5, 212.4, 190.7, 194.0, 191.7]
speed_up = [1.54, 1.39, 1.52, 1.53, 1.54]

# Plot Execution Time
p1 = plot(time_single_thread, label="Single Thread", marker=:o, ylabel="Execution Time (s)", title="Execution Time Comparison", legend=:topleft)
plot!(p1, time_multi_thread, label="Multi Thread", marker=:o)

# Plot Speedup
p2 = plot(speed_up, label="Speedup", marker=:o, ylabel="Speedup", title="Speedup Achieved", legend=:topleft)

# Display the plots
plot(p1, p2, layout=(1, 2), size=(1000, 400))
