t = Template(;
           user="amandanicotina",
           authors=["Amanda Nicotina", "David Dodel"],
           plugins=[
               License(name="MIT"),
               Git(),
               GitHubActions(),
           ],
       )



# Define the file name and path
output_file = "output.txt"

# Open the file in write mode and write content
open(output_file, "w") do file
    # Write content to the file
    write(file, "This is an example of writing to a text file in Julia.\n")
    write(file, "Julia is fun!!\n")
    write(file, "PhD, on the other hand, not so fun.")
end

println("File written successfully to $output_file")

using GrapeMR

spin = Spin([0.0,0.0,1.0], [0.1], [0.01], [0.0], [1.0], ["bla"], ["bla"]);
control_field = @time spline_RF(1000, 0.5, 1.0) 
iso = dynamics(control_field, spin[1])

@time begin
    adj_ini = [1.0, 2.0, 3.0, 4.0];
    adj = backward_propagation(adj_ini, control_field, iso);
end

@time begin
    χ0 = [1.0, 2.0, 3.0, 4.0] 
    χ  = hcat(zeros(eltype(χ0), 4, size(control_field.B1x)[2]), χ0);
    χ  = backward_propagation!(χ, control_field, iso);
end
