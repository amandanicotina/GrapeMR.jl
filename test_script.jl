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
