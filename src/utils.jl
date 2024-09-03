function CreateConfigFile(path)
    io = open(path*"/config_file.jl","w")
    println(io,"
        using RobinHood\n\n
        filename=\"quiver_of_choice\"\n\n
        Quiver_choice(n) #use a quiver label among 1,2,3\n\n
        Run_Multiple_Strings(filename)
    ")
    close(io)
end

function CreateConfigFile()
    path=pwd()
    io = open(path*"/config_file.jl","w")
    println(io,"
        using RobinHood\n\n
        filename=\"quiver_of_choice\"\n\n
        Run_Multiple_Strings(filename)
    ")
    close(io)
end