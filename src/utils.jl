function CreateConfigFile(path)
    io = open(path*"/config_file.jl","w")
    println(io,"using RobinHood

global zstar_list=[] #list of values of z*
global P_list=[] #list of values of P

filename=\"quiver_of_choice\" #the name for the .h5 data file that will be created
Quiver_choice(n) #use a quiver label among 1,2,3

Run_Multiple_Strings(filename,P_list,zstar_list)")
    close(io)
    mkpath("./results")
    return
end

function CreateConfigFile()
    path=pwd()
    io = open(path*"/config_file.jl","w")
    println(io,"using RobinHood

global zstar_list=[] #list of values of z*
global P_list=[] #list of values of P

filename=\"quiver_of_choice\" #the name for the .h5 data file that will be created
Quiver_choice(n) #use a quiver label among 1,2,3

Run_Multiple_Strings(filename,P_list,zstar_list)")
    close(io)
    mkpath("./results")
    return
end

export CreateConfigFile