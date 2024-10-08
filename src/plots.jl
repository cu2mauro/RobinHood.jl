function Z0_Plot(filename,P,zstar)
    file = h5open("results/"*filename,"r")
    Eint=read(file["P$(P)_z$(zstar)/Eint"])
    sols=read(file["P$(P)_z$(zstar)/sols"])
    Lint=read(file["P$(P)_z$(zstar)/Lint"])

    z₀ = sols[:,Int(end*3/4+0.5)]

    plot(plot(z₀, Lint, ylabel="L", lc=:red, xlabel="z₀", label=""),
    plot(z₀, Eint, ylabel="E", lc=:red, xlabel="z₀", label=""),dpi=500)
end
export Z0_Plot

function R0_Plot(filename,P,zstar)
    file = h5open("results/"*filename,"r")
    Eint=read(file["P$(P)_z$(zstar)/Eint"])
    sols=read(file["P$(P)_z$(zstar)/sols"])
    Lint=read(file["P$(P)_z$(zstar)/Lint"])

    r₀ = sols[:,Int(end/4+.5)]

    plot(plot(r₀, Lint, ylabel="L", lc=:blue, xlabel="r₀", label=""),
    plot(r₀, Eint, ylabel="E", lc=:blue, xlabel="r₀", label=""),dpi=500)
end
export R0_Plot

function LE_Plot(filename,P,zstar)
    p=plot(dpi=500)
    file = h5open("results/"*filename,"r")
    Eint=read(file["P$(P)_z$(zstar)/Eint"])
    Lint=read(file["P$(P)_z$(zstar)/Lint"])
    
    p=plot!(Lint, Eint, label="P=$(P), z*=$(zstar)", xlabel="L", ylabel="E", lw=3,dpi=500)
    close(file)
    plot!(dpi=500)
end
function LE_Plot(filename)
    p=plot(dpi=500)
    file = h5open("results/"*filename,"r")

    P_list=read(file["P_list"])
    zstar_list=read(file["zstar_list"])

    for P in P_list
        for zstar in zstar_list
            Eint=read(file["P$(P)_z$(zstar)/Eint"])
            Lint=read(file["P$(P)_z$(zstar)/Lint"])
            p=plot!(Lint, Eint, label="P=$(P), z*=$(zstar)", xlabel="L", ylabel="E", lw=3,dpi=500)
        end
    end
    close(file)
    plot!(dpi=500)
end
export LE_Plot

function String_Plot(filename,P,zstar,i)
    p=plot()
    file = h5open("results/"*filename,"r")

    Eint=read(file["P$(P)_z$(zstar)/Eint"])
    sols=read(file["P$(P)_z$(zstar)/sols"])
    Lint=read(file["P$(P)_z$(zstar)/Lint"])
    p=plot!(sols[i,Int(end/2)+1:end],interval(Int(size(sols,2)/4+0.5),Lint[i]), ylabel="x",xlabel="z",label="",dpi=500)
    close(file)
    plot!(title="P=$(P), z*=$(zstar)",dpi=500)
end
function String_Plot(filename,P,zstar)
    p=plot()
    file = h5open("results/"*filename,"r")

    Eint=read(file["P$(P)_z$(zstar)/Eint"])
    sols=read(file["P$(P)_z$(zstar)/sols"])
    Lint=read(file["P$(P)_z$(zstar)/Lint"])
    for i in 1:1:length(Lint)
        p=plot!(sols[i,Int(end/2)+1:end],interval(Int(size(sols,2)/4+0.5),Lint[i]), ylabel="x",xlabel="z",label="",dpi=500)
    end
    close(file)
    plot!(title="P=$(P), z*=$(zstar)",dpi=500)
end
function String_Plot(filename)
    p=plot(dpi=500)
    file = h5open("results/"*filename,"r")

    P_list=read(file["P_list"])
    zstar_list=read(file["zstar_list"])

    for P in P_list
        for zstar in zstar_list
            Eint=read(file["P$(P)_z$(zstar)/Eint"])
            sols=read(file["P$(P)_z$(zstar)/sols"])
            Lint=read(file["P$(P)_z$(zstar)/Lint"])
            for i in 1:1:length(Lint)
                p=plot!(sols[i,Int(end/2)+1:end],interval(Int(size(sols,2)/4+0.5),Lint[i]), ylabel="x",xlabel="z",label="",dpi=500)
            end
        end
    end
    close(file)
    plot!(dpi=500)
end
export String_Plot

function String_3DPlot(filename,P,zstar,i)
    p=plot()
    file = h5open("results/"*filename,"r")

    Eint=read(file["P$(P)_z$(zstar)/Eint"])
    sols=read(file["P$(P)_z$(zstar)/sols"])
    Lint=read(file["P$(P)_z$(zstar)/Lint"])
    p=plot!(sols[i,Int(end/2)+1:end],sols[i,1:Int(end/2)],interval(Int(size(sols,2)/4+0.5),Lint[i]), ylabel="r",xlabel="z",zlabel="x",label="",dpi=500)
    close(file)
    plot!(title="P=$(P), z*=$(zstar)",dpi=500)
end
function String_3DPlot(filename,P,zstar)
    p=plot()
    file = h5open("results/"*filename,"r")

    Eint=read(file["P$(P)_z$(zstar)/Eint"])
    sols=read(file["P$(P)_z$(zstar)/sols"])
    Lint=read(file["P$(P)_z$(zstar)/Lint"])
    for i in 1:1:length(Lint)
        p=plot!(sols[i,Int(end/2)+1:end],sols[i,1:Int(end/2)],interval(Int(size(sols,2)/4+0.5),Lint[i]), ylabel="r",xlabel="z",zlabel="x",label="",dpi=500)
    end
    close(file)
    plot!(title="P=$(P), z*=$(zstar)",dpi=500)
end
export String_3DPlot