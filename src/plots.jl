using Plots
using HDF5

#include("backgrounds/background_simple.jl")

function Z0_Plot(filename,P,etast)
    file = h5open("results/"*filename,"r")
    Eint=read(file["P$(P)_eta$(etast)/Eint"])
    sols=read(file["P$(P)_eta$(etast)/sols"])
    Lint=read(file["P$(P)_eta$(etast)/Lint"])

    z₀ = sols[:,Int(end*3/4+0.5)]

    plot(plot(z₀, Lint, ylabel="L", lc=:red, xlabel="z₀", label=""),
    plot(z₀, Eint, ylabel="E", lc=:red, xlabel="z₀", label=""),dpi=500)
end
export Z0_Plot

function R0_Plot(filename,P,etast)
    file = h5open("results/"*filename,"r")
    Eint=read(file["P$(P)_eta$(etast)/Eint"])
    sols=read(file["P$(P)_eta$(etast)/sols"])
    Lint=read(file["P$(P)_eta$(etast)/Lint"])

    r₀ = sols[:,Int(end/4+.5)]

    plot(plot(r₀, Lint, ylabel="L", lc=:blue, xlabel="r₀", label=""),
    plot(r₀, Eint, ylabel="E", lc=:blue, xlabel="r₀", label=""),dpi=500)
end
export R0_Plot

function LE_Plot(filename,P,etast)
    p=plot(dpi=500)
    file = h5open("results/"*filename,"r")
    Eint=read(file["P$(P)_eta$(etast)/Eint"])
    Lint=read(file["P$(P)_eta$(etast)/Lint"])
    
    p=plot!(Lint, Eint, label="P=$(P), z*=$(etast)", xlabel="L", ylabel="E", lw=3,dpi=500)
    close(file)
    plot!(dpi=500)
end
function LE_Plot(filename)
    p=plot(dpi=500)
    file = h5open("results/"*filename,"r")

    P_list=read(file["P_list"])
    etast_list=read(file["etast_list"])

    for P in P_list
        for etast in etast_list
            Eint=read(file["P$(P)_eta$(etast)/Eint"])
            Lint=read(file["P$(P)_eta$(etast)/Lint"])
            p=plot!(Lint, Eint, label="P=$(P), z*=$(etast)", xlabel="L", ylabel="E", lw=3,dpi=500)
        end
    end
    close(file)
    plot!(dpi=500)
end
export LE_Plot

function String_Plot(filename,P,etast,i)
    p=plot()
    file = h5open("results/"*filename,"r")

    Eint=read(file["P$(P)_eta$(etast)/Eint"])
    sols=read(file["P$(P)_eta$(etast)/sols"])
    Lint=read(file["P$(P)_eta$(etast)/Lint"])
    p=plot!(sols[i,Int(end/2)+1:end],interval(Int(size(sols,2)/4+0.5),Lint[i]), ylabel="x",xlabel="z",label="",dpi=500)
    close(file)
    plot!(title="P=$(P), z*=$(etast)",dpi=500)
end
function String_Plot(filename,P,etast)
    p=plot()
    file = h5open("results/"*filename,"r")

    Eint=read(file["P$(P)_eta$(etast)/Eint"])
    sols=read(file["P$(P)_eta$(etast)/sols"])
    Lint=read(file["P$(P)_eta$(etast)/Lint"])
    for i in 1:1:length(Lint)
        p=plot!(sols[i,Int(end/2)+1:end],interval(Int(size(sols,2)/4+0.5),Lint[i]), ylabel="x",xlabel="z",label="",dpi=500)
    end
    close(file)
    plot!(title="P=$(P), z*=$(etast)",dpi=500)
end
function String_Plot(filename)
    p=plot(dpi=500)
    file = h5open("results/"*filename,"r")

    P_list=read(file["P_list"])
    etast_list=read(file["etast_list"])

    for P in P_list
        for etast in etast_list
            Eint=read(file["P$(P)_eta$(etast)/Eint"])
            sols=read(file["P$(P)_eta$(etast)/sols"])
            Lint=read(file["P$(P)_eta$(etast)/Lint"])
            for i in 1:1:length(Lint)
                p=plot!(sols[i,Int(end/2)+1:end],interval(Int(size(sols,2)/4+0.5),Lint[i]), ylabel="x",xlabel="z",label="",dpi=500)
            end
        end
    end
    close(file)
    plot!(dpi=500)
end
export String_Plot

function String_3DPlot(filename,P,etast,i)
    p=plot()
    file = h5open("results/"*filename,"r")

    Eint=read(file["P$(P)_eta$(etast)/Eint"])
    sols=read(file["P$(P)_eta$(etast)/sols"])
    Lint=read(file["P$(P)_eta$(etast)/Lint"])
    p=plot!(sols[i,Int(end/2)+1:end],sols[i,1:Int(end/2)],interval(Int(size(sols,2)/4+0.5),Lint[i]), ylabel="r",xlabel="z",zlabel="x",label="",dpi=500)
    close(file)
    plot!(title="P=$(P), z*=$(etast)",dpi=500)
end
function String_3DPlot(filename,P,etast)
    p=plot()
    file = h5open("results/"*filename,"r")

    Eint=read(file["P$(P)_eta$(etast)/Eint"])
    sols=read(file["P$(P)_eta$(etast)/sols"])
    Lint=read(file["P$(P)_eta$(etast)/Lint"])
    for i in 1:1:length(Lint)
        p=plot!(sols[i,Int(end/2)+1:end],sols[i,1:Int(end/2)],interval(Int(size(sols,2)/4+0.5),Lint[i]), ylabel="r",xlabel="z",zlabel="x",label="",dpi=500)
    end
    close(file)
    plot!(title="P=$(P), z*=$(etast)",dpi=500)
end
export String_3DPlot