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
function LE_Plot!(filename,P,zstar)
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
function LE_Plot!(filename)
    file = h5open("results/"*filename,"r")

    P_list=read(file["P_list"])
    zstar_list=read(file["zstar_list"])

    for P in P_list
        for zstar in zstar_list
            Eint=read(file["P$(P)_z$(zstar)/Eint"])
            Lint=read(file["P$(P)_z$(zstar)/Lint"])
            plot!(Lint, Eint, label="P=$(P), z*=$(zstar)", xlabel="L", ylabel="E", lw=3,dpi=500)
        end
    end
    close(file)
    plot!(dpi=500)
end
export LE_Plot!

function String_Plot(filename,P,zstar,i)
    p=plot()
    file = h5open("results/"*filename,"r")

    sols=read(file["P$(P)_z$(zstar)/sols"])
    Lint=read(file["P$(P)_z$(zstar)/Lint"])
    p=plot!(sols[i,Int(end/2)+1:end],interval(Int(size(sols,2)/4+0.5),Lint[i]), ylabel="x",xlabel="z",label="",dpi=500)
    close(file)
    plot!(title="P=$(P), z*=$(zstar)",dpi=500)
end
function String_Plot(filename,P,zstar)
    p=plot()
    file = h5open("results/"*filename,"r")

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

    sols=read(file["P$(P)_z$(zstar)/sols"])
    Lint=read(file["P$(P)_z$(zstar)/Lint"])
    p=plot!(interval(Int(size(sols,2)/4+0.5),Lint[i]),sols[i,1:Int(end/2)], sols[i,Int(end/2)+1:end], ylabel="r",xlabel="x",zlabel="z",label="",dpi=500)
    close(file)
    plot!(title="P=$(P), z*=$(zstar)",dpi=500)
end
function String_3DPlot!(filename,P,zstar,i)
    file = h5open("results/"*filename,"r")

    sols=read(file["P$(P)_z$(zstar)/sols"])
    Lint=read(file["P$(P)_z$(zstar)/Lint"])
    plot!(interval(Int(size(sols,2)/4+0.5),Lint[i]),sols[i,1:Int(end/2)], sols[i,Int(end/2)+1:end], ylabel="r",xlabel="x",zlabel="z",label="",dpi=500)
    close(file)
    plot!(title="P=$(P), z*=$(zstar)",dpi=500)
end
function String_3DPlot(filename,P,zstar)
    p=plot()
    file = h5open("results/"*filename,"r")

    sols=read(file["P$(P)_z$(zstar)/sols"])
    Lint=read(file["P$(P)_z$(zstar)/Lint"])
    for i in 1:1:length(Lint)
        p=plot!(interval(Int(size(sols,2)/4+0.5),Lint[i]),sols[i,1:Int(end/2)], sols[i,Int(end/2)+1:end], ylabel="r",xlabel="x",zlabel="z",label="",dpi=500)
    end
    close(file)
    plot!(title="P=$(P), z*=$(zstar)",dpi=500)
end
function String_3DPlot!(filename,P,zstar)
    file = h5open("results/"*filename,"r")

    sols=read(file["P$(P)_z$(zstar)/sols"])
    Lint=read(file["P$(P)_z$(zstar)/Lint"])
    for i in 1:1:length(Lint)
        plot!(interval(Int(size(sols,2)/4+0.5),Lint[i]),sols[i,1:Int(end/2)], sols[i,Int(end/2)+1:end], ylabel="r",xlabel="x",zlabel="z",label="",dpi=500)
    end
    close(file)
    plot!(title="P=$(P), z*=$(zstar)",dpi=500)
end
export String_3DPlot
export String_3DPlot!

function String_3DPlot_curve(filename,P,zstar)
    p=plot()
    file = h5open("results/"*filename,"r")

    sols=read(file["P$(P)_z$(zstar)/sols"])
    Lint=read(file["P$(P)_z$(zstar)/Lint"])
    for i in 1:1:length(Lint)
        p=plot!(BSplineManifold(vec2man(sols[i,:]),KV),size=(1500,750),xlabel='x',ylabel='r',zlabel='z')
    end
    close(file)
    plot!(title="P=$(P), z*=$(zstar)",dpi=500)
end
export String_3DPlot_curve
function String_3DPlot_curve!(filename,P,zstar)
    file = h5open("results/"*filename,"r")

    sols=read(file["P$(P)_z$(zstar)/sols"])
    Lint=read(file["P$(P)_z$(zstar)/Lint"])
    for i in 1:1:length(Lint)
        plot!(BSplineManifold(vec2man(sols[i,:]),KV),size=(1500,750),xlabel='x',ylabel='r',zlabel='z')
    end
    close(file)
    plot!(title="P=$(P), z*=$(zstar)",dpi=500)
end
export String_3DPlot_curve


function String_Plot_curve(filename,P,zstar,i,charx,chary)
    p=plot()
    file = h5open("results/"*filename,"r")

    sols=read(file["P$(P)_z$(zstar)/sols"])
    Lint=read(file["P$(P)_z$(zstar)/Lint"])

    M=BSplineManifold(vec2man(sols[i,:]),KV)

    ss=Array(range(0,pi,100))
    x=[M(i)[1] for i in ss]
    r=[M(i)[2] for i in ss]
    z=[M(i)[3] for i in ss]

    if charx=='x'
        xx=x
    elseif charx=='r'
        xx=r
    elseif charx=='z'
        xx=z
    elseif charx=='s'
        charx='σ'
        xx=ss
    else
        error("\nChoose a good label.")
        return
    end

    if chary=='x'
        yy=x
    elseif chary=='r'
        yy=r
    elseif chary=='z'
        yy=z
    elseif chary=='s'
        chary='σ'
        yy=ss
    else
        error("\nChoose a good label.")
        return
    end

    p=plot!(xx,yy,xlabel=charx,ylabel=chary)
    close(file)
    plot!(title="P=$(P), z*=$(zstar)",dpi=500)
end