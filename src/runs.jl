using OptimizationOptimJL, ReverseDiff
using Interpolations
using Base.Threads
using HDF5

include("backgrounds/background_simple.jl")
include("plots.jl")
include("actions.jl")

function Run_Single_String(filename) #to be made, still old version
    rmax = 10e1 # is cutoff

    # physical parameters
    L = 0.13

    # initialization
    I = interval(ll,L)

    r0 = 10 .* ones(length(I))
    z0 = 8 .* ones(length(I))
    c0 = [r0;z0]

    # action
    SNG(c, I) = action(c, I)

    # optimization
    eqconst = [rmax, rmax, etast, etast]
    lbounds = [fill(rst,length(I));fill(etast,length(I))]
    ubounds = [fill(rmax,length(I));fill(P,length(I))]
    cons(res, c, I) = (res .= [c[1], c[Int(end/2)], c[Int(end/2)+1], c[end]])
    optprob = OptimizationFunction(SNG, Optimization.AutoReverseDiff(true), cons = cons)
    prob = OptimizationProblem(optprob, c0, I,; lcons = eqconst, ucons = eqconst, lb = lbounds, ub = ubounds)
    sol = solve(prob, IPNewton());

    # plots and comparisons
    E=SNG(sol.u,I)
    r₀=interpolate((I,), sol.u[1:Int(end/2)], Gridded(Linear()))(0)
    z₀=interpolate((I,), sol.u[Int(end/2)+1:end], Gridded(Linear()))(0)
    println("The value of r₀ is ", r₀,".")
    plot(plot(I,sol.u[1:Int(end/2)], ylabel="r", lc=:blue, xlabel="x", label=""),plot!(I,fill(r₀,length(I)), label="r₀", lc=:blue, xlabel="x", ls=:dot),
    plot(I,sol.u[Int(end/2)+1:end], ylabel="z", lc=:red,xlabel="x", label=""),plot!(I,fill(z₀,length(I)), label="z₀", lc=:red, xlabel="x", ls=:dot),plot_title="L=$(L)")
end
export Run_Single_String

function Run_Multiple_Strings(filename)
    file = h5open("results/"*filename,"w")
    println("\nData file named ",filename,".h5 was created.")
    rmax = 8e1 # is cutoff
    SNG(c, I) = action(c, I)
    cons(res, c, I) = (res .= [c[1], c[Int(end/2)], c[Int(end/2)+1], c[end]])

    #preinitialization
    ll=60 #only even to make it work faster
    Lint=(Array(range(cbrt(0.01),cbrt(0.10),length=40))).^3

    II=1:1:length(Lint)
    Eint=similar(Lint)
    I = interval(ll,1)
    sols=Array{Float64}(undef,length(Lint),2*length(I))

    println("\nGetting started with the loop now.")
    file["P_list"]=P_list[:]
    file["etast_list"]=etast_list[:]
    for PP in P_list
        for ee in etast_list
            global P=PP
            global etast=ee
            create_group(file, "P$(P)_eta$(etast)")
            #Threads.@threads for i in II
            for i in II
                # initialization
                L = Lint[i]
                I = interval(ll,L)
                r0 = rmax/10 .* ones(length(I))
                z0 = etast .* ones(length(I)) #+ 0.1 .* ones(length(I)) - rand(length(I))./10
                c0 = [r0;z0]
                # optimization
                eqconst = [rmax, rmax, etast, etast]
                lbounds = [fill(rst,length(I));fill(0,length(I))]
                ubounds = [fill(rmax,length(I));fill(P,length(I))]
                optprob = OptimizationFunction(SNG, Optimization.AutoReverseDiff(true), cons = cons)
                prob = OptimizationProblem(optprob, c0, I,; lcons = eqconst, ucons = eqconst, lb = lbounds, ub = ubounds)
                sol = solve(prob, IPNewton(),g_tol=1e-12,x_tol=1e-4)
                sols[i,:] = sol.u
                Eint[i] = SNG(sol.u,I)
                #=  # STRING TO BRANE
                    eta_ext=[etast;fill(P-1,length(I)-2);etast]
                    # second optimization
                    eqconst2 = [rmax, rmax]
                    lbounds2 = fill(rst,length(I))
                    ubounds2 = fill(rmax,length(I))
                    cons2(res, r, I) = (res .= [r[1], r[end]])
                    SNG2(r, I) = action([r;eta_ext], I)
                    optprob2 = OptimizationFunction(SNG2, Optimization.AutoReverseDiff(true), cons = cons2)
                    prob2 = OptimizationProblem(optprob2, r0, I,; lcons = eqconst2, ucons = eqconst2, lb = lbounds2, ub = ubounds2)
                    sol2 = solve(prob2, IPNewton());
                    Eint2[i] = SNG2(sol2.u,I) =#
            end
        file["P$(P)_eta$(etast)/Eint"]=Eint
        file["P$(P)_eta$(etast)/sols"]=sols
        file["P$(P)_eta$(etast)/Lint"]=Array(Lint)
        println("\nFinished iteration with P=$(P) and z*=$(etast).")
        end
    end
    close(file)

    println("\nRun finished! Time for plotting!")
end
export Run_Multiple_Strings