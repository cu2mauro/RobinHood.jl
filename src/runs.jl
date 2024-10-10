function Run_Single_String(filename) #to be made, still old version
    #=rmax = 10e1 # is cutoff

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
    eqconst = [rmax, rmax, zstar, zstar]
    lbounds = [fill(rst,length(I));fill(zstar,length(I))]
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
    plot(I,sol.u[Int(end/2)+1:end], ylabel="z", lc=:red,xlabel="x", label=""),plot!(I,fill(z₀,length(I)), label="z₀", lc=:red, xlabel="x", ls=:dot),plot_title="L=$(L)") =#
end
export Run_Single_String

function Run_Multiple_Strings(filename,P_list,zstar_list)

    snapping=false

    file = h5open("results/"*filename,"w")
    println("\nData file named ",filename,".h5 was created.")

    rmax = 8e1 # is cutoff
    NN=50 #number of points ALONG the curve for integration
    ncp=7 #number of CONTROL POINTS of the curve
    Nstrings=20 #number of different strings to optimize
    deg=2 #degree of spline

    global KV=BSplineSpace{deg}(KnotVector(knotvec(ncp,deg))) #knot vector for control points
    ss=Array(range(0,pi,NN)) #range for integration
    SNG(c, ss) = action(c, ss, KV)
    cons(res, c, I) = (res .= [c[1],c[2],c[3],c[end-2],c[end-1],c[end]])

    Lint=(Array(range(sqrt(0.001),cbrt(0.2),length=Nstrings))).^2 #separations
    II=1:1:Nstrings #to iterate all the strings 
    Eint=similar(Lint) #energies
    sols=Array{Float64}(undef,Nstrings,3*ncp) #solutions

    c0 = [SizedVector(i, rmax/2 , 1) for i in range(-Lint[1]/3,Lint[1]/3,ncp)]
    cc=fill([0.0,0.0,0.0],length(c0))
    for i in eachindex(c0) cc[i]=c0[i] end
    ccc=vcat(cc...)

    println("\nGetting started with the loop now.")
    file["P_list"]=P_list[:]
    file["zstar_list"]=zstar_list[:]
    file["KV"]=KV.knotvector.vector
    for PP in P_list
        for zz in zstar_list
            global P=PP
            global zstar=zz
            ccc=vcat(cc...)
            create_group(file, "P$(P)_z$(zstar)")
            snap_flag=false #
            #Threads.@threads for i in II
            for i in II
                # initialization
                L = Lint[i]
                #ccc = try sols[i-1,:] catch temp ccc end #use latest results as starting poitn if available
                # optimization
                if snap_flag==false
                    eqconst = [-L/2,rmax,zstar,L/2, rmax, zstar]
                    lbounds=Float64[]; for i in 1:ncp lbounds=vcat(lbounds,[-L/2,rst,0]) end
                    ubounds=Float64[]; for i in 1:ncp ubounds=vcat(ubounds,[L/2,rmax,P]) end
                    optprob = OptimizationFunction(SNG, Optimization.AutoFiniteDiff(), cons = cons)
                    prob = OptimizationProblem(optprob, ccc, ss,; lcons = eqconst, ucons = eqconst, lb=lbounds, ub=ubounds)
                    sol = solve(prob, Optim.IPNewton(),g_tol=1e-12,x_tol=1e-8)
                    sols[i,:] = sol.u
                    Eint[i] = SNG(sol.u,ss)
                    zvals=[sols[i,j] for j in 3:3:3*ncp]
                    if any(abs.(zvals.-P/2).<0.15) && snapping==true
                        snap_flag=true
                    end
                end
                if snap_flag==true
                    z_ext=[zstar;fill(P/2,length(I)-2);zstar]
                    # second optimization
                    eqconst2 = [rmax, rmax]
                    lbounds2 = fill(rst,length(I))
                    ubounds2 = fill(rmax,length(I))
                    cons2(res, r, I) = (res .= [r[1], r[end]])
                    SNG2(r, I) = action([r;z_ext], I)
                    optprob2 = OptimizationFunction(SNG2, Optimization.AutoReverseDiff(true), cons = cons2)
                    prob2 = OptimizationProblem(optprob2, r0, I,; lcons = eqconst2, ucons = eqconst2, lb = lbounds2, ub = ubounds2)
                    sol2 = solve(prob2, IPNewton())
                    sols[i,:] = [sol2.u;z_ext]
                    Eint[i] = SNG2(sol2.u,I)
                end
            end
        file["P$(P)_z$(zstar)/Eint"]=Eint
        file["P$(P)_z$(zstar)/sols"]=sols
        file["P$(P)_z$(zstar)/Lint"]=Array(Lint)
        println("\nFinished iteration with P=$(P) and z*=$(zstar).")
        end
    end
    close(file)

    println("\nRun finished! Time for plotting!")
end
export Run_Multiple_Strings