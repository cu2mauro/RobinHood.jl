function lagrangian(s,x,xs,r,rs,z,zs)
    L = @. NaNMath.sqrt(F2(r,z) * xs^2 + G2(r,z) * rs^2 + S2(r,z) * zs^2)                                          
	return L
end
export lagrangian

function action(c,ss,KV)
    #now c is control points
    cmat=reshape(c,(3,length(c)/3))
    cm=Array{MVector{3,Float64}}(undef,length(c))
    for i in 1:length(c)
        cm[i]=MVector{3}(cmat[:,i])
    end
    M=BSplineManifold(cm,KV)
    hh=[ss[2:1:end];0]-ss
    pop!(hh)
    x=[M(i)[1] for i in ss]
    r=[M(i)[2] for i in ss]
    z=[M(i)[3] for i in ss]
    s=ss[1:length(hh)]+hh./2
    xx=interpolate((ss,), r, Gridded(Linear()))(s)
    xs=diff(r)./hh
    rr=interpolate((ss,), r, Gridded(Linear()))(s)
    rs=diff(r)./hh
    zz=interpolate((ss,), z, Gridded(Linear()))(s)
    zs=diff(z)./hh
    Lag=lagrangian(s,xx,xs,rr,rs,zz,zs)
    S=sum(hh .* Lag)
    return S
end
export action

function interval(N::Int)
    I=Array(range(0,π,length=N))
    I=cat([0,0,0],I,[π,π,π];dims=1)
    return I
end
export interval