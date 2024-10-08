function lagrangian(s,x,xs,r,rs,z,zs)
    L = @. NaNMath.sqrt(F2(r,z) * xs^2 + G2(r,z) * rs^2 + S2(r,z) * zs^2)                                          
	return L
end
export lagrangian

function action(c,ss,KV)
    cmat=transpose(reshape(c,3,Int(length(c)/3)))
    cm=Array{Vector{Float64}}(undef,Int(length(c)/3))
    for i in 1:Int(length(c)/3)
        cm[i]=cmat[i,:]
    end
    M=BSplineManifold(cm,KV)
    hh=[ss[2:1:end];0]-ss
    pop!(hh)
    x=[M(i)[1] for i in ss]
    r=[M(i)[2] for i in ss]
    z=[M(i)[3] for i in ss]
    s=ss[1:length(hh)]+hh./2
    xx=interpolate((ss,), x, Gridded(Linear()))(s)
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
    I=cat(I;dims=1)
    return I
end
export interval

function vec2man(c)
    cmat=transpose(reshape(c,3,Int(length(c)/3)))
    cm = [SizedVector(i, 54 , 2.3) for i in range(-L/2,L/2,Int(length(c)/3))]
    for i in 1:Int(length(c)/3)
        cm[i]=cmat[i,:]
    end
    return cm
end