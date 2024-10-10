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
    xs=diff(x)./hh
    rr=interpolate((ss,), r, Gridded(Linear()))(s)
    rs=diff(r)./hh
    zz=interpolate((ss,), z, Gridded(Linear()))(s)
    zs=diff(z)./hh
    Lag=lagrangian(s,xx,xs,rr,rs,zz,zs)
    S=sum(hh .* Lag)
    return S
end
export action

function knotvec(N::Int,deg::Int)
    I=Array(range(0,π,length=N-deg+1))
    for i in 1:deg 
        I=cat([0],I,[π];dims=1) 
    end
    return I
end
export knotvec

function interval(N::Int,L)
    I=Vector{Float64}(undef, 2N-1)
    I=sqrt.(Array(range(0,(L/2)^2,length=N)))
    I=[-I[end:-1:1];I]
    filter!(e->hash(e)!=hash(-0.0),I)
    return I
end
export interval

function vec2man(c)
    cmat=transpose(reshape(c,3,Int(length(c)/3)))
    cm = [SizedVector(i, 54 , 2.3) for i in range(-1,1,Int(length(c)/3))]
    for i in 1:Int(length(c)/3)
        cm[i]=cmat[i,:]
    end
    return cm
end