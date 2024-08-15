using Interpolations
using NaNMath

include("backgrounds/background_simple.jl")

function lagrangian(x,r,rx,z,zx)
    L = @. NaNMath.sqrt(F2(r,z) + G2(r,z) * rx^2 + S2(r,z) * zx^2)                                          
	return L
end
export lagrangian

function action(c,I)
    hh=[I[2:1:end];0]-I
    pop!(hh)
    r=c[1:Int(length(c)/2)]
    z=c[Int(length(c)/2+1):end]
    bx=I[1:length(hh)]+hh./2
    rs=interpolate((I,), r, Gridded(Linear()))(bx)
    rx=diff(r)./hh
    zs=interpolate((I,), z, Gridded(Linear()))(bx)
    zx=diff(z)./hh
    Lag=lagrangian(bx,rs,rx,zs,zx)
    S=sum(hh .* Lag)
    return S
end
export action

function interval(N::Int,L)
    I=Vector{Float64}(undef, 2N-1)
    I=sqrt.(Array(range(0,(L/2)^2,length=N)))
    I=[-I[end:-1:1];I]
    filter!(e->hash(e)!=hash(-0.0),I)
    return I
end
export interval