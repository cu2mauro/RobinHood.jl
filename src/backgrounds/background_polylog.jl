using Polylogarithms
using QuadGK
using Richardson, Plots

# constants
global const c=1
global const μ=1
global const gg=sqrt(9/2)
global const N=1 #8?
global const P=10
global const sst=eps()
global const rol=1.11070395
global const etast=1

# functions
H(r) = @. 1- (c^2/r^3)
f(r) = @. -μ/r^3 + 2/9 * gg^2 * r^2 * H(r)^2
X(r) = @. H(r)^(-1/4)
V(s,eta) = @. Vt(s,eta)/s #never used alone
Vt(s,eta) = @. N*P^3/(2*pi^3) * real(polylog(3,-ℯ^(-pi/P*(abs(s)+im+im*eta)))-polylog(3,-ℯ^(-pi/P*(abs(s)-im+im*eta)))) # factor maybe N*P^2/(pi^2)
f1(r,s,eta) = @. (3*pi)/(2*X(r)^2) * sqrt(s^2 + (3X(r)^4*s*dsV(s,eta))/de2V(s,eta))
f3(r,s,eta) = @. (X(r)^2 * de2V(s,eta))/(3*s*dsV(s,eta))
F2(r,s,eta) = @. (2/9 * gg^2)^2 *f1(r,s,eta)^2 *H(r)*r^4
G2(r,s,eta) = @. (2/9 * gg^2)^2 *f1(r,s,eta)^2 *H(r)* r^2/f(r)
S2(r,s,eta) = @. (2/9 * gg^2) *f1(r,s,eta)^2 *f3(r,s,eta)*sqrt(H(r))* r^2
F(r,s,eta) = @. sqrt(F2(r,s,eta))
G(r,s,eta) = @. sqrt(G2(r,s,eta))
S(r,s,eta) = @. sqrt(S2(r,s,eta))

de2Vt(s,eta) = @. N*P/(2*pi) * real(polylog(1,-ℯ^(-pi/P*(abs(s)-im+im*eta)))-polylog(1,-ℯ^(-pi/P*(abs(s)+im+im*eta))))
dsVt(s,eta) = @. N*P^2/(2*pi^2) * real(polylog(2,-ℯ^(-pi/P*(abs(s)-im+im*eta)))-polylog(2,-ℯ^(-pi/P*(abs(s)+im+im*eta))))
de2V(s,eta) = @. de2Vt(s,eta)/s
dsV(s,eta) = @. (s*dsVt(s,eta) - Vt(s,eta))/ s^2

# special case 3

Veff(r,s,eta,e0)= @. F(r,s,eta)/(F(r,s,e0)*S(r,s,eta)) * sqrt(F2(r,s,eta)-F2(r,s,e0))

Lqq(r,s,e0)= @. 2 * first(quadgk(eta -> 1/Veff(r,s,eta,e0),etast,e0,rtol=1e-3))
Lqqb_limzero(e0)=first(Richardson.extrapolate(1.0, rtol=1e-5) do x                                                                                                                            
    @. Lqq(rol,x,e0) /(27/4 * pi^2 * H(rol) * rol^4)                                                 
end)
Lqqb(e0) = @. Lqq(rol,sst,e0) *sqrt((27/4 * pi^2 * H(rol) * rol^4)) / S(rol,sst,4.2)

Eqq(r,s,e0)=@. F(r,s,e0) * Lqq(r,s,e0) + 2 * first(quadgk(eta -> S(r,s,eta)/F(r,s,eta) * sqrt(F2(r,s,eta)-F2(r,s,e0)),etast,e0,rtol=1e-3))
Eqqb_limzero(e0)=first(Richardson.extrapolate(10.0, rtol=1e-5) do x                                                                                                                            
    @. Eqq(rol,x,e0) / S(rol,x,e0)                                                 
end)
Eqqb(e0) = @. Eqq(rol,sst,e0) / S(rol,sst,e0)

#= # values for plots
etarange=1.02:0.2:9.0
Lbrange=Lqqb(etarange)
Ebrange=Eqqb(etarange) =#