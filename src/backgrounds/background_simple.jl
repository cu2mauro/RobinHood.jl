#constants
μ=0 #for susy
q=1 #can be changed
N=5

#interval_functions
ht(t) = 0.5 * (sign(t) + 1)
iv(t, a, b) = ht(t-a) - ht(t-b)

function Quiver_choice(quiverlabel)

    if quiverlabel == 1
        #EX1 IN PAPER ◯-◯-◯-◯-◯-☐
        @eval α(z,P) = @. -81 * π^2 * N/6 * (((1-P^2)*z + z^3)*iv(z,0,P-1)+((2P^2-3P+1)*(z-P) + (P-1)*(P-z)^3)*iv(z,P-1,P)) 
        @eval α2(z,P) = @. -81 * π^2 * N/6 * (6z*iv(z,0,P-1)+(6*(P-1)*(P-z))*iv(z,P-1,P)) 
        @eval α(z,P) = @. -81 * π^2 * N/6 * (((1-P^2)*z + z^3)*iv(z,0,P-1)+((2P^2-3P+1)*(z-P) + (P-1)*(P-z)^3)*iv(z,P-1,P)) 
        @eval α2(z,P) = @. -81 * π^2 * N/6 * (6z*iv(z,0,P-1)+(6*(P-1)*(P-z))*iv(z,P-1,P)) 

    elseif quiverlabel == 2
        #EX2 IN PAPER ◯-◯-◯-☐-◯-◯-◯
        @eval α(z,P) = @. -81 * π^2 * N * (((-P^2)/8*z + z^3/6)*iv(z,0,Int(P/2))+((-P^2)/8*(P-z) + 1/6*(P-z)^3)*iv(z,Int(P/2),P));
        @eval α2(z,P) = @. -81 * π^2 * N * (z*iv(z,0,Int(P/2))+(P-z)*iv(z,Int(P/2),P));
        @eval α(z,P) = @. -81 * π^2 * N * (((-P^2)/8*z + z^3/6)*iv(z,0,Int(P/2))+((-P^2)/8*(P-z) + 1/6*(P-z)^3)*iv(z,Int(P/2),P));
        @eval α2(z,P) = @. -81 * π^2 * N * (z*iv(z,0,Int(P/2))+(P-z)*iv(z,Int(P/2),P));

    elseif quiverlabel == 3
        #EX3 IN PAPER ☐-◯-◯-◯-◯-◯-☐ 
        @eval α(z,P) = @. -81 * π^2 * N * (((1-P)/2*z + z^3/6) * iv(z,0,1)+(1/6 - P/2*z + 1/2*z^2) * iv(z,1,P-1)+((1-P)/2*(P-z) + 1/6*(P-z)^3) * iv(z,P-1,P)) 
        @eval α2(z,P) = @. -81 * π^2 * N * (z * iv(z,0,1)+1  * iv(z,1,P-1)+(P-z) * iv(z,P-1,P))
        @eval α(z,P) = @. -81 * π^2 * N * (((1-P)/2*z + z^3/6) * iv(z,0,1)+(1/6 - P/2*z + 1/2*z^2) * iv(z,1,P-1)+((1-P)/2*(P-z) + 1/6*(P-z)^3) * iv(z,P-1,P)) 
        @eval α2(z,P) = @. -81 * π^2 * N * (z * iv(z,0,1)+1  * iv(z,1,P-1)+(P-z) * iv(z,P-1,P))

    else
        error("Please choose a good quiver label.")
        return
    end

    
    #ex not in paper ◯-◯-◯-◯-◯-◯-◯-◯
    #α(z) = @. -27/2 * π^2 * N * (z^3 - P*(P+2)*z) 
    #α2(z) = @. -81 * π^2 * N * z  

    #ex not in paper ◯-◯-☐-◯-◯-◯-☐-◯-◯ 
    #global kk=4
    #global qq=5
    #a1 = kk/P * (-2*P^2-2*P*qq+qq^2+kk*P+2*kk*qq)
    #a2 = kk/P * (-2*kk*P-2*P^2+2*kk*qq-2*P*qq+qq^2)
    #a3 = -kk/P *(P^2+P*qq+qq^2+kk*P+2*kk*qq)
    #a4 = kk^3
    #α(z) = @. -81 * π^2 * N/6 * ((a1*z+z^3) * iv(z,0,kk)+(a4+a2*z+3*kk*z^2) * iv(z,kk,kk+qq)+(a3*(P-z)+kk/(P-kk-qq)*(P-z)^3) * iv(z,kk+qq,P)) 
    #α2(z) = @. -81 * π^2 * N * (z * iv(z,0,kk)+ kk * iv(z,kk,kk+qq)+kk/(P-kk-qq)*(P-z) * iv(z,kk+qq,P))

    #action_functions
    @eval f(r) = @. (1 - μ/r^4 - q^2/r^6)
    @eval F2(r,z,P) = @. - (α(z,P)/α2(z,P)) * r^4
    @eval G2(r,z,P) = @. - α(z,P) / (α2(z,P)*f(r))
    @eval S2(r,z) = @. r^2 / 6
    @eval F(r,z,P) = @. sqrt(F2(r,z,P))
    @eval G(r,z,P) = @. sqrt(G2(r,z,P))
    @eval S(r,z) = @. sqrt(S2(r,z))
    return
end
export Quiver_choice