
"""
    logistic(x::Cdouble,thX::Cdouble=33.275)

    returns the value of the logistic function at x, 0 if x<-thX, 1 if x>thX
        1/(1+e^{-x})

    logistic(x::Cdouble,a::Cdouble,b::Cdouble,α::Cdouble=1.0,thX::Cdouble=33.275)

    returns the value of the logistic function at x, 0 if αx<-thX, 1 if αx>thX
        a + (b-a)/(1+(1/α)*e^{-αx})
"""
function logistic(x::Cdouble,thX::Cdouble=33.275)
    if x<-thX
        y = 0.0
    elseif x>thX
        y = 1.0
    else
        y = 1.0/(1.0+exp(-x))
    end
    y
end
function logistic(x_::Cdouble,a_::Cdouble,b_::Cdouble,alpha_::Cdouble=1.0,thX::Cdouble=33.275)
    if (alpha_*x_)<-thX
        y = a_
    elseif (alpha_*x_)>thX
        y = b_
    else
        y = a_ + (b_-a_)/(1.0+(1.0/alpha_)*exp(-alpha_*x_))
    end
    y
end

"""
    function used for tests of the integral of composition of functions

        f_logistic(r::Array{Cdouble,1},r0::Cdouble,Δr::Cdouble;A::Cdouble=1.0)

    input:

      - r: radial distance
      - r0: critical distance 
      - Δr: transition distance
      - A: amplitude (default at 1)

    output: 

      - ρ0/(1.0+exp((r-μ0)/Δr))
"""
function f_logistic(r::Array{Cdouble,1},r0::Cdouble,Δr::Cdouble;A::Cdouble=1.0)
    A*logistic.(-(r.-r0)/Δr)
end

"""
    g_sqrtquad(τ::Cdouble,a::Cdouble,b::Cdouble,c::Cdouble)

    inner function of the composition in the integral ∫_τ1^τ2 f(g(τ)) dτ

    Inspired by the distance from a central point to a point on a parametric line 

    computes √(a*τ^2+b*τ+c) and should be called with a>0, and b^2<4ac so that a*τ^2+b*τ+c⩾0
"""
function g_sqrtquad(τ::Cdouble,a::Cdouble,b::Cdouble,c::Cdouble)
    sqrt(a*τ^2+b*τ+c)
end

"""
    g_sqrtquad_inv_pos(r::Cdouble,a::Cdouble,b::Cdouble,c::Cdouble)

    increasing inverse of sqrt(a*τ^2+b*τ+c), i.e. τ⩾-b/(2a)

    should be called with a>0 and b^2<4ac

    returns √((r^2 + (b^2/(4a)) - c)/a) - b/(2a)
"""
function g_sqrtquad_inv_pos(r::Cdouble,a::Cdouble,b::Cdouble,c::Cdouble)
    if (a<1.0e-14) # isapprox(a,0.0,atol=1.0e-14))
        throw("MINOTAUR: g_sqrtquad_inv_pos recieved a degenerate quadratic, cannot proceed (a=0)")
    end
    if ((b^2)<(4ac)) 
        throw("MINOTAUR: g_sqrtquad_inv_pos recieved a degenerate case, cannot proceed (Δ>0)")
    end
    sqrt((1.0/a)*(r^2 + ((b^2)/(4a)) - c)) - (b/(2a))
end

"""
    g_sqrtquad_inv_neg(r::Cdouble,a::Cdouble,b::Cdouble,c::Cdouble)

    decreasing inverse of sqrt(a*τ^2+b*τ+c), i.e. τ<-b/(2a)

    should be called with a>0 and b^2<4ac

    returns -√((r^2 + (b^2/(4a)) - c)/a) - b/(2a)
"""
function g_sqrtquad_inv_neg(r::Cdouble,a::Cdouble,b::Cdouble,c::Cdouble)
    if (a<1.0e-14) # isapprox(a,0.0,atol=1.0e-14))
        throw("MINOTAUR: g_sqrtquad_inv_neg recieved a degenerate quadratic, cannot proceed (a=0)")
    end
    if ((b^2)<(4ac)) 
        throw("MINOTAUR: g_sqrtquad_inv_neg recieved a degenerate case, cannot proceed (Δ>0)")
    end
    -sqrt((1.0/a)*(r^2 + ((b^2)/(4a)) - c)) - (b/(2a))
end


"""
    interval(τ_min_0::Cdouble,τ_max_0::Cdouble,a::Cdouble,b::Cdouble,c::Cdouble)

    Let τ_min=min(τ_min_0,τ_max_0) ⩽ τ_max=max(τ_min_0,τ_max_0)

    This function determines the interval boundaries corresponding to g([τ_min,τ_max]) with g(τ) = √(a*τ^2+b*τ+c) (where a>0 and b^2<4ac so that a*τ^2+b*τ+c⩾0)

    The determination of the boundaries depends on the values τ_min and τ_max relative to τ0 = -b/(2a) (corrseponding to the argmin of g)

    if τ_min and τ_max ⩾ τ0
        g is increasing, so rmin = g(τ_min), and rmax = g(τ_max). rmid = Inf to indicate the lack of intermediate point
    if τ_min and τ_max ⩽ τ0
        g is decreasing, so rmin = g(τ_max), and rmax = g(τ_min). rmid = Inf to indicate the lack of intermediate point
    if τ_min⩽τ0⩽τ_max
        g is decreasing over [τ_min,τ0], and increasing over [τ0,τ_max], rmid = g(τ0), rmin = g(τ_min), and rmax = g(τ_max)

"""
function interval(τ_min_0::Cdouble,τ_max_0::Cdouble,a::Cdouble,b::Cdouble,c::Cdouble)
    if (τ_max_0<=τ_min_0)
        τ_min = τ_max_0
        τ_max = τ_min_0
    else
        τ_min = τ_min_0
        τ_max = τ_max_0
    end
    rlow = 0.0
    rmid = Inf
    rup  = 0.0
    τ0 = -b/(2a)
    if ((τ_min>=τ0) & (τ_max>=τ0)) # g_sqrtquad_inv_pos
        rlow = g_sqrtquad(τ_min,a,b,c)
        rup = g_sqrtquad(τ_max,a,b,c)
    elseif ((τ_min<τ0) & (τ_max<τ0)) # g_sqrtquad_inv_neg
        rlow = g_sqrtquad(τ_max,a,b,c)
        rup = g_sqrtquad(τ_min,a,b,c)
    else # if ((τ_min<τ0) & (τ_max>=τ0)) # g_sqrtquad_inv_neg then g_sqrtquad_inv_neg
        rmid = g_sqrtquad(τ0,a,b,c)
        rlow = g_sqrtquad(τ_min,a,b,c) # g_sqrtquad_inv_neg
        rup = g_sqrtquad(τ_max,a,b,c)  # g_sqrtquad_inv_neg
    # else # should never happen
    #     rmid = g_sqrtquad(τ0,a,b,c)
    #     rlow = g_sqrtquad(τ_max,a,b,c)
    #     rup = g_sqrtquad(τ_min,a,b,c)
    end
    rlow,rmid,rup # if isinf(rmid), integrate over [rlow,rmax], otherwise, integrate of [rmid,rlow] and [rmid,rup]
end


"""
    quadrature_fg(f::Function,τ_min::Cdouble,τ_max::Cdouble,a::Cdouble,b::Cdouble,c::Cdouble,Δr0::Cdouble;Nr_min::Int64=20)

    computes the trapezoid quadrature for the integral

        ∫_{τ_min}^{τ_max} f(g(τ)) dτ

    with τ_min<τ_max, and g(τ) = √(a*τ^2+b*τ+c) (where a>0 and b^2<4ac so that a*τ^2+b*τ+c⩾0)

    The minumim resolution for the variable r for f is Δr0, and the discretization of r should have at least Nr_min nodes.





    quadrature_fg(fgs::Array{Cdouble,1},τs::Array{Cdouble,1})

    computes the quadrature ∑_n f(g(τ_n)) ∫ e_n(τ) dτ ≃ ∫_{τ_min}^{τ_max} f(g(τ)) dτ

    with fgs[n] = f(g(τs[n])) and e_n piecewise linear basis functions 

"""
function quadrature_fg(f::Function,τ_min::Cdouble,τ_max::Cdouble,a::Cdouble,b::Cdouble,c::Cdouble,Δr0::Cdouble;Nr_min::Int64=20)
    val = 0.0
    if (!isapprox(τ_min,τ_max,atol=1.0e-14))
        rlow,rmid,rup = interval(τ_min,τ_max,a,b,c)
        if (!isinf(rmid))
            Nr = max(Nr_min,ceil(Int64,(rup-rlow)/Δr0))
            r_disc = collect(LinRange(rlow,rup,Nr));
            if (τ_min>=τ0)
                τ_disc = [g_sqrtquad_inv_pos(r,a,b,c) for r in r_disc]
            else
                τ_disc = [g_sqrtquad_inv_neg(r,a,b,c) for r in r_disc[Nr:-1:1]]
            end
            fs = [f(τ) for τ in τ_disc]

            # trapezoid rule 
            # sum(0.5*[τ_disc[2]-τ_disc[1]; τ_disc[3:end]-τ_disc[1:end-2]; τ_disc[end]-τ_disc[end-1]].*fs)
            for i in 1:(Nr-1)
                val = val + (1.0/2.0)*(fs[i]+fs[i+1])*(τ_disc[i+1]-τ_disc[i])
            end
        else
            # decreasing g
            Nr1 = max(ceil(Int64,Nr_min/2),ceil(Int64,(rlow-rmid)/Δr0))
            r_disc_1 = collect(LinRange(rmid,rlow,Nr1));
            τ_disc_1 = [g_sqrtquad_inv_neg(r,a,b,c) for r in r_disc_1[Nr1:-1:1]]
            fs_1 = [f(τ) for τ in τ_disc_1]

            # increasing g 
            Nr2 = max(ceil(Int64,Nr_min/2),ceil(Int64,(rup-rmid)/Δr0))
            r_disc_2 = collect(LinRange(rmid,rlup,Nr2));
            τ_disc_2 = [g_sqrtquad_inv_pos(r,a,b,c) for r in r_disc_2]
            fs_2 = [f(τ) for τ in τ_disc_2]

            # trapezoid rule 
            # sum(0.5*[τ_disc[2]-τ_disc[1]; τ_disc[3:end]-τ_disc[1:end-2]; τ_disc[end]-τ_disc[end-1]].*fs)
            for i in 1:(Nr1-1)
                val = val + (1.0/2.0)*(fs_1[i]+fs_1[i+1])*(τ_disc_1[i+1]-τ_disc_1[i])
            end
            for i in 1:(Nr2-1)
                val = val + (1.0/2.0)*(fs_2[i]+fs_2[i+1])*(τ_disc_2[i+1]-τ_disc_2[i])
            end
        end
    end
    val 
end
function quadrature_fg(fgs::Array{Cdouble,1},τs::Array{Cdouble,1})
    sum(0.5*[τs[2]-τs[1]; τs[3:end]-τs[1:end-2]; τs[end]-τs[end-1]].*fgs)
end






