"""
    quadrature_fg(f::Function,g::Function,g_inv::Function,τ_min::Cdouble,τ_max::Cdouble,Δr0::Cdouble;Nr_min::Int64=20)

    g(τ) and g^{-1}(r) are bijections over the interval Ω = [min(τ_min,τ_max),min(τ_min,τ_max)]

    f∘g is integrable and finite over Ω

    ∫_{τ_min}^{τ_max} f(g(τ)) dτ

"""
function quadrature_fg(f::Function,g::Function,g_inv::Function,τ_min::Cdouble,τ_max::Cdouble,Δr0::Cdouble;Nr_min::Int64=20)
    val = 0.0
    r_τ_min = g(τ_min)
    r_τ_max = g(τ_max)
    Nr = max(Nr_min,ceil(Int64,abs(r_τ_max-r_τ_min)/Δr0))
    r_disc = collect(LinRange(r_τ_min,r_τ_max,Nr));
    τ_disc = [g_inv(r) for r in r_disc]
    fs = [f(g(τ))  for τ in τ_disc]
    # trapezoid rule 
    for i in 1:(Nr-1)
        val = val + (1.0/2.0)*(fs[i]+fs[i+1])*(τ_disc[i+1]-τ_disc[i])
    end
    val 
end

#TODO: deal with the case lim_{τ->τ_min+} g(τ) = g(τ_min), and lim_{r->g(τ_min)} f(r) = ±∞. Idea ?