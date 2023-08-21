"""
    determines the distance on the parametric line starting from (xm,ym,zm) in the direction (x0,y0,z0)

    τ_root_sphere(x0::Cdouble,y0::Cdouble,z0::Cdouble,xm::Cdouble,ym::Cdouble,zm::Cdouble,R2::Cdouble,rm2::Cdouble)

    input: 

      - xm, ym and zm: origin of the parametric line (τ=0)
      - x0, y0 and z0: coordinates of a point on the parametric line
      - R2: square radius of the intersecting sphere 
      - rm2: square distance between the center of the sphere and (xm,ym,zm), rm^2=xm^2+ym^2+zm^2  (R2-rm2⩾0)

    output: 

      - return the distance between the (xm,ym,zm) and the intersection point (sphere of radius R and straight line between (xm,ym,zm) and (x0,y0,z0))
"""
function τ_root_sphere(x0::Cdouble,y0::Cdouble,z0::Cdouble,xm::Cdouble,ym::Cdouble,zm::Cdouble,R2::Cdouble,rm2::Cdouble)
    # angular direction (xm,ym,zm) to (x0,y0,z0)
    τ_max = sqrt((x0-xm)^2 + (y0-ym)^2 + (z0-zm)^2)
    α = (x0-xm)/τ_max
    β = (y0-ym)/τ_max
    γ = (z0-zm)/τ_max
    # discriminant
    Δ = 4*((xm*α + ym*β+ zm*γ)^2 + (R2-rm2))
    # first root 
    τ1 = -(xm*α + ym*β + zm*γ + 0.5sqrt(Δ))
    # second root
    τ2 = -(xm*α + ym*β + zm*γ - 0.5sqrt(Δ))
    # projected square distances
    d1 = ((xm+α*τ1)-x0)^2 + ((ym+β*τ1)-y0)^2 + ((zm+γ*τ1)-z0)^2
    d2 = ((xm+α*τ2)-x0)^2 + ((ym+β*τ2)-y0)^2 + ((zm+γ*τ2)-z0)^2
    # select the root such that the distance to the point (x0,y0,z0) is shortest
    if (d1<d2)
        τ_bar = τ1
    else
        τ_bar = τ2
    end

    # return the distance and the direction cos 
    τ_bar,α,β,γ
end


"""
    quadrature evaluating the integral ∫_{τmin}^{τmax} f(g(τ)) dτ for a sphere of radius μ0

    with g(τ)=√(a*τ^2+b*τ+c) with a, b and c defined for each position in the sphere (depend on the location in the sphere and the direction of the parametric line)

    f(r) is a mass function depending depending on the distance from the center of the sphere 

    the integral boundaries are such that the parametric curve satifies τmin->(r[n],φ[j],θ[k]) and in τmax->(r(τmax),φ(τmax),θ(τmax)) with r(τmax)=μ0+κΔr 
    τmin is set to 0 without loss of generality and τmax is determined for each position in the sphere 

    quadrature_fg_sphere(r::Array{Cdouble,1},φ::Array{Cdouble,1},θ::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,Δr::Cdouble;κ::Cdouble=5.0,Nτ::Int64=20)

    input:

      - r,φ,θ: spherical coordinates of the discretization of the sphere (polar axis Oz)
      - x0,y0,z0: cartesian coordinate of a point belonging to the parametric straight line of integration
      - μ0: radius of the sphere 
      - Δr: transition distance of the mass function f
    
    optional input:

      - κ:  maximum radius for the integration is r=μ0+κΔr
      - Nτ: number of discretization nodes 

    output:

      - FGQ_rφθ: 3D array containing the quadrature for each discretization nodes (r,φ,θ)
"""
function quadrature_fg_sphere(r::Array{Cdouble,1},φ::Array{Cdouble,1},θ::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,Δr::Cdouble;κ::Cdouble=5.0,Nτ::Int64=20)
    # square radius of the sphere 
    R2 = (μ0 + κ*Δr)^2

    # compute quadrature for each discretization nodes (r,φ,θ)
    FGQ_rφθ = zeros(Cdouble,length(r),length(φ),length(θ));
    for n in eachindex(r)
        rm2 = r[n]^2
        if (R2>=rm2)
            for j in eachindex(φ)
                zm   = r[n]*cos(φ[j])
                rxym = r[n]*sin(φ[j])
                for k in eachindex(θ)
                    # cartesian coordinates
                    xm = rxym*cos(θ[k])
                    ym = rxym*sin(θ[k])

                    # quadrature ∫_0^{τ_root} f(g(τ)) dτ
                    τ_root,α,β,γ = τ_root_sphere(x0,y0,z0,xm,ym,zm,R2,rm2)
                    # this is not optimal, the details are at the transition, so, the discretization of τ should reflect that
                    τ_disc = collect(LinRange(0.0,τ_root,Nτ)) 
                    r_param = sqrt.((xm.+α*τ_disc).^2 .+ (ym.+β*τ_disc).^2 .+ (zm.+γ*τ_disc).^2)
                    # a = α^2 + β^2 + γ^2
                    # b = 2.0*(α*xm + β*ym + γ*zm)
                    # c = xm^2 + ym^2 + zm^2;
                    # r_param = g_sqrtquad.(τ_disc,a,b,c)
                    ρ_disc = f_logistic(r_param,μ0,Δr;A=1.0);
                    # quadrature f(g(τ))
                    FGQ_rφθ[n,j,k] = quadrature_fg(ρ_disc,τ_disc)
                end
            end
        end
    end

    # return the quadrature for each discretization nodes of the sphere 
    FGQ_rφθ
end

"""
    quadrature_fg_sphere_g_opt_f(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,Δr::Cdouble;κ::Cdouble=5.0,Nτ::Int64=20)

    quadrature evaluating the integral ∫_{τmin}^{τmax} f(g(τ)) dτ for a sphere of radius μ0 

    Contrary to quadrature_fg_sphere the quadrature is adjusted to f instead of g

    see quadrature_fg_sphere and quadrature_fg(f::Function,τ_min::Cdouble,τ_max::Cdouble,a::Cdouble,b::Cdouble,c::Cdouble,Δr0::Cdouble;Nr_min::Int64=20)
"""
function quadrature_fg_sphere_g_opt_f(r::Array{Cdouble,1},φ::Array{Cdouble,1},θ::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,Δr::Cdouble;κ::Cdouble=5.0,Nτ::Int64=20)
    # square radius of the sphere 
    R2 = (μ0 + κ*Δr)^2

    # compute quadrature for each discretization nodes (r,φ,θ)
    FGQ_rφθ = zeros(Cdouble,length(r),length(φ),length(θ));
    fint = (x::Cdouble->f_logistic(x,μ0,Δr;A=1.0))
    for n in eachindex(r)
        rm2 = r[n]^2
        if (R2>=rm2)
            for j in eachindex(φ)
                zm   = r[n]*cos(φ[j])
                rxym = r[n]*sin(φ[j])
                for k in eachindex(θ)
                    # cartesian coordinates
                    xm = rxym*cos(θ[k])
                    ym = rxym*sin(θ[k])

                    # quadrature ∫_0^{τ_root} f(g(τ)) dτ
                    τ_root,α,β,γ = τ_root_sphere(x0,y0,z0,xm,ym,zm,R2,rm2) # if b^2≃4*a*c, then the center of the cylinder is on the parametric line 
                    a = α^2 + β^2 + γ^2
                    b = 2.0*(α*xm + β*ym + γ*zm)
                    c = xm^2 + ym^2 + zm^2;
                    FGQ_rφθ[n,j,k] = quadrature_fg(fint,0.0,τ_root,a,b,c,Δr/5.0;Nr_min=Nτ) 
                end
            end
        end
    end

    # return the quadrature for each discretization nodes of the sphere 
    FGQ_rφθ
end