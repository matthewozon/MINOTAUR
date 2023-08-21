"""
    determines the distance on the parametric line starting from (xm,ym,zm) in the direction (x0,y0,z0) intersecting with the cylinder centered on the Oy axis and of radius R

    τ_root_cylinder(x0::Cdouble,y0::Cdouble,z0::Cdouble,xm::Cdouble,ym::Cdouble,zm::Cdouble,R2::Cdouble,rm2::Cdouble)

    input: 

      - xm, ym and zm: origin of the parametric line (τ=0)
      - x0, y0 and z0: coordinates of a point on the parametric line
      - R2: square radius of the intersecting cylinder 
      - rm2: square distance between the center of the cylinder and (xm,ym,zm), rm^2=xm^2+ym^2  (R2-rm2⩾0)

    output: 

      - return the distance between the (xm,ym,zm) and the intersection point (cylinder of radius R along the Oy axis and straight line between (xm,ym,zm) and (x0,y0,z0))
"""
function τ_root_cylinder(x0::Cdouble,y0::Cdouble,z0::Cdouble,xm::Cdouble,ym::Cdouble,zm::Cdouble,R2::Cdouble,rm2::Cdouble)
    # angular direction (xm,ym,zm) to (x0,y0,z0)
    τ_max = sqrt((x0-xm)^2 + (y0-ym)^2 + (z0-zm)^2)
    α = (x0-xm)/τ_max
    γ = (z0-zm)/τ_max
    # discriminant
    Δ = 4*((zm*γ + xm*α)^2 + (γ^2+α^2)*(R2-rm2))
    # first root 
    τ1 = -(zm*γ + xm*α + 0.5sqrt(Δ))/(γ^2+α^2)
    # second root
    τ2 = -(zm*γ + xm*α - 0.5sqrt(Δ))/(γ^2+α^2)
    # projected square distances
    d1 = ((zm+γ*τ1)-z0)^2 + ((xm+α*τ1)-x0)^2
    d2 = ((zm+γ*τ2)-z0)^2 + ((xm+α*τ2)-x0)^2
    # select the root such that the distance to the apperture is shortest
    if (d1<d2)
        τ_bar = τ1
    else
        τ_bar = τ2
    end

    τ_bar,α,γ
end


"""
    quadrature evaluating the integral ∫_{τmin}^{τmax} f(g(τ)) dτ for a cylinder of radius μ0 (along the axis Oy)


    with g(τ)=√(a*τ^2+b*τ+c) with a, b and c defined for each position in the cylinder (depend on the location in the cylinder and the direction of the parametric line)

    f(r) is a mass function depending depending on the distance from the center of the cylinder 

    the integral boundaries are such that the parametric curve satifies τmin->(r[n],θ[j],y[k]) and in τmax->(r(τmax),θ(τmax),y(τmax)) with r(τmax)=μ0+κΔr 
    τmin is set to 0 without loss of generality and τmax is determined for each position in the cylinder 


    quadrature_fg_cylinder(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,Δr::Cdouble;κ::Cdouble=5.0,Nτ::Int64=20)

    input:

      - r,θ,y: cylindrical coordinates of the discretization of the cylinder (symmetry axis Oy)
      - x0,y0,z0: cartesian coordinate of a point belonging to the parametric straight line of integration
      - μ0: radius of the cylinder
      - Δr: transition distance of the mass function f
    
    optional input:

      - κ:  maximum radius for the integration is r=μ0+κΔr
      - Nτ: number of discretization nodes 

    output:

      - FGQ_rθy: 3D array containing the quadrature for each discretization nodes (r,θ,y)
      
"""
function quadrature_fg_cylinder(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,Δr::Cdouble;κ::Cdouble=5.0,Nτ::Int64=20)
    # square radius of the cylinder
    R2 = (μ0 + κ*Δr)^2

    # compute quadrature for each discretization nodes (r,θ,y)
    FGQ_rθy = zeros(Cdouble,length(r),length(θ),length(y));
    for n in eachindex(r)
        rm2 = r[n]^2
        if (R2>=rm2)
            for j in eachindex(θ)
                # cartesian coordinates
                xm  = r[n]*sin(θ[j])
                zm  = r[n]*cos(θ[j])
                for k in eachindex(y)
                    # quadrature ∫_0^{τ_root} f(g(τ)) dτ
                    τ_root,α,γ = τ_root_cylinder(x0,y0,z0,xm,y[k],zm,R2,rm2)
                    # this is not optimal, the details are at the transition, so, the discretization of τ should reflect that
                    τ_disc = collect(LinRange(0.0,τ_root,Nτ)) 
                    # r_param = sqrt.((xm.+α*τ_disc).^2 .+ (zm.+γ*τ_disc).^2)
                    a = α^2 + γ^2
                    b = 2.0*(α*xm + γ*zm)
                    c = xm^2 + zm^2;
                    r_param = g_sqrtquad.(τ_disc,a,b,c)
                    ρ_disc = f_logistic(r_param,μ0,Δr;A=1.0);
                    FGQ_rθy[n,j,k] = quadrature_fg(ρ_disc,τ_disc)
                end
            end
        end
    end

    # return 
    FGQ_rθy
end


"""
    quadrature_fg_cylinder_g_opt_f(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,Δr::Cdouble;κ::Cdouble=5.0,Nτ::Int64=20)

    quadrature evaluating the integral ∫_{τmin}^{τmax} f(g(τ)) dτ for a cylinder of radius μ0 (along the axis Oy)

    Contrary to quadrature_fg_cylinder the quadrature is adjusted to f instead of g

    see quadrature_fg_cylinder and quadrature_fg(f::Function,τ_min::Cdouble,τ_max::Cdouble,a::Cdouble,b::Cdouble,c::Cdouble,Δr0::Cdouble;Nr_min::Int64=20)
"""
function quadrature_fg_cylinder_g_opt_f(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,Δr::Cdouble;κ::Cdouble=5.0,Nτ::Int64=20)
    # square radius of the cylinder
    R2 = (μ0 + κ*Δr)^2

    # compute quadrature for each discretization nodes (r,θ,y)
    FGQ_rθy = zeros(Cdouble,length(r),length(θ),length(y));
    fint = (x::Cdouble->f_logistic(x,μ0,Δr;A=1.0))
    for n in eachindex(r)
        rm2 = r[n]^2
        if (R2>=rm2) # (R2>=rm2)
            for j in eachindex(θ)
                # cartesian coordinates
                xm  = r[n]*sin(θ[j])
                zm  = r[n]*cos(θ[j])
                for k in eachindex(y)
                    # quadrature ∫_0^{τ_root} f(g(τ)) dτ
                    τ_root,α,γ = τ_root_cylinder(x0,y0,z0,xm,y[k],zm,R2,rm2) # if b^2≃4*a*c, then the center of the cylinder is on the parametric line 
                    a = α^2 + γ^2
                    b = 2.0*(α*xm + γ*zm)
                    c = xm^2 + zm^2;
                    FGQ_rθy[n,j,k] = quadrature_fg(fint,0.0,τ_root,a,b,c,Δr/5.0;Nr_min=Nτ) 
                end
            end
        end
    end

    # return 
    FGQ_rθy
end
