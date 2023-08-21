"""
This is the [`MINOTAUR`](@ref), it contains 
* [`MINOTAUR.MINOTAUR`](@ref)
* [`MINOTAUR.logistic`](@ref)
* [`MINOTAUR.f_logistic`](@ref)
* [`MINOTAUR.g_sqrtquad`](@ref)
* [`MINOTAUR.g_sqrtquad_inv_pos`](@ref)
* [`MINOTAUR.g_sqrtquad_inv_neg`](@ref)
* [`MINOTAUR.interval`](@ref)
* [`MINOTAUR.quadrature_fg`](@ref)
* [`MINOTAUR.τ_root_sphere`](@ref)
* [`MINOTAUR.quadrature_fg_sphere`](@ref)
* [`MINOTAUR.τ_root_cylinder`](@ref)
* [`MINOTAUR.quadrature_fg_cylinder`](@ref)



    Attempt to optimize the discretization for the quadrature to compute the integral

        ∫_τ1^τ2 f(g(τ)) dτ
    
    assuming that g is a bijection, and that g^{-1} assumes non-infinite values

    We choose to discretize τ so that enough details of f are represented

    We choose a subdivision r_i that captures the details of f over the interval [min(g(τ1),g(τ2)),max(g(τ1),g(τ2))]

"""
module MINOTAUR

export logistic, f_logistic
export g_sqrtquad, g_sqrtquad_inv_pos, g_sqrtquad_inv_neg
export interval, quadrature_fg
include("discretization.jl")

export τ_root_sphere, quadrature_fg_sphere
include("sphere.jl")

export τ_root_cylinder,quadrature_fg_cylinder
include("cylinder.jl")
end
