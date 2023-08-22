[![MINAUTOR CI](https://github.com/matthewozon/MINOTAUR/actions/workflows/CI.yml/badge.svg)](https://github.com/matthewozon/MINOTAUR/actions/workflows/CI.yml)
[![Documentation](https://github.com/matthewozon/MINOTAUR/actions/workflows/documentation.yml/badge.svg)](https://github.com/matthewozon/MINOTAUR/actions/workflows/documentation.yml)

# MINOTAUR
Julia package for computing quadrature for integrals of composed functions

```math
\int_{\tau_1}^{\tau_2} f(g(\tau)) \mathrm{d}\tau
```

name: coMposition INtegral cOmpuTing quAdratURe

# Examples

## Error w.r.t. the numebr of discretization nodes

The computation of the integral

```math
\int_{e^{-\frac{11\pi}{4}}}^{e^{\frac{\pi}{4}}} \sin\left(\ln \tau \right) \mathrm{d}\tau = \left[\frac{\tau}{2}\left(\sin\left(\ln x\right) - \cos\left(\ln x\right)\right)\right]_{e^{-\frac{11\pi}{4}}}^{e^{\frac{\pi}{4}}}
```

can be carried out using MINOTAUR and compared to the result of another quadrature with the same number of discretization nodes not adjusted to the $\sin$ function (linearly spaced nodes $\tau_i = \tau_{i-1} + \Delta\tau$) with the code:

```
f     = (r::Cdouble->sin(r))
g     = (τ::Cdouble->log(τ))
g_inv = (r::Cdouble->exp(r))
FG    = (τ::Cdouble->0.5τ*(sin(log(τ)) - cos(log(τ)))) 

τ_min = exp(-3π/4-2π) # 0.001 
τ_max = exp(π/4) # 2.0
Δr0 = π/2.0 # big enough so that the optional argument Nr_min is the number of discretization nodes to compute the quadrature

# integral value
val = FG(τ_max) - FG(τ_min)

# quadrature 
Nr_min_array     = [10; 20; 50; 100; 200; 500; 1000; 2000; 5000; 10000; 20000]
int_fg_array     = zeros(Cdouble,length(Nr_min_array))
int_fg_lin_array = zeros(Cdouble,length(Nr_min_array))
for j in eachindex(Nr_min_array)
  int_fg_array[j] = quadrature_fg(f,g,g_inv,τ_min,τ_max,Δr0;Nr_min=Nr_min_array[j])
  
  # linear discretization
  τ_disc_lin = collect(LinRange(τ_min,τ_max,Nr_min_array[j]))
  fs_lin = [f(g(τ))  for τ in τ_disc_lin]
  # trapezoid rule 
  for i in 1:(Nr_min_array[j]-1)
    int_fg_lin_array[j] = int_fg_lin_array[j] + (1.0/2.0)*(fs_lin[i]+fs_lin[i+1])*(τ_disc_lin[i+1]-τ_disc_lin[i])
  end
end
```

In this figure: left) the discretization adjusted for $f$ is plotted against the linear discretization , and right) the absolute value of the error between the quadrature and the true value of the integral (in blue, using the adjusted discretization, in orange, using the linear discretization of $\tau$)

![quadrature_wrt_nb_nodes](https://github.com/matthewozon/MINOTAUR/assets/7929598/6fc90605-a205-4ae6-94f7-fc4059b6c0e8)


# Install

In the Julia REPL:

```
] add https://github.com/matthewozon/MINOTAUR
```

