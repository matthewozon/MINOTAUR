using Test
using MINOTAUR

function dummy_test()
  true
end

@testset "dummy test" begin
  @test dummy_test()
end



function test_logistic()
  x = [-100.0; -33.276; -33.274; 0.0; 33.274; 33.276; 100.0]
  r = logistic.(x,33.275)
  cond1 = (isapprox(r[1],0.0,atol=1.0e-14)) & (isapprox(r[1],r[2],atol=1.0e-14)) & (isapprox(r[2],r[3],atol=1.0e-14))
  cond2 = (isapprox(r[4],0.5,atol=1.0e-14))
  cond3 = (isapprox(r[5],1.0,atol=1.0e-14)) & (isapprox(r[5],r[6],atol=1.0e-14)) & (isapprox(r[6],r[7],atol=1.0e-14))

  r = logistic.(x,0.0,1.0,1.0,33.275)
  cond4 = (isapprox(r[1],0.0,atol=1.0e-14)) & (isapprox(r[1],r[2],atol=1.0e-14)) & (isapprox(r[2],r[3],atol=1.0e-14))
  cond5 = (isapprox(r[4],0.5,atol=1.0e-14))
  cond6 = (isapprox(r[5],1.0,atol=1.0e-14)) & (isapprox(r[5],r[6],atol=1.0e-14)) & (isapprox(r[6],r[7],atol=1.0e-14))

  r = f_logistic(-x,0.0,1.0;A=1.0)
  cond7 = (isapprox(r[1],0.0,atol=1.0e-14)) & (isapprox(r[1],r[2],atol=1.0e-14)) & (isapprox(r[2],r[3],atol=1.0e-14))
  cond8 = (isapprox(r[4],0.5,atol=1.0e-14))
  cond9 = (isapprox(r[5],1.0,atol=1.0e-14)) & (isapprox(r[5],r[6],atol=1.0e-14)) & (isapprox(r[6],r[7],atol=1.0e-14))

  cond1 & cond2 & cond3 & cond4 & cond5 & cond6 & cond7 & cond8 & cond9 
end


function test_g()
  # a>0.0
  a = 1.0;
  c = 1.0
  # b^2<4ac
  b = 1.0
  τ = [0.0; 1.0; -b/(2a)]
  r = g_sqrtquad.(τ,a,b,c)

  cond1 = (isapprox(r[1],sqrt(c),atol=1.0e-14)) & (isapprox(r[2],sqrt(a+b+c),atol=1.0e-14)) & (isapprox(r[3],sqrt(c-((b^2)/(4a))),atol=1.0e-14))

  # a>0.0
  a = 1.0;
  c = 1.0
  # b^2=4ac
  b = 2.0
  τ = [0.0; 1.0; -b/(2a)]
  r = g_sqrtquad.(τ,a,b,c)
  cond2 = (isapprox(r[1],sqrt(c),atol=1.0e-14)) & (isapprox(r[2],sqrt(a+b+c),atol=1.0e-14)) & (isapprox(r[3],sqrt(c-((b^2)/(4a))),atol=1.0e-14))

  cond1 & cond2
end


function test_g_sqrtquad_inv()
  # a>0.0
  a = 1.0;
  c = 1.0
  # b^2<4ac
  b = 1.0

  # τ ⩾ -b/(2a)
  τ = [0.0; 1.0; -b/(2a)]
  r = g_sqrtquad.(τ,a,b,c)
  τ_inv_pos =  g_sqrtquad_inv_pos.(r,a,b,c)
  cond1 = all(isapprox.(τ,τ_inv_pos,atol=1.0e-14))

  # τ ⩽ -b/(2a)
  τ = [-2.0; -1.0; -b/(2a)]
  r = g_sqrtquad.(τ,a,b,c)
  τ_inv_neg =  g_sqrtquad_inv_neg.(r,a,b,c)
  cond2 = all(isapprox.(τ,τ_inv_neg,atol=1.0e-14))

  cond1 & cond2
end

function test_interval()
  # a>0.0
  a = 1.0;
  c = 1.0
  # b^2<4ac
  b = 1.0
  # limit in τ
  τ0 = -b/(2a)
  
  # pos
  τ_min_0 = τ0 + 1.0
  τ_max_0 = τ_min_0 + 1.0
  rlow,rmid,rup = interval(τ_min_0,τ_max_0,a,b,c)
  cond1 = (isinf(rmid)) & (rlow<rup) & (isapprox(rlow,g_sqrtquad(τ_min_0,a,b,c),atol=1.0e-14)) & (isapprox(rup,g_sqrtquad(τ_max_0,a,b,c),atol=1.0e-14))

  # neg
  τ_max_0 = τ0 - 1.0
  τ_min_0 = τ_max_0 - 1.0
  rlow,rmid,rup = interval(τ_min_0,τ_max_0,a,b,c)
  cond2 = (isinf(rmid)) & (rlow<rup) & (isapprox(rlow,g_sqrtquad(τ_max_0,a,b,c),atol=1.0e-14)) & (isapprox(rup,g_sqrtquad(τ_min_0,a,b,c),atol=1.0e-14))

  # neg and pos 
  τ_max_0 = τ0 + 1.0
  τ_min_0 = τ0 - 0.5
  rlow,rmid,rup = interval(τ_min_0,τ_max_0,a,b,c)
  cond3 = (!isinf(rmid)) & (rmid<rup) & (rmid<rlow) & (isapprox(rlow,g_sqrtquad(τ_min_0,a,b,c),atol=1.0e-14)) & (isapprox(rup,g_sqrtquad(τ_max_0,a,b,c),atol=1.0e-14)) & (isapprox(rmid,g_sqrtquad(τ0,a,b,c),atol=1.0e-14))

  cond1 & cond2 & cond3
end

function test_quadrature()
  # a>0.0
  a = 1.0;
  c = 1.0
  # b^2<4ac
  b = 1.0
  # limit in τ
  τ0 = -b/(2a)

  # integration interval is a singleton, i.e. ∫=0
  f = (x::Cdouble->x)
  Δr0 = 0.1;
  Nr_min = 20000
  val = quadrature_fg(f,τ0,τ0,a,b,c,Δr0;Nr_min=Nr_min)
  cond1 = isapprox(val,0.0,atol=1.0e-14)

  # cases: τmin and τmax ⩾ τ0. f(x)=x and f(x)=x^2
  # analytical primitive
  f = (x::Cdouble->x)
  τmax = τ0 + 1.0;
  τmin = τ0 # 0.0
  y02 = (1/(4*a^2))*(4a*c-b^2)
  primitive_fg = (y::Cdouble->sqrt(a)*( (0.5y*sqrt(y02+y^2)) + (0.5y02)*log(abs(y+sqrt(y02+y^2))) ))
  integral_fg = primitive_fg(τmax-τ0) - primitive_fg(τmin-τ0)
  val = quadrature_fg(f,τmin,τmax,a,b,c,Δr0;Nr_min=Nr_min)

  cond2 = isapprox(integral_fg,val,atol=(τmax-τmin)/Nr_min)

  # another analytical case 
  f = (x::Cdouble->x^2)
  primitive_fg = (y::Cdouble->((a/3.0)*y^3) + ((b/2)*y^2) + c*y)
  integral_fg = primitive_fg(τmax) - primitive_fg(τmin);
  val = quadrature_fg(f,τmin,τmax,a,b,c,Δr0;Nr_min=Nr_min)

  cond3 = isapprox(integral_fg,val,atol=(τmax-τmin)/Nr_min)


  # cases: τmin ⩽ τ0  and τmax ⩾ τ0. f(x)=x and f(x)=x^2
  # analytical primitive
  f = (x::Cdouble->x)
  τmax = τ0 + 1.0;
  τmin = τ0 - 0.5
  y02 = (1/(4*a^2))*(4a*c-b^2)
  primitive_fg = (y::Cdouble->sqrt(a)*( (0.5y*sqrt(y02+y^2)) + (0.5y02)*log(abs(y+sqrt(y02+y^2))) ))
  integral_fg = primitive_fg(τmax-τ0) - primitive_fg(τmin-τ0)
  val = quadrature_fg(f,τmin,τmax,a,b,c,Δr0;Nr_min=Nr_min)

  cond4 = isapprox(integral_fg,val,atol=(τmax-τmin)/Nr_min)

  # another analytical case 
  f = (x::Cdouble->x^2)
  primitive_fg = (y::Cdouble->((a/3.0)*y^3) + ((b/2)*y^2) + c*y)
  integral_fg = primitive_fg(τmax) - primitive_fg(τmin);
  val = quadrature_fg(f,τmin,τmax,a,b,c,Δr0;Nr_min=Nr_min)

  cond5 = isapprox(integral_fg,val,atol=(τmax-τmin)/Nr_min)


  # cases: τmin and τmax ⩽ τ0. f(x)=x and f(x)=x^2
  # analytical primitive
  f = (x::Cdouble->x)
  τmax = τ0;
  τmin = τ0-1 # 0.0
  y02 = (1/(4*a^2))*(4a*c-b^2)
  primitive_fg = (y::Cdouble->sqrt(a)*( (0.5y*sqrt(y02+y^2)) + (0.5y02)*log(abs(y+sqrt(y02+y^2))) ))
  integral_fg = primitive_fg(τmax-τ0) - primitive_fg(τmin-τ0)
  val = quadrature_fg(f,τmin,τmax,a,b,c,Δr0;Nr_min=Nr_min)

  cond6 = isapprox(integral_fg,val,atol=(τmax-τmin)/Nr_min)

  # another analytical case 
  f = (x::Cdouble->x^2)
  primitive_fg = (y::Cdouble->((a/3.0)*y^3) + ((b/2)*y^2) + c*y)
  integral_fg = primitive_fg(τmax) - primitive_fg(τmin);
  val = quadrature_fg(f,τmin,τmax,a,b,c,Δr0;Nr_min=Nr_min)

  cond7 = isapprox(integral_fg,val,atol=(τmax-τmin)/Nr_min)


  cond1 & cond2 & cond3 & cond4 & cond5 & cond6 & cond7
end


function test_quadrature_bijection()
  # sin(log(x))
  f     = (r::Cdouble->sin(r))
  g     = (τ::Cdouble->log(τ))
  g_inv = (r::Cdouble->exp(r))

  f_lin     = (r::Cdouble->sin(log(r)))
  g_lin     = (τ::Cdouble->τ)
  g_lin_inv = (r::Cdouble->r)

  FG    = (τ::Cdouble->0.5τ*(sin(log(τ)) - cos(log(τ)))) 

  τ_min = exp(-3π/4-2π) # 0.001 
  τ_max = exp(π/4) # 2.0
  Δr0 = π/2.0 # big enough so that the optional argument Nr_min is the number of discretization nodes to compute the quadrature
  Nr_min = 100;

  # adjusted quadrature
  int_fg = quadrature_fg(f,g,g_inv,τ_min,τ_max,Δr0;Nr_min=Nr_min)

  # linear quadrature
  int_fg_lin = quadrature_fg(f_lin,g_lin,g_lin_inv,τ_min,τ_max,Δr0;Nr_min=Nr_min)

  # integral value
  val = FG(τ_max) - FG(τ_min)

  err_fg  = int_fg-val
  err_lin = int_fg_lin-val

  cond1 = (isapprox(err_fg,0.0,atol=2.0*(τ_max-τ_min)/Nr_min)) & (isapprox(err_lin,0.0,atol=2.0*(τ_max-τ_min)/Nr_min))


  # cos(log(x))
  f     = (r::Cdouble->cos(r))
  g     = (τ::Cdouble->log(τ))
  g_inv = (r::Cdouble->exp(r))

  f_lin     = (r::Cdouble->cos(log(r)))
  g_lin     = (τ::Cdouble->τ)
  g_lin_inv = (r::Cdouble->r)

  FG    = (τ::Cdouble->0.5τ*(sin(log(τ)) + cos(log(τ)))) 

  τ_min = exp(-3π/4-2π) # 0.001 
  τ_max = exp(π/4) # 2.0
  Δr0 = π/2.0 # big enough so that the optional argument Nr_min is the number of discretization nodes to compute the quadrature
  Nr_min = 100;

  # adjusted quadrature
  int_fg = quadrature_fg(f,g,g_inv,τ_min,τ_max,Δr0;Nr_min=Nr_min)

  # linear quadrature
  int_fg_lin = quadrature_fg(f_lin,g_lin,g_lin_inv,τ_min,τ_max,Δr0;Nr_min=Nr_min)

  # integral value
  val = FG(τ_max) - FG(τ_min)

  err_fg  = int_fg-val
  err_lin = int_fg_lin-val

  cond2 = (isapprox(err_fg,0.0,atol=2.0*(τ_max-τ_min)/Nr_min)) & (isapprox(err_lin,0.0,atol=2.0*(τ_max-τ_min)/Nr_min))


  # (log(x))^2
  f     = (r::Cdouble->r^2)
  g     = (τ::Cdouble->log(τ))
  g_inv = (r::Cdouble->exp(r))

  f_lin     = (r::Cdouble->(log(r))^2)
  g_lin     = (τ::Cdouble->τ)
  g_lin_inv = (r::Cdouble->r)

  FG    = (τ::Cdouble->τ*((log(τ))^2) - 2τ*log(τ) + 2τ) 

  τ_min = 0.25 # 0.001 
  τ_max = 2.0 # 2.0
  Δr0   = 0.5 # big enough so that the optional argument Nr_min is the number of discretization nodes to compute the quadrature
  Nr_min = 100;

  # adjusted quadrature
  int_fg = quadrature_fg(f,g,g_inv,τ_min,τ_max,Δr0;Nr_min=Nr_min)

  # linear quadrature
  int_fg_lin = quadrature_fg(f_lin,g_lin,g_lin_inv,τ_min,τ_max,Δr0;Nr_min=Nr_min)

  # integral value
  val = FG(τ_max) - FG(τ_min)

  err_fg  = int_fg-val
  err_lin = int_fg_lin-val

  cond3 = (isapprox(err_fg,0.0,atol=2.0*(τ_max-τ_min)/Nr_min)) & (isapprox(err_lin,0.0,atol=2.0*(τ_max-τ_min)/Nr_min))



  # 1/sqrt(x), x>0.0
  f     = (r::Cdouble->1.0/r)
  g     = (τ::Cdouble->sqrt(τ))
  g_inv = (r::Cdouble->r^2)

  f_lin     = (r::Cdouble->1.0/sqrt(r))
  g_lin     = (τ::Cdouble->τ)
  g_lin_inv = (r::Cdouble->r)

  FG    = (τ::Cdouble->2sqrt(τ)) 

  τ_min = 0.01
  τ_max = 1.0
  Δr0   = 0.5 # big enough so that the optional argument Nr_min is the number of discretization nodes to compute the quadrature
  Nr_min = 100;

  # adjusted quadrature
  int_fg = quadrature_fg(f,g,g_inv,τ_min,τ_max,Δr0;Nr_min=Nr_min)

  # linear quadrature
  int_fg_lin = quadrature_fg(f_lin,g_lin,g_lin_inv,τ_min,τ_max,Δr0;Nr_min=Nr_min)

  # integral value
  val = FG(τ_max) - FG(τ_min)

  err_fg  = int_fg-val
  err_lin = int_fg_lin-val

  cond4 = (isapprox(err_fg,0.0,atol=2.0*(τ_max-τ_min)/Nr_min)) & (isapprox(err_lin,0.0,atol=2.0*(τ_max-τ_min)/Nr_min))




  # 1/tan(x), x>0.0
  f     = (r::Cdouble->1.0/r)
  g     = (τ::Cdouble->tan(τ))
  g_inv = (r::Cdouble->atan(r))

  f_lin     = (r::Cdouble->1.0/tan(r))
  g_lin     = (τ::Cdouble->τ)
  g_lin_inv = (r::Cdouble->r)

  FG    = (τ::Cdouble->log(abs(sin(τ)))) 

  τ_min = π/100.0
  τ_max = π/4.0

  Δr0   = 0.5 # big enough so that the optional argument Nr_min is the number of discretization nodes to compute the quadrature
  Nr_min = 100;

  # adjusted quadrature
  int_fg = quadrature_fg(f,g,g_inv,τ_min,τ_max,Δr0;Nr_min=Nr_min)

  # linear quadrature
  int_fg_lin = quadrature_fg(f_lin,g_lin,g_lin_inv,τ_min,τ_max,Δr0;Nr_min=Nr_min)

  # integral value
  val = FG(τ_max) - FG(τ_min)

  err_fg  = int_fg-val
  err_lin = int_fg_lin-val

  cond5 = (isapprox(err_fg,0.0,atol=2.0*(τ_max-τ_min)/Nr_min)) & (isapprox(err_lin,0.0,atol=2.0*(τ_max-τ_min)/Nr_min))

  cond1 & cond2 & cond3 & cond4 & cond5
end

@testset "Discretization" begin
  @test test_logistic()
  @test test_g()
  @test test_interval()
  @test test_quadrature()
  @test test_quadrature_bijection()
end


# dependence of the error against the number of discfretization node
# f     = (r::Cdouble->sin(r))
# g     = (τ::Cdouble->log(τ))
# g_inv = (r::Cdouble->exp(r))

# f_lin     = (r::Cdouble->sin(log(r)))
# g_lin     = (τ::Cdouble->τ)
# g_lin_inv = (r::Cdouble->r)

# FG    = (τ::Cdouble->0.5τ*(sin(log(τ)) - cos(log(τ)))) 

# τ_min = exp(-3π/4-2π) # 0.001 
# τ_max = exp(π/4) # 2.0
# Δr0 = π/2.0 # big enough so that the optional argument Nr_min is the number of discretization nodes to compute the quadrature

# # # adjusted quadrature
# # int_fg_array[j] = quadrature_fg(f,g,g_inv,τ_min,τ_max,Δr0;Nr_min=Nr_min_array[j])

# # # linear quadrature
# # int_fg_lin_array[j] = quadrature_fg(f_lin,g_lin,g_lin_inv,τ_min,τ_max,Δr0;Nr_min=Nr_min_array[j])
# # integral value
# val = FG(τ_max) - FG(τ_min)

# # quadrature 
# Nr_min_array     = [10; 20; 50; 100; 200; 500; 1000; 2000; 5000; 10000; 20000]
# int_fg_array     = zeros(Cdouble,length(Nr_min_array))
# int_fg_lin_array = zeros(Cdouble,length(Nr_min_array))
# for j in eachindex(Nr_min_array)
#   int_fg_array[j]     = quadrature_fg(f,g,g_inv,τ_min,τ_max,Δr0;Nr_min=Nr_min_array[j])
#   int_fg_lin_array[j] = quadrature_fg(f_lin,g_lin,g_lin_inv,τ_min,τ_max,Δr0;Nr_min=Nr_min_array[j])
# end

# figure(figsize=[10,5]); 
# ax0 = subplot(121)
# τ_disc = [g_inv(r) for r in collect(LinRange(g(τ_min),g(τ_max),20))]
# scatter(collect(LinRange(τ_min,τ_max,20)),τ_disc)
# xlim(τ_min-0.05,τ_max+0.05)
# xlabel("Linear discretization",fontsize=12)
# ylabel("\$g^{-1}\$-discretization",fontsize=12)
# xticks(fontsize=12)
# yticks(fontsize=12)
# ax = subplot(122)
# loglog(Nr_min_array,abs.(int_fg_array.-val));
# loglog(Nr_min_array,abs.(int_fg_lin_array.-val));
# xlim(Nr_min_array[1],Nr_min_array[end])
# xlabel("Number of discretization nodes",fontsize=12)
# ylabel("quadrature error",fontsize=12)
# xticks(fontsize=12)
# yticks(fontsize=12)
# legend(["discretization adjusted to f", "linear discretization"],fontsize=12)
# tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)

# # s = @sprintf "\$\\int_{\\tau_1}^{\\tau_2}\\sin(\\log(\\tau)) d\\tau\$ = %1.10e" val
# # s = @sprintf "\$\\int_{\\tau_1}^{\\tau_2}\\sin(\\log(\\tau)) d\\tau\$ = %1.2f" val
# s = @sprintf "\$\\int_{\\tau_1}^{\\tau_2}\\sin(\\log(\\tau)) d\\tau\$ = 0"
# ax.annotate(s, xy=(3, 1),  xycoords="axes fraction", xytext=(0.1, 0.3), textcoords="axes fraction", color="black",fontsize=14)

# savefig("quadrature_wrt_nb_nodes.png")
# savefig("quadrature_wrt_nb_nodes.pdf")


function test_τ_root_cylinder()
  # case 1
  x0 = 1.0; y0 = 0.0; z0 = 1.0
  xm = 0.0; ym = 0.0; zm = 0.0; rm2 = xm^2+zm^2
  R2 = 1.0;
  dist_c,α,γ = τ_root_cylinder(x0,y0,z0,xm,ym,zm,R2,rm2)
  cond1 = (isapprox(dist_c,1.0,atol=1.0e-14)) & (isapprox(α,sqrt(2)/2,atol=1.0e-14)) & (isapprox(γ,sqrt(2)/2,atol=1.0e-14))

  # case 2
  x0 = 1.0; y0 = 0.0; z0 = 0.0
  xm = 0.0; ym = 0.0; zm = 0.0; rm2 = xm^2+zm^2
  R2 = 1.0;
  dist_c,α,γ = τ_root_cylinder(x0,y0,z0,xm,ym,zm,R2,rm2)
  cond2 = (isapprox(dist_c,1.0,atol=1.0e-14)) & (isapprox(α,1.0,atol=1.0e-14)) & (isapprox(γ,0.0,atol=1.0e-14))

  # case 3
  x0 = 1.0; y0 = 0.0; z0 = 1.0
  xm = 0.0; ym = 0.0; zm = 1.0; rm2 = xm^2+zm^2
  R2 = 1.0;
  dist_c,α,γ = τ_root_cylinder(x0,y0,z0,xm,ym,zm,R2,rm2)
  cond3 = (isapprox(dist_c,0.0,atol=1.0e-14)) & (isapprox(α,1.0,atol=1.0e-14)) & (isapprox(γ,0.0,atol=1.0e-14))

  # case 4
  x0 = 1.0; y0 = 0.0; z0 = 1.0
  xm = 1.0; ym = 0.0; zm = 0.0; rm2 = xm^2+zm^2
  R2 = 1.0;
  dist_c,α,γ = τ_root_cylinder(x0,y0,z0,xm,ym,zm,R2,rm2)
  cond4 = (isapprox(dist_c,0.0,atol=1.0e-14)) & (isapprox(α,0.0,atol=1.0e-14)) & (isapprox(γ,1.0,atol=1.0e-14))

  # case 5
  x0 = 1.0; y0 = 0.0; z0 = 1.0
  xm = 0.5; ym = 0.0; zm = 0.5; rm2 = xm^2+zm^2
  R2 = 1.0;
  dist_c,α,γ = τ_root_cylinder(x0,y0,z0,xm,ym,zm,R2,rm2)
  cond5 = (isapprox(dist_c,sqrt(R2)-sqrt(rm2),atol=1.0e-14)) & (isapprox(α,sqrt(2)/2,atol=1.0e-14)) & (isapprox(γ,sqrt(2)/2,atol=1.0e-14))

  cond1 & cond2 & cond3 & cond4 & cond5
end


function test_τ_root_sphere()
  # case 1
  x0 = 1.0; y0 = 0.0; z0 = 1.0
  xm = 0.0; ym = 0.0; zm = 0.0; rm2 = xm^2+ym^2+zm^2
  R2 = 1.0;
  dist_c,α,β,γ = τ_root_sphere(x0,y0,z0,xm,ym,zm,R2,rm2)
  cond1 = (isapprox(dist_c,1.0,atol=1.0e-14)) & (isapprox(α,sqrt(2)/2,atol=1.0e-14)) & (isapprox(β,0.0,atol=1.0e-14)) & (isapprox(γ,sqrt(2)/2,atol=1.0e-14))

  # case 2
  x0 = 1.0; y0 = 0.0; z0 = 0.0
  xm = 0.0; ym = 0.0; zm = 0.0; rm2 = xm^2+ym^2+zm^2
  R2 = 1.0;
  dist_c,α,β,γ = τ_root_sphere(x0,y0,z0,xm,ym,zm,R2,rm2)
  cond2 = (isapprox(dist_c,1.0,atol=1.0e-14)) & (isapprox(α,1.0,atol=1.0e-14)) & (isapprox(β,0.0,atol=1.0e-14)) & (isapprox(γ,0.0,atol=1.0e-14))

  # case 3
  x0 = 1.0; y0 = 0.0; z0 = 1.0
  xm = 0.0; ym = 0.0; zm = 1.0; rm2 = xm^2+ym^2+zm^2
  R2 = 1.0;
  dist_c,α,β,γ = τ_root_sphere(x0,y0,z0,xm,ym,zm,R2,rm2)
  cond3 = (isapprox(dist_c,0.0,atol=1.0e-14)) & (isapprox(α,1.0,atol=1.0e-14)) & (isapprox(β,0.0,atol=1.0e-14)) & (isapprox(γ,0.0,atol=1.0e-14))

  # case 4
  x0 = 1.0; y0 = 0.0; z0 = 1.0
  xm = 1.0; ym = 0.0; zm = 0.0; rm2 = xm^2+ym^2+zm^2
  R2 = 1.0;
  dist_c,α,β,γ = τ_root_sphere(x0,y0,z0,xm,ym,zm,R2,rm2)
  cond4 = (isapprox(dist_c,0.0,atol=1.0e-14)) & (isapprox(α,0.0,atol=1.0e-14)) & (isapprox(β,0.0,atol=1.0e-14)) & (isapprox(γ,1.0,atol=1.0e-14))

  # case 5
  x0 = 1.0; y0 = 0.0; z0 = 1.0
  xm = 0.5; ym = 0.0; zm = 0.5; rm2 = xm^2+ym^2+zm^2
  R2 = 1.0;
  dist_c,α,β,γ = τ_root_sphere(x0,y0,z0,xm,ym,zm,R2,rm2)
  cond5 = (isapprox(dist_c,sqrt(R2)-sqrt(rm2),atol=1.0e-14)) & (isapprox(α,sqrt(2)/2,atol=1.0e-14)) & (isapprox(β,0.0,atol=1.0e-14)) & (isapprox(γ,sqrt(2)/2,atol=1.0e-14))

  # case 6
  x0 = 1.0; y0 = 1.0; z0 = 1.0
  xm = 0.5; ym = 0.5; zm = 0.5; rm2 = xm^2+ym^2+zm^2
  R2 = 1.0;
  dist_c,α,β,γ = τ_root_sphere(x0,y0,z0,xm,ym,zm,R2,rm2)
  cond6 = (isapprox(dist_c,sqrt(R2)-sqrt(rm2),atol=1.0e-14)) & (isapprox(α,sqrt(3)/3,atol=1.0e-14)) & (isapprox(β,sqrt(3)/3,atol=1.0e-14)) & (isapprox(γ,sqrt(3)/3,atol=1.0e-14))
  
  cond1 & cond2 & cond3 & cond4 & cond5 & cond6
end


@testset "Distance" begin
  @test test_τ_root_cylinder()
  @test test_τ_root_sphere()
end

function test_quadrature_fg_cylinder()
  x0 = 2.0; y0 = 2.0; z0 = 2.0;
  μ0 = 1.0;
  Δr = 0.5;
  rmin = 0.0;
  rmax = 1.75;
  Nr = 50;
  r = collect(LinRange(rmin,rmax,Nr));
  Nθ = 101;
  θ = collect(LinRange(0.0,2.0π,Nθ));
  Ny = 50;
  ymin = 0.0;
  ymax = 4.0;
  y = collect(LinRange(ymin,ymax,Ny));

  # compare the two integration methods
  FGQ_rθy = quadrature_fg_cylinder(r,θ,y,x0,y0,z0,μ0,Δr;κ=1.0,Nτ=200);
  FGQ_rθy_g_opt_f = MINOTAUR.quadrature_fg_cylinder_g_opt_f(r,θ,y,x0,y0,z0,μ0,Δr;κ=1.0,Nτ=200);

  cond1 = (!any(isinf.(FGQ_rθy))) & (!any(isnan.(FGQ_rθy))) & (all(FGQ_rθy.>=0.0))
  cond2 = (!any(isinf.(FGQ_rθy_g_opt_f))) & (!any(isnan.(FGQ_rθy_g_opt_f))) & (all(FGQ_rθy_g_opt_f.>=0.0))
  cond3 = ((100*sum(abs.(FGQ_rθy-FGQ_rθy_g_opt_f))/sum(abs.(FGQ_rθy)))<0.1) # if less than 0.1% relative absolute value difference, it's just numerical discrepancies

  cond1 & cond2 & cond3
end

function test_quadrature_fg_sphere()
  x0 = 2.0; y0 = 2.0; z0 = 2.0;
  μ0 = 1.0;
  Δr = 0.5;
  rmin = 0.0;
  rmax = 1.75;
  Nr = 50;
  r = collect(LinRange(rmin,rmax,Nr));
  Nφ = 50;
  φ = collect(LinRange(0.0,1.0π,Nφ));
  Nθ = 101;
  θ = collect(LinRange(0.0,2.0π,Nθ));
  

  # compare the two integration methods
  FGQ_rφθ = quadrature_fg_sphere(r,φ,θ,x0,y0,z0,μ0,Δr;κ=1.0,Nτ=200);
  FGQ_rφθ_g_opt_f = MINOTAUR.quadrature_fg_sphere_g_opt_f(r,φ,θ,x0,y0,z0,μ0,Δr;κ=1.0,Nτ=200);

  cond1 = (!any(isinf.(FGQ_rφθ))) & (!any(isnan.(FGQ_rφθ))) & (all(FGQ_rφθ.>=0.0))
  cond2 = (!any(isinf.(FGQ_rφθ_g_opt_f))) & (!any(isnan.(FGQ_rφθ_g_opt_f))) & (all(FGQ_rφθ_g_opt_f.>=0.0))
  cond3 = ((100*sum(abs.(FGQ_rφθ-FGQ_rφθ_g_opt_f))/sum(abs.(FGQ_rφθ)))<0.1) # if less than 0.1% relative absolute value difference, it's just numerical discrepancies

  cond1 & cond2 & cond3
end


@testset "Cylinder and sphere" begin
  @test test_quadrature_fg_cylinder()
  @test test_quadrature_fg_sphere()
end


# Some plots might be interesting for offline debugging
# using PyPlot
# fig = figure() #figsize=[10,5]
# ax1 = subplot(111,polar=true)
# ax1.set_rticks([μ0/4, 2μ0/4, 3μ0/4, μ0])
# yticks(fontsize=12)
# xticks(fontsize=12)
# ax1.set_rlabel_position(35.0)
# ax1.plot(atan(x0,z0)*ones(Cdouble,2),[0.0; μ0], color="red",label="\$\\theta=45\$")
# # pcm1 = ax1.pcolormesh(θ,r,100abs.(FGQ_rθy_g_opt_f[:,:,3]-FGQ_rθy[:,:,3])./FGQ_rθy[:,:,3],edgecolors="face")
# pcm1 = ax1.pcolormesh(θ,r,100abs.(FGQ_rθy_g_opt_f[:,:,3]-FGQ_rθy[:,:,3])./FGQ_rθy_g_opt_f[:,:,3],edgecolors="face")
# cax1 = fig.add_axes([0.03, .1, 0.02, 0.3])
# cb1 = fig.colorbar(pcm1, orientation="vertical", cax=cax1, shrink=0.6)
# cb1.set_label("relative error [%]", fontsize=12)
# cb1.ax.tick_params(labelsize=12)
# ax1.legend(loc="lower left", bbox_to_anchor=(.5 + cos(atan(x0,z0))/2, .5 + sin(atan(x0,z0))/2),fontsize=12)
# tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)


# figure(); plot(r,FGQ_rθy[:,7,3]); plot(r,FGQ_rθy_g[:,7,3]); plot(r,FGQ_rθy_g_opt_f[:,7,3]); 
# figure(); plot(r,FGQ_rθy[:,17,3]); plot(r,FGQ_rθy_g[:,17,3]); plot(r,FGQ_rθy_g_opt_f[:,17,3]); 
# figure(); plot(r,FGQ_rθy[:,7,3]-FGQ_rθy_g[:,7,3]); plot(r,FGQ_rθy[:,7,3]-FGQ_rθy_g_opt_f[:,7,3]); plot(r,FGQ_rθy_g[:,7,3]-FGQ_rθy_g_opt_f[:,7,3])
# figure(); plot(r,FGQ_rθy[:,17,3]-FGQ_rθy_g[:,17,3]); plot(r,FGQ_rθy[:,17,3]-FGQ_rθy_g_opt_f[:,17,3]); plot(r,FGQ_rθy_g[:,17,3]-FGQ_rθy_g_opt_f[:,17,3])
# figure(); plot(r,100.0*(FGQ_rθy[:,17,3]-FGQ_rθy_g_opt_f[:,17,3])./FGQ_rθy[:,17,3]); plot(r,100.0*(FGQ_rθy_g[:,17,3]-FGQ_rθy_g_opt_f[:,17,3])./FGQ_rθy[:,17,3])