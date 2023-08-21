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


@testset "Discretization" begin
  @test test_logistic()
  @test test_g()
  @test test_interval()
  @test test_quadrature()
end