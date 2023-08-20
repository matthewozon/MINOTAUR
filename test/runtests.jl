using Test
using MINOTAUR

function dummy_test()
  true
end

@testset "dummy test" begin
  @test dummy_test()
end
