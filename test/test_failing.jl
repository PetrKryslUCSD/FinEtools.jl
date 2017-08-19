module mphunm1
using FinEtools
using Base.Test
function test()
    u(expression) = phun(expression)
    t = 0.333*u("s")
    v = 2.0*u("m/s")
    display(v*t)
end
end
using mphunm1
mphunm1.test()

module mphunm3
using FinEtools
using Base.Test
function test()
    for sou in [:SI :US :IMPERIAL :CGS :SIMM]
        t = 0.333*phun(system_of_units = sou, "s")
        v = 2.0*phun(system_of_units = sou, "m/s")
        display(v*t/phun(system_of_units = sou, "m"))
        @test abs(v*t/phun(system_of_units = sou, "m") - 0.666) < 1.0e-6
    end
end
end
using mphunm3
mphunm3.test()
