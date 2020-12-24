# code for testing our package
using PSVrelief, Test

@testset "steam" begin
    # test case steam
    # From API 520, Part I, 2000 (section 3.7.2)
    # W = 69615 kg/h
    # P = 11032*1.1 + 101.3 = 12236 kPa
    # saturated
    # API gets 1100 mm2, and this is with a Napier factor of 51.5, we are using 51.45
    @test (abs(PSVrelief.PSVsteamSize(69615.0, 12236.0, "Sat") - 1100.0) < 2.0)
end

@testset "vapour" begin
    # test case vapour
    # From API 520, Part I, 2013 (section 5.6.3.2.3)
    # W = 24270.0 kg/h
    # P = 670.0 kPa
    # MW = 51
    # Z = 0.9
    # k = Cp/Cv = 1.11, or C = 328.
    # T = 348 - 273
    # API gets 3698 mm2
    # W, P, Tcelcius, MW, k, Z
    @test (abs(PSVrelief.PSVvaporSize(24270.0, 670.0, 75.0, 51.0, 1.11, 0.9) - 3698) < 2.0)
end

@testset "liquid1" begin
    # test case liquid1
    # from API 520, Part I, 2013 (section 5.8.1.3)
    # flow = 6814 litres/min
    # density = 900 kg/m3
    # P = 1896 kpag (incudes 10%)
    # Pback = 345 kPag
    # viscosity 2000 SSU (API uses archaic units, the table gives 431.7 cSt)
    # API gets 3180 mm2 with a backpressure correction of 0.97, or 3085 without the backpressure correciton
    # PSVliquidSize(W, P, Pback, d, mu)
    # PSVliquidSize(900.0*6814*60.0/1000.0, 1896.0, 345.0, 900.0, 431.7*0.900)
    @test (abs(PSVrelief.PSVliquidSize(900.0*6814*60.0/1000.0, 1896.0, 345.0, 900.0, 431.7*0.900) - 3085.0) < 10.0)
end

@testset "liquid2" begin
    # test case liquid2
    # from API 520, Part I, 2013 (section 5.8.1.3)
    # flow = 6814 litres/min
    # density = 900 kg/m3
    # P = 1896 kpag (incudes 10%)
    # Pback = 345 kPag
    # viscosity 2000 SSU (API uses archaic units, the table gives 431.7 cSt)
    # we will flip around the sizing case to get the rating case
    # API gets 3180 mm2 with a backpressure correction of 0.97, or 3085 without the backpressure correciton
    # PSVliquidSize(W, P, Pback, d, mu)
    targetFlow = 900.0*6814*60.0/1000.0
    targetArea = PSVrelief.PSVliquidSize(targetFlow, 1896.0, 345.0, 900.0, 431.7*0.900)
    @test (abs(PSVrelief.PSVliquidRate(targetArea, 1896.0, 345.0, 900.0, 431.7*0.900) - targetFlow) < 1.0)
end