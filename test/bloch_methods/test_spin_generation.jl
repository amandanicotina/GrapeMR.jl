using Test
using GrapeMR  # or the specific module name where `Spin` and `generate_spins` are defined

@testset "Spin struct and generate_spins function" begin
    # Define input values
    m_init = [0.0, 0.0, 1.0]
    t1s = [1.0, 2.0]
    t2s = [0.1, 0.2]
    b0s = [-0.1, 0.0, 0.1]
    b1s = [0.8, 1.0]
    targets = ["mz", "zero"]
    labels = ["spin1", "spin2"]

    # Generate spins
    spins = generate_spins(m_init, t1s, t2s, b0s, b1s, targets, labels)

    # Total number of spins = length(t1s) * length(b0s) * length(b1s)
    expected_total = length(t1s) * length(b0s) * length(b1s)
    @test length(spins) == expected_total

    # Check types and values
    for s in spins
        @test s isa Spin
        @test s.m_init == m_init
        @test s.T1 in t1s
        @test s.T2 in t2s
        @test s.b0_inho in b0s
        @test s.b1_inho in b1s
        @test s.target in string.(targets)
        @test s.label in string.(labels)
        @test s.n_spins == expected_total
    end

    # Test edge case: empty input
    empty_spins = generate_spins(m_init, Float64[], Float64[], b0s, b1s, Symbol[], Symbol[])
    @test isempty(empty_spins)
end
