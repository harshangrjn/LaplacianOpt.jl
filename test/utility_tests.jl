@testset "utility tests: get_rounded_zeros_and_ones!" begin
    A1 = [0 -0.99999; 0.999999 2.000000001; 1.0000001 -1.000001; -0.000001 0.0000001]
    A2 = [0 -0.99999; 0.999999 2.000000001; 1.0000001 -1.000001; -0.000001 0.0000001]

    LO.get_rounded_zeros_and_ones!(A1, 1E-5)
    @test A1[1,2] == -1
    @test A1[2,1] ==  1
    @test A1[2,2] >  2
    @test A1[3,1] == 1
    @test A1[3,2] == -1
    @test A1[4,1] == 0
    @test A1[4,2] == 0

    LO.get_rounded_zeros_and_ones!(A2, 1E-6)
    @test A2[1,2] > -1

end

@testset "utility tests: get_optimal_graph_edges" begin

    num_nodes = 3 
    z_val = [0.0 1 1; 1 0 0; 1 0 0]
    edges = LO.get_optimal_graph_edges(num_nodes, z_val)
    @test edges[1] == (1,2)
    @test edges[2] == (1,3)

end