@testset "utility tests: round_zeros_ones!" begin
    A1 = [0 -0.99999; 0.999999 2.000000001; 1.0000001 -1.000001; -0.000001 0.0000001]
    A2 = [0 -0.99999; 0.999999 2.000000001; 1.0000001 -1.000001; -0.000001 0.0000001]

    LOpt.round_zeros_ones!(A1, tol = 1E-5)
    @test A1[1, 2] == -1
    @test A1[2, 1] == 1
    @test A1[2, 2] > 2
    @test A1[3, 1] == 1
    @test A1[3, 2] == -1
    @test A1[4, 1] == 0
    @test A1[4, 2] == 0

    LOpt.round_zeros_ones!(A2, tol = 1E-6)
    @test A2[1, 2] > -1
end

@testset "utility tests: optimal_graph_edges" begin
    z_val = [0.0 1 1; 1 0 0; 1 0 0]
    edges = LOpt.optimal_graph_edges(z_val)
    @test edges[1] == (1, 2)
    @test edges[2] == (1, 3)
end

@testset "utility tests: algebraic_connectivity and fiedler_vector" begin
    adj_mat_1 = [0 1 2; 1 0 3; 2 3 0.0]
    graph_data = LOpt.GraphData(adj_mat_1)
    @test isapprox(graph_data.ac, 4.267949192, atol = 1E-6)
    v = [0.7886751345948128, -0.5773502691896261, -0.21132486540518727]
    @test isapprox(v, graph_data.fiedler, atol = 1E-6)

    adj_mat_2 = [0 1 2; 1 0 0; 2 0 0.0]
    graph_data = LOpt.GraphData(adj_mat_2)
    @test isapprox(graph_data.ac, 1.26794919243, atol = 1E-6)
    v = [-0.2113248654051879, 0.7886751345948118, -0.5773502691896268]
    @test isapprox(v, graph_data.fiedler, atol = 1E-6)

    # Disconnected graph
    adj_mat_3 = [0 1 1 0; 1 0 0 0; 1 0 0 0; 0 0 0 0.0]
    graph_data = LOpt.GraphData(adj_mat_3)
    @test isapprox(graph_data.ac, 0, atol = 1E-6)
    v = [0.5773502691896258, 0.5773502691896273, 0.5773502691896241, 0.0]
    @test isapprox(v, graph_data.fiedler, atol = 1E-6)
end

@testset "utility tests: cheeger constant evaluation" begin

    # Complete bipartite graph
    n = 6
    n1 = Int(ceil(n / 2))
    G = zeros(n, n)
    for i in 1:n1, j in (n1+1):n
        G[i, j] = 1.0
        G[j, i] = G[i, j]
    end

    cheeger =
        LOpt.cheeger_constant(G, glpk_optimizer, formulation_type = "set_cardinality")
    @test isapprox(cheeger["cheeger_constant"], 5 / 3, atol = 1E-6)

    cheeger = LOpt.cheeger_constant(G, glpk_optimizer, formulation_type = "set_volume")
    @test isapprox(cheeger["cheeger_constant"], 5 / 9, atol = 1E-6)
end
