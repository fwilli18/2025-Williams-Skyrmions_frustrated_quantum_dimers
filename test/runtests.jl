using DrWatson, Test
@quickactivate "2025-Williams-Skyrmions_frustrated_quantum_dimers"

# Here you include files using `srcdir`
# include(srcdir("file.jl"))

# Run test suite
println("Starting tests")
ti = time()

@testset "2025-Williams-Skyrmions_frustrated_quantum_dimers tests" begin
    @test 1 == 1
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
