#!/usr/bin/env julia
#
# To run the Manifolds.jl Tests as a script with Command line options,
# You have to have TestEvn in your environment.
#
# Command Line Options (stil WIP)
#
include("runtests_cmd_utils.jl")

if (abspath(PROGRAM_FILE) == @__FILE__) # Run as script
    if "--help" ∈ lowercase.(ARGS)
        println(
            """
    # Manifolds.jl Test Suite Command Line Arguments

    * `--exclude-manifold` - Specify manifold(s) to exclude (either from all or from the `--manifold` argument)
    * `--exclude-function` - Specify manifold(s) to exclude (either from all or from the `--manifold` argument)
    * `--help`             - print this help
    * `--manifold`         - Specify manifold(s) to test (only these are tested)
    * `--function`         - Specify manifold(s) to test (only these are tested)
    * `--element-type`     - (planned feature) Specify element type(s) to use in default arrays (default: `Float64`)

    ## Examples

    * `./runtests.jl` runs the default (full suite)
    * `./runtests.jl --manifold Sphere` only run the tests on the sphere(s)
    * `./runtests.jl --function exp` only run the `exp`onential map tests
    * `./runtests.jl --exclude-function exp` run all tests _except `exp`onential map tests
    """,
        )
        exit()
    end
    using TestEnv
    TestEnv.activate("Manifolds")
    using Manifolds
    M = Sphere(2)
    i = 1
    only_manifolds = Type[]
    exclude_manifolds = Type[]
    only_functions = Function[]
    exclude_functions = Function[]
    while i <= length(ARGS)
        if startswith(ARGS[i], "--")
            cmd = ARGS[i][3:end]
            if lowercase(cmd) == "manifold"
                (j, candidates) = collect_args(i + 1, ARGS)
                global i = j #advance to next field
                global only_manifolds = process_manifolds(candidates)
            end
            if lowercase(cmd) == "exclude-manifold"
                (j, candidates) = collect_args(i + 1, ARGS)
                global i = j #advance to next field
                global exclude_manifolds = process_manifolds(candidates)
            end
            if lowercase(cmd) == "function"
                (j, candidates) = collect_args(i + 1, ARGS)
                global i = j #advance to next field
                global only_functions = process_functions(candidates)
            end
            if lowercase(cmd) == "exclude-function"
                (j, candidates) = collect_args(i + 1, ARGS)
                global i = j #advance to next field
                global exclude_functions = process_manifolds(candidates)
            end
        else
            global i += 1
        end
    end
    # Set these in ENV
    ENV["MANIFOLDS_TESTSUITE_ONLY_MANIFOLDS"] = only_manifolds
    ENV["MANIFOLDS_TESTSUITE_EXCLUDE_MANIFOLDS"] = only_manifolds
    ENV["MANIFOLDS_TESTSUITE_ONLY_FUNCTIONS"] = only_functions
    ENV["MANIFOLDS_TESTSUITE_EXCLUDE_FUNCTIONS"] = exclude_functions
    print(ENV)
    # Setup test variants
    #...then continue with the test run in normal mode below
end
include("header.jl")
include("aqua.jl")

# new tests here

# run old / legacy tests
include("run_legacy_tests.jl")
