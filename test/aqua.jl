using Aqua, Manifolds, Test, StaticArrays # Last package temporary

@testset "Aqua.jl" begin
    Aqua.test_all(
        Manifolds;
        ambiguities=(
            exclude=[
                *,
                ==,
                allocate,
                Manifolds.allocate_coordinates,
                allocate_result,
                check_point,
                copyto!,
                compose,
                compose!,
                get_coordinates,
                get_coordinates!,
                get_vector,
                get_vector!,
                getindex,
                inv,
                inv!,
                mean,
                Manifolds.mul!,
                reshape,
                similar,
                setindex!,
                view,
                Manifolds.AbstractGroupOperation,
            ],
            broken=false,
        ),
        piracies=(
            treat_as_own=[
                AbstractDecoratorManifold, # MAybe fix?
                AbstractManifold, # Maybe fix?
                AbstractNumbers,  # Maybe fix?
                ProductManifold,  # Maybe fix?
                ExtrinsicEstimation, # already deprecated
                Manifolds.EmptyTrait, # Maybe fix?
                Manifolds.TraitList, # Maybe fix?
                Manifolds.GeodesicInterpolationWithinRadius, # Probably fix
                # StaticArray, # Definetly fix!
                allocate, # Maybe fix
            ],
        ),
    )
end
