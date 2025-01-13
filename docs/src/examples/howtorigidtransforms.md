# [How to do Rigid Transforms](@id how_to_do_rigid_transforms)

This tutorial follows the previous tutorial [How to work with Rotations](@ref how_to_work_with_rotations).  The examples aim to demonstrate how familiar rigid body transformations in 2D can be achieved with the JuliaManifolds packages.  The previous tutorial showed how to use the `SpecialOrthogonal(2)` manifold by itself.  This tutorial shows how the `ProductManifold` of `Translation` with `Rotations` is combined into a predefined [`Manifolds.SpecialEuclidean`](@ref).

The rigid transforms examples discussed here are in the same spirit as packages such as [CoordinateTransformations.jl](https://github.com/JuliaGeometry/CoordinateTransformations.jl).  Manifolds.jl instead builds coordinate transformations from the ground up using fundamental abstractions about generalized manifolds.  The `SpecialEuclidean(2)` group manifold used in this tutorial is one specific manifold, which is the product of two underlying manifolds.

!!! note
    While many of the operations listed below at first seem excessive for basic rigid transform operations, this is not the case.  Once the data formats and representations are setup, the actual `compose` and `apply` operations are super fast and syntactically easy.  Significant efforts have gone into defining a general framework for manifold operations that can go beyond operations such as Riemannian geometry.

## Rigid Transforms Setup

Let's load load some packages first
```julia
using Manifolds
using LinearAlgebra
using StaticArrays
```

Consider a combination of translation with a rotation of some avatar or axis-dyad on the xy-plane, where each location (i.e. the combo or position and rotation) which will be represented as some point `p` on a group `G` manifold `M`.  This is indicated by:
```julia
# xy dimension in this case is 2
n = 2
# the group manifold also known as SE(2)
G = SpecialEuclidean(n)
M = base_manifold(G)
# can separately take the Translation and Rotations manifolds within M
M_T = M[1]
M_R = M[2]
```

Next, define the default basis to use later and also a starting reference rotation and translation
```julia
B = DefaultOrthogonalBasis()
# choose the identity reference translation and rotation
Tr0, R0 = SA[0.0; 0], SA[1.0 0; 0 1]
```

To structure the tutorial, let's "walk" the perimeter of a rectangle, passing through various points on the manifold.  To do so, we need to understand how coordinates, vectors in the Lie algebra, and group elements are used.

## From Coordinates, SE(2)

The coordinates are easily readible / sharable numbers that we humans frequently stack together as a vector.  In this case, the three degrees of freedom representing a transformation, namely `[Δx; Δy; Δθ]`.  To reduce the notational load, this example will drop the `Δ` and simply state `[x;y;θ]` as the transform coordinates (i.e. a point on the `SpecialEuclidean(2)` manifold, but not yet in the group data representation format).  

First, define the individual manifold tangent vectors (translation, rotation) of the identity element coordinates
```julia
t0, w0 = [0.0; 0],  get_vector(M_R, R0, 0, B)
```

!!! note
    The `hat` notation could also be used in this case:
    ```julia
    _t0, _w0 = [0.0; 0],  hat(M_R, R0, 0)
    @assert isapprox(t0, _t0); @assert isapprox(w0, _w0)
    ```

Next, define delta elements that each represent a certain "rigid body transform" so that the rectangular path can be completed.  Starting from the "origin" (i.e. reference point) above, we define relative segments between points using the coordinates and we will `compose` these segments together in an upcoming section below:
```julia
# walk forward on side x and turn left
t1a, w1a = [1.0; 0],  get_vector(M_R, R0, 0  , B)
t1b, w1b = [1.0; 0],  get_vector(M_R, R0, π/2, B)

# walk positive y and turn left again
t2, w2 = [1.0; 0],  get_vector(M_R, R0, π/2, B)

# walk negative x and turn left
t3, w3 = [2.0; 0],  get_vector(M_R, R0, π/2, B)

# walk negative y back onto origin and turn left
t4, w4 = [1.0; 0],  get_vector(M_R, R0, π/2, B)
```

Notice that each of these coordinates are **not** defined against the origin, but are fully relative translations and rotations.

Before we can compose these separate translation and rotation pairs, we need to leverage some pairing data structure by which Manifolds.jl can identify a Product manifold operation is needed.  There are currently two ways to do so.  The `compose` hereafter can work with either the `ProductRepr` or `ProductArray` representations.

### ProductRepr for SpecialEuclidean(2)

Let's first define the identity group element on the manifold directly:
```julia
p0 = ProductRepr(Tr0, R0)
```

!!! note
    Using the language or reference frames, this is the identity reference to which other points will be defined.  For example, consider measuring various points relative to a common "reference frame", then `p0` might represent that reference point in-and-about the other points on the same `SpecialEuclidean(2)` manifold.

Above, we defined individual tranlate and rotate vectors such as `t2,w2`.  The next step towards using `compose` for rigid body transformations is to properly pair (i.e. data structuring) the Lie algebra elements and then map to the equivalent Lie group element on the group manifold `G` or `M`.

The trivial case is a zero vector from the point `p0`
```julia
# make Lie Algebra element 
x0 = ProductRepr(t0, w0)

p0_ = exp(G, p0, x0)  # or,  exp(M, p0, x0)
```

The exponential maps the Lie algebra to the associated Lie group, but in this trivial zero-vector case should be exactly the same point on the manifold
```julia
@assert isapprox(p0_.parts[1], p0.parts[1]); @assert isapprox(p0_.parts[2], p0.parts[2])
```

All the other relative vector segments are mapped onto the manifold, and carefully not relative to the point `p0` in each case -- we will soon show why
```julia
# calculate the exponential mapping from point p0
x1a = ProductRepr(t1a, w1a)   # Lie algebra element
p1a = exp(G, p0, x1a)         # Lie group element

x1b = ProductRepr(t1b, w1b)
p1b = exp(G, p0, x1b)

x2 = ProductRepr(t2, w2)
p2 = exp(G, p0, x2)

x3 = ProductRepr(t3, w3)
p3 = exp(G, p0, x3)

x4 = ProductRepr(t4, w4)
p4 = exp(G, p0, x4)
```

If you want to continue using `ProductRepr` then skip ahead to [Compose Rigid Transforms, 2D](@ref compose_rigid_transforms_2d).

### ProductArray for SpecialEuclidean(2)

As mentioned earlier, there are two data representation formats available for product manifolds.  This paragraph will repeat the same steps as the `ProductRepr` paragraph above, but using the `ProductArray` representation instead.

!!! note
    `ProductArray` might seem a bit more laborious than `ProductRepr` at first, but the increased structure allows for faster implementation details within the Manifolds.jl ecosystem.

The underlying data is represented according to a shape specification
```julia
reshaper = Manifolds.StaticReshaper()
shape_G = Manifolds.ShapeSpecification(reshaper, M)
shape_se = Manifolds.ShapeSpecification(reshaper, M.manifolds...)
```

And the associated Lie group elements (still `SpecialEuclidean(2)`) are mapped from the vector elements above
```julia
p0  = Manifolds.prod_point(shape_se, (t0,exp(M_R,  R0, w0))... )
p1a = Manifolds.prod_point(shape_se, (t1a,exp(M_R, R0, w1a))... )
p1b = Manifolds.prod_point(shape_se, (t1b,exp(M_R, R0, w1b))... )
p2  = Manifolds.prod_point(shape_se, (t2,exp(M_R,  R0, w2))... )
p3  = Manifolds.prod_point(shape_se, (t3,exp(M_R,  R0, w3))... )
p4  = Manifolds.prod_point(shape_se, (t4,exp(M_R,  R0, w4))... )
```

Just as with the `ProductRepr` approach, each of the points, e.g. `p2`, represent a relative rigid transform -- note again that the exponential map was done relative to the choosen identity `R0`.  The identity point for translations `M_T` are trivial, i.e. `[0;0]`, and therefore not elaborated as is done with the Rotations under `M_R`.

## [Compose Rigid Transforms, 2D](@id compose_rigid_transforms_2d)

In this section all the details from above are brought together via the action `compose`.  Recall the objective for this tutorial is to walk around a rectangle using the `Manifolds.SpecialEuclidean(2)` mechanizations and abstractions.  Note that either the `ProductRepr` or `ProductArray` representations will dispatch the same for the code below, but the internal data representation formats will differ.

The fine grain details in this section are perhaps not as straight forward as they first appear, specifically how relative transforms (which are points on the manifold `SpecialEuclidean(2)`) relate to the consecutive positions of the walked rectangle (which are also points on the manifold `SpecialEuclidean(2)`).

```julia
## Consecutively compose the "rigid transforms" (i.e. points) together around a rectangle

p_01a   = compose(G, p0,     p1a)
p_01b   = compose(G, p_01a,  p1b)
p_012   = compose(G, p_01b,  p2)
p_0123  = compose(G, p_012,  p3)
p_01234 = compose(G, p_0123, p4)
```

Staring from the choosen origin `p0` and looking in the direction along the choosen x-axis (i.e. `θ = 0`), walk forwards 1.0 units towards a new point `p_01a`.  This relative transform from `p0` on the manifold to the new point `p_01a` on the manifold is captured in the group element `p1a` (also `SpecialEuclidean(2)`).

Next, walk from `p_01a` to `p_01b` according to the relative transformation captured by the group element `p1b` which is again forwards 1.0 units, but also followed by a positive rotation of `π/2` radians.  Taking the default basis `B` as an xy-tangent plane, we understand a positive rotation according to the right-hand-rule as going from the local x-axis towards the y-axis.

!!! note
    As a sneak peak, we expect this new point `p_01b` in human readible coordinates to be at `[2;0;π/2]` relative to the choosen reference point `p0`.  We will show and confirm this computation below.

It is worth reiterating that the relative group actions `p1a` and `p1b` are Lie group elements that were constructed from local Lie algebra vectors using `p0` in each case.  Here, during the `compose` operation those "relative transforms" are chained together (i.e. composed) into new (totally separate) points `p_01a, p_01b` on the `SpecialEuclidean(2)` manifold.  Our cartoon interpretation of walking around a rectangle mentally easiest when operating relative to the starting point `p0`.  The full use of Manifolds.jl goes well beyond this "rigid transform" construct.

To complete the walk-about, we compose `p2` onto the previous to get a new point `p_012` which is now the opposite corner of the rectangle.  Two more consecutive "transforms" `p3, p4` complete the rectangle through points `p_0123` and `p_02134`, respectively, and each time turning `π/2` radians.  We therefore expect point `p_01234` to exactly the same as the reference starting point `p0`.

## To Coordinates, SE(2)

Thus far, we converted user inputs from a convenient coordinate represenation into two computational representations (either `ProductRepr` or `ProductArray`).  The next step is to be able to convert back from the manifold data representations back to coordinates.  In this example, we first take the logmap of a group element in `SpecialEuclidean(2)` to get the associated Lie algebra element, and then extract the coordinates from that vector, and lets do so for the point `p_01` as elluded to above
```julia
x_01b_ = get_coordinates(G, p0, log(G, p0, p_01b), B)

# check the coordinates are as expected
@assert isapprox(     x_01b_[1:2], [2;0], atol = 1e-10 )
@assert isapprox( abs(x_01b_[3]),   π/2 , atol = 1e-10 )
```

Notice that both the `log` and `get_coordinates` operations are performed relative to our choosen reference point `p0`.  These operations therefore are extracting the tangent space vector (and its coordinates) relative to the origin.  If we were to extract the vector coordinates of `p_01b` relative to `p_01a`, we would expect to get the same numerical coordinate values used to construction the relative transform `p1b`.

Similarly, we check that the opposite corner of the walked rectangle is as we expect
```julia
x_012_ = get_coordinates(G, p0, log(G, p0, p_012), B)
@assert isapprox(     x_012_[1:2], [2;1], atol = 1e-10)
@assert isapprox( abs(x_012_[3]),      π, atol = 1e-10 )
```

And also the past point on the path `p_01234` lands back on the reference point `p0` as expected
```julia
x_01234_ = get_coordinates(G, p0, log(G, p0, p_01234), B)
@assert isapprox( x_01234_, [0;0;0], atol = 1e-10 )
```

This concludes the tutorial with emphasis on functions:
- `get_vector` / `hat`,
- `exp` vs. `log`,
- `get_coordinates` / `vee`, and
- `compose`.
