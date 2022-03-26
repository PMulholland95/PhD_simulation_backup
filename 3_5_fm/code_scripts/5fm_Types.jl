"""
    abstract type AbstractFluid

The abstract type to represent a fluid model
"""
abstract type AbstractFluidModel end

"""
    abstract type AbstractGeometryParameters

The abstract type to represent the geometry parameters that are input
for a fluid model calcluation
"""
abstract type AbstractGeometryParameters end

"""
    abstract type AbstractFieldlineGeometry

The abstract type to represent the geometry coefficients along a
fieldline for use in a fluid model calculation
"""
abstract type AbstractFieldlineGeometry end

"""
    abstract type SolverMethod

      Abstract type to repesent an iterative method algorithm

      # Available Algorithms
      # - `JacobiDavidson()` : Complex Jacobi-Davidson method
"""
abstract type SolverMethod end

"""
    NullParameters

Singleton struct representing no geometry parameters
"""
struct NullParameters <: AbstractGeometryParameters end


"""
    PlasmaParameters

Composite type holding the plasma parameters for a fluid model calculations,
current fields are `∇T`, `∇n`, `τ` and `β`

"""
struct PlasmaParameters
  # Normalized temperature gradient: ∇T = -L_ref/T dT/dr"
  ∇Tᵢ::AbstractFloat

  ∇Tₑ::AbstractFloat

  # Normalized density gradient: ∇n = -L_ref/n dn/dr
  ∇n::AbstractFloat

  # Ion-to-electron temperature ration
  τ::AbstractFloat

  # Normalized ratio of plasma pressure to magnetic pressure
  β::AbstractFloat
end

"""
    PlasmaParameters(∇T::Real,∇n::Real,τ::Real=0,β::Real=0)

Constructor for `PlasmaParameters` type taking real arguments
"""
function PlasmaParameters(∇Tᵢ::Real,
			  ∇Tₑ::Real,
                          ∇n::Real;
                          τ::Real=1.0,
                          β::Real=0.0,
                         )
  return PlasmaParameters(convert(Float64,∇Tᵢ),
			  convert(Float64,∇Tₑ),
                          convert(Float64,∇n),
                          convert(Float64,τ),
                          convert(Float64,β))
end



"""
    LinearProblem

Composite type holding the properties for solving the linear drift wave problem
"""
struct LinearProblem
  solver::SolverMethod
  solnInterval::AbstractFloat
  ε::AbstractFloat
  derivativeErrorOrder::Integer
  hyperdiffusionOrder::Integer
  quadNodes::AbstractVector{AbstractFloat}
  quadWeights::AbstractVector{AbstractFloat}
end


"""
    LinearProblem(; kwargs...)

# Arguments
- `solver::SolverMethod`=JDQZ() : Solver method used for linear solver
- `solnInterval::AbstractFloat`=8π : Interval in radians over which to compute linear solutions
- `ε::AbstractFloat`=200.0 : Numerical hyperdiffusion parameter
- `errorOrder::Integer`=4 : Error order for finite difference discretization
- `hyperdiffusionOrder::Integer`=2 : Exponent of hyperdiffusion operator ``(iε)ⁿ∂ⁿ/∂ηⁿ``
- `quadPoints::Integer`=1000 : Number of points using in Gauss-Legendre quadrature

# See also: (`SolverMethod`)[@ref]
"""
function LinearProblem(;
                       solver=JDQZ(),
                       solnInterval=8π,
                       ε=200.0,
                       errorOrder=4,
                       hyperdiffusionOrder=2,
                       quadPoints=1000,
                      )
  nodes, weights = gausslegendre(quadPoints)
  return LinearProblem(solver,solnInterval,
                       ε,errorOrder,hyperdiffusionOrder,
                       nodes,weights)
end
