struct JDQZ <: SolverMethod end
struct direct <: SolverMethod end

"""
    linearSolver(kx::AbstractFloat,
                 ky::AbstractFloat,
                 model::FluidModel{T,NF,NV,FG};
                 target=0.0+im,
                 nev=3,
                 kwargs...
                ) where {T,NF,NV,FG <: AbstractFieldlineGeometry}

Apply the linear solution method specified by `solver` to compute `nev` eigenvalues/eigenmodes
for the spectral point specified by `(kx,ky)` using a target point specified by `target`,
the meaning of `target` is dependent on algorithm type. Solver specific arguments can be
passed through `kwargs`.
"""
function linearSolver(kx::AbstractFloat,
                      ky::AbstractFloat,
                      model::FluidModel{T,NF,NV,FG};
                      target::Complex{T}=0.0+1.0im,
                      kwargs...
                     ) where {T <: AbstractFloat,NF,NV,FG <: AbstractFieldlineGeometry}
  # Iterative methods typically require sparse matrices
  # Currently just apply sparse() to the result of denseLinearOperators
  A, B = sparse.(denseLinearOperators(kx,ky,model))
  evals, evecs = _linearSolver(model.problem.solver,A,B,target;kwargs...)

  # This is a simple rule to try to select the correct eigenmode
  # and should possibly be moved out of this function
  #evIndex = findmin(norm.(evals .- target,2))[2]
  #return evals[evIndex], evecs[:,evIndex]
  return evals, evecs
end

function linearSolver!(lhs::AbstractArray,
                       rhs::AbstractArray,
                       kx::AbstractFloat,
                       ky::AbstractFloat,
                       model::FluidModel{T,NF,NV,FG};
                       target::Complex{T}=one(T)*im,
                       kwargs...
                      ) where {T <: AbstractFloat,
                               NF, NV,
                               FG <: AbstractFieldlineGeometry}
  return _linearSolver(model.problem.solver,lhs,rhs,target;kwargs...)
end

function _linearSolver(::JDQZ,
                       lhs::AbstractSparseMatrix{Complex{T},I},
                       rhs::AbstractSparseMatrix{Complex{T},I},
                       target::Complex{T}=0.0+1.0im;
                       kwargs...
                      ) where {T <: AbstractFloat,I}
  probSize = size(lhs,1)
  pschur, res = jdqz(lhs,rhs,
                     solver=haskey(kwargs,:solver) ? kwargs[:solver] : GMRES(probSize,iterations=7),
                     target= Near(real(target)+10im*imag(target)),
                     pairs=haskey(kwargs,:nev) ? kwargs[:nev] : 3,
                     subspace_dimensions = haskey(kwargs,:subspace_dimensions) ? kwargs[:subspace_dimensions] : 50:200,
                     max_iter = haskey(kwargs,:max_iter) ? kwargs[:max_iter] : 1000,
                     tolerance = haskey(kwargs,:tolerance) ? convert(Float64,kwargs[:tolerance]) : convert(Float64,sqrt(eps(real(eltype(lhs))))),
                     numT = Complex{T},
                     verbosity=1);

  Q = pschur.Q.basis
  Z = pschur.Z.basis
  evals, evecs = eigen(Z'*lhs*Q,Z'*rhs*Q)
  return evals, Q*evecs
end

#=
function computeConditionNumber(kx::Float64,ky::Float64,geom::FieldlineGeometry;errorOrder=2)
  A, B = denseLinearOperators(kx,ky,geom,errorOrder=errorOrder)
  return cond(inv(B)*A)
end

function computeConditionNumber(kx::AbstractRange,ky::AbstractRange,geom::FieldlineGeometry;errorOrder=2)
  conditionNumbers = Array{Float64}(undef,length(kx),length(ky))
  for (iy,iky) in enumerate(ky)
    for (ix,ikx) in enumerate(kx)
      conditionNumbers[ix,iy] = computeConditionNumber(ikx,iky,geom,errorOrder=errorOrder)
    end
  end
  return conditionNumbers
end

=#
#function initialValueSolver(kx::Float64,ky::Float64,geom::FieldlineGeometry)
#
#end

function computeLinearRoots!(kx::AbstractFloat,
                             ky::AbstractFloat,
                             model::FluidModel{T,6,1,FG};
                            ) where {T, FG <: AbstractFieldlineGeometry}
  kxIndex, kyIndex = index(kx,ky,model.spectralGrid)
  indexRange = model.indices[kxIndex,kyIndex]
  pointRange = model.points[kxIndex,kyIndex]
  roots = computeLinearRoots(kx,ky,pointRange,model.eigenvectors[1,kxIndex,kyIndex],
                             @view(model.geometry[indexRange]),model.plasma,model.problem)
  model.eigenvalues[kxIndex,kyIndex] = roots
end

"""
    computeLinearRoots(kx::AbstractFloat,
                       ky::AbstractFloat,
                       pointRange::AbstractRange,
                       vec::Interpolations.Extrapolation,
                       plasma::PlasmaParameters,
                       geom::StructArray{FG},
                       problem::LinearProblem,
                      ) where {FG <: AbstractFieldlineGeometry}

For a given eigenfunction specfied by `vec` and geometry interval, compute the unstable,
stable and marginally stable solutions, this can provide a check on the direct
eigenvalue calculation.
"""
function computeLinearRoots(kx::AbstractFloat,
                            ky::AbstractFloat,
                            pointRange::AbstractRange,
                            vec::Interpolations.Extrapolation,
                            geom::StructArray{FG},
                            plasma::PlasmaParameters,
                            problem::LinearProblem,
                           ) where {FG <: AbstractFieldlineGeometry}

  #indexRange = intervalRangeIndices(kx,ky,model)
  #geomView = @view(geom[indexRange])
  Bk_vec, Dk_vec = buildBkDk(kx,ky,geom)
  #pointRange = model.geometry.points[indexRange]
  firstPoint = first(pointRange)
  lastPoint = last(pointRange)

  Bk = CubicSplineInterpolation(pointRange,Bk_vec)
  Dk = CubicSplineInterpolation(pointRange,Dk_vec)
  sqrtg = CubicSplineInterpolation(pointRange,geom.jac)
  B = CubicSplineInterpolation(pointRange,geom.Bhat)
  #ϕ = CubicSplineInterpolation(pointRange,vec)

  # Gauss-Legendre quadrature with predefined knots and weights is
  # fast, easiest to do sequentially rather than parallel
  L3 = quadgl1D(L3_integrand,problem.quadNodes,problem.quadWeights,
                firstPoint,lastPoint,vec,B,sqrtg,Bk,plasma;quadT=Float64)

  L2 = quadgl1D(L2_integrand,problem.quadNodes,problem.quadWeights,
                firstPoint,lastPoint,ky,vec,B,sqrtg,Bk,Dk,plasma;quadT=Float64)

  L1 = quadgl1D(L1_integrand,problem.quadNodes,problem.quadWeights,
                firstPoint,lastPoint,ky,vec,B,sqrtg,Bk,Dk,plasma;quadT=Float64)

  L0 = quadgl1D(L0_integrand,problem.quadNodes,problem.quadWeights,
                firstPoint,lastPoint,ky,vec,B,sqrtg,plasma;quadT=Float64)

  # Use the Polynomials package to find the roots of the cubic equation
  # ω³ + (L₂/L₃)ω² + (L₁/L₃)ω + L₀/L₁ = 0
  roots = Polynomials.roots(Polynomial([L0/L3,L1/L3,L2/L3,1.0]))

  # Sort the modes, it may be more efficient to move this outside the function
  unstableMode = findfirst(x->imag(x)>0,roots)
  stableMode = findfirst(x->imag(x)<0,roots)
  marginalMode = findfirst(x->isapprox(imag(x),0.0),roots)
  return SVector{3,ComplexF64}(roots[!isnothing(unstableMode) ? unstableMode : 1],
                               roots[!isnothing(stableMode) ? stableMode : 2],
                               roots[!isnothing(marginalMode) ? marginalMode : 3])
end

function marginalMode(kx::AbstractFloat,
                      ky::AbstractFloat,
                      model::FluidModel{T,6,1,FG};
                     ) where {T, FG <: AbstractFieldlineGeometry}
  kxIndex, kyIndex = index(kx,ky,model.spectralGrid)
  indexRange = model.indices[kxIndex,kyIndex]
  geomView = @view(model.geometry[indexRange])
  pointRange = model.points[kxIndex,kyIndex]

  return marginalMode(kx,ky,model.eigenvalues[kxIndex,kyIndex][1],
                      pointRange,model.eigenvectors[1,kxIndex,kyIndex],
                      geomView,model.plasma,model.problem,quadT=T)
end

"""
    marginalMode(kx::AbstractFloat,
                 ky::AbstractFloat,
                 ω::Complex{T},
                 pointRange::AbstractRange,
                 vec::Interpolations.Extrapolation,
                 geom::StructArray{FG},
                 plasma::PlasmaParameters
                 problem::LinearProblem,
                ) where {T, FG <: AbstractFieldlineGeometry}

For a given eigenfunction with unstable eigenvalue `ω` and inteprolated eigenmode `vec`
over the geometry interval `geom`, compute only the marginally stable solution,
which can also provide a check on `computeLinearRoots()`

# See also: [`computeLinearRoots`](@ref)
"""
function marginalMode(kx::AbstractFloat,
                      ky::AbstractFloat,
                      ω::Complex{T},
                      pointRange::AbstractRange,
                      vec::Interpolations.Extrapolation,
                      geom::StructArray{FG},
                      plasma::PlasmaParameters,
                      problem::LinearProblem;
                      quadT::Type=Float64,
                     ) where {T,FG <: AbstractFieldlineGeometry}
  Bk_vec, Dk_vec = buildBkDk(kx,ky,geom)
  firstPoint = first(pointRange)
  lastPoint = last(pointRange)

  Bk = CubicSplineInterpolation(pointRange,Bk_vec)
  Dk = CubicSplineInterpolation(pointRange,Dk_vec)
  sqrtg = CubicSplineInterpolation(pointRange,geom.jac)
  B = CubicSplineInterpolation(pointRange,geom.Bhat)

  L3 = quadgl1D(L3_integrand,problem.quadNodes,problem.quadWeights,
                firstPoint,lastPoint,vec,B,sqrtg,Bk,plasma;quadT=quadT)
  L2 = quadgl1D(L2_integrand,problem.quadNodes,problem.quadWeights,
                firstPoint,lastPoint,ky,vec,B,sqrtg,Bk,Dk,plasma;quadT=quadT)

  return -L2/L3 - 2*real(ω)
end

function L3_integrand(x::T,
                      ϕ::Interpolations.Extrapolation,
                      B::Interpolations.Extrapolation,
                      sqrtg::Interpolations.Extrapolation,
                      Bk::Interpolations.Extrapolation,
                      plasma::PlasmaParameters,
                    ) where T
  return sqrtg(x)*B(x)*(1+Bk(x)*(1+5/3*plasma.τ))*abs2(ϕ(x))
end

function L2_integrand(x::T,
                      ky::T,
                      ϕ::Interpolations.Extrapolation,
                      B::Interpolations.Extrapolation,
                      sqrtg::Interpolations.Extrapolation,
                      Bk::Interpolations.Extrapolation,
                      Dk::Interpolations.Extrapolation,
                      plasma::PlasmaParameters,
                    ) where T
  return abs2(ϕ(x))*sqrtg(x)*B(x)*(ky*
                                    ((plasma.∇Tᵢ-2/3*plasma.∇n)*plasma.τ*Bk(x) - plasma.∇n)/B(x)
                                    -Dk(x)*(1+5/3*plasma.τ))
end

function L1_integrand(x::T,
                      ky::T,
                      ϕ::Interpolations.Extrapolation,
                      B::Interpolations.Extrapolation,
                      sqrtg::Interpolations.Extrapolation,
                      Bk::Interpolations.Extrapolation,
                      Dk::Interpolations.Extrapolation,
                      plasma::PlasmaParameters,
                    ) where T
  return -sqrtg(x)*B(x)*(ky*Dk(x)*(plasma.∇Tᵢ-2/3*plasma.∇n)/B(x)*plasma.τ)*abs2(ϕ(x)) - (1+5/3*plasma.τ)*abs2(Interpolations.gradient(ϕ,x)[])/(sqrtg(x)*B(x))
end

function L0_integrand(x::T,
                      ky::T,
                      ϕ::Interpolations.Extrapolation,
                      B::Interpolations.Extrapolation,
                      sqrtg::Interpolations.Extrapolation,
                      plasma::PlasmaParameters,
                    ) where T
  return -ky*(plasma.∇Tᵢ-2/3*plasma.∇n)/B(x)*plasma.τ*abs2(Interpolations.gradient(ϕ,x)[])/(sqrtg(x)*B(x))
end

function kparallel2(kx::AbstractFloat,
                    ky::AbstractFloat,
                    model::FluidModel{T,6,1,FG};
                   ) where {T, FG <: AbstractFieldlineGeometry}
  kxIndex, kyIndex = index(kx,ky,model.spectralGrid)
  indexRange = model.indices[kxIndex,kyIndex]
  geomView = @view(model.geometry[indexRange])
  pointRange = model.points[kxIndex,kyIndex]
  return kparallel2(pointRange,model.eigenvectors[kxIndex,kyIndex],
                    geomView,model.problem;quadT=T)
end

"""
    kparallel2(pointRange::AbstractRange,
               vec::Interpolations.Extrapolation,
               geom::StructArray{FG},
               problem::LinearProblem;
               quadT::Type=Float64
              ) where {FG <: AbstractFieldlineGeometry}

For a given eigenfunction specfied by `vec` and geometry interval, compute the
effective paralle wavenumber squared: ⟨k∥²⟩
"""
function kparallel2(pointRange::AbstractRange,
                    vec::Interpolations.Extrapolation,
                    geom::StructArray{FG},
                    problem::LinearProblem;
                    quadT::Type=Float64,
                   ) where FG <: AbstractFieldlineGeometry
  sqrtg = CubicSplineInterpolation(pointRange,geom.jac)
  B = CubicSplineInterpolation(pointRange,geom.Bhat)

  return kparallel2(pointRange,vec,B,sqrtg,problem;quadT=quadT)
end

function kparallel2(pointRange::AbstractRange,
                    vec::Interpolations.Extrapolation,
                    B::Interpolations.Extrapolation,
                    sqrtg::Interpolations.Extrapolation,
                    problem::LinearProblem;
                    quadT::Type=Float64,
                   )
  ϕ_norm = quadgl1D(fieldlineNorm,problem.quadNodes,problem.quadWeights,
                    first(pointRange),last(pointRange),vec,B,sqrtg;quadT=quadT)
  dϕ_avg = quadgl1D(kparallel2_integrand,problem.quadNodes,problem.quadWeights,
                    first(pointRange),last(pointRange),vec,B,sqrtg;quadT=quadT)
  return dϕ_avg/ϕ_norm
end

function fieldlineNorm(x::T,
                       f::Interpolations.Extrapolation,
                       B::Interpolations.Extrapolation,
                       sqrtg::Interpolations.Extrapolation,
                      ) where T
  return sqrtg(x)*B(x)*abs2(f(x))
end

function kparallel2_integrand(x::T,
                              ϕ::Interpolations.Extrapolation,
                              B::Interpolations.Extrapolation,
                              sqrtg::Interpolations.Extrapolation,
                             ) where T
  return abs2(Interpolations.gradient(ϕ,x)[])/(sqrtg(x)*B(x))
end
