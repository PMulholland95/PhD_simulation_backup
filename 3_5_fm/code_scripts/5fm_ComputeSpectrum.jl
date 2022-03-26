function computeLinearSpectrum!(model::FluidModel{T, NF, NV, FG};
                               ) where {T <: AbstractFloat, NF, NV, FG <: AbstractFieldlineGeometry}
  kxRange = model.spectralGrid.kxRange
  kyRange = model.spectralGrid.kyRange
  kyStartIndex, ω_target = initialGuess(model)

  #rank = NF*length(first(model.indices))
  #lhs = spzeros(Complex{T},rank,rank)
  #rhs = spzeros(Complex{T},rank,rank)+I

  init = true
  for kyIndex in kyStartIndex:-1:1
    ky = kyRange[kyIndex]
    for kx in kxRange
      println("Computing linear solution at (kx,ky) = (",kx,",",ky,")")
      computeLinearMode!(kx,ky,model,init=init,scanDir=:down)
      if init
        init = false
      end
    end
  end

  scanDir = :up
  for kyIndex = kyStartIndex+1:1:nKyModes(model)
    ky = kyRange[kyIndex]
    for kx in kxRange
      println("Computing linear solution at (kx,ky) = (",kx,",",ky,")")
      computeLinearMode!(kx,ky,model,init=init,scanDir=:up)
    end
  end
end


function computeLinearMode!(kx::AbstractFloat,
                            ky::AbstractFloat,
                            model::FluidModel{T, 6, 1, FG};
                            init::Bool=false,
                            scanDir::Symbol=:none
                           ) where {T, FG <: AbstractFieldlineGeometry}

  kxIndex, kyIndex = index(kx,ky,model.spectralGrid)
  indexRange = model.indices[kxIndex,kyIndex]
  pointRange = model.points[kxIndex,kyIndex]
  geomView = @view(model.geometry[indexRange])
  npoints = length(pointRange)

  ω_target = getTarget(kxIndex,kyIndex,model,scanDir,init=init)

  if imag(ω_target) > zero(T)
    # Solve for the eigenmodes
    ωs, βs = linearSolver(kx,ky,model,target=ω_target,nev=5)

    # Sort the modes
    ω_index = sortModes(kx,ky,pointRange,ωs,βs,geomView,model.plasma,model.problem)
    #β_vec = βs[:,ω_index]./maximum(abs.(βs[:,ω_index]))

    # Update the spectrum
    updateSpectrum!(model,kxIndex,kyIndex,ωs[ω_index],βs[:,ω_index])
  else
    updateSpectrum!(model,kxIndex,kyIndex,zero(Complex{T}),nothing)
  end
end

function computeLinearMode!(lhs::AbstractSparseArray,
                            rhs::AbstractSparseArray,
                            kxIndex::Integer,
                            kyIndex::Integer,
                            model::FluidModel{T, 6, 1, FG};
                            init::Bool=false,
                            scanDir::Symbol=:none,
                           ) where {T <: AbstractFloat,
                                    VT <: AbstractVector{Complex{T}},
                                    FG <: AbstractFieldlineGeometry}
  indexRange = model.indices[kxIndex,kyIndex]
  pointRange = model.points[kxIndex,kyIndex]
  geomView = @view(model.geometry[indexRange])
  npoints = length(pointRange)
  kx = model.spectralGrid.kxRange[kxIndex]
  ky = model.spectralGrid.kyRange[kyIndex]

  ω_target = getTarget(kxIndex,kyIndex,model,scanDir,init=init)

  if imag(ω_target) > zero(T)
    # Solve for the eigenmodes
    ωs, βs = linearSolver!(lhs,rhs,kx,ky,model,target=ω_target,nev=5)

    # Sort the modes
    ω_index = sortModes(kx,ky,pointRange,ωs,βs,geomView,model.plasma,model.problem)
    #evec .= βs[:,ω_index]./maximum(abs.(βs[:,ω_index]))

    return ωs[ω_index], βs[:,ω_index]
  else
    return zero(Complex{T}), nothing
  end
end


function updateSpectrum!(model::FluidModel{T,NF,NV,FG},
                         kxIndex::Integer,
                         kyIndex::Integer,
                         ω::Complex{T},
                         β::Union{AbstractVector{Complex{T}},Nothing};
                        ) where {T <: AbstractFloat,
                                 NF, NV,
                                 FG <: AbstractFieldlineGeometry}

  indexRange = model.indices[kxIndex,kyIndex]
  pointRange = model.points[kxIndex,kyIndex]
  geomView = @view(model.geometry[indexRange])
  npoints = length(pointRange)
  kx = model.spectralGrid.kxRange[kxIndex]
  ky = model.spectralGrid.kyRange[kyIndex]

  if imag(ω) > zero(T)
    # This block could be used to make the interval smaller
    #=
    β_sum = Vector{T}(undef,length(β_view))
    cumsum!(β_sum,abs2.(β_view)/maximum(abs2.(β_view)))
    lowerIndex = findfirst(x->x>1e-8,β_sum)
    upperIndex = findlast(x->x<(1.0-1e-8),β_sum)
    model.indices[kxIndex,kyIndex] = model.indices[kxIndex,kyIndex][lowerIndex:upperIndex]
    model.points[kxIndex,kyIndex] = model.points[kxIndex,kyIndex][lowerIndex:upperIndex]
    indexRange = indexRange[lowerIndex:upperIndex]
    pointRange = pointRange[lowerIndex:upperIndex]
    geomView = @view(model.geometry[indexRange])
    npoints = length(pointRange)
    =#

    # Fit the spline
    evec = CubicSplineInterpolation(pointRange,β[1:npoints]./maximum(abs.(β[1:npoints])))
    normVal = vecNorm(pointRange,evec,geomView,model.problem)
    ω0 = marginalMode(kx,ky,ω,pointRange,evec,
                      geomView,model.plasma,model.problem)
    roots = SVector(ω,conj(ω),ω0)
    cMat, cMatInv = couplingMatrix(ky,roots,model.plasma)

    # Update the model
    model.eigenvalues[kxIndex,kyIndex] = roots
    #model.eigenvectors[kxIndex,kyIndex] = evec
    model.eigenvectors[kxIndex,kyIndex] = β[1:npoints]./maximum(abs.(β[1:npoints]))
    model.vecNorms[kxIndex,kyIndex] = normVal
    model.couplingMatrix[kxIndex,kyIndex] = cMat
    model.invCouplingMatrix[kxIndex,kyIndex] = cMatInv
    model.instability[kxIndex,kyIndex] = true
    if kxIndex > 1
      nkx = nKxModes(model.spectralGrid)
      kxIndex = 2*nkx - kxIndex + 1
      tempPointRange = model.points[kxIndex,kyIndex]
      pointRange = length(tempPointRange) == npoints ? tempPointRange : range(first(tempPointRange),step=step(tempPointRange),length=npoints)
      # For the symmetric case, the eigenvector is just reversed
      #model.eigenvectors[kxIndex,kyIndex] = CubicSplineInterpolation(pointRange,reverse(β[1:npoints]./maximum(abs.(β[1:npoints]))))
      model.eigenvectors[kxIndex,kyIndex] = reverse(β[1:npoints]./maximum(abs.(β[1:npoints])))
      model.vecNorms[kxIndex,kyIndex] = normVal
      model.eigenvalues[kxIndex,kyIndex] = roots
      model.couplingMatrix[kxIndex,kyIndex] = cMat
      model.invCouplingMatrix[kxIndex,kyIndex] = cMatInv
      model.instability[kxIndex,kyIndex] = true
    end
  else
    model.eigenvalues[kxIndex,kyIndex] = SVector{3,Complex{T}}(zeros(Complex{T},3))
    model.eigenvalues[2*nKxModes(model.spectralGrid)-kxIndex+1,kyIndex] = SVector{3,Complex{T}}(zeros(Complex{T},3))
  end
end

function getTarget(kxIndex::Integer,
                   kyIndex::Integer,
                   model::FluidModel{T,NF,NV,FG},
                   scanDir::Symbol=:down;
                   init::Bool=false
                  ) where {T,NF,NV,FG}

  if init
    indexRange = model.indices[kxIndex,kyIndex]
    pointRange = model.points[kxIndex,kyIndex]
    ι = model.geomParameters.ι
    nfp = model.geomParameters.nfp
    σ = nfp > 1 ? π/abs(ι-nfp) : 0.4π

    ϕ = Vector{Float64}(undef,length(pointRange))
    geomView = @view(model.geometry[indexRange])
    map!(x->gaussianMode(x,σ),ϕ,pointRange)
    vec = CubicSplineInterpolation(pointRange,ϕ)

    ω = computeLinearRoots(0.0,model.spectralGrid.kyRange[kyIndex],pointRange,
                           vec,geomView,model.plasma,model.problem)
    return convert(Complex{T},ω[1])
  else
    if scanDir === :none
      return zero(Complex{T})
    else
      index = CartesianIndex(kxIndex,kyIndex)
      prevIndex = index - CartesianIndex(1,0)
      if prevIndex.I[1] < 1
        prevIndex = scanDir == :down ? index + CartesianIndex(0,1) : index - CartesianIndex(0,1)
      end
      return model.eigenvalues[prevIndex][1]
    end
  end
end

function initialGuess(model::FluidModel{T, NF, NV, FG},
                     ) where {T <: AbstractFloat, NF, NV, FG <: AbstractFieldlineGeometry}
  kyRange = model.spectralGrid.kyRange
  medianKyIndex = convert(Int,iseven(length(kyRange)) ? median(1:length(kyRange)+1) : median(1:length(kyRange)))
  ι = model.geomParameters.ι
  nfp = model.geomParameters.nfp

  kxIndex, kyIndex = index(0.0,first(kyRange),model.spectralGrid)
  indexRange = model.indices[kxIndex,kyIndex]
  pointRange = model.points[kxIndex,kyIndex]
  geomView = @view(model.geometry[indexRange])
  ϕ = Vector{Float64}(undef,length(pointRange))

  σ = nfp > 1 ? π/abs(ι-nfp) : 0.4π
  map!(x->gaussianMode(x,σ),ϕ,pointRange)
  vec = CubicSplineInterpolation(pointRange,ϕ)

  for kyIndex = medianKyIndex:length(kyRange)
    ω = computeLinearRoots(0.0,kyRange[kyIndex],pointRange,
                           vec,geomView,model.plasma,model.problem)
    if imag(ω[1]) > 0
      return kyIndex, ω[1]
    end
  end

  for kyIndex = medianKyIndex-1:-1:1
    ω = computeLinearRoots(0.0,kyRange[kyIndex],pointRange,
                           vec,geomView,model.plasma,model.problem)
    if imag(ω[1]) > 0
      return kyIndex, ω[1]
    end
  end

  # If we are here, something has probably gone wrong, but try it anyway
  return medianKyIndex, one(T)*im
end

"""
    sortModes(kx,ky,target,evals,evecs,geom,model)

Sorts the modes given by `evals` and `evecs` according to a combintation of growth rate,
k∥² and possibly phase
"""
function sortModes(kx::AbstractFloat,
                   ky::AbstractFloat,
                   pointRange::AbstractRange,
                   evals::AbstractVector,
                   evecs::AbstractMatrix,
                   geom::StructArray{FG},
                   plasma::PlasmaParameters,
                   problem::LinearProblem;
                   weights::Vector{T}=[0.5,0.25,0.25],
                 ) where {T,FG <: AbstractFieldlineGeometry}
  nev = length(evals)
  size(evecs,2) == nev ||
    throw(DimensionMismatch("Number of eigenvalues does not match number of eigenvectors"))


  gammas = imag(evals)./maximum(imag(evals))
  kparallel = zeros(T,nev)
  phaseScore = zeros(T,nev)

  #targetDist = zeros(Float64,nev)

  #indexRange = intervalRangeIndices(kx,ky,model)
  #pointRange = model.geometry.points[indexRange]
  #geomView = @view(geom[indexRange])

  sqrtg = CubicSplineInterpolation(pointRange,geom.jac)
  Bhat = CubicSplineInterpolation(pointRange,geom.Bhat)

  for evIndex in eachindex(evals)
    # For the three-field model ϕ will always span the first N points
    ϕ_view = view(evecs,1:length(pointRange),evIndex)
    ϕ = CubicSplineInterpolation(pointRange,ϕ_view)
    kparallel[evIndex] = kparallel2(pointRange,ϕ,Bhat,sqrtg,problem)
    phase = crossPhase(view(evecs,:,evIndex),3,1,3)
    phaseScore[evIndex] = sign(phase)*(1.0 - abs(0.5π - abs(phase))/(0.5π))
  end
  kparallel .= 1.0 ./ (kparallel ./ minimum(kparallel))
  #for evIndex in eachindex(evals)
  #  println(evIndex,' ',weights[1]*gammas[evIndex],' ',weights[2]*kparallel[evIndex],' ',weights[3]*phaseScore[3])
  #end
  modeScore = (weights[1] .* gammas .+  weights[2] .* kparallel
               .+ weights[3] .* phaseScore)/sum(weights)
  return findmax(modeScore)[2]
end

#=
function is_boundary_mode(evec::AbstractVector{T};
                          atol::Float64=1e-6,
                         ) where T
  max_evec = maximum(abs.(evec))
  return abs(first(evec))/max_evec <= atol && abs(last(evec))/max_evec <= atol
end

function is_boundary_mode(evec::Interpolations.Extrapolation;
                          atol::Float64=1e-6,
                         )
  return is_boundary_mode(evec.itp.itp.coefs,atol=atol)
end
=#
