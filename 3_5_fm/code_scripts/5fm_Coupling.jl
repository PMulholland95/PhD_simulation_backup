
function overlapInterval(range1::AbstractRange,range2::AbstractRange)
  if first(range1) >= last(range2)
    return nothing
  elseif last(range1) <= first(range2)
    return nothing
  else
    return range(max(first(range1),first(range2)),
                 step=step(range1),
                 stop=min(last(range1),last(range2)))
  end
end

function kxTriplet(kx1::Integer,
                   kx2::Integer,
                   nkx::Integer;
                  )
  return abs(kx1 - kx2) <= nkx-1 ? kx2 > kx1 ? 2*nkx - (kx2 - kx1) : kx2 < kx1 ? kx1 - kx2 + 1 : 0 : 0
end

function kxTriplet(kx1::AbstractFloat,
                   kx2::AbstractFloat,
                   maxKx::AbstractFloat;
                  )
  return abs(kx1 - kx2) <= maxKx ? kx1 - kx2 : nothing
end

function kyTriplet(ky1::Integer,
                   ky2::Integer,
                   nky::Integer;
                  )
  return abs(ky1-ky2) <= nky ? abs(ky1-ky2) : 0
end

function kyTriplet(ky1::AbstractFloat,
                   ky2::AbstractFloat,
                   maxKy::AbstractFloat;
                  )
  return abs(ky1-ky2) <= maxKy ? ky1 - ky2 : nothing
end

function overlapIntegrand(x::AbstractFloat,
                          u::Interpolations.Extrapolation,
                          v::Interpolations.Extrapolation,
                          sqrtg::Interpolations.Extrapolation,
                          modB::Interpolations.Extrapolation;
                         )
  return sqrtg(x)*modB(x)*abs(u(x))*abs(v(x))
end

function overlapIntegrand(x::AbstractFloat,
                          u::Interpolations.Extrapolation,
                          v::Interpolations.Extrapolation,
                          w::Interpolations.Extrapolation,
                          sqrtg::Interpolations.Extrapolation,
                          modB::Interpolations.Extrapolation;
                         )
  return sqrtg(x)*modB(x)*abs(u(x))*abs(v(x))*abs(w(x))
end

function overlapIntegral(kx1::AbstractFloat,
                         ky1::AbstractFloat,
                         kx2::AbstractFloat,
                         ky2::AbstractFloat,
                         kx3::AbstractFloat,
                         ky3::AbstractFloat,
                         model::FluidModel{T,NF,NV,FG};
                        ) where {T <: AbstractFloat,
                                 NF, NV,
                                 FG <: AbstractFieldlineGeometry}
  ikx1 = getKxIndex(kx1,model.spectralGrid)
  iky1 = getKyIndex(ky1,model.spectralGrid)
  ikx2 = getKxIndex(kx2,model.spectralGrid)
  iky2 = getKyIndex(ky2,model.spectralGrid)
  ikx3 = getKxIndex(kx3,model.spectralGrid)
  iky3 = getKyIndex(ky3,model.spectralGrid)
  return overlapIntegral(ikx1,iky1,ikx2,iky2,ikx3,iky3,model)
end

function overlapIntegral(kx1::Integer,
                         ky1::Integer,
                         kx2::Integer,
                         ky2::Integer,
                         kx3::Integer,
                         ky3::Integer,
                         model::FluidModel{T,NF,NV,FG};
                        ) where {T <: AbstractFloat,
                                 NF, NV,
                                 FG <: AbstractFieldlineGeometry}
  overlapIndices12 = overlapInterval(model.indices[kx1,ky1],model.indices[kx2,ky2])
  overlapIndices23 = overlapInterval(model.indices[kx2,ky2],model.indices[kx3,ky3])
  overlapIndices = (!isnothing(overlapIndices12) && !isnothing(overlapIndices23)) ?
                   overlapInterval(overlapIndices12,overlapIndices23) : nothing

  if !isnothing(overlapIndices)
    #overlapPoints = overlapInterval(model.points[kx1,ky1],model.points[kx2,ky2])
    overlapPoints = model.geomParameters.points[overlapIndices]
    geom = @view(model.geometry[overlapIndices])

    u = CubicSplineInterpolation(model.points[kx1,ky1],conj.(model.eigenvectors[kx1,ky1]))
    v = CubicSplineInterpolation(model.points[kx2,ky2],model.eigenvectors[kx2,ky2])
    w = CubicSplineInterpolation(model.points[kx3,ky3],model.eigenvectors[kx3,ky3])

    return overlapIntegral(overlapPoints,u,v,w,geom,model.problem;quadT=T)/((model.vecNorms[kx1,ky1]*model.vecNorms[kx2,ky2]*model.vecNorms[kx3,ky3])^(1/3))

  else
    return zero(Complex{T})
  end
end

function overlapIntegral(overlapPoints::AbstractRange,
                         u::Interpolations.Extrapolation,
                         v::Interpolations.Extrapolation,
                         w::Interpolations.Extrapolation,
                         geom::StructArray{FG},
                         problem::LinearProblem;
                         quadT::Type=T
                        ) where {T, FG <: AbstractFieldlineGeometry}

  modB = CubicSplineInterpolation(overlapPoints,geom.Bhat)
  sqrtg = CubicSplineInterpolation(overlapPoints,geom.jac)
  return quadgl1D(overlapIntegrand,problem.quadNodes,problem.quadWeights,
                  first(overlapPoints),last(overlapPoints),u,v,w,sqrtg,modB,quadT=quadT)
end

function overlapIntegral(kx1::AbstractFloat,
                         ky1::AbstractFloat,
                         kx2::AbstractFloat,
                         ky2::AbstractFloat,
                         model::FluidModel{T,NF,NV,FG};
                        ) where {T <: AbstractFloat,
                                 NF, NV,
                                 FG <: AbstractFieldlineGeometry}
  ikx1 = getKxIndex(kx1)
  iky1 = getKyIndex(ky1)
  ikx2 = getKxIndex(kx2)
  iky2 = getKyIndex(ky2)
  return overlapIntegral(ikx1,iky1,ikx2,iky2,model)
end

function overlapIntegral(kx1::Integer,
                         ky1::Integer,
                         kx2::Integer,
                         ky2::Integer,
                         model::FluidModel{T,NF,NV,FG};
                        ) where {T <: AbstractFloat,
                                 NF, NV,
                                 FG <: AbstractFieldlineGeometry}
  overlapIndices = overlapInterval(model.indices[kx1,ky1],model.indices[kx2,ky2])

  if !isnothing(overlapIndices)
    #overlapPoints = overlapInterval(model.points[kx1,ky1],model.points[kx2,ky2])
    overlapPoints = model.geomParameters.points[overlapIndices]
    geom = @view(model.geometry[overlapIndices])

    u = CubicSplineInterpolation(model.points[kx1,ky1],model.eigenvectors[kx1,ky1])
    v = CubicSplineInterpolation(model.points[kx2,ky2],model.eigenvectors[kx2,ky2])

    return overlapIntegral(overlapPoints,u,v,geom,model.problem;quadT=T)/(sqrt(model.vecNorms[kx1,ky1])*sqrt(model.vecNorms[kx2,ky2]))

  else
    return zero(Complex{T})
  end
end

function overlapIntegral(overlapPoints::AbstractRange,
                         u::Interpolations.Extrapolation,
                         v::Interpolations.Extrapolation,
                         geom::StructArray{FG},
                         problem::LinearProblem;
                         quadT::Type=T
                        ) where {T, FG <: AbstractFieldlineGeometry}

  modB = CubicSplineInterpolation(overlapPoints,geom.Bhat)
  sqrtg = CubicSplineInterpolation(overlapPoints,geom.jac)
  return quadgl1D(overlapIntegrand,problem.quadNodes,problem.quadWeights,
                  first(overlapPoints),last(overlapPoints),u,v,sqrtg,modB,quadT=quadT)
end

function couplingMatrix(ky::T,
                        eigenvalues::SVector{3,Complex{T}},
                        plasma::PlasmaParameters,
                       ) where {T <: AbstractFloat}
  ηi = plasma.∇Tᵢ - 2/3*plasma.∇n
  τ = plasma.τ
  M11 = 1.0
  M21 = ky*ηi/eigenvalues[1]
  M22 = ky*ηi/eigenvalues[2]
  M23 = ky*ηi/eigenvalues[3]
  M31 = im/eigenvalues[1]*(1+5/3*τ+τ*ky*ηi/eigenvalues[1])
  M32 = im/eigenvalues[2]*(1+5/3*τ+τ*ky*ηi/eigenvalues[2])
  M33 = im/eigenvalues[3]*(1+5/3*τ+τ*ky*ηi/eigenvalues[3])
  massMat = @SMatrix [M11 M11 M11;
                      M21 M22 M23;
                      M31 M32 M33]

  return massMat, inv(massMat)
end

function computeEigenvectors(invCouplingMatrix::AbstractMatrix,
                             fieldVector::AbstractVector{Complex{T}};
                            ) where {T <: AbstractFloat}
  size(invCouplingMatrix,1) == size(invCouplingMatrix,2) || throw(DimensionMismatch("Inverse coupling matrix must be square!"))
  nFields = size(invCouplingMatrix,1)
  nPoints = div(length(fieldVector),nFields)
  eigenvectors = Array{Complex{T},2}(undef,nPoints,nFields)

  for i = 1:nPoints
    for j = 1:nFields
      for k = 1:nFields
        eigenvectors[i,j] = invCouplingMatrix[j,k]*fieldVector[nPoints*(k-1)+i] 
      end
    end
  end

  return eigenvectors
end

function computeCouplings(model::FluidModel{T,6,1,FG};
                         ) where {T <: AbstractFloat,
                                  FG <: AbstractFieldlineGeometry}
  nkx = nKxModes(model.spectralGrid)
  nky = nKyModes(model.spectralGrid)
  kxRange = model.spectralGrid.kxRange
  kyRange = model.spectralGrid.kyRange
  tauMatrix_ZF = zeros(Complex{T},nkx,nky)
  tauMatrix_NZF = similar(tauMatrix_ZF)
  coeffMatrix_ZF = similar(tauMatrix_ZF)
  coeffMatrix_NZF = similar(tauMatrix_ZF)
  maxKMatrix = Matrix{Tuple{Int,Int}}(undef,nkx,nky)

  ω = model.eigenvalues

  for ky1 in kyRange
    iky1 = getKyIndex(ky1,model.spectralGrid)
    for ky2 in kyRange
      ky3 = kyTriplet(ky1,ky2,last(kyRange))
      iky2 = getKyIndex(ky2,model.spectralGrid)
      iky3 = getKyIndex(ky3,model.spectralGrid)
      for kx1 in kxRange
        ikx1 = getKxIndex(kx1,model.spectralGrid)
        tau_zf = zero(Complex{T})
        tau_nzf = zero(Complex{T})
        for kx2 in kxRange
          if !isapprox((kx1*ky2 - kx2*ky1),zero(T))
            kx3 = kxTriplet(kx1,kx2,last(kxRange))
            ikx2 = getKxIndex(kx2,model.spectralGrid)
            ikx3 = getKxIndex(kx3,model.spectralGrid)
            if !isnothing(kx3) && !isnothing(ky3) && !isnothing(ikx3) && !isnothing(iky3)
              println("kx: ",kx1,' ',ikx1,' ',kx2,' ',ikx2,' ',kx3,' ',ikx3)
              println("ky: ",ky1,' ',iky1,' ',ky2,' ',iky2,' ',ky3,' ',iky3)
              τ = -im/(ω[ikx3,iky3][3] + ω[ikx2,iky2][2] - conj(ω[ikx1,iky1][1]))
              println(τ,' ',tau_nzf)
              tau_nzf = abs(τ) > abs(tau_nzf) ? τ : tau_nzf
            end
          end
        end
        if (iky1 != iky2)
          println(ikx1,' ',iky1,' ',tau_nzf)
          tauMatrix_NZF[ikx1,iky1] = tau_nzf
        end
      end
    end
  end
  return tauMatrix_NZF
end



function threeWaveCoupling(model::FluidModel{T,6,1,FG};
                           only_max::Bool=true,
                          ) where {T <: AbstractFloat,
                                   FG <: AbstractFieldlineGeometry}
  nkx = nKxModes(model.spectralGrid)
  nky = nKyModes(model.spectralGrid)
  kxRange = model.spectralGrid.kxRange
  kyRange = model.spectralGrid.kyRange
  maxKx = last(kxRange)
  maxKy = last(kyRange)
  negKxRange = reverse(-kxRange)


  tauMatrix_ZF = zeros(Complex{T},nkx,nky)
  tauMatrix_NZF = zeros(Complex{T},nkx,nky)
  coeffMatrix_ZF = zeros(Complex{T},nkx,nky)
  coeffMatrix_NZF = zeros(Complex{T},nkx,nky)

  workTau_ZF = zeros(Complex{T},size(model.eigenvalues,1))
  workTau_NZF = zeros(Complex{T},size(model.eigenvalues))
  workCoeffs_ZF = zeros(Complex{T},size(model.eigenvalues,1))
  workCoeffs_NZF = zeros(Complex{T},size(model.eigenvalues))

  @batch for I in CartesianIndices(tauMatrix_ZF)
    τ_ZF, C_ZF, τ_NZF, C_NZF = threeWaveCoupling!(workTau_ZF,workTau_NZF,
                                                  workCoeffs_ZF,workCoeffs_NZF,
                                                  I[1],I[2],model)
    tauMatrix_ZF[I] = τ_ZF
    coeffMatrix_ZF[I] = C_ZF
    tauMatrix_NZF[I] = τ_NZF
    coeffMatrix_NZF[I] = C_NZF
    fill!(workTau_ZF,zero(Complex{T}))
    fill!(workCoeffs_ZF,zero(Complex{T}))
    fill!(workTau_NZF,zero(Complex{T}))
    fill!(workCoeffs_NZF,zero(Complex{T}))
  end

  return tauMatrix_ZF, coeffMatrix_ZF, tauMatrix_NZF, coeffMatrix_NZF
end

function threeWaveCoupling!(tauVector_ZF::AbstractVector,
                            tauMatrix_NZF::AbstractArray,
                            coeffVector_ZF::AbstractVector,
                            coeffMatrix_NZF::AbstractArray,
                            ikx1::Integer,
                            iky1::Integer,
                            model::FluidModel{T,NF,NV,FG};
                           ) where {T <: AbstractFloat,
                                    NF, NV,
                                    FG <: AbstractFieldlineGeometry}

  kxRange = model.spectralGrid.kxRange
  kyRange = model.spectralGrid.kyRange
  maxKx = last(kxRange)
  maxKy = last(kyRange)
  kx1 = kxRange[ikx1]
  ky1 = kyRange[iky1]

  ω = model.eigenvalues
  couplingMatrix = model.couplingMatrix
  invCouplingMatrix = model.invCouplingMatrix
  sgrid = model.spectralGrid

  for ky2 in kyRange
    if ky2 == ky1
      iky2 = getKyIndex(ky2,model.spectralGrid)
      for kx2 in -maxKx+kx1:step(kxRange):maxKx
        kx3 = kxTriplet(kx1,kx2,maxKx)
        ikx2 = getKxIndex(kx2,model.spectralGrid)
        if !isnothing(kx3) && !isnothing(ikx2) && kx2 != kx1 && kx2 != -kx1
          τ = -im/(ω[ikx2,iky2][2] - conj(ω[ikx1,iky1][1]))
          C_pqF = ky1*(kx1-kx2) *(invCouplingMatrix[ikx1,iky1][1,2]*
                                  couplingMatrix[ikx2,iky2][2,2] +
                                  invCouplingMatrix[ikx1,iky1][1,3]*
                                  couplingMatrix[ikx2,iky2][3,2])
          overlap = overlapIntegral(ikx1,iky1,ikx2,iky2,model)
          tauVector_ZF[ikx2] = τ
          coeffVector_ZF[ikx2] = C_pqF*kxSpec(kx1,sgrid)*kxSpec(kx2,sgrid)*kxSpec(kx2,sgrid)*kySpec(ky1,sgrid)*kySpec(ky2,sgrid)
        end
      end
    else
      ky3 = kyTriplet(ky1,ky2,maxKy)
      iky2 = getKyIndex(ky2,model.spectralGrid)
      iky3 = getKyIndex(ky3,model.spectralGrid)
      for kx2 in -maxKx+kx1:step(kxRange):maxKx
        if !isapprox((kx1*ky2 - kx2*ky1),zero(T))
          kx3 = kxTriplet(kx1,kx2,maxKx)
          ikx2 = getKxIndex(kx2,model.spectralGrid)
          ikx3 = getKxIndex(kx3,model.spectralGrid)
          if !isnothing(ikx2) && !isnothing(kx3) && !isnothing(ky3) && !isnothing(ikx3) && !isnothing(iky3)
            #println("kx: ",kx1,' ',ikx1,' ',kx2,' ',ikx2,' ',kx3,' ',ikx3)
            #println("ky: ",ky1,' ',iky1,' ',ky2,' ',iky2,' ',ky3,' ',iky3)
            τ = -im/(ω[ikx3,iky3][3] + ω[ikx2,iky2][2] - conj(ω[ikx1,iky1][1]))
            C_pqr = (kx1*ky2 - ky1*kx2)*(invCouplingMatrix[ikx1,iky1][1,2]*
                                         couplingMatrix[ikx2,iky2][2,2] +
                                         0.5*invCouplingMatrix[ikx1,iky1][1,3]*
                                         couplingMatrix[ikx2,iky2][3,2])
            overlap = overlapIntegral(ikx1,iky1,ikx2,iky2,ikx3,iky3,model)
            #C_ZF, C_NZF = couplingCoefficients(kx1,ky1,1,kx2,ky2,2,kx3,ky3,3,model)
            #println(τ,' ',tau_nzf)
            tauMatrix_NZF[ikx2,iky2] = τ
            coeffMatrix_NZF[ikx2,iky2] = C_pqr*kxSpec(kx1,sgrid)*kxSpec(kx2,sgrid)*kxSpec(kx3,sgrid)*kySpec(ky1,sgrid)*kySpec(ky2,sgrid)*kySpec(ky3,sgrid)*overlap
          end
        end
      end
    end
  end
  maxIndex_ZF = findmax(abs.(tauVector_ZF.*coeffVector_ZF))[2]
  maxIndex_NZF = findmax(abs.(tauMatrix_NZF.*coeffMatrix_NZF))[2]
  return (tauVector_ZF[maxIndex_ZF], coeffVector_ZF[maxIndex_ZF],
          tauMatrix_NZF[maxIndex_NZF], coeffMatrix_NZF[maxIndex_NZF])
end

function couplingCoefficients(kx1::AbstractFloat,
                              ky1::AbstractFloat,
                              p::Integer,
                              kx2::AbstractFloat,
                              ky2::AbstractFloat,
                              q::Integer,
                              kx3::AbstractFloat,
                              ky3::AbstractFloat,
                              r::Integer,
                              model::FluidModel{T,6,1,FG};
                             ) where{T <: AbstractFloat,
                                     FG <: AbstractFieldlineGeometry}
  ikx1, iky1 = index(kx1,ky2,model.spectralGrid)
  ikx2, iky2 = index(kx2,ky2,model.spectralGrid)
  ikx3, iky3 = index(kx3,ky3,model.spectralGrid)
  C_pqF = ky1*(kx1-kx2) *(model.invCouplingMatrix[ikx1,iky1][p,2]*
                          model.couplingMatrix[ikx2,iky2][2,q] +
                          model.invCouplingMatrix[ikx1,iky1][p,3]*
                          model.couplingMatrix[ikx2,iky2][3,q])

  # For the three wave model, any entry in row 1 is 1
  C_pqr = (kx1*ky2 - ky1*kx2)*(model.invCouplingMatrix[ikx1,iky1][p,2]*
                               #model.couplingMatrix[ikx3,iky3][1,r]*
                               model.couplingMatrix[ikx2,iky2][2,q] +
                               0.5*model.invCouplingMatrix[ikx1,iky1][p,3]*
                               #model.couplingMatrix[ikx3,iky3][1,r]*
                               model.couplingMatrix[ikx2,iky2][3,q])
  return C_pqF, C_pqr
end
