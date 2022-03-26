## My attempt at implementing the 5FM into Julia

## Here, using 3FM LinearOperators.jl script as a starting point

## Changed: "model::FluidModel{T,3,1,GeneGeometry} -> model::FluidModel{T,6,1,GeneGeometry}"

#="""
    denseLinearOperators(kx::AbstractFloat,ky::AbstractFloat,
                         model::FluidModel{T,3,1,GeneGeometry}) where T

Create the dense matrix operators `(A,B)` for the generalized eigenvalue problem 
A⋅x = ω B⋅x for the three-field fluid model at a given `(kx,ky)`
from the fieldline geometry specified by `geom` and fluid model parameters given by `model`.
"""
=#

function denseLinearOperators(kx::AbstractFloat,
                              ky::AbstractFloat,
                              model::FluidModel{T,6,1,GeneGeometry},
                             ) where T
  # Linear operators using GENE normalizations
  # Find the range of indices for a given (kx,ky) using the intervalRange
  # defined in the model
  # indexRange = intervalRangeIndices(kx,ky,model)
  kxIndex, kyIndex = index(kx,ky,model.spectralGrid)
  indexRange = model.indices[kxIndex,kyIndex]
  geomView = @view(model.geometry[indexRange])

#=
  nPoints = length(indexRange)
  rank = 6*nPoints
  lhs = zeros(Complex{T},(rank,rank))
  rhs = Matrix{Complex{T}}(I,(rank,rank))
=#

  # For convenience extract the relevant plasma parameters
  ## Added electron temperature gradient
  ∇n = model.plasma.∇n
  ∇Tᵢ= model.plasma.∇Tᵢ
  ∇Tₑ= model.plasma.∇Tₑ
  τ = model.plasma.τ 
  β = model.plasma.β
  npol = model.problem.solnInterval/2π ### Check that this is correct
  Bhat = geomView.Bhat
  jac = geomView.jac

  # Set the FLR and curvature vectors spanning the indices given by indexRange
  Bk, Dk = buildBkDk(kx,ky,geomView) 

  ### May not need this given details from Ian's script below - using mu1, mu2, etc.
  # Set the numerical hyperdiffusion parameters
  dθ = step(model.geomParameters.points)
  hyperDiffConstant = model.problem.ε
  hypz = hyperDiffConstant*(im*0.5*dθ)^model.problem.hyperdiffusionOrder

  ## Added from Ians' script
  probSize = length(geomView.Bhat) 
  omegap = ∇n+∇Tᵢ
  z_domain = -(npol):(2*npol)/(probSize-1):(npol); # setting the z domain correctly
  z_domain = transpose(z_domain); ### Not sure what this line is doing in Ian's script - assuming he's taking the transpose   
  jacBinv = zeros(1,probSize);
  jacBinv = jacBinv + 1/(jac.*Bhat); # inverse of the quantity J*B along the field line, important for calculating parallel derivatives correctly
  dz = z_domain[2] - z_domain[1];

  ## Added from Ian's script
  ### Not sure how this hyperdiffusion method works in comparison to 3FM approach?
  mu1 = 2.0; # hyperdiffusion coefficient Eq 1
  mu2 = 2.0; # hyperdiffusion coefficient Eq 2
  mu3 = 2.0; # hyperdiffusion coefficient Eq 3
  nu = 0; # k_perp hyperdiffusion coefficient
  dBkdz=zeros(probSize,1); # preallocating for the parallel derivative
  
	### Not sure about this dBkdz term, or the following loops

	for i=1
		dBkdz[i]= jacBinv[i]/π*((0/2)*Bk[i]+(1/2)*Bk[i+1])/dz; # calculates dBkdz from Bk(z) profile
	end

	for i=probSize
		dBkdz[i]= jacBinv[i]/π*((0/2)*Bk[i]-(1/2)*Bk[i-1])/dz; # calculates dBkdz from Bk(z) profile
	end

	for i=2:probSize-1
		dBkdz[i]= jacBinv[i]/π*(-(1/2)*Bk[i-1]+(1/2)*Bk[i+1])/dz;
	end

	## Following is 3FM procedure for obtaining matrices A and B (from eigenvalue problem A⋅x = ω B⋅x)
	## I'll start by using Ian's more verbose version of calculating A and B - can then maybe make a more concise version once we confirm this works with the rest of 5FM-Julia-implementation

#=

  for index = 1:nPoints
    Bk_θ = Bk[index]
    Dk_θ = Dk[index]
    Bhat = geomView.Bhat[index]
    jac = geomView.jac[index]

    # Compute the finite difference splines for the first order derivative operator
    # e.g. ∂/∂z → [c₁,c₂,…,cₙ₋₁,cₙ] where n = model.grid.derivativeErrorOrder + 1
    fdm = parallelDerivative(index,nPoints,1,model.problem.derivativeErrorOrder)  
    # Compute the finite difference spline for the n-th order derivative operator
    # where n is specified by model.grid.hyperdiffusionOrder
    hyperFdm = parallelDerivative(index,nPoints,
                                  model.problem.hyperdiffusionOrder,
                                  model.problem.derivativeErrorOrder)

    # Deal with boundaries where the point beyond the range is assumed to be zero
    derivativeIndices = findall(x->((x+index)>0 && (x+index)<=nPoints),fdm.grid)
    hyperDiffIndices = findall(x->((x+index)>0 && (x+index)<=nPoints),hyperFdm.grid)

    fdmGrid = fdm.grid[derivativeIndices]
    # Apply the prefactor -i/√gB to each derivative term
    fdmCoefs = convert.(Complex{T},fdm.coefs[derivativeIndices])*(-im/(dθ*jac*Bhat))

    hyperDiffGrid = hyperFdm.grid[hyperDiffIndices]
    # Approximate the n-th order hyperdiffusion operator by taking only the
    # leading order term
    hyperDiffCoefs = hyperFdm.coefs[hyperDiffIndices]*
                     (-im*hypz/(dθ*jac*Bhat)^model.problem.hyperdiffusionOrder)

    # Use array views to efficiently access the underlying data to apply the
    # finite differnce splines. For the three-field model with operator rank
    # 3N × 3N, the system can be broken into N × N blocks:
    #     || Φₖ ||    || A1 | A2 | A3 ||   || Φₖ ||
    # A ⋅ || Uₖ || ⇒  || A4 | A5 | A6 || ⋅ || Uₖ ||
    #     || Tₖ ||    || A7 | A8 | A9 ||   || Tₖ ||

    # Left-hand side operator
    # Block 1
    lhsView = view(lhs,index,index .+ hyperDiffGrid)
    lhsView .+= hyperDiffCoefs
    lhs[index,index] += ky*∇n/Bhat+ Dk_θ*(1+5*τ/3)

    # Block 2
    lhsView = view(lhs,index,(nPoints+index) .+ fdmGrid)
    lhsView .+= fdmCoefs

    # Block 3
    lhs[index,2*nPoints+index] += Dk_θ*τ

    # Block 4
    lhsView = view(lhs,nPoints+index,index .+ fdmGrid)
    lhsView .+= fdmCoefs*(1+5*τ/3)

    # Block 5
    lhsView = view(lhs,nPoints+index,(nPoints+index) .+ hyperDiffGrid)
    lhsView .+= hyperDiffCoefs

    # Block 6
    lhsView = view(lhs,nPoints+index,(2*nPoints+index) .+ fdmGrid)
    lhsView .+= fdmCoefs*τ

    # Block 7
    lhs[2*nPoints+index,index] += ky*(∇T - 2*∇n/3)/Bhat #+ 10*τ*Dk_θ/9  
    #lhs[2*nPoints+index,2*nPoints+index] += 5*τ*Dk_θ/3
    
    # Block 8
    #
    # Block 9
    lhsView = view(lhs,2*nPoints+index,(2*nPoints+index) .+ hyperDiffGrid)
    lhsView .+= hyperDiffCoefs

    # Right hand side operator
    # Block 1
    rhs[index,index] += Bk_θ*(1+5*τ/3)

    #Block 3
    rhs[index,2*nPoints+index] += Bk_θ*τ

  end 
  return lhs, rhs
end

=#

## From Ian's script - to the end of his document


# Constructing the matrices A and B

nFields = 6  ##I assume this can be taken directly from the FluidModel when using "setupLinearModel" in Interfaces.jl

	B = zeros(nFields*probSize,nFields*probSize); ### B = rhs in 3FM
	A = zeros(nFields*probSize,nFields*probSize); ### A = lhs in 3FM

	for i=1:probSize
	    B[i,i]=1.0+(5/(3*τ))*Bk[i];
	    B[i,i+probSize]= Bk[i];
	    B[i,i+3*probSize]= Bk[i];
	    B[i+probSize,i]= (5/(3*τ))*Bk[i];
	    B[i+probSize,i+probSize]= Bk[i];
	    B[i+probSize,i+3*probSize]= Bk[i];
	    B[i+2*probSize,i+2*probSize]= 1;
	    B[i+2*probSize,i+4*probSize]= 1; # the A_parallel term in the ion parallel momentum equation, generally doesn't make a big difference in the calculation
	    B[i+3*probSize,i+3*probSize]= 1;
	    B[i+4*probSize,i+4*probSize]= 1;
	end

	for i=1:probSize
		A[i,i]=Dk[i]*5/(3*τ);
		A[i,i+probSize] = ky*∇n+Dk[i];
	    A[i,i+2*probSize] = 0; # parallel derivative
	    A[i,i+3*probSize] = Dk[i];
	    A[probSize+i,i] = Dk[i]*(1+5/(3*τ));
	    A[probSize+i,i+3*probSize] = Dk[i];
	    A[probSize+i,i+4*probSize] = (-2*1i/β)*dBkdz[i]; # parallel derivative
	    A[probSize+i,i+5*probSize] = Dk[i]; # electron temperature term
	    A[2*probSize+i,i] =  5/(3*τ)*0; # parallel derivative
	    A[2*probSize+i,i+probSize] = 0; # parallel derivative
	    A[2*probSize+i,i+3*probSize] = 0; # parallel derivative
	    A[2*probSize+i,i+4*probSize] = -ky/τ*(∇n+∇Tᵢ);
	    A[3*probSize+i,i] = 0; #10/(9*τ^2)*Dk[i); # ADDITIONAL TERM
	    A[3*probSize+i,i+probSize] = ky/τ*(∇Tᵢ-(2/3)*∇n);
	    A[3*probSize+i,i+3*probSize] = 0; #5/(3*τ)*Dk[i); # ADDITIONAL TERM
	    A[4*probSize+i,i] = (-1)*0; # parallel derivative
	    A[4*probSize+i,i+probSize] = 0; # minus parallel derivative of above line
	    A[4*probSize+i,i+4*probSize] = ky*∇n;
	    A[5*probSize+i,i+4*probSize] = ∇Tₑ*ky;
	    A[5*probSize+i,i+5*probSize] = 0; # parallel derivative shenanigans
	end

	for j=1
	    A[j,j+2*probSize] = (-1im*-jacBinv[j]/π*(0/2))/dz; # parallel derivative
	    A[j,j+2*probSize+1] = (-1im*-jacBinv[j]/π)*(1/2)/dz;
	    
	    A[probSize+j,j+4*probSize] = A[probSize+j,j+4*probSize] + 2*Bk[j]/β*(-1im*-jacBinv[j]/π*(0/2))/dz; # parallel derivative
	    A[probSize+j,j+4*probSize+1] = A[probSize+j,j+4*probSize+1] + 2*Bk[j]/β*(-1im*-jacBinv[j]/π)*(1/2)/dz;
	    
	    A[2*probSize+j,j] = (-1im*-jacBinv[j]/π*(5/(3*τ))*(0/2))/dz; # parallel derivative
	    A[2*probSize+j,j+1] = (-1im*-jacBinv[j]/π)*(5/(3*τ))*(1/2)/dz;
	    
	    A[2*probSize+j,j+probSize] = (-1im*-jacBinv[j]/π*(0/2))/dz; # parallel derivative
	    A[2*probSize+j,j+probSize+1] = (-1im*-jacBinv[j]/π)*(1/2)/dz;
	    
	    A[2*probSize+j,j+3*probSize] = (-1im*-jacBinv[j]/π*(0/2))/dz; # parallel derivative
	    A[2*probSize+j,j+3*probSize+1] = (-1im*-jacBinv[j]/π)*(1/2)/dz;
	    
	    A[4*probSize+j,j] = (-1)*(-1im*-jacBinv[j]/π*(0/2))/dz; # parallel derivative
	    A[4*probSize+j,j+1] = (-1)*(-1im*-jacBinv[j]/π)*(1/2)/dz;
	    
	    A[4*probSize+j,j+probSize] = (-1im*-jacBinv[j]/π*(0/2))/dz; # parallel derivative
	    A[4*probSize+j,j+probSize+1] = (-1im*-jacBinv[j]/π)*(1/2)/dz;
	    
	    A[5*probSize+j,j+5*probSize] = (-1)*(-1im*-jacBinv[j]/π*(0/2))/dz; # parallel derivative
	    A[5*probSize+j,j+5*probSize+1] = (-1)*(-1im*-jacBinv[j]/π)*(1/2)/dz;
	end

	for j=probSize
	    A[j,j+2*probSize] = (-1im*-jacBinv[j]/π*(0/2))/dz; # parallel derivative
	    A[j,j+2*probSize-1] = (-1im*-jacBinv[j]/π)*(-1/2)/dz;
	    
	    A[probSize+j,j+4*probSize] = A[probSize+j,j+4*probSize] + 2*Bk[j]/β*(-1im*-jacBinv[j]/π*(0/2))/dz; # parallel derivative
	    A[probSize+j,j+4*probSize-1] = A[probSize+j,j+4*probSize-1] + 2*Bk[j]/β*(-1im*-jacBinv[j]/π)*(-1/2)/dz;
	    
	    A[2*probSize+j,j] = (-1im*-jacBinv[j]/π*(5/(3*τ))*(0/2))/dz; # parallel derivative
	    A[2*probSize+j,j-1] = (-1im*-jacBinv[j]/π)*(5/(3*τ))*(-1/2)/dz;
	    
	    A[2*probSize+j,j+probSize] = (-1im*-jacBinv[j]/π*(0/2))/dz; # parallel derivative
	    A[2*probSize+j,j+probSize-1] = (-1im*-jacBinv[j]/π)*(-1/2)/dz;
	    
	    A[2*probSize+j,j+3*probSize] = (-1im*-jacBinv[j]/π*(0/2))/dz; # parallel derivative
	    A[2*probSize+j,j+3*probSize-1] = (-1im*-jacBinv[j]/π)*(-1/2)/dz;
	    
	    A[4*probSize+j,j] = (-1)*(-1im*-jacBinv[j]/π*(0/2))/dz; # parallel derivative
	    A[4*probSize+j,j-1] = (-1)*(-1im*-jacBinv[j]/π)*(-1/2)/dz;
	    
	    A[4*probSize+j,j+probSize] = (-1im*-jacBinv[j]/π*(0/2))/dz; # parallel derivative
	    A[4*probSize+j,j+probSize-1] = (-1im*-jacBinv[j]/π)*(-1/2)/dz;
	    
	    A[5*probSize+j,j+5*probSize] = (-1)*(-1im*-jacBinv[j]/π*(0/2))/dz; # parallel derivative
	    A[5*probSize+j,j+5*probSize-1] = (-1)*(-1im*-jacBinv[j]/π)*(-1/2)/dz;
	end

	for j=2:probSize-1
	    A[j,j+2*probSize-1] = (-1im*-jacBinv[j]/π*(-1/2))/dz; # parallel derivative
	    A[j,j+2*probSize+1] = (-1im*-jacBinv[j]/π)*(1/2)/dz;
	    
	    A[probSize+j,j+4*probSize-1] = A[probSize+j,j+4*probSize-1] + 2*Bk[j]/β*(-1im*-jacBinv[j]/π*(-1/2))/dz; # parallel derivative
	    A[probSize+j,j+4*probSize+1] = A[probSize+j,j+4*probSize+1] + 2*Bk[j]/β*(-1im*-jacBinv[j]/π)*(1/2)/dz;
	    
	    A[2*probSize+j,j-1] = (-1im*-jacBinv[j]/π*(5/(3*τ))*(-1/2))/dz; # parallel derivative
	    A[2*probSize+j,j+1] = (-1im*-jacBinv[j]/π)*(5/(3*τ))*(1/2)/dz;
	    
	    A[2*probSize+j,j+probSize-1] = (-1im*-jacBinv[j]/π*(-1/2))/dz; # parallel derivative
	    A[2*probSize+j,j+probSize+1] = (-1im*-jacBinv[j]/π)*(1/2)/dz;
	    
	    A[2*probSize+j,j+3*probSize-1] = (-1im*-jacBinv[j]/π*(-1/2))/dz; # parallel derivative
	    A[2*probSize+j,j+3*probSize+1] = (-1im*-jacBinv[j]/π)*(1/2)/dz;
	    
	    A[4*probSize+j,j-1] = (-1)*(-1im*-jacBinv[j]/π*(-1/2))/dz; # parallel derivative
	    A[4*probSize+j,j+1] = (-1)*(-1im*-jacBinv[j]/π)*(1/2)/dz;
	    
	    A[4*probSize+j,j+probSize-1] = (-1im*-jacBinv[j]/π*(-1/2))/dz; # parallel derivative
	    A[4*probSize+j,j+probSize+1] = (-1im*-jacBinv[j]/π)*(1/2)/dz;
	    
	    A[5*probSize+j,j+5*probSize-1] = (-1)*(-1im*-jacBinv[j]/π*(-1/2))/dz; # parallel derivative
	    A[5*probSize+j,j+5*probSize+1] = (-1)*(-1im*-jacBinv[j]/π)*(1/2)/dz;
	end

	for k=3:probSize-2
	    A[k,k-2]=A[k,k-2]-1im*(-jacBinv[k]/π)^4*mu1*(1/1)/dz^4;
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k]/π)^4*mu1*(-4/1)/dz^4;
	    A[k,k]=A[k,k]-1im*(-jacBinv[k]/π)^4*mu1*(6/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k]/π)^4*mu1*(-4/1)/dz^4;
	    A[k,k+2]=A[k,k+2]-1im*(-jacBinv[k]/π)^4*mu1*(1/1)/dz^4;
	end

	for k=3+probSize:2*probSize-2
	    A[k,k-2]=A[k,k-2]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(1/1)/dz^4;
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(-4/1)/dz^4;
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(6/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(-4/1)/dz^4;
	    A[k,k+2]=A[k,k+2]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(1/1)/dz^4;
	end

	for k=3+2*probSize:3*probSize-2
	    A[k,k-2]=A[k,k-2]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(1/1)/dz^4;
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(6/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k+2]=A[k,k+2]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(1/1)/dz^4;
	end

	for k=3+3*probSize:4*probSize-2
	    A[k,k-2]=A[k,k-2]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(1/1)/dz^4;
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(6/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k+2]=A[k,k+2]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(1/1)/dz^4;
	end

	for k=3+4*probSize:5*probSize-2
	    A[k,k-2]=A[k,k-2]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(1/1)/dz^4;
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(6/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k+2]=A[k,k+2]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(1/1)/dz^4;
	end

	for k=3+5*probSize:6*probSize-2
	    A[k,k-2]=A[k,k-2]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(1/1)/dz^4;
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(6/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k+2]=A[k,k+2]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(1/1)/dz^4;
	end

	for k=1
	    A[k,k]=A[k,k]-1im*(-jacBinv[k]/π)^4*mu1*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k]/π)^4*mu1*(6/1)/dz^4;
	    A[k,k+2]=A[k,k+2]-1im*(-jacBinv[k]/π)^4*mu1*(-4/1)/dz^4;
	    A[k,k+3]=A[k,k+3]-1im*(-jacBinv[k]/π)^4*mu1*(1/1)/dz^4;
	end

	for k=probSize+1
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(6/1)/dz^4;
	    A[k,k+2]=A[k,k+2]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(-4/1)/dz^4;
	    A[k,k+3]=A[k,k+3]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(1/1)/dz^4;
	end

	for k=2*probSize+1
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(6/1)/dz^4;
	    A[k,k+2]=A[k,k+2]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k+3]=A[k,k+3]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(1/1)/dz^4;
	end

	for k=3*probSize+1
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(6/1)/dz^4;
	    A[k,k+2]=A[k,k+2]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k+3]=A[k,k+3]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(1/1)/dz^4;
	end

	for k=4*probSize+1
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(6/1)/dz^4;
	    A[k,k+2]=A[k,k+2]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k+3]=A[k,k+3]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(1/1)/dz^4;
	end

	for k=5*probSize+1
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(6/1)/dz^4;
	    A[k,k+2]=A[k,k+2]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k+3]=A[k,k+3]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(1/1)/dz^4;
	end

	for k=2
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k]/π)^4*mu2*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k]=A[k,k]-1im*(-jacBinv[k]/π)^4*mu2*(6/1)/dz^4;
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k]/π)^4*mu2*(-4/1)/dz^4;
	    A[k,k+2]=A[k,k+2]-1im*(-jacBinv[k]/π)^4*mu2*(1/1)/dz^4;
	end

	for k=probSize+2
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(6/1)/dz^4;
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(-4/1)/dz^4;
	    A[k,k+2]=A[k,k+2]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(1/1)/dz^4;
	end

	for k=2*probSize+2
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(6/1)/dz^4;
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k+2]=A[k,k+2]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(1/1)/dz^4;
	end

	for k=3*probSize+2
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(6/1)/dz^4;
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k+2]=A[k,k+2]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(1/1)/dz^4;
	end

	for k=4*probSize+2
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(6/1)/dz^4;
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k+2]=A[k,k+2]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(1/1)/dz^4;
	end

	for k=5*probSize+2
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(6/1)/dz^4;
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k+2]=A[k,k+2]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(1/1)/dz^4;
	end

	for k=probSize
	    A[k,k]=A[k,k]-1im*(-jacBinv[k]/π)^4*mu1*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k]/π)^4*mu1*(6/1)/dz^4;
	    A[k,k-2]=A[k,k-2]-1im*(-jacBinv[k]/π)^4*mu1*(-4/1)/dz^4;
	    A[k,k-3]=A[k,k-3]-1im*(-jacBinv[k]/π)^4*mu1*(1/1)/dz^4;
	end

	for k=2*probSize
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(6/1)/dz^4;
	    A[k,k-2]=A[k,k-2]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(-4/1)/dz^4;
	    A[k,k-3]=A[k,k-3]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(1/1)/dz^4;
	end

	for k=3*probSize
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(6/1)/dz^4;
	    A[k,k-2]=A[k,k-2]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k-3]=A[k,k-3]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(1/1)/dz^4;
	end

	for k=4*probSize
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(6/1)/dz^4;
	    A[k,k-2]=A[k,k-2]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k-3]=A[k,k-3]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(1/1)/dz^4;
	end

	for k=5*probSize
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(6/1)/dz^4;
	    A[k,k-2]=A[k,k-2]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k-3]=A[k,k-3]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(1/1)/dz^4;
	end

	for k=6*probSize
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(6/1)/dz^4;
	    A[k,k-2]=A[k,k-2]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k-3]=A[k,k-3]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(1/1)/dz^4;
	end

	for k=probSize-1
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k]/π)^4*mu1*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k]=A[k,k]-1im*(-jacBinv[k]/π)^4*mu1*(6/1)/dz^4;
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k]/π)^4*mu1*(-4/1)/dz^4;
	    A[k,k-2]=A[k,k-2]-1im*(-jacBinv[k]/π)^4*mu1*(1/1)/dz^4;
	end

	for k=2*probSize-1
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(6/1)/dz^4;
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(-4/1)/dz^4;
	    A[k,k-2]=A[k,k-2]-1im*(-jacBinv[k-probSize]/π)^4*mu2*(1/1)/dz^4;
	end

	for k=3*probSize-1
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(6/1)/dz^4;
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k-2]=A[k,k-2]-1im*(-jacBinv[k-2*probSize]/π)^4*mu3*(1/1)/dz^4;
	end

	for k=4*probSize-1
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(6/1)/dz^4;
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k-2]=A[k,k-2]-1im*(-jacBinv[k-3*probSize]/π)^4*mu3*(1/1)/dz^4;
	end

	for k=5*probSize-1
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(6/1)/dz^4;
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k-2]=A[k,k-2]-1im*(-jacBinv[k-4*probSize]/π)^4*mu3*(1/1)/dz^4;
	end

	for k=6*probSize-1
	    A[k,k+1]=A[k,k+1]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    # A[k,k]=A[k,k]-Bk[k]*Bk[k]*nu;
	    A[k,k]=A[k,k]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(6/1)/dz^4;
	    A[k,k-1]=A[k,k-1]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(-4/1)/dz^4;
	    A[k,k-2]=A[k,k-2]-1im*(-jacBinv[k-5*probSize]/π)^4*mu3*(1/1)/dz^4;
	end
	
	lhs = A
	rhs = B

	return lhs, rhs
end
