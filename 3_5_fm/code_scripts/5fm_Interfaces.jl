function setupLinearModel(surface::VmecSurface;
                          nFields=6,
                          nVectors=1,
                          prec=Float64,
                          nkx=11,
                          nky=10,
                          Δkx = 0.05,
                          Δky = 0.05,
                          kxStart = 0.0,
                          kyStart = 0.0,
                          kxSpecLimit = 0.2,
                          kySpecLimit = 0.6,
                          kxSpecPower = 2.0,
                          kySpecPower = 4/3,
                          ∇Tᵢ = 3.0,
                          ∇Tₑ = 0.0,
                          ∇n = 0.0,
                          τ = 1.0,
                          β = 1.0,
                          solver=JDQZ(),
                          solnInterval=8π,
                          ε=200.0,
                          errorOrder=4,
                          hyperdiffusionOrder=2,
                          quadPoints=1000,
                          points_per_2pi=128,
                         ) 
  iota = surface.iota[1];
  s_hat = VMEC.shat(surface);
  nfp = surface.nfp;

  plasma = PlasmaParameters(∇Tᵢ,∇Tₑ,∇n,τ,β)
  specGrid = SpectralGrid(nkx=nkx,nky=nky,Δkx=Δkx,Δky=Δky,
                          kxStart=kxStart,kyStart=kyStart,
                          kxSpecLimit=kxSpecLimit,kySpecLimit=kySpecLimit,
                          kxSpecPower=kxSpecPower,kySpecPower=kySpecPower)
  problem = LinearProblem(solver=solver,solnInterval=solnInterval,
                          ε=ε,errorOrder=errorOrder,
                          hyperdiffusionOrder=hyperdiffusionOrder,
                          quadPoints=quadPoints)

  maxAngle = last(specGrid.kxRange)/(abs(s_hat)*first(specGrid.kyRange)) + solnInterval

  # Defined in θ
  angleRange = -maxAngle:2π/points_per_2pi:maxAngle

  p = MagneticCoordinateCurve(PestCoordinates,surface.s*surface.phi[end]/(2π)*surface.signgs,0.0,-angleRange/iota[1])

  metric, modB, sqrtg, K1, K2, ∂B∂θ = VMEC.geneGeometryCoefficients(VMEC.GeneFromPest(),p,surface)
  θ = map(i->-i.ζ*iota,p)
  geom = FieldlineGeometry(metric,modB,sqrtg,K1,K2,∂B∂θ,θ)
  geomPars = GeneParameters(iota,s_hat,nfp,vector2range(θ))

  return FluidModel(nFields,nVectors,plasma,specGrid,problem,geom,geomPars,prec=prec)
  
end
