#Creating and saving a fluid model

PTSM = PlasmaTurbulenceSaturationModel

function woutToModelSpectrum3fm(woutfile::String, output::String;
                     s0=0.5,
		     nF=3,
                     gradT=3.0,
		     gradn=0.0,
		     beta=0.0,
		     solnInterval=8π
		     )

  wout = NetCDF.open(woutfile)
  vmec, vmec_data = readVmecWout(wout)
  vsurf = VmecSurface(s0,vmec)
  model = PlasmaTurbulenceSaturationModel.setupLinearModel(vsurf; nFields=nF, ∇T=gradT, ∇n=gradn, β=beta, solnInterval=solnInterval);

  PTSM.computeLinearSpectrum!(model);

  PTSM.save(model,filename=output)

end

function woutToModelSingleMode3fm(woutfile::String,
		     kx::Float64,
		     ky::Float64,
                     s0::Float64,
		     nF::Int,
                     gradT::Float64,
		     gradn::Float64,
		     beta::Float64,
		     output::String)
  wout = NetCDF.open(woutfile)
  vmec, vmec_data = readVmecWout(wout)
  vsurf = VmecSurface(s0,vmec)
  model = PlasmaTurbulenceSaturationModel.setupLinearModel(vsurf; nFields=nF, ∇Tᵢ=gradT, ∇n=gradn, β=beta, solnInterval=solnInterval);

  PTSM.computeLinearMode!(kx=kx,ky=ky,model=model);

  PTSM.save(model,filename=output)

end

function woutToModelSpectrum5fm(woutfile::String, output::String;
                     s0=0.5,
		     nF=6,
                     gradTi=3.0,
                     gradTe=0.0,
		     gradn=0.0,
		     beta=0.0,
		     solnInterval=8π
		     )

  wout = NetCDF.open(woutfile)
  vmec, vmec_data = readVmecWout(wout)
  vsurf = VmecSurface(s0,vmec)
  model = PlasmaTurbulenceSaturationModel.setupLinearModel(vsurf; nFields=nF, ∇Tᵢ=gradTi, ∇Tₑ=gradTe, ∇n=gradn, β=beta, solnInterval=solnInterval);

  PTSM.computeLinearSpectrum!(model);

  PTSM.save(model,filename=output)

end

function woutToModelSingleMode5fm(woutfile::String,
		     kx::Float64,
		     ky::Float64,
                     s0::Float64,
		     nF::Int,
                     gradTi::Float64,
                     gradTe::Float64,
		     gradn::Float64,
		     beta::Float64,
		     npol::Int,
		     output::String)
  wout = NetCDF.open(woutfile)
  vmec, vmec_data = readVmecWout(wout)
  vsurf = VmecSurface(s0,vmec)
  model = PlasmaTurbulenceSaturationModel.setupLinearModel(vsurf; nFields=nF, ∇Tᵢ=gradTi, ∇Tₑ=gradTe, ∇n=gradn, β=beta, solnInterval=npol*2π);

  PTSM.computeLinearMode!(kx=kx,ky=ky,model=model);

  PTSM.save(model,filename=output)

end
