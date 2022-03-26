#Creating gist file from VMEC out (wout file)

function woutToGist(woutfile::String, output::String;
                     s0=0.5,
                     α=0.0,
		     nz0=128,
		     npol=1
		     )

	VMEC.geneGeometryFile(woutfile,s0,α,nz0,npol,output)
end

