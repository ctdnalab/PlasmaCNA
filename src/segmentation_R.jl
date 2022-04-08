
function read_segs_dnacopy(inFile::String)
	df = CSV.File(inFile,delim='\t') |> DataFrame
	allS = Vector{Tuple{Any,Int,Int}}()
	for r in eachrow(df)
		s = (r[Symbol("output.chrom")],r[Symbol("output.loc.start")],r[Symbol("output.loc.end")])
		push!(allS,s)
	end
	return allS
end

function annotate_segments!(df::DataFrame,segData::Vector{Tuple{Any,Int,Int}})
	df[!,:segment] .= 0
	for r in range(1,stop=size(df,1))
		for (idx,s) in enumerate(segData)
			if s[1] != df[r,:chrom]
				continue
			else
				if hasintersect([s[2],s[3]],[df[r,:startPos],df[r,:endPos]]) == true
					df[r,:segment] = idx
				end
			end
		end
	end
end

function update_bin_segs(binFolder::String,segFolder::String,outFolder::String;segType::String="dnacopy")
	run(`mkdir -p $outFolder`)
	binFiles = readdir(binFolder)
	
	for f in binFiles
		fPath = joinpath(binFolder,f)
		segFile = joinpath(segFolder,replace(f,".bins."=>".seg."))
		df = CSV.File(fPath,delim='\t') |> DataFrame
		
		if segType == "dnacopy"
			segs = read_segs_dnacopy(segFile)
		end
		annotate_segments!(df,segs)
		CSV.write(joinpath(outFolder,f),df,delim='\t')
	end
end