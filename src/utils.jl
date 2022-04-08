
function hasintersect(r1::Vector{A},r2::Vector{B}) where {A<:Real,B<:Real}
	@assert length(r1) == 2 && length(r2) == 2
	@assert r1[1] <= r1[2] && r2[1] <= r2[2]
	if r1[1] <= r2[2] && r2[1] <= r1[2]
		return true
	end
	return false
end

function getoverlap(r1::Vector{A},r2::Vector{B}) where {A<:Real,B<:Real}
	return length(range(maximum([r1[1],r2[1]]),stop=minimum([r1[2],r2[2]])))
end

"Generates windows of a fixed size along a chromosome"
function get_fixed_windows(chr::String,rSize::Int,chromEnd::Int;stepSize::Int=-1)
	allWindows = Vector{Window}()
	if stepSize == -1
		stepSize = rSize
	end
	bins = range(1,step=stepSize,stop=chromEnd)

	for (idx,b) in enumerate(bins)
		binEnd = b+rSize-1
		if binEnd > chromEnd
			binEnd = chromEnd
		end
		push!(allWindows,Window(chr,idx,b,binEnd,chromEnd,Dict{String,Float64}()))
	end
	return allWindows
end

function tbl2wig(df::DataFrame,returnCol::Symbol,outFile::String;chromName::Symbol=:chrom,startPos::Symbol=:startPos,endPos::Symbol=:endPos,returnType::Type=Int)
	allChrom  = unique(df[chromName])
	out = open(outFile,"w")
	for c in allChrom
		subD = filter(row -> row[chromName] == c,df)
		step = subD[1,endPos]
		if size(subD,1) == 1
			span = subD[end,endPos]
		else
			span = step
		end
		write(out,"fixedStep chrom=$c start=1 step=$step span=$span\n")
		for r in eachrow(subD)
			if returnType == Int
				x = Int(floor(r[returnCol]))
				write(out,"$(x)\n")
			else
				write(out,"$(r[returnCol])\n")
			end
		end
	end
	close(out)
end