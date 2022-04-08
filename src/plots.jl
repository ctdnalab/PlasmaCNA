function getStdTheme(;mainFont::String="ArialMT",pointSize::Float64=0.8,
	lineWidth::Float64=0.35,majorFontSize::Int=10,minorFontSize::Int=8,
	background=colorant"#ffffff")

	if background == colorant"#ffffff"
		gridColor = colorant"#c9c9c9"
	else
		gridColor = colorant"#ffffff"
	end

	t = Theme(
		panel_fill=background,
		grid_color=gridColor,
		grid_line_style=:solid,
		background_color=background,
		major_label_color=colorant"#282728",
		minor_label_color=colorant"#282728",
		key_title_color=colorant"#282728",
		key_label_color=colorant"#282728",
		major_label_font_size=(majorFontSize)pt, # font size for axes label and title
		minor_label_font_size=(minorFontSize)pt, #font size for numbers on axes
		highlight_width=0pt,
		point_size=(pointSize)mm,
		line_width=(lineWidth)mm,
		default_color=colorant"#3b76d6",
		plot_padding=[0.5mm,0.5mm,0.5mm,0.5mm]
		)

	if isnothing(mainFont) == false
		t.major_label_font=mainFont
		t.minor_label_font=mainFont
		t.key_title_font=mainFont
		t.key_label_font=mainFont
	end
	return t
end


function _bin_based_labels(df::DataFrame)
	posRange = collect(range(1,stop=size(df,1)))
	chroms = unique(df[!,:chrom])
	chromStarts = [1]
	for (idx,c) in enumerate(chroms[1:end])
		dfSize = size(df |> @filter(_.chrom == c) |> DataFrame,1)
		push!(chromStarts,chromStarts[idx]+dfSize)
	end
	tickLabels = Dict(zip(chromStarts[1:end-1],[string(x) for x in chroms]))
	return (chromStarts[1:end-1],tickLabels,posRange)
end

"""
    plot_plasmacna_segs(df,outFile;drawsegs=true,sampleName="",
	cnvCutoff=0.3,outType="pdf",sortChroms=false,
	yMin=0,yMax=0,plotwidth=300,plotheight=65)

Generates a plot of bins and segements across the genome from a bin DataFrame. 
Segment lines are drawn at either the median depth of all bins in the segment, or 
at the value provided in the segAvg column if found in the dataframe.\n\n
If yMin and yMax are set to 0, the min/max y values used for the plot will be
-2/2 if all bin depths are between -2/2. If any bin depths exceed these limits, 
the min/max y values will be set 0.1 less/more than the min/max bin depth respectively.\n\n
Segment gain/loss coloring is determined by the cnvCutoff argument. A segment's copy number is
calculated as binDepth / cnvCutoff + 2, and segments are colored as gain/loss if their
copy number is > 2.04 or < 1.96.

# Arguments
*  df::DataFrame Bin DataFrame to plot
*  outFile::String Path to output file. Do not include file extension
*  drawsegs::Bool Draw colored horizontal lines across bins in the same segment (default true)
*  sampleName::String Sample name, used as x axis label (default "")
*  cnvCutoff::Float64 Value to determine gain/loss coloring (default 0.3)
*  outType::String pdf or svg. Appropriate file extension will be added to outFile (default pdf)
*  sortChroms::Bool Sort chrom column by name before generating plot (default false)
*  yMin::Real Minimum y axis value (default 0)
*  yMax::Real Maximum y axis value (default 0)
*  plotwidth::Real Plot width in mm (default 300)
*  plotheight::Real Plot height in mm (default 65)
"""
function plot_plasmacna_segs(df::DataFrame,outFile::String;
	drawsegs::Bool=true,sampleName::String="",cnvCutoff::Float64=0.3,
	outType::String="pdf",sortChroms::Bool=false,yMin::Real=0,yMax::Real=0,
	plotwidth::Real=300,plotheight::Real=65)

	COL = Dict{String,Any}()
	COL["points"] = colorant"#878484"
	COL["gain"] = colorant"#db2929"
	COL["loss"] = colorant"#30a836"
	COL["base"] = colorant"#3b76d6"
	COL["background"] = colorant"#ffffff"

	t = getStdTheme(pointSize=0.2,lineWidth=0.35,background=COL["background"])
	Gadfly.push_theme(t)

	if sampleName == ""
		sampleName = string([x for x in split(basename(outFile),".")][1])
		sampleName = replace(sampleName,".tsv"=>"")
	end


	chromStarts,tickLabels,posRange = _bin_based_labels(df)
	colNames = Dict{Symbol,Bool}()
	for s in names(df)
		colNames[Symbol(s)] = true
	end

	if sortChroms == true
		sort!(df,:chrom)
	end

	if yMin == 0
		if minimum(filter(!isnan,df[!,:depth])) > -2.0
			yMin = -2.0
		else
			yMin = minimum(filter(!isnan,df[!,:depth])) - .1
		end
	end

	if yMax == 0
		if maximum(filter(!isnan,df[!,:depth])) < 2.0
			yMax = 2.0
		else
			yMax = maximum(filter(!isnan,df[!,:depth])) + .1
		end
	end

	layerCount = 1
	lBase = layer(x=posRange,y=df[!,:depth],Geom.point,order=layerCount,color=[COL["points"]])
	layerCount += 1
	allL = [lBase]

	if drawsegs == true
		segs = unique(df[!,:segment])
		for s in segs
			segBins = [x for x in range(1,stop=size(df,1)) if df[x,:segment] == s]
			sStart = posRange[segBins[1]]
			sEnd = posRange[segBins[end]]

			if haskey(colNames,:segAvg) == true
				yVal = filter(!isnan,df[segBins,:segAvg])[1]
			else
				yVal = median(filter(!isnan,df[segBins,:depth]))
			end
			copyNum = yVal / cnvCutoff + 2
			if haskey(colNames,:copynum) == true
				if df[segBins[1],:copynum] != 0
					copyNum = df[segBins[1],:copynum]
				end
			end

			if copyNum >= 2.04
				l = layer(x=[sStart],xend=[sEnd],y=[yVal],yend=[yVal],
					Geom.segment,color=[COL["gain"]],order=layerCount)
			elseif copyNum < 1.96
				l = layer(x=[sStart],xend=[sEnd],y=[yVal],yend=[yVal],
					Geom.segment,color=[COL["loss"]],order=layerCount)
			else
				l = layer(x=[sStart],xend=[sEnd],y=[yVal],yend=[yVal],
					Geom.segment,color=[COL["base"]],order=layerCount)
			end
			layerCount += 1
			push!(allL,l)
		end
	end

	p = plot(allL...,Guide.xticks(ticks=chromStarts),
		Scale.x_continuous(labels = x -> tickLabels[x]),
		Guide.xlabel(sampleName),Guide.ylabel("Log2 Ratio"),
		Coord.cartesian(xmin=1,xmax=posRange[end],ymin=yMin,ymax=yMax));

	if outType == "pdf"
		p |> PDF("$outFile.pdf",(plotwidth)mm,(plotheight)mm)
	elseif outType == "svg"
		p |> SVG("$outFile.svg",(plotwidth)mm,(plotheight)mm)
	end
	Gadfly.pop_theme()
	return p
end

"""
    plot_plasmacna_folder(inFolder,outFolder;yMin=-2.0,yMax=2.0,plotwidth=300,plotheight=65)

Takes a folder of bin tsv files and generates a folder of seg plots for each file.
# Arguments
*  inFolder::String Path to folder of bin tsvs
*  outFolder::String Path to output folder
*  yMin::Real Minimum y axis value (default -2.0)
*  yMax::Real Maximum y axis value (default 2.0)
*  plotwidth::Real Plot width in mm (default 300)
*  plotheight::Real Plot height in mm (default 65)
"""
function plot_plasmacna_folder(inFolder::String,outFolder::String;
	yMin::Real=-2.0,yMax::Real=2.0,plotwidth::Real=300,plotheight::Real=65)
	files = readdir(inFolder)
	run(`mkdir -p $outFolder`)
	for f in files
		fPath = joinpath(inFolder,f)
		n = replace(f,".bins.tsv"=>"")
		n = replace(f,".tsv"=>"")
		outN = joinpath(outFolder,n)
		df = CSV.File(fPath,delim='\t') |> DataFrame
		if typeof(df[1,:chrom]) == Int
			sort!(df,:chrom)
		end
		p = plot_plasmacna_segs(df,joinpath(outFolder,outN),yMin=yMin,yMax=yMax,
			plotwidth=plotwidth,plotheight=plotheight)
	end
end

function plot_matched_samples(inFolders::Vector{String},methodNames::Vector{String},
	outFolder::String; yMin::Real=0,yMax::Real=0,cnvCutoff::Float64=0.3,
	outType::String="PDF")
	run(`mkdir -p $outFolder`)
	D = Dict{String,Vector{String}}()
	for fol in inFolders
		allF = readdir(fol)
		for f in allF
			s = string(split(f,".")[1])
			if haskey(D,s) == false
				D[s] = String[]
			end
			push!(D[s],f)
		end
	end

	for s in keys(D)
		if length(D[s]) < length(inFolders)
			continue
		end
		allP = Any[]
		for (idx,fol) in enumerate(inFolders)
			df = CSV.File(joinpath(fol,D[s][idx]),delim='\t') |> DataFrame
			sort!(df,:chrom)
			p1 = plot_plasmacna_segs(df,"",sampleName=methodNames[idx],outType="none",
				yMin=yMin,yMax=yMax,cnvCutoff=cnvCutoff)
			push!(allP,p1)
		end
		pF = title(vstack(allP...),s)
		if outType == "PDF"
			pF |> PDF("$(joinpath(outFolder,s)).pdf",300mm,(65*length(inFolders))mm);
		elseif outType == "SVG"
			pF |> SVG("$(joinpath(outFolder,s)).svg",300mm,(65*length(inFolders))mm);
		end
	end
end

#function coord_based_labels(df::DataFrame)
#	posRange = [df[1,:endPos] - df[1,:startPos]/2]
#	total = 0
#	for r in eachrow(df)
#		dist = r[:endPos] - r[:startPos]
#		push!(posRange,(dist/2)+total)
#		total += dist
#	end
#	posRange = posRange[1:end-1]
#
#	chroms = sort(unique(df[!,:chrom]))
#	chromSizes = []
#	for c in chroms
#		dfS = df |> @filter(_.chrom == c) |> DataFrame
#		push!(chromSizes,maximum(dfS[!,:endPos]))
#	end
#
#	chromStarts = Int[1]
#	total = 1
#	for (idx,c) in enumerate(chroms[2:end])
#		total += chromSizes[idx]
#		push!(chromStarts,total)
#	end
#
#	tickLabels = Dict(zip(chromStarts,[string(x) for x in chroms]))
#	return (chromStarts,tickLabels,posRange)
#end

