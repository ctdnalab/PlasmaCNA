using Distributed
using QuickArgParse
@everywhere using PlasmaCNA
@everywhere using DataFrames,CSV,Query

@everywhere function runCmd(s::Cmd)
	run(s)
end

@everywhere mutable struct PSInputs
	bamFile::String
	normalBins::DataFrame
	bins::DataFrame
	toolPath::String
	outDir::String
end

@everywhere function plasmaSeqProcessing(p::PSInputs)
	bamF = basename(p.bamFile)
	n = replace(bamF,".bam"=>"")
	toolPath = p.toolPath
	binBed = joinpath(p.outDir,"binCoords.bed")
	genomeFile = joinpath(p.outDir,"chromOrder.tsv")

	get_bedtools_depth!(p.bamFile,p.bins,binBed,genomeFile,bedtools=toolPath)

	p.bins[!,:normDepth] .= p.normalBins[!,:depth]

	#Only include bins that have reads in normal and sample
	bins = p.bins |> @filter(_.normDepth > 0) |> DataFrame
	bins = bins |> @filter(_.depth > 0) |> DataFrame

	#Normalize sample and control by mean read depth and gc content
	normalize_depth!(bins,"mean")
	normalize_depth!(bins,"mean",depthCol=:normDepth)
	normalize_depth!(bins,"log")
	normalize_depth!(bins,"log",depthCol=:normDepth)
	normalize_depth!(bins,"gc")
	normalize_depth!(bins,"gc",depthCol=:normDepth)

	#Normalize sample by normal depth
	normalize_depth!(bins,"col",normCol=:normDepth)
	bins = bins |> @filter(_.map < Inf) |> DataFrame
	CSV.write("$(p.outDir)/chgTables/$n.bins.tsv",bins,delim='\t')
	return
end

##############################################################################################################
function main()
	#Setup
	req = ["BinTable","BamDir","NormalBam","OutDir","DNACopyPath"]
	reqhelp = ["Table of genomic bins, generated by flexbin_generation.jl",
		"Directory of WGS bam files to analyze","Single bam file of healthy WGS data",
		"Directory for output files",
		"Path to DNAcopy R script, found in PlasmaCNA/bin"]
	opts= ["toolpath","RscriptPath"]
	optsHelp = ["Path to bedtools executable calc",
	"Path to RScript executable for running DNACopy"]
	optsDefault = ["bedtools","Rscript"]
	title = "PlasmaSeq"
	desc = "A Julia implementation of PlasmaSeq (DOI: 10.1186/gm439). 
Can be multithreaded using julia -p [cores]. Bams must be coordinate sorted."

	R = process_reqs(req,reqHelp=reqhelp,optTag=opts,optHelp=optsHelp,
		optDefault=optsDefault,title=title,desc=desc)
	build_usage!(R)
	params = parse_args(R)

	inFolder = abspath(params["BamDir"])
	normal = abspath(params["NormalBam"])
	outFolder = abspath(params["OutDir"])
	dnacopyPath = abspath(params["DNACopyPath"])
	rscriptPath = params["RscriptPath"]
	bedtoolsPath = params["toolpath"]
	run(`mkdir -p $(outFolder)`)
	run(`mkdir -p $(outFolder)/chgTables`)
	run(`mkdir -p $(outFolder)/segTables`)
	run(`mkdir -p $(outFolder)/finalTables`)
	run(`mkdir -p $(outFolder)/plots`)

	#Reading bins and findings bams
	bins = CSV.File(params["BinTable"],delim='\t') |> DataFrame
	bams = readdir(inFolder)
	bams = [x for x in bams if occursin("bai",x) == false]
	write_bedtools_tempfiles(bins,outFolder)
	binBed = joinpath(outFolder,"binCoords.bed")
	genomeFile = joinpath(outFolder,"chromOrder.tsv")

	#Processing normal reference sample data
	normalBins = copy(bins)
	get_bedtools_depth!(normal,normalBins,binBed,genomeFile,bedtools=bedtoolsPath)

	#Processing tumor sample bam data and writing normalized depth tables
	allParams = PSInputs[]
	for b in bams
		sampleBins = copy(bins)
		p = joinpath(inFolder,b)
		bParams = PSInputs(p,normalBins,sampleBins,bedtoolsPath,outFolder)
		push!(allParams,bParams)
	end
	pmap(plasmaSeqProcessing,allParams)

	#Generating commands to run DNAcopy segmentation
	cmds = Cmd[]
	for b in bams
		n = replace(b,".bam"=>"")
		segCmd = `$(rscriptPath) $(dnacopyPath) --IN=$(outFolder)/chgTables/$n.bins.tsv --OUT=$(outFolder)/segTables/$n.seg.tsv`
		push!(cmds,segCmd)
	end

	#Running DNAcopy segmentation and updating tables with segment data
	pmap(runCmd,cmds)

	update_bin_segs("$(outFolder)/chgTables","$(outFolder)/segTables","$(outFolder)/finalTables")
	plot_plasmacna_folder("$(outFolder)/finalTables","$(outFolder)/plots")
end

main()
