using QuickArgParse
using PlasmaCNA
using CSV
using DataFrames

function fixbin_generation(chromSize::String,binSize::Int,mapFile::String,refGenome::String,
	outName::String,outDir::String;wig::Bool=false)

	bins = get_fixedbin_coords(chromSize,binSize*1000)
	CSV.write("$(outDir)/$(outName).fixedbin.$(binSize)kb.raw.tsv",bins,delim='\t')
	annotate_bin_composition!(refGenome,bins)
	CSV.write("$(outDir)/$(outName).fixedbin.$(binSize)kb.comp.tsv",bins,delim='\t')
	annotate_bin_mappability!(mapFile,bins)
	CSV.write("$(outDir)/$(outName).fixedbin.$(binSize)kb.final.tsv",bins,delim='\t')

	if wig == true
		df = CSV.File("$(outDir)/$(outName).fixedbin.$(binSize)kb.final.tsv") |> DataFrame
		tbl2wig(df,:gc,"$(outDir)/$(outName).fixedbin.$(binSize)kb.gc.wig",returnType=Float64)
		tbl2wig(df,:map,"$(outDir)/$(outName).fixedbin.$(binSize)kb.mappability.wig",returnType=Float64)
	end
end

function main()
	req = ["ChromSize","BinSize","MapFile","RefGenome","OutputName","OutputDir"]
	reqhelp = ["Tab delimited file of chromosomes to use and their sizes in bp","Size of bins in kb",
		"Mappability wig file from GenMap","Fasta file of the reference genome","Name for output files",
		"Directory for output files"]
	flags = ["w"]
	flagHelp = ["Generate separate wig files for GC and mappability (ie for ichorCNA)"]
	title = "Fixed Size Bin Generation"
	desc = "Generates fixed size bins across the genome and\nannotates with GC content and mappability"

	R = process_reqs(req,reqHelp=reqhelp,flags=flags,flagHelp=flagHelp,title=title,desc=desc)
	build_usage!(R)
	A = parse_args(R)

	fixbin_generation(A["ChromSize"],parse(Int,A["BinSize"]),A["MapFile"],
		A["RefGenome"],A["OutputName"],A["OutputDir"],wig=A["w"])
end

main()