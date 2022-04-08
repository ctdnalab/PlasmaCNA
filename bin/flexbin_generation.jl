using PlasmaCNA
using QuickArgParse
using CSV

function process_chroms(chromString::String;chrtag="")
	chromData = [string(x) for x in split(chromString,",")]
	finalChroms = String[]
	for c in chromData
		if occursin("-",c) == true
			cRange = [parse(Int,x) for x in split(c,'-')]
			cNums = sort(collect(range(cRange[1],stop=cRange[2])))
			cNames = [string(x) for x in cNums]
			cNames = ["$(chrtag)$(x)" for x in cNames]
			for n in cNames
				push!(finalChroms,n)
			end
		else
			c = string(c)
			push!(finalChroms,c)
		end
	end
	return finalChroms
end


function main()
	req = ["BinNumber","MapFile","RefGenome","OutputName","OutputDir"]
	reqhelp = ["Number of bins to generate","Mappability wig file from GenMap",
		"Fasta file of the reference genome","Name for output files","Directory for output files"]
	opts= ["chroms","chrtag","minmap"]
	optsHelp = ["Specify chroms to build bins on, in the same order as the ref genome fasta file.",
		"String to add to numbered chromosomes. Can be NONE",
		"Minimum mappability for a good base"]
	optsDefault = ["1-22,chrX,chrY","chr","0.9"]
	title = "Flexible Size Bin Generation"
	desc = "Generates flexible size bins across the genome\nwith approximately equal number mappable bases\nand annotates with GC content and mappability"

	R = process_reqs(req,reqHelp=reqhelp,optTag=opts,optHelp=optsHelp,
		optDefault=optsDefault,title=title,desc=desc)
	build_usage!(R)
	A = parse_args(R)

	binNum = parse(Int,A["BinNumber"])
	minMap = parse(Float64,A["minmap"])
	refGenome = A["RefGenome"]
	mapFile = A["MapFile"]
	addTag = A["chrtag"]
	outName = A["OutputName"]
	outDir = A["OutputDir"]
	if uppercase(addTag) == "NONE"
		addTag = ""
	end

	goodChroms = process_chroms(A["chroms"],chrtag=addTag)

	run(`mkdir -p $outDir`)
	bins = get_flexbin_coords(binNum,minMap,mapFile,goodChrom=goodChroms)
	CSV.write("$(outDir)/$(outName).flexbin.$(binNum).raw.tsv",bins,delim='\t')
	annotate_bin_composition!(refGenome,bins)
	CSV.write("$(outDir)/$(outName).flexbin.$(binNum).comp.tsv",bins,delim='\t')
	annotate_bin_mappability!(mapFile,bins)
	CSV.write("$(outDir)/$(outName).flexbin.$(binNum).final.tsv",bins,delim='\t')
end


main()