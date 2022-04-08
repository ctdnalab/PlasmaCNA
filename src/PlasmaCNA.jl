module PlasmaCNA
	using FASTX
	using BioSequences
	using BioAlignments
	using GZip

	using CSV
	using DataFrames
	using Query

	using Statistics
	using StatsBase
	using Loess

	using Gadfly
	using Cairo
	using Colors

	include("utils.jl")
	include("bam_processing.jl")
	include("reference_proccessing.jl")
	include("segmentation_R.jl")
	include("plots.jl")

	#Generating and processing reference genome bins
	export get_fixed_windows
	export get_fixedbin_coords
	export get_flexbin_coords
	export annotate_bin_mappability!
	export annotate_bin_composition!

	#Reading and normalizing bin depths from bam/cram
	export write_bedtools_tempfiles
	export get_bedtools_depth!
	export normalize_depth!

	export read_segs_dnacopy
	export annotate_segments!
	export update_bin_segs

	export hasintersect
	export tbl2wig
	export plot_plasmacna_folder
	export plot_plasmacna_segs
end
