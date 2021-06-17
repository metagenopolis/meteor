#!/usr/bin/env ruby

require 'optparse'
require 'ostruct'
require 'pp'
require 'fileutils'

class OptParse
	# Return a structure describing the options.
	def self.parse(args)
		# The options specified on the command line will be collected in *options*.
		# We set default values here.
		options = OpenStruct.new
		options.fastqDir = ""
		options.projectName = ""
		options.maskSampleName = ""
		options.techno = ""
		options.isCompressed = false
		options.isDispatched = false
		
		opts = OptionParser.new
		opts.banner = "Usage: MeteorImportFastq.rb [options]"

		# Mandatory arguments.	   
		opts.separator ""
		opts.separator "Mandatory arguments:"
		opts.on("-i", "--fastqDir fastqDir", String, "fastq directory") { |v| options.fastqDir = v.to_s }		
		opts.on("-p", "--projectName projectName", String, "project name (ansi-string without space)") { |v| options.projectName = v.to_s }
		opts.on("-m", "--maskSampleName maskSampleName", String, "regular expression for extracting sample name") { |v| options.maskSampleName = v.to_s }
		opts.on("-t", "--techno techno", String, "sequencing technology (SOLiD, Illumina, Proton)") { |v| options.techno = v.to_s }
		
		# Optional arguments
		opts.separator ""
		opts.separator "Optional arguments:"
		opts.on("-c", "fastq files are compressed") { options.isCompressed = true }
		opts.on("-d", "fastq files are already dispatched in directory ([-m] option not necessary)") { options.isDispatched = true }

		# Help argument
		opts.separator ""
		opts.on_tail("-h", "Show this message") do
			puts opts
			exit
		end

		# Display usage if no arguments given
		if args.empty?
			puts opts
			exit
		end

		begin
			opts.parse!(args)
		rescue OptionParser::ParseError => e
			puts e
		end

		if options.fastqDir == ""
			STDERR.puts "Missing argument : fastq directory [-i]"
			exit
		end
		
		if options.projectName==""
			STDERR.puts "Missing argument : project name [-p]"
			exit
		end
		
		if !options.isDispatched
			if options.maskSampleName==""
				STDERR.puts "Missing argument : mask for sample name [-m]"
				exit
			end
		end
		
		if !File.exists?(File.expand_path(options.fastqDir))
			STDERR.puts "Error : fastq directory " + options.fastqDir + " not found"
			exit
		end
				
		return options
	end  # parse()
end  # class OptParse

### MAIN ###

options = OptParse.parse(ARGV)

# process all fastq files
if options.isDispatched
	fastq_file_array = Dir.glob(File.join(options.fastqDir,"*","*.{fq,fastq}*")) # get fastq in subdirectories
else
	fastq_file_array = Dir.glob(File.join(options.fastqDir, "*.{fq,fastq}*"))
end

fastq_file_array.each do |fastq_file|
	STDOUT.puts ["Import", fastq_file].join(' ')
	full_sample_name = File.basename(fastq_file)

	if options.isCompressed
		full_sample_name = full_sample_name.gsub(/\.(bz|bz2|bzip2|gz|gzip|zip|tar.gz|tar.gzip|tar.bz|tar.bz2|tar.bzip2|tar.zip)\z/,'')
	end
	full_sample_name = full_sample_name.gsub(/\.(fq|fastq)\z/,'')
	
	# extract paired-end info
	if ( full_sample_name =~ /.*_single$/ )
		tag = "single"
	elsif ( full_sample_name =~ /.*_1$/ )
		tag = "1"
	elsif ( full_sample_name =~ /.*_2$/ )
		tag = "2"
	else
		tag = ""
	end
	
	if options.isDispatched
		sample_name = File.basename(File.dirname(fastq_file))
	else
		# split full sample name (in fact library/run name) in order to extract sample_name according to regex mask 
		full_sample_name_array = /#{options.maskSampleName}/.match(full_sample_name)
		sample_name = full_sample_name_array[0]
	end
	
	if options.isDispatched
		sample_dir = File.dirname(fastq_file)
	else
		# create directory for the sample and move fastq file into
		sample_dir = File.join(File.dirname(fastq_file), sample_name)
		FileUtils.mkdir_p(sample_dir)
		FileUtils.mv(fastq_file,sample_dir)
	end
	
	# create and edit census_stage_0.ini INI file (for Meteor)  
	Dir.chdir(sample_dir)

	census_stage_file_name = full_sample_name + "_census_stage_0.ini"
	census_stage_file = File.new(census_stage_file_name,"w")

	# census_ini_file description with one line for each used xsq file
	census_stage_file.puts "[sample_info]"
	census_stage_file.puts "sample_name="+sample_name
	census_stage_file.puts "condition_name=NA"
	census_stage_file.puts "project_name="+options.projectName
	census_stage_file.puts "sequencing_date=1900-01-01"
	census_stage_file.puts "sequencing_device="+options.techno
	census_stage_file.puts "census_status=0"
	census_stage_file.puts "read_length=-1"
	census_stage_file.puts "tag="+tag
	census_stage_file.puts "full_sample_name="+full_sample_name
	census_stage_file.puts "[sample_file]"
	census_stage_file.puts "fastq_file="+File.basename(fastq_file)
	if options.isCompressed
		census_stage_file.puts "is_compressed=1"
	end
	census_stage_file.close

	Dir.chdir(options.fastqDir)
end

STDOUT.puts "Importation done!"
STDOUT.flush

