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
		options.CSFastaQualDir = ""
		options.projectName = ""
		options.maskSampleName = ""
		options.techno = "SOLiD"
		options.isDispatched = false
		
		opts = OptionParser.new
		opts.banner = "Usage: MeteorImportCSFastaQual.rb [options]"

		# Mandatory arguments.	   
		opts.separator ""
		opts.separator "Mandatory arguments:"
		opts.on("-i", "--csfastaDir csfastaDir", String, "CSFasta directory") { |v| options.CSFastaQualDir = v.to_s }		
		opts.on("-p", "--projectName projectName", String, "project name (ansi-string without space)") { |v| options.projectName = v.to_s }
		opts.on("-m", "--maskSampleName maskSampleName", String, "regular expression for extracting sample name") { |v| options.maskSampleName = v.to_s }
		
		# Optional arguments
		opts.separator ""
		opts.separator "Optional arguments:"
		opts.on("-d", "files are already dispatched in directory ([-m] option not necessary)") { options.isDispatched = true }

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

		if options.CSFastaQualDir == ""
			STDERR.puts "Missing argument : CSFasta directory [-i]"
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
		
		if !File.exists?(File.expand_path(options.CSFastaQualDir))
			STDERR.puts "Error : CSFasta directory " + options.CSFastaQualDir + " not found"
			exit
		end
				
		return options
	end  # parse()
end  # class OptParse

### MAIN ###

options = OptParse.parse(ARGV)

# process all CSFasta files
if options.isDispatched
	csfasta_file_array = Dir.glob(File.join(options.CSFastaQualDir,"*","*.csfasta")) # get csfasta in subdirectories
else
	csfasta_file_array = Dir.glob(File.join(options.CSFastaQualDir, "*.csfasta"))
end

csfasta_file_array.each do |csfasta_file|
	STDOUT.puts ["Import", csfasta_file].join(' ')
	
	# retrieve quality file
	qual_file = csfasta_file.gsub('csfasta','QV.qual')
	# test if qual file exists, else add only qual extension
	if !File.exists?(qual_file)
		qual_file = csfasta_file.gsub('csfasta','qual')
		# retest if new file exists, else qual_file is blank
		if !File.exists?(qual_file)
			qual_file = ''
		end
	end
	
	full_sample_name = File.basename(csfasta_file)
	full_sample_name = full_sample_name.gsub('.csfasta','')
	
	tag = ''
	
	if options.isDispatched
		sample_name = File.basename(File.dirname(csfasta_file))
	else
		# split full sample name (in fact library/run name) in order to extract sample_name according to regex mask 
		full_sample_name_array = /#{options.maskSampleName}/.match(full_sample_name)
		sample_name = full_sample_name_array[0]
	end
	
	if options.isDispatched
		sample_dir = File.dirname(csfasta_file)
	else
		# create directory for the sample and move fastq file into
		sample_dir = File.join(File.dirname(csfasta_file), sample_name)
		FileUtils.mkdir_p(sample_dir)
		FileUtils.mv(csfasta_file,sample_dir)
		if qual_file != ''
			FileUtils.mv(qual_file,sample_dir)
		end
	end
	
	# open csfasta file and calculate read length
	f = File.open(csfasta_file, "r")
	while true
		s = f.readline
		if s.chars.first == '>'	
			s = f.readline # read sequence
			s = s.strip
			read_length = s.length - 1
			break
		end
	end
	f.close
	
	# create and edit census_stage_0.ini INI file (for Meteor)  
	Dir.chdir(sample_dir)

	census_stage_file_name = full_sample_name + "_census_stage_0.ini"
	census_stage_file = File.new(census_stage_file_name,"w")

	census_stage_file.puts "[sample_info]"
	census_stage_file.puts "sample_name="+sample_name
	census_stage_file.puts "condition_name=NA"
	census_stage_file.puts "project_name="+options.projectName
	census_stage_file.puts "sequencing_date=1900-01-01"
	census_stage_file.puts "sequencing_device="+options.techno
	census_stage_file.puts "census_status=0"
	census_stage_file.puts "read_length="+read_length.to_s
	census_stage_file.puts "tag="+tag
	census_stage_file.puts "full_sample_name="+full_sample_name
	census_stage_file.puts "[sample_file]"
	census_stage_file.puts "csfasta_file="+File.basename(csfasta_file)
	if qual_file != ''
		census_stage_file.puts "qual_file="+File.basename(qual_file)
	else
		census_stage_file.puts "qual_file="
	end
	census_stage_file.close

	Dir.chdir(options.CSFastaQualDir)
end

STDOUT.puts "Importation done!"
STDOUT.flush

