#!/usr/bin/env ruby
#
require 'rubygems'
require 'optparse'
require 'ostruct'
require 'pp'
require 'fileutils'
require 'inifile'

FASTA_DIR = 'fasta'
DATABASE_DIR = 'database'
BOWTIE_BUILD_1_EXE = 'bowtie-build-l'
BOWTIE_BUILD_2_EXE = 'bowtie2-build-l'
COLORSPACE_INDEX = '_colorspace_index'
DNASPACE_INDEX = '_dnaspace_index'
BOWTIE_THREADS = 4
REFERENCE_INI_FILE_PREFIX = '_reference.ini'

REFERENCE_INFO_SECTION = 'reference_info'
REFERENCE_NAME_STR = 'reference_name'
REFERENCE_ENTRY_TYPE_STR = 'entry_type'
REFERENCE_DATE_STR = 'reference_date'
DATABASE_TYPE_STR = 'database_type'
HAS_LITE_INFO = 'has_lite_info'

REFERENCE_FILE_SECTION = 'reference_file'
IS_LARGE_REFERENCE_STR = 'is_large_reference'
REFERENCE_DATABASE_DIR_STR = 'database_dir'
REFERENCE_FASTA_DIR_STR = 'fasta_dir'
REFERENCE_FASTA_FILE_COUNT_STR = 'fasta_file_count'
REFERENCE_FASTA_FILENAME_STR = 'fasta_filename'

REFERENCE_BOWTIE_INDEX_SECTION = 'bowtie_index'
REFERENCE_BOWTIE2_INDEX_SECTION = 'bowtie2_index'
REFERENCE_INDEX_COUNT = 'index_count'
REFERENCE_IS_COLOR_SPACE_BOWTIE_INDEXED = 'is_color_space_indexed'
REFERENCE_IS_DNA_SPACE_BOWTIE_INDEXED = 'is_DNA_space_indexed'
REFERENCE_COLOR_SPACE_BOWTIE_INDEX_PREFIX_NAME_STR = 'color_space_bowtie_index_prefix_name'
REFERENCE_DNA_SPACE_BOWTIE_INDEX_PREFIX_NAME_STR = 'dna_space_bowtie_index_prefix_name'
	
class OptParse
	# Return a structure describing the options.
	def self.parse(args)
		# The options specified on the command line will be collected in *options*.
		# We set default values here.
		options = OpenStruct.new
		options.inputFastaFileName = ""
		options.referenceRootDir = ""
		options.referenceName = ""
		options.indexBowtie1 = true
		options.indexBowtie2 = true
		options.threads = 1
		
		opts = OptionParser.new
		opts.banner = "Usage: MeteorReferenceBuilder.rb [options]"

		# Mandatory arguments.	   
		opts.separator ""
		opts.separator "Mandatory arguments:"
		opts.on("-i", "--input input", String, "input fasta filename") { |v| options.inputFastaFileName = v.to_s }		
		opts.on("-p", "--refRootDir refRootDir", String, "root-path of the reference repository") { |v| options.referenceRootDir = v.to_s }
		opts.on("-n", "--refName refName", String, "name of the reference (ansi-string without space)") { |v| options.referenceName = v.to_s }
		
		# Optional arguments
		opts.separator ""
		opts.separator "Optional arguments:"
		opts.on("-1", "--no-bowtie1 no-bowtie1", String, "no index for bowtie1") { |v| options.indexBowtie1 = false }
		opts.on("-2", "--no-bowtie2 no-bowtie2", String, "no index for bowtie2") { |v| options.indexBowtie2 = false }
		opts.on("-t", "--threads threads", Integer, "thread count for bowtie2 (if available)") { |v| options.threads = v.to_i }
		
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

		if options.inputFastaFileName == ""
			STDERR.puts "Missing argument : input fasta filename [-i]"
			exit
		end
		
		if options.referenceRootDir==""
			STDERR.puts "Missing argument : reference repository root path [-p]"
			exit
		end

		if options.referenceName==""
			STDERR.puts "Missing argument : reference name [-n]"
			exit
		end
		
		if !File.exists?(File.expand_path(options.inputFastaFileName))
			STDERR.puts "Error : input fasta file " + options.inputFastaFileName + " not found"
			exit
		end
		
		if !File.exists?(File.expand_path(options.referenceRootDir))
			STDERR.puts "Error : reference root directory " + options.referenceRootDir + " not found"
			exit
		end
				
		return options
	end  # parse()
end  # class OptParse

### MAIN ###

options = OptParse.parse(ARGV)

# 1- Prepare repository

referenceDir = File.join(options.referenceRootDir, options.referenceName)
FileUtils.mkdir_p(referenceDir)
FileUtils.mkdir_p(File.join(referenceDir, FASTA_DIR))
FileUtils.mkdir_p(File.join(referenceDir, DATABASE_DIR))

# 2- Read input fasta file

# fasta file with headear as numeric index
outputFastaFileName = File.join(referenceDir, FASTA_DIR, [options.referenceName, '.fa'].join(''))
# tabulated text file with numeric index and gene length in columns
outputAnnotationFileName = File.join(referenceDir, DATABASE_DIR, [options.referenceName, '_lite_annotation'].join(''))

### FG 2011-11-12, avoid overwriting already done numeric indexing
if options.inputFastaFileName != outputFastaFileName
	
	inputFastaFile = File.open(options.inputFastaFileName)
	outputFastaFile = File.new(outputFastaFileName,"w")
	outputAnnotationFile = File.new(outputAnnotationFileName,"w")

	geneID = 1 # initialize numeric gene index

	$/ = ">" 
	inputFastaFile.gets 
	while rec = inputFastaFile.gets
		rec.chomp!
		nl = rec.index("\n") 
		header = rec[0..nl-1] 
		seq = rec[nl+1..-1] 
		seq.gsub!(/\n/,'') 
	
		outputFastaFile.puts [">", geneID.to_s].join("")
		outputFastaFile.puts seq
	
		outputAnnotationFile.puts [geneID.to_s, seq.length.to_s].join("\t")
	
		geneID = geneID + 1
	end 

	inputFastaFile.close
	outputFastaFile.close
	outputAnnotationFile.close
end

# 3- prepare index for bowtie
Dir.chdir(File.join(referenceDir, FASTA_DIR))

if options.indexBowtie1
	# bowtie1 - colorspace
	cmd = [BOWTIE_BUILD_1_EXE, '-f', '-C']
	cmd.push(File.join(referenceDir, FASTA_DIR, [options.referenceName, '.fa'].join('')))
	cmd.push(File.join(referenceDir, FASTA_DIR, [options.referenceName, COLORSPACE_INDEX].join('')))
	puts "Executing command:"
	puts cmd.join(' ')
	system(cmd.join(' '))
	#puts cmd.join(' ')

	# bowtie1 - dnaspace
	cmd = [BOWTIE_BUILD_1_EXE, '-f']
	cmd.push(File.join(referenceDir, FASTA_DIR, [options.referenceName, '.fa'].join('')))
	cmd.push(File.join(referenceDir, FASTA_DIR, [options.referenceName, DNASPACE_INDEX].join('')))
	puts "Executing command:"
	puts cmd.join(' ')
	system(cmd.join(' '))
end

if options.indexBowtie2
	# bowtie2 - dnaspace
	cmd = [BOWTIE_BUILD_2_EXE, '-f']
	if options.threads > 1
		cmd.push('--threads', options.threads.to_s)
	end
	cmd.push(File.join(referenceDir, FASTA_DIR, [options.referenceName, '.fa'].join('')))
	cmd.push(File.join(referenceDir, FASTA_DIR, [options.referenceName, DNASPACE_INDEX].join('')))
	puts "Executing command:"
	puts cmd.join(' ')
	system(cmd.join(' '))
end

# 4- Edit INI file
referenceINIFileName = File.join(referenceDir, options.referenceName + REFERENCE_INI_FILE_PREFIX)
if !File.exist?(referenceINIFileName)
	referenceINIFile = IniFile.new(:filename => referenceINIFileName)
else
	referenceINIFile = IniFile.load(referenceINIFileName)
end

referenceINIFile[REFERENCE_INFO_SECTION][REFERENCE_NAME_STR] = options.referenceName
referenceINIFile[REFERENCE_INFO_SECTION][REFERENCE_ENTRY_TYPE_STR] = 'fragment'
referenceINIFile[REFERENCE_INFO_SECTION][DATABASE_TYPE_STR] = 'text'
referenceINIFile[REFERENCE_INFO_SECTION][HAS_LITE_INFO] = 1
referenceINIFile[REFERENCE_INFO_SECTION][REFERENCE_DATE_STR] = Time.now.strftime("%Y-%m-%d")

referenceINIFile[REFERENCE_FILE_SECTION][REFERENCE_DATABASE_DIR_STR] = DATABASE_DIR
referenceINIFile[REFERENCE_FILE_SECTION][REFERENCE_FASTA_DIR_STR] = FASTA_DIR
referenceINIFile[REFERENCE_FILE_SECTION][REFERENCE_FASTA_FILE_COUNT_STR] = 1
referenceINIFile[REFERENCE_FILE_SECTION][[REFERENCE_FASTA_FILENAME_STR,'1'].join('_')] = File.basename(outputFastaFileName)

if options.indexBowtie1
	referenceINIFile[REFERENCE_BOWTIE_INDEX_SECTION][IS_LARGE_REFERENCE_STR] = 1
	referenceINIFile[REFERENCE_BOWTIE_INDEX_SECTION][REFERENCE_IS_COLOR_SPACE_BOWTIE_INDEXED] = 1
	referenceINIFile[REFERENCE_BOWTIE_INDEX_SECTION][[REFERENCE_COLOR_SPACE_BOWTIE_INDEX_PREFIX_NAME_STR,'1'].join('_')] = [options.referenceName, COLORSPACE_INDEX].join('')
	referenceINIFile[REFERENCE_BOWTIE_INDEX_SECTION][REFERENCE_IS_DNA_SPACE_BOWTIE_INDEXED] = 1
	referenceINIFile[REFERENCE_BOWTIE_INDEX_SECTION][[REFERENCE_DNA_SPACE_BOWTIE_INDEX_PREFIX_NAME_STR,'1'].join('_')] = [options.referenceName, DNASPACE_INDEX].join('')
end

if options.indexBowtie2
	referenceINIFile[REFERENCE_BOWTIE2_INDEX_SECTION][IS_LARGE_REFERENCE_STR] = 1
	referenceINIFile[REFERENCE_BOWTIE2_INDEX_SECTION][REFERENCE_IS_DNA_SPACE_BOWTIE_INDEXED] = 1
	referenceINIFile[REFERENCE_BOWTIE2_INDEX_SECTION][[REFERENCE_DNA_SPACE_BOWTIE_INDEX_PREFIX_NAME_STR,'1'].join('_')] = [options.referenceName, DNASPACE_INDEX].join('')
end

referenceINIFile.save

STDOUT.puts "Meteor reference building done!"      
STDOUT.flush
