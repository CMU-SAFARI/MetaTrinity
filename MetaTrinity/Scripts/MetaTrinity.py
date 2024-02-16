#! /usr/bin/env python
import argparse, math, os, subprocess, sys, tempfile, time
# Import metalign modules
import containment_search as select
import read_mapping as mapper


def metalign_parseargs():  # handle user arguments
	parser = argparse.ArgumentParser(description='Runs full metalign pipeline on input reads file(s).')
	parser.add_argument('reads', help='Path to reads file.')
	parser.add_argument('data', help='Path to data/ directory with the files from setup_data.sh')
	parser.add_argument('--cutoff', type=int, default=0.0001, help='Minimap cutoff value. Default is 0.0001.')
	parser.add_argument('--db_dir', default = 'AUTO', help='Directory with all organism files in the full database.')
	parser.add_argument('--dbinfo_in', default='AUTO', help='Location of db_info file. Default: data/db_info.txt')
	parser.add_argument('--keep_temp_files', action='store_true', help='Retain KMC files after this script finishes.')
	parser.add_argument('--input_type', default='AUTO', choices=['fastq', 'fasta', 'AUTO'],help='Type of input file (fastq/fasta). Default: try to auto-determine')
	parser.add_argument('--length_normalize', action='store_true', help='Normalize abundances by genome length.')
	parser.add_argument('--low_mem', action='store_true',help='Run in low memory mode, with inexact multimapped processing.')
	parser.add_argument('--metalign_results', default='NONE', help='Give location of Metalign-seeding results if already done.')
	parser.add_argument('--min_abundance', type=float, default=10**-4,help='Minimum abundance for a taxa to be included in the results. Default: 10^(-4).')
	parser.add_argument('--no_quantify_unmapped', action='store_true',help='Do not factor in unmapped reads in abundance estimation.')
	parser.add_argument('--output', default='abundances.tsv', help='Output abundances file. Default: abundances.tsv')
	parser.add_argument('--pct_id', type=float, default=0.5,help='Minimum percent identity from reference to count a hit.')
	parser.add_argument('--precise', action='store_true',help='Run in precise mode. Overwrites --read_cutoff and --min_abundance to 100 and 0.1.')
	parser.add_argument('--rank_renormalize', action='store_true',help='Renormalize abundances to 100 pct. at each rank, e.g if an organism has a species but not genus label.')
	parser.add_argument('--read_cutoff', type=int, default=1, help='Number of reads to count an organism as present.')
	parser.add_argument('--sampleID', default='NONE', help='Sample ID for output. Defaults to input file name(s).')
	parser.add_argument('--sensitive', action='store_true', help='Run in sensitive mode. Sets --cutoff value to 0.0.')
	parser.add_argument('--strain_level', action='store_true', help='Profile strains (off by default).')
	parser.add_argument('--temp_dir', default='AUTO/', help='Directory to write temporary files to.')
	parser.add_argument('--threads', type=int, default=4, help='Number of compute threads for Minimap2/KMC. Default: 4')
	parser.add_argument('--verbose', action='store_true', help='Print verbose output.')
	parser.add_argument('--minimap_n', type=int, default=3)
	parser.add_argument('--mmi_dir',  default = 'AUTO', help='Minimap-Threader: Directory containing all mmi-files')
	parser.add_argument('--translation',  default = 'AUTO', help='Accession to taxid for subset DB generation')
	parser.add_argument('--filter', default = 'base-counting', choices=['adjacency-filter', 'base-counting', 'edlib', 'grim_original', 'grim_original_tweak', 'hd', 'magnet', 'qgram', 'shd', 'shouji', 'sneakysnake'], help='algorithm for read mapping')
	parser.add_argument('--edit_dist_threshold', type=int, default=15, help='-r edit distance threshold for minimap2.')
	args = parser.parse_args()
	return args


def main():
	args = metalign_parseargs()
	if not args.data.endswith('/'):
		args.data += '/'
	if args.temp_dir == 'AUTO/':
		args.temp_dir = tempfile.mkdtemp(prefix=args.data)
	if not args.temp_dir.endswith('/'):
		args.temp_dir += '/'
	# Set arguments that default to AUTO
	if args.dbinfo_in == 'AUTO':
		args.dbinfo_in = args.data + 'db_info.txt'
	if args.db_dir == 'AUTO':
		args.db_dir = args.data + 'organism_files/'
	if args.input_type == 'AUTO':
		splits = args.reads.split('.')
		if splits[-1] == 'gz':  # gz doesn't help determine file type
			splits = splits[:-1]
		if splits[-1] in ['fq', 'fastq']:
			args.input_type = 'fastq'
		elif splits[-1] in ['fa', 'fna', 'fasta']:
			args.input_type = 'fasta'
		else:
			sys.exit('Could not auto-determine file type. Use --input_type.')

	# handle precise and sensitive modes
	if args.sensitive and args.precise:
		sys.exit('You cannot use both --sensitive and --precise.')
	if args.sensitive:
		args.cutoff = 0.0
	elif args.precise:
		args.read_cutoff = 100
		args.min_abundance = 0.1

	# Ensure that arguments agree between scripts
	args.db = args.temp_dir + 'subset_db.fna'
	args.dbinfo = args.temp_dir + 'subset_db_info.txt'
	args.dbinfo_out = args.dbinfo
	args.infiles = [args.reads]  # read_mapping expects a list

	# Run the database selection and map/profile routines
	select.select_main(args)  # runs containment_search routine
	mapper.map_main(args)  # runs read_mapping routine
	if not args.keep_temp_files:  # clean up
		subprocess.Popen(['rm', '-r', args.temp_dir]).wait()


if __name__ == '__main__':
	main()
#
