#! /usr/bin/env python
import argparse, math, os, subprocess, sys, tempfile, shutil


def select_parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description='Run CMash and select a subset of the whole database to align to.')
	parser.add_argument('reads', help='Path to reads file.')
	parser.add_argument('data', help='Path to data/ directory with the files from setup_data.sh')
	parser.add_argument('--metalign_results', default='NONE', help='Give location of Metalign-seeding results if already done.')
	parser.add_argument('--cutoff', type=int, default=0.0001, help='Seed count cutoff value. Default is 0.0001.')
	parser.add_argument('--db', default='AUTO', help='Where to write subset database. Default: temp_dir/subset_db.fna')
	parser.add_argument('--db_dir', default='AUTO', help='Directory with all organism files in the full database.')
	parser.add_argument('--dbinfo_in', default='AUTO', help='Specify location of db_info file. Default is data/db_info.txt')
	parser.add_argument('--dbinfo_out', default='AUTO',
		help='Where to write subset db_info. Default: temp_dir/subset_db_info.txt')
	parser.add_argument('--input_type', default='AUTO', choices=['fastq', 'fasta', 'AUTO'],
		help='Type of input file (fastq/fasta). Default: try to auto-determine')
	parser.add_argument('--keep_temp_files', action='store_true', help='Retain KMC files after this script finishes.')
	parser.add_argument('--strain_level', action='store_true',
		help='Include all strains above cutoff. Default: 1 strain per species.')
	parser.add_argument('--temp_dir', default='AUTO/', help='Directory to write temporary files to.')
	parser.add_argument('--threads', type=int, default=4, help='How many compute threads for KMC to use. Default: 4')
	parser.add_argument('--minimap_n', type=int, default=3, help='Minimap: Discard chains consisting of <INT> number of minimizers')
	parser.add_argument('--mmi_dir',  default = 'AUTO', help='Minimap-Threader: Directory containing all mmi-files')
	parser.add_argument('--translation',  default = 'AUTO', help='Accession to taxid for subset DB generation')
	parser.add_argument('--filter', default='base-counting', choices=['adjacency-filter', 'base-counting', 'edlib', 'grim_original', 'grim_original_tweak', 'hd', 'magnet', 'qgram', 'shd', 'shouji', 'sneakysnake'])
	parser.add_argument('--edit_dist_threshold', type=int, default=15, help='-r edit distance threshold for minimap2.')
	args = parser.parse_args()
	return args


def read_dbinfo(args):
	taxid2info = {}
	with(open(args.dbinfo_in, 'r')) as infile:
		infile.readline()  # skip header line
		for line in infile:
			splits = line.strip().split('\t')
			acc, taxid = splits[0], splits[2]
			if taxid not in taxid2info:
				# first element stores all accessions for this taxid
				taxid2info[taxid] = [[splits[0]], splits[1]]
				taxid2info[taxid].extend(splits[3:])
			else:
				taxid2info[taxid][0].append(acc)
	return taxid2info

def run_minimap_and_cutoff(args, taxid2info):
	if args.metalign_results == 'NONE':
		seed_count = subprocess.check_output(["../MetaFast/ContainmentSearch/cs", "-n", str(args.minimap_n), args.mmi_dir, 
		args.reads, args.translation, args.temp_dir + "ContainmentResults.csv"]).decode('UTF-8').splitlines()[-1]
		print(seed_count)

	else:
		seed_count = args.metalign_results

	organisms_to_include, species_included = [], {}
	with(open(seed_count, 'r')) as metalign_results:
		#metalign_results.readline()  # skip header line
		for line in metalign_results:
			splits = line.strip().split(',')
			organism, containment_index = splits[0], float(splits[-1])
			if containment_index >= args.cutoff:
				if not args.strain_level:
					taxid = organism.split('taxid_')[1].split(
						'_genomic.fna')[0].replace('_', '.')
					species = taxid2info[taxid][3].split('|')[-2]
					if species not in species_included or species == '':
						species_included[species] = 1
					else:
						continue
				organisms_to_include.append(organism)
	return organisms_to_include

def make_db_and_dbinfo(args, organisms_to_include, taxid2info):
	open(args.db, 'w').close()  # clear cmash results; no longer needed
	with(open(args.db, 'a')) as outfile:
		for organism in organisms_to_include:
			organism_fname = args.db_dir + organism
			# write organisms to full db via cat to append-mode file handler
			subprocess.Popen(['zcat', organism_fname], stdout=outfile).wait()

	with(open(args.dbinfo_out, 'w')) as outfile:
		# write header lines
		outfile.write('Accesion\tLength\tTaxID\tLineage\tTaxID_Lineage\n')
		outfile.write('Unmapped\t0\tUnmapped\t|||||||Unmapped\t|||||||Unmapped\n')
		for organism in organisms_to_include:
			taxid = organism.split('taxid_')[1].split(
				'_genomic.fna')[0].replace('_', '.')
			length = taxid2info[taxid][1]
			namelin, taxlin = taxid2info[taxid][2], taxid2info[taxid][3]
			for acc in taxid2info[taxid][0]:
				outfile.write('\t'.join([acc,length,taxid,namelin,taxlin]) + '\n')


def select_main(args = None):
	if args == None:
		args = select_parseargs()
	if not args.data.endswith('/'):
		args.data += '/'
	if args.db_dir == 'AUTO':
		args.db_dir = args.data + 'organism_files/'
	if not args.db_dir.endswith('/'):
		args.db_dir += '/'
	if args.temp_dir == 'AUTO/':
		args.temp_dir = tempfile.mkdtemp(prefix=args.data)
	if not args.temp_dir.endswith('/'):
		args.temp_dir += '/'
	if not os.path.exists(args.temp_dir):
		os.makedirs(args.temp_dir)
	if args.dbinfo_in == 'AUTO':
		args.dbinfo_in = args.data + 'db_info.txt'
	if args.dbinfo_out == 'AUTO':
		args.dbinfo_out = args.temp_dir + 'subset_db_info.txt'
	if args.db == 'AUTO':
		args.db = args.temp_dir + 'subset_db.fna'
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

	taxid2info = read_dbinfo(args)
	organisms_to_include = run_minimap_and_cutoff(args, taxid2info)
	make_db_and_dbinfo(args, organisms_to_include, taxid2info)


if __name__ == '__main__':
	args = select_parseargs()
	select_main(args)
#
