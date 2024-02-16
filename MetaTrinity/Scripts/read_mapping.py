#! /usr/bin/env python
import argparse, os, subprocess, sys, time


start = time.time()  # start a program timer
RANKS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']


def echo(msg, verbose):
	if not verbose:  # only print if args.verbose flag was used
		return
	global start
	seconds = time.time() - start
	m, s = divmod(seconds, 60)
	h, m = divmod(m, 60)
	hms = "%02d:%02d:%02d" % (h, m, s)
	print('['+hms+'] ' + msg)


def profile_parseargs():  # handle user arguments
	parser = argparse.ArgumentParser(description='Compute abundance estimations for species in a sample.')
	parser.add_argument('infiles', nargs='+', help='sam or reads file(s) (space-delimited if multiple). Required.')
	parser.add_argument('data', help='Path to data/ directory with the files from setup_data.sh')
	parser.add_argument('--db', default='NONE', help='Path to database from containment_search. Required if read files given')
	parser.add_argument('--dbinfo', default='AUTO', help='Location of db_info file. Default: data/db_info.txt')
	parser.add_argument('--input_type', default='AUTO', choices=['fastq', 'fasta', 'sam', 'AUTO'],
		help='Type of input file (fastq/fasta/sam). Default: try to automatically determine')
	parser.add_argument('--length_normalize', action='store_true', help='Normalize abundances by genome length.')
	parser.add_argument('--low_mem', action='store_true',
		help='Run in low memory mode, with inexact multimapped processing.')
	parser.add_argument('--min_abundance', type=float, default=10**-4,
		help='Minimum abundance for a taxa to be included in the results. Default: 10^(-4).')
	parser.add_argument('--rank_renormalize', action='store_true',
		help='Renormalize abundances to 100 pct. at each rank, e.g if an organism has a species but not genus label.')
	parser.add_argument('--output', default='abundances.tsv', help='Output abundances file. Default: abundances.txt')
	parser.add_argument('--pct_id', type=float, default=0.5,
		help='Minimum percent identity from reference to count a hit.')
	parser.add_argument('--no_quantify_unmapped', action='store_true',
		help='Do not factor in unmapped reads in abundance estimation.')
	parser.add_argument('--read_cutoff', type=int, default=1, help='Number of reads to count an organism as present.')
	parser.add_argument('--sampleID', default='NONE', help='Sample ID for output. Defaults to input file name(s).')
	parser.add_argument('--threads', type=int, default=4, help='Number of compute threads for Minimap2. Default: 4')
	parser.add_argument('--verbose', action='store_true', help='Print verbose output.')
	parser.add_argument('--filter', default='base-counting', choices=['adjacency-filter', 'base-counting', 'edlib', 'grim_original', 'grim_original_tweak', 'hd', 'magnet', 'qgram', 'shd', 'shouji', 'sneakysnake'])
	parser.add_argument('--edit_dist_threshold', type=int, default=15, help='-r edit distance threshold for minimap2.')
	args = parser.parse_args()
	return args


# Given a taxonomic lineage, get the rank of the taxid
def get_taxid_rank(taxlin):
	# get number of empty tax ranks at end of taxlin
	end_empty, splits = 0, taxlin.split('|')
	for i in range(1, len(taxlin)+1):
		if splits[-i] == '':
			end_empty += 1
			continue
		break
	return RANKS[-(end_empty+1)]


# Parses information in dbinfo file, which maps NCBI accession to
# 	accession length, taxid, name lineage, and taxid lineage.
# Also maps all lowest-level taxids to organism length, which is the sum of
#  	accession lengths for accessions with that taxid, and lineage info
def get_acc2info(args):
	echo('Reading dbinfo file...', args.verbose)
	acc2info, taxid2info = {}, {}
	with(open(args.dbinfo, 'r')) as infofile:
		infofile.readline()  # skip header line
		for line in infofile:
			acc, acclen, taxid, namelin, taxlin = line.strip().split('\t')
			rank = get_taxid_rank(taxlin)
			if rank == 'strain' and acc != 'Unmapped':
				taxid += '.1'  # CAMI formatting specification
				taxlin += '.1'
			acclen = int(acclen)
			acc2info[acc] = [acclen, taxid, namelin, taxlin]
			if taxid in taxid2info:
				taxid2info[taxid][0] += acclen
			else:
				taxid2info[taxid] = [acclen, rank, namelin, taxlin]
	return acc2info, taxid2info


# Determine whether to count a mapping hit based on user specifications
# Return True if hit should be filtered out, False otherwise
def filter_line(args, splits):
	cigar = splits[5]  # quality of mapping determined by CIGAR string
	matched_len, total_len, cur = 0, 0, 0
	for ch in cigar:
		if not ch.isalpha():  # ch is a number
			cur = (cur * 10) + int(ch)
		else:  # ch is a letter
			if ch == 'M' or ch == '=':
				matched_len += cur
			total_len += cur
			cur = 0
	edit_distance = int(splits[11][5:])
	if float(matched_len) / float(total_len) < args.pct_id:
		return True
	return False  # if read passes quality checks, don't filter


# Parses FLAG field in sam line
def parse_flag(flag, cigar):
	# 1st or 2nd read in a pair, both false if single end
	pair1 = (flag & 1 != 0) and (flag & 64 != 0)
	pair2 = (flag & 1 != 0) and (flag & 128 != 0)
	chimeric = (flag & 2048 != 0)
	# "bad" here means unmapped or cigar string unavailable
	is_bad = (flag & 4 != 0) or (cigar == '*')  # or (flag & 2048 != 0)
	return pair1, pair2, chimeric, is_bad


# for paired end reads, return read hits to references that both paired ends hit
def intersect_read_hits(read_hits, pair1maps, pair2maps):
	if pair1maps == 0 or pair2maps == 0:  # one end unmapped means no intersect
		return [], []
	# gather all ref hits then partition into hits for each pair
	all_ref_hits = [hit[2] for hit in read_hits]
	pair1refs, pair2refs = all_ref_hits[:pair1maps], all_ref_hits[pair1maps:]
	# now intersect the lists and return hits to references in the intersect
	intersect = set([ref for ref in pair1refs if ref in pair2refs])
	#intersect = set([ref for ref in all_ref_hits])
	intersect_hits = [hit for hit in read_hits if hit[2] in intersect]
	return [intersect, intersect_hits]


# Given hits for a read, extract total hit length and quality scores for both
#  	ends of a read pair (if paired), and filter reads using user specification
def clean_read_hits(args, read_hits, pair1maps, pair2maps):
	hitlen, readquals = 0, ''
	filtered_hits = []
	for hit in range(len(read_hits)):
		# remove user-filtered & chimeric reads, update # of pair end map counts
		if filter_line(args, read_hits[hit]) or read_hits[hit][1][2]:
			filtered_hits.append(hit)
			if read_hits[hit][1][0]:
				pair1maps -= 1
			elif read_hits[hit][1][1]:
				pair2maps -= 1

		if read_hits[hit][9] != '*':  # first hit for a read / paired end
			readquals += read_hits[hit][10]
			hitlen += len(read_hits[hit][9])
	read_hits = [read_hits[i] for i in range(len(read_hits))
					if i not in filtered_hits]
	return read_hits, [pair1maps, pair2maps, hitlen, readquals]


# Given all hits for a read, decide whether to use it and whether multimapped
# Returns multimapped reads, and if not, taxid and hit length of uniq map
def process_read(args, read_hits, pair1, pair2, pair1maps, pair2maps):
	read_hits, properties = clean_read_hits(args,read_hits,pair1maps,pair2maps)
	pair1maps, pair2maps, hitlen, readquals = properties
	if len(read_hits) == 0:  # all lines filtered
		return [], 'Ambiguous', '', -1
	if pair1 or pair2:  # paired read
		if pair1maps + pair2maps == 1:
			# if uniq mapped to one end and unmapped to other, count mapped end
			return [], read_hits[0][2], readquals, hitlen

		intersect, intersect_hits = intersect_read_hits(
			read_hits, pair1maps, pair2maps)
		if len(intersect) == 0:  # one end unmapped, other multimapped
			return [], 'Ambiguous', '', -1  #we consider this case too ambiguous
		elif len(intersect) == 1:  # read pairs only agree on one taxid
			return [], read_hits[0][2], readquals, hitlen  # considered uniq map
		else:  # both ends multimapped, handle later
			return intersect_hits, '', readquals, hitlen

	else:  # single end
		if pair1maps > 1:  # multimapped
			return read_hits, '', readquals, hitlen  # multimapped
		else:
			#return False, read_hits[0][2], readquals, hitlen
			return [], read_hits[0][2], readquals, hitlen


# Remove hits to taxids not in taxids2abs (no unique mappings to that taxid)
def preprocess_multimapped(args, multimapped, taxids2abs):
	for i in range(len(multimapped)):
		hitlen = multimapped[i][-1]
		multimapped[i] = [hit for hit in multimapped[i][:-1]
			if hit in taxids2abs]
		if len(multimapped[i]) > 0:
			multimapped[i].append(hitlen)  # ensure hitlen is kept
	multimapped = [read for read in multimapped if len(read) > 0]
	return multimapped


# Processes minimap2 output and returns a dict of taxids mapped to number of
#  	uniquely mapped reads & bases for that taxid, & a list of multimapped reads.
def map_and_process(args, instream, acc2info, taxid2info):
	taxids2abs, multimapped = {}, []  # taxids to abundances, multimapped reads
	low_mem_mmap = {}  # tracks multimapped read hits per taxon in low_mem mode
	taxids2abs['Unmapped'] = ([0.0, 0.0] + taxid2info['Unmapped'])  #placeholder
	prev_read, read_hits = '', [] # read tracker, all hits for read (full lines)
	pair1maps, pair2maps = 0, 0  # reads mapped to each pair (single = pair1)
	tot_rds = 0  # total number of reads in the file

	for line in instream:
		if not args.input_type == 'sam':
			line = line.decode('utf-8')
			if not line:  # process finished
				break
		if line.startswith('@'):
			continue
		splits = line.strip().split()
		if len(splits) < 6:
			continue
		pair1, pair2, chimer, is_bad = parse_flag(int(splits[1]), splits[5])
		if is_bad:  # unmapped or cigar string unavailable
			continue
		splits[1] = [pair1, pair2, chimer, is_bad]  # store flag fields
		# here we change accession to taxid since we want to
		#  	profile by organisms, not accessions
		splits[2] = acc2info[splits[2]][1]

		read, ref = splits[0], splits[2]
		if read != prev_read:
			tot_rds += 1
			if tot_rds % 100000 == 0:
				echo('Processed ' + str(tot_rds) + ' read hits.', args.verbose)
			# get uniq hit taxid and hitlen, or multimapped hits intersect
			intersect_hits, taxid, readquals, hitlen = process_read(
				args, read_hits, pair1, pair2, pair1maps, pair2maps)
			# reset these read-specific variables
			prev_read, read_hits, pair1maps, pair2maps = read, [], 0, 0
			if taxid == 'Ambiguous':  # ambiguous mapping (see process_read)
				if not args.no_quantify_unmapped:
					taxids2abs['Unmapped'][0] += 1.0  # reads hit
				continue
			if taxid != '' and args.length_normalize:
				hitlen /= taxid2info[taxid][0]  # normalize by genome length
			if intersect_hits == []:  # unique hit
				if taxid in taxids2abs:
					taxids2abs[taxid][0] += 1  # reads hit
					taxids2abs[taxid][1] += hitlen  # bases hit
				else:  # also store lineage for taxid
					taxids2abs[taxid] = [1, hitlen] + taxid2info[taxid]
			else:  # multimapped hit
				# store taxids hit and length (in bases) of hits
				#intersect_hits = [[hit[2], len(hit[9])]
				#	for hit in intersect_hits]
				if not args.low_mem:  # store set of taxids per multimapped read
					intersect_hits = [hit[2] for hit in intersect_hits]
					intersect_hits.append(hitlen)  # total hit length
					multimapped.append(intersect_hits)
				else:  # low_mem just stores overall num. of hit bases per taxid
					for hit in intersect_hits:
						taxid = hit[2]
						if taxid in low_mem_mmap:
							low_mem_mmap[taxid] += len(hitlen)
						else:
							low_mem_mmap[taxid] = len(hitlen)
		#else:
		pair1maps += pair1 or not(pair1 or pair2)  # pair1 or single
		pair2maps += pair2  # unchanged if pair2 false
		read_hits.append(splits)
	if not args.no_quantify_unmapped:
		if tot_rds == 0:
			sys.exit('No reads mapped. Aborting...')
		taxids2abs['Unmapped'][1] = taxids2abs['Unmapped'][0] / float(tot_rds)
	return taxids2abs, multimapped, low_mem_mmap


# Assign multimapped reads to specific organisms via proportional method,
#  	e.g. proportional to uniquely mapped reads; method used by MiCoP1
def resolve_multi_prop(args, taxids2abs, multimapped, low_mem_mmap, taxid2info):
	echo('Assigning multimapped reads...', args.verbose)
	# handle low_mem case: simply add hits weighted by unique abundance
	if args.low_mem:
		sum_abs = float(sum([taxids2abs[tax][1] for tax in taxids2abs]))
		for taxid in low_mem_mmap:
			if taxid not in taxids2abs:  # no uniquely-mapped reads
				continue
			proportion = taxids2abs[taxid][1] / sum_abs
			weighted_hits = low_mem_mmap[taxid] * proportion
			if args.length_normalize:  # normalize by genome length if requested
				weighted_hits /= taxid2info[taxid][0]
			taxids2abs[taxid][1] += weighted_hits
		return taxids2abs

	# ensures early read assignments don't affect proportions of later ones
	to_add = {}
	for read_hits in multimapped:
		# get abundances of all taxids in read_hits
		all_taxids = list(set([hit for hit in read_hits[:-1]
			if hit in taxids2abs]))
		if len(all_taxids) == 0:  # all hits were to taxids with no unique hits
			continue
		taxid_abs = [taxids2abs[tax][1] for tax in all_taxids]

		## now get cumulative proportional abs. of taxids relative to each other
		sumabs = sum(taxid_abs)
		if sumabs == 0.0:
			continue
		proportions = [ab / sumabs for ab in taxid_abs]
		hitlen = read_hits[-1]

		# divide hit length proportionally among hit taxids; divided assignment
		for i in range(len(all_taxids)):
			this_hitlen = proportions[i] * hitlen
			if args.length_normalize:
				this_hitlen /= taxid2info[all_taxids[i]][0]
			if all_taxids[i] in to_add:
				to_add[all_taxids[i]] += this_hitlen
			else:
				to_add[all_taxids[i]] = this_hitlen
	for taxid in to_add:  # add in the multimapped portions
		taxids2abs[taxid][1] += to_add[taxid]
	return taxids2abs


# Renormalize each taxonomic rank so each rank sums to 100% abundance
def rank_renormalize(args, clades2abs, only_strains=False):
	rank_totals = {}
	for i in RANKS:
		rank_totals[i] = 0.0
	mapped_pct = 100.0
	if not args.no_quantify_unmapped:  # normalize using pct of mapped reads
		if 'Unmapped' in clades2abs:
			mapped_pct = 100.0 - (100.0 * clades2abs['Unmapped'][-1])
	for clade in clades2abs:
		if clade == 'Unmapped':
			continue
		rank, ab = clades2abs[clade][1], clades2abs[clade][-1]
		if only_strains and rank != 'strain':
			continue
		rank_totals[rank] += ab  # add this to the rank sum total

	for clade in clades2abs:  # here's the normalization
		if clade == 'Unmapped':
			continue
		rank = clades2abs[clade][1]
		if only_strains and rank != 'strain':
			continue
		clades2abs[clade][-1]/= (rank_totals[clades2abs[clade][1]] / mapped_pct)
	return clades2abs


# given initital taxids2abs, some taxids will be higher than strain level;
#  	we insert "unknown" taxa down to strain level for normalization purposes
def gen_lower_taxa(taxids2abs):
	to_add = {}  # lower taxa to add after iteration
	for taxid in taxids2abs:
		taxid, rank, taxlin, namelin, ab = taxids2abs[taxid]
		if rank == 'strain':  # already lowest level
			continue
		# make a new taxid and name for this strain indicating that it is a
		#  	placeholder for a higher taxon that we cannot further narrow down
		rankpos = RANKS.index(rank)
		lowest_name = namelin.split('|')[rankpos]
		new_name = lowest_name + ' unknown strain'
		new_taxid = taxid+'.0'
		to_add[new_taxid] = [new_taxid, 'strain', taxlin + new_taxid,
			namelin + new_name, ab]

	# add in the new taxa
	for taxa in to_add:
		taxids2abs[taxa] = to_add[taxa]
	# clean out higher taxa listings -- they will return when we fill the tree
	taxids2abs = {k:v for k,v in taxids2abs.items() if v[1] == 'strain'}
	return taxids2abs


# Get abundances for all clades in the tree and put in CAMI format
def tree_results_cami(args, taxids2abs):
	# rearrange fields to be in CAMI format
	for taxid in taxids2abs:
		old = taxids2abs[taxid]
		taxids2abs[taxid] = [taxid, old[3], old[5], old[4], old[1]]
	taxids2abs = gen_lower_taxa(taxids2abs)  # ensures everything strain level
	# always renormalize strains, to ensure legitimate profile
	taxids2abs = rank_renormalize(args, taxids2abs, only_strains=True)

	# Now compute higher clade abundances
	clades2abs = {k:v for k,v in taxids2abs.items()}
	for taxid in taxids2abs:
		taxlin = taxids2abs[taxid][2].split('|')
		namelin = taxids2abs[taxid][3].split('|')
		for i in range(len(taxlin)-1):
			clade = taxlin[i]
			if clade == '':  # unspecified at this level
				continue
			if clade in clades2abs:  # already have clade entry, add abundance
				clades2abs[clade][-1] += taxids2abs[taxid][-1]
			else:
				# determine CAMI fields for this clade
				clade_taxid, clade_rank = clade, RANKS[i]
				clade_taxlin = '|'.join(taxlin[:i+1])
				clade_namelin = '|'.join(namelin[:i+1])
				clade_ab = taxids2abs[taxid][-1]  # currently just lower tax ab
				clades2abs[clade_taxid] = [clade_taxid, clade_rank,
					clade_taxlin, clade_namelin, clade_ab]

	if args.rank_renormalize:
		clades2abs = rank_renormalize(args, clades2abs)
	return clades2abs


# Processes uniquely-mapped reads, then estimates abundances using
#  	uniquely-mapped abundances and multimapped read information
def compute_abundances(args, infile, acc2info, tax2info):
	# taxids and higher clades to abundances, and multimapped reads dict
	taxids2abs, clades2abs, multimapped, low_mem_mmap = {}, [], {}, {}
	# run mapping and process to get uniq map abundances & multimapped reads
	#taxids2abs, multimapped = map_and_process(args, infile, acc2info, tax2info)

	if args.input_type == 'sam': # input stream from sam file
		instream = open(infile, 'r')
	else:  # run minimap2 and stream its output as input
		mapper = subprocess.Popen(['../MetaFast/ReadMapping/rm', '-ax', 'sr', '-t', '1', '-2', '-n', '3', '-r', str(args.edit_dist_threshold), '--filter='+str(args.filter), '--secondary=yes', args.db, infile], stdout=subprocess.PIPE, bufsize=1)
		instream = iter(mapper.stdout.readline, "")
	taxids2abs, multimapped, low_mem_mmap = map_and_process(args,
		instream, acc2info, tax2info)
	if args.input_type == 'sam':
		instream.close()
	else:
		mapper.stdout.close()
		mapper.wait()
	if len(multimapped) > 0:
		multimapped = preprocess_multimapped(args, multimapped, taxids2abs)

	# filter out organisms below the read cutoff set by the user, if applicable
	taxids2abs = {k:v for k,v in taxids2abs.items() if v[0] > args.read_cutoff}
	if len(multimapped) > 0 or len(low_mem_mmap) > 0:
		taxids2abs = resolve_multi_prop(args, taxids2abs,
			multimapped, low_mem_mmap, tax2info)
	results = tree_results_cami(args, taxids2abs)
	return results


# Combines and averages results across all input files,
#  	and packs information into easy-to-write form
def gather_results(args, acc2info, taxid2info):
	results = {}
	for infile in args.infiles:
		echo('Computing abundances for input file: ' + infile, args.verbose)
		file_res = compute_abundances(args, infile, acc2info, taxid2info)
		for clade in file_res:
			if clade not in results:
				results[clade] = file_res[clade]
			else:
				results[clade][-1] += file_res[clade][-1]

	if 'Unmapped' in results:
		del results['Unmapped']
	echo('Compiling and writing results...', args.verbose)
	rank_results = {}
	for i in range(len(RANKS)):
		rank_results[i] = []
	for clade in results:
		results[clade][4] = results[clade][4] / len(args.infiles)  #avg of files
		rank = RANKS.index(results[clade][1])
		if rank == 7:  # strain; add extra CAMI genomeID and OTU fields
			taxid = results[clade][0]
			cami_genid, cami_otu = taxid, taxid.split('.')[0]
			results[clade].extend([cami_genid, cami_otu])
		rank_results[rank].append(results[clade])  # results per rank
	return rank_results


# Writes results out in CAMI format
def write_results(args, rank_results):
	with(open(args.output, 'w')) as outfile:
		# Print some CAMI format header lines
		if args.sampleID == 'NONE':
			outfile.write('@SampleID:' + ','.join(args.infiles) + '\n')
		else:
			outfile.write('@SampleID:' + args.sampleID + '\n')
		outfile.write('@Version:Metalign\n')
		outfile.write('@Ranks: ' +
			'superkingdom|phylum|class|order|family|genus|species|strain\n\n')
		outfile.write('@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\t' +
			'PERCENTAGE\t_CAMI_genomeID\t_CAMI_OTU\n')

		for i in range(len(RANKS)):
			lines = rank_results[i]  # all lines to write for this tax level
			# now sort clades in rank by descending abundance
			lines.sort(key=lambda x: 100.0-x[4])
			if lines == None:
				continue
			for line in lines:
				if line[4] < args.min_abundance:
					continue
				if line[4] < 0.00001:
					line[4] = 0.00001
				else:
					line[4] = float('%.5f' % line[4])
				line = [str(i) for i in line]
				outfile.write('\t'.join(line)+'\n')


def map_main(args = None):
	if args == None:
		args = profile_parseargs()
	if args.pct_id > 1.0 or args.pct_id < 0.0:
		sys.exit('Error: --pct_id must be between 0.0 and 1.0, inclusive.')
	if args.db == 'NONE' and not args.infiles[0].endswith('sam'):
		sys.exit('Error: --db must be specified unless sam files are provided.')
	if not args.data.endswith('/'):
		args.data += '/'
	if args.dbinfo == 'AUTO':
		args.dbinfo = args.data + 'db_info.txt'
	if args.input_type == 'AUTO':
		splits = args.infiles[0].split('.')
		if splits[-1] == 'gz':  # gz doesn't help determine file type
			splits = splits[:-1]
		if splits[-1] in ['fq', 'fastq']:
			args.input_type = 'fastq'
		elif splits[-1] in ['fa', 'fna', 'fasta']:
			args.input_type = 'fasta'
		elif splits[-1] == 'sam':
			args.input_type = 'sam'
		else:
			sys.exit('Could not auto-determine file type. Use --input_type.')
	open(args.output, 'w').close()  # test to see if writeable

	# maps NCBI accession to length, taxid, name lineage, taxid lineage
	acc2info, taxid2info = get_acc2info(args)
	# gathers results for all infiles, combines, and organizes into tax levels
	rank_results = gather_results(args, acc2info, taxid2info)
	write_results(args, rank_results)


if __name__ == '__main__':
	args = profile_parseargs()
	map_main(args)
#
