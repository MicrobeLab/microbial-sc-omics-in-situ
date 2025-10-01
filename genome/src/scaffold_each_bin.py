'''
usage: python scaffold_each_bin.py [input_file] [output_file_prefix] [b/g]
'''

import sys

in_file = sys.argv[1]
out_prefix = sys.argv[2]
fmt = sys.argv[3]
assert fmt == 'b' or fmt == 'g'

# ---- hard code ----- 
min_scaffold_per_bin = 2
min_length_per_scaffold = 500
min_geno_size = 10000
# --------------------


def parse_input_line(line, fmt = 'b'):
	line = line.strip()
	if fmt == 'b':
		lls = line.split('\t')
	elif fmt == 'g':
		lls = line.split(',')
	else:
		lls = None
	scaff = lls[0]
	bin_info = lls[1]
	scaff_len = int(scaff.split('_')[3])
	return (scaff, bin_info, scaff_len)



# creare three dicts: 
# 1. bin -> genome_zie
# 2. bin -> #scaffolds
# 3. bin -> scaffold_ids


bin2gsz = {}
bin2num = {}
bin2ids = {}

with open(in_file) as fh:
	for l in fh:
		node_id, bin_name, node_len = parse_input_line(l, fmt)
		if node_len < min_length_per_scaffold:
			continue
		if not bin_name in bin2gsz:
			bin2gsz[bin_name] = node_len
			bin2num[bin_name] = 1
			bin2ids[bin_name] = [node_id]
		else:
			bin2gsz[bin_name] += node_len
			bin2num[bin_name] += 1
			bin2ids[bin_name].append(node_id)


for each_bin in bin2gsz:
	if bin2gsz[each_bin] < min_geno_size or bin2num[each_bin] < min_scaffold_per_bin:
		continue
	else:
		lst_nodes = bin2ids[each_bin]
		if ' ' in each_bin:
			each_bin = each_bin.replace(' ', '-')
		out_filename = out_prefix + '_' + each_bin + '.txt'
		with open(out_filename, 'w') as fh_out:
			for node in lst_nodes:
				fh_out.write(node + '\n')

				