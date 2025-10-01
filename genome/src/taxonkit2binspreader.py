'''
usage: python taxonkit2binspreader.py [input_file] [output_file] [b/g]
'''

import sys

min_contig_per_bin = 1
fmt = sys.argv[3]
assert fmt == 'b' or fmt == 'g'


def clade2family(lineage):
	pos_f = lineage.find("f__")
	pos_semicol = lineage.find(";", pos_f)
	assigned = True
	if lineage[pos_f:pos_semicol] == "f__":
		assigned = False
		f_name = None
	else:
		f_name = lineage[(pos_f+3):pos_semicol]
	return (assigned, f_name)



dict_family2count = {}

fhi = open(sys.argv[1])

for l in fhi:
	lls = l.strip().split('\t')
	node = lls[0]
	lng = lls[1]
	ass, nam = clade2family(lng)
	if not ass:
		continue
	if nam in dict_family2count:
		dict_family2count[nam] += 1
	else:
		dict_family2count[nam] = 1

fhi.close()


good_fam = []

for f in dict_family2count:
	cou = dict_family2count[f]
	if cou >= min_contig_per_bin:
		good_fam.append(f)

dict_fam2idx = {}
idx = 0
for f in good_fam:
	idx += 1
	dict_fam2idx[f] = idx


fhin = open(sys.argv[1])
fho = open(sys.argv[2], 'w')

for l in fhin:
	lls = l.strip().split('\t')
	node = lls[0]
	lng = lls[1]
	ass, nam = clade2family(lng)
	if not ass:
		continue
	if not nam in good_fam:
		continue
	if fmt == 'b':
		outl = node + '\t' + 'bin_' + nam + '\n'
	if fmt == 'g':
		outl = node + ',' + str(dict_fam2idx[nam]) + '\n'
	fho.write(outl)

fhin.close()
fho.close()
