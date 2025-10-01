import sys

r1=sys.argv[1]
r2=sys.argv[2]
out_prefix = sys.argv[3]

with open(out_prefix + ".yaml", "w") as f:
    f.write('['+'\n')
    f.write('{'+'\n')
    f.write('orientation: "fr",'+'\n')
    f.write('type: "paired-end",'+'\n')
    f.write(f'right reads: ["{r1}"], left reads: ["{r2}"]' + '\n')
    f.write('}'+'\n')
    f.write(']'+'\n')
