import sys
import pyfaidx

if len(sys.argv) != 4:
    print('Find uncovered sequence from alignment')
    print('Usage : paf_gap_seq.py fasta_loc paf_loc output_loc')
    sys.exit(1)
else:
    fasta_loc, paf_loc, output_file = sys.argv[1:]

asm_seq = pyfaidx.Fasta(fasta_loc)

paf_data = []
paf_line_data = []
qry_ctg = ''
with open(paf_loc, 'r') as f:
    for l in f:
        paf_line = l.split('\t')

        if qry_ctg == '':
            qry_ctg = paf_line[0]
            paf_line_data = [qry_ctg, int(paf_line[1]), []]

        if qry_ctg != paf_line[0]:
            qry_ctg = paf_line[0]
            paf_data.append(paf_line_data)
            paf_line_data = [qry_ctg, int(paf_line[1]), []]

        paf_line_data[2].append([int(paf_line[2]), int(paf_line[3]) - 1])

paf_data.append(paf_line_data)

paf_gap_data = []
for ctg_data in paf_data:
    paf_gap_data.append((ctg_data[0], ctg_data[1], []))
    ctg_data[2] = [[-1, -1]] + sorted(ctg_data[2], key=lambda t : (t[0], t[1])) + [[ctg_data[1], ctg_data[1]]]

    bnd = -1
    for st, nd in ctg_data[2]:
        if bnd + 1 <= st - 1:
            paf_gap_data[-1][2].append([bnd + 1, st - 1])
        
        bnd = max(nd, bnd)

baseline = 0
with open(output_file, 'w') as f:
    for ctg_name, ctg_len, gap_list in paf_gap_data:
        for st, nd in gap_list:
            if st == 0 or nd == ctg_len - 1 or nd - st + 1 > baseline:
                gap_seq = asm_seq[ctg_name][st:nd + 1]
                f.write('>' + gap_seq.fancy_name + '\n')
                f.write(str(gap_seq) + '\n')