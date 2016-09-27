import sys
import io
import os.path as op

with io.open(sys.argv[1]) as fd:
    lines = [l.strip() for l in fd.readlines()]

new_lines = []
for i, line in enumerate(lines):
    if i == 0:
        start_pad = 200
        stop_pad = 20
    elif i == len(lines) - 1:
        start_pad = 20
        stop_pad = 200
    else:
        start_pad = 20
        stop_pad = 20
    rec = line.split('\t')
    rec[1] = str(int(rec[1]) - start_pad)
    rec[2] = str(int(rec[2]) - stop_pad)
    new_lines.append('\t'.join(rec) + '\n')

o_fn = op.join(op.dirname(sys.argv[1]), 'pad_' + op.basename(sys.argv[1]))
with io.open(o_fn, 'w') as fd:
    fd.write(''.join(new_lines))
