from __future__ import print_function
import os
import sys

if __name__ == '__main__':
    if len(sys.argv) == 0:
        print("usage: par_cc_3d_perf.py file_dir")
    perf_file_dir = sys.argv[1]
    n_perf_file = len(os.listdir(perf_file_dir))
    block_dims = ''
    read_time, local_label_time, global_resolve_time, \
        reassign_time, count_time = 0, 0, 0, 0, 0
    for perf_file in os.listdir(perf_file_dir):
        with open(os.path.join(perf_file_dir, perf_file), 'r') as pf:
            for l in pf:
                if l.find('block dimension') >= 0:
                    block_dims = l.strip()\
                        .replace('block dimension (layers, height, width) = ', '')
                if l.find('time to read input: ') >= 0:
                    read_time += float(
                        l.strip()
                         .replace('time to read input: ', '')
                         .replace(' ms', ''))
                if l.find('time to label local image ') >= 0:
                    local_label_time += float(
                        l.strip()
                         .replace('time to label local image ', '')
                         .replace(' ms', ''))
                if l.find('time to communiate and resolve global label ') >= 0:
                    global_resolve_time += float(
                        l.strip()
                         .replace('time to communiate and resolve global label ', '')
                         .replace(' ms', ''))
                if l.find('time to reassign to resolved global label ') >= 0:
                    reassign_time += float(
                        l.strip()
                         .replace('time to reassign to resolved global label ', '')
                         .replace(' ms', ''))
                if l.find('time to count number of components ') >= 0:
                    count_time += float(
                        l.strip()
                         .replace('time to count number of components ', '')
                         .replace(' ms', ''))
    with open(os.path.join(perf_file_dir, 'perf_{}.txt'.format(n_perf_file)), 'w') as pf:
        pf.write('number of processors = {}\n'.format(n_perf_file))
        pf.write('block dimension = {}\n'.format(block_dims))
        pf.write('read input time = {} ms\n'.format(read_time / n_perf_file))
        pf.write('local label time = {} ms\n'.format(local_label_time / n_perf_file))
        pf.write('global resolve time = {} ms\n'.format(global_resolve_time / n_perf_file))
        pf.write('reassign label time = {} ms\n'.format(reassign_time / n_perf_file))
        pf.write('count component time = {} ms\n'.format(count_time / n_perf_file))
