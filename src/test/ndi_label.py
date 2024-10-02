from __future__ import print_function
import os
import sys
import numpy as np
import scipy.ndimage as ndi


def gen_random_data(i_, l0=1024, l1=None, l2=None):
    outdir = os.path.join(os.path.dirname(os.path.relpath(__file__)),
                          'connected_component')
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    if l1 is None:
        l1 = l0
    if l2 is None:
        l2 = l0
    volume = np.random.randint(0, high=2, size=(l0, l1, l2), dtype=np.uint8)
    np.save(os.path.join(outdir, 'random{}.npy'.format(i_)), volume)
    label, n = ndi.label(volume)
    np.save(os.path.join(outdir, 'ndi_label{}.npy'.format(i_)), label)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("usage: python ndi_label n l0 [l1] [l2]\n")
    n = int(sys.argv[1])
    for i in range(0, n):
        if len(sys.argv) == 3:
            gen_random_data(i, int(sys.argv[2]))
        if len(sys.argv) == 4:
            gen_random_data(i, int(sys.argv[2]), int(sys.argv[3]))
        if len(sys.argv) == 5:
            gen_random_data(i, int(sys.argv[2]),
                            int(sys.argv[3]), int(sys.argv[4]))