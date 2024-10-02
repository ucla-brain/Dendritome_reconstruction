from __future__ import division, print_function
import os
import sys
import math
import common_dirs

dev_nodes = {"c2001", "c2002", "c2003", "c2004"}


def cmp_host(host1, host2):
    if len(host1) != 3 or len(host2) != 3:
        raise ValueError("each argument should contain hostname, ncpu and load")
    # load is first key
    if host1[2] < host2[2]:
        return 1
    elif host1[2] > host2[2]:
        return -1
    elif host1[1] > host2[1]:
        return 1
    elif host1[1] > host2[1]:
        return -1
    else:
        return 0


def rank_host(qhost_out_path):
    mpi_hosts = []
    with open(qhost_out_path, 'r') as qf:
        for l in qf:
            if l[0:2] != 'c2':
                continue
            hostname, arch, ncpu, nsoc, ncor, nthr, load, \
            memtot, memuse, swapto, swapus = l.split()
            if hostname not in dev_nodes:
                try:
                    ncpu = int(ncpu)
                except ValueError:
                    print("discarding hostname = {} due to error parsing ncpu"
                          .format(hostname))
                try:
                    load = float(load)
                except ValueError:
                    print("discarding hostname = {}, load - ".format(hostname))
                    continue
                mpi_hosts.append((hostname, ncpu, load))
    mpi_hosts.sort(cmp=cmp_host, reverse=True)
    with open(os.path.join(common_dirs.proj_dir(), "host_rank"), 'w') as f:
        for mpi_host in mpi_hosts:
            f.write("{} {} {}\n".format(mpi_host[0], mpi_host[1], mpi_host[2]))
    return mpi_hosts


def multi_thread_hostfile(ranked_hosts, n_proc, n_thread_per_proc):
    n_jobs = 0
    with open(os.path.join(common_dirs.proj_dir(), "hostfile"), 'w') as f:
        for mpi_host in ranked_hosts:
            n_jobs_on_host = (mpi_host[1] - math.ceil(mpi_host[2])) \
                             / n_thread_per_proc
            if n_jobs_on_host > 0:
                f.write("{}.c.ini.usc.edu:{}\n".format(mpi_host[0], int(n_jobs_on_host)))
                n_jobs += n_jobs_on_host
    if n_jobs < n_proc:
        print("insufficient cpu resources\n")
        exit(1)


def round_robin(ranked_hosts):
    total_cluster_cores = sum([host[1] for host in ranked_hosts])
    with open(os.path.join(common_dirs.proj_dir(), "hostfile"), 'w') as f:
        for n in range(0, total_cluster_cores):
            f.write("{}.c.ini.usc.edu:{}\n".
                    format(ranked_hosts[n % len(ranked_hosts)][0], 1))


def main():
    if len(sys.argv) < 3:
        print("usage: python gen_hostfile n_procs n_thread_per_proc\n")
        exit(0)
    n_proc = int(sys.argv[1])
    n_thread_per_proc = int(sys.argv[2])

    qhost_out_path = os.path.join(common_dirs.proj_dir(), 'qhost.out')
    if not os.path.isfile(qhost_out_path):
        raise ValueError(qhost_out_path + " does not exist")
    ranked_hosts = rank_host(qhost_out_path)

    if n_thread_per_proc > 1:
        multi_thread_hostfile(ranked_hosts, n_proc, n_thread_per_proc)

    else:
        round_robin(ranked_hosts)


# give 2 arguments: n_procs, nthreads / proc
if __name__ == "__main__":
    main()
