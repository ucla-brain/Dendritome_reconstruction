from __future__ import print_function
import os
import socket
import shutil
import argparse
import tarfile
import subprocess
import multiprocessing
import common_dirs


def hostname():
    return socket.gethostname()


def cmake_bin_path():
    if hostname() == 'c2001':
        cmake_path = os.path.join(common_dirs.mcp_util_dir(),
                                  'cmake3.8/cmake-3.8.2-Linux-x86_64/bin',
                                  'cmake')
        assert(os.path.exists(cmake_path))
        return cmake_path
    else:
        return '/usr/local/cmake-3.8.2/bin/cmake'


def cluster_cc_path():
    cc = '/usr/local/gcc-7.2.0/bin/gcc'
    assert(os.path.exists(cc))
    return cc


def cluster_cxx_path():
    cxx = '/usr/local/gcc-7.2.0/bin/g++'
    assert(os.path.exists(cxx))
    return cxx


def cluster_mpicc_path():
    mpicc = '/usr/local/mpich-3.1.3/bin/mpicc'
    assert(os.path.exists(mpicc))
    return mpicc


def cluster_mpicxx_path():
    mpicxx = '/usr/local/mpich-3.1.3/bin/mpicxx'
    assert(os.path.exists(mpicxx))
    return mpicxx


def cluster_mpiexec_path():
    mpiexec = '/usr/local/mpich-3.1.3/bin/mpiexec'
    assert(os.path.exists(mpiexec))
    return mpiexec


def cluster_cmake_linker_flags():
    # cluster libc version is 2.12, too low for some libraries with versioned
    # symbols. glibc-2.23 is compiled. linker flags should define
    # specify rpath (run time), rpath-link (link time) and
    # dynamic-linker (path to the dynamic linker ld-linux-*.so)
    glibc_dir = '/ifs/loni/faculty/dong/mcp/utils/glibc-2.23'
    glibc_lib_dir = os.path.join(glibc_dir, 'lib')
    linker_path = os.path.join(glibc_lib_dir, 'ld-linux-x86-64.so.2')
    linker_flags = '-Wl,-rpath,{0},-rpath-link,{0},-dynamic-linker,{1}' \
        .format(glibc_lib_dir, linker_path)
    return ['-DCMAKE_SHARED_LINKER_FLAGS:STRING={}'.format(linker_flags),
            '-DCMAKE_EXE_LINKER_FLAGS:STRING={}'.format(linker_flags)]

def clone_or_update_git_repository(repo_dir, repository_url):
    if os.path.exists(repo_dir):
        print('git', 'pull in ', repo_dir)
        subprocess.call(['git', 'pull'], cwd=repo_dir, shell=False)
    else:
        print('git', 'clone', repository_url, repo_dir)
        subprocess.check_call(['git', 'clone', repository_url,
                               repo_dir], shell=False, )


def export_git_repository(repo_dir, target_dir, branch='', tag=''):
    if not branch:
        branch = 'master'
    shutil.rmtree(target_dir, ignore_errors=True)
    print('cloning repo {0} branch {1}'
          .format(os.path.basename(repo_dir), branch))
    subprocess.check_call(['git', 'clone', '--shared', '--branch',
                           branch, repo_dir, target_dir], shell=False)
    if tag:
        print('checking out tag', tag)
        subprocess.check_call(['git', 'checkout', tag],
                              cwd=target_dir, shell=False)

    print(os.path.join(target_dir, '.git'))
    shutil.rmtree(os.path.join(target_dir, '.git'), ignore_errors=False)


def unpack_package(package_path, unpack_dir):
    print('unpacking', package_path, 'to', unpack_dir, '...')
    package_name = os.path.basename(package_path)
    if package_name.lower().endswith('.tar.gz') or \
       package_name.lower().endswith('.tar.bz2') or \
       package_name.lower().endswith('.tar.xz') or \
       package_name.lower().endswith('.tgz'):
        with tarfile.open(package_path, mode='r|*') as tf:
            tf.extractall(path=unpack_dir)


def make_build_dir(src_dir):
    build_dir = os.path.normpath(os.path.join(src_dir,
                                              '__' + os.path.basename(src_dir) +
                                              'build'))
    shutil.rmtree(build_dir, ignore_errors=True)
    if not os.path.exists(build_dir):
        os.makedirs(build_dir)
    return build_dir


def cmake_common_cmd(install_dir, build_type='Release'):
    cmake = cmake_bin_path()
    if not os.path.exists(cmake):
        cmake = 'cmake'
    cmd = [cmake]
    # currently assuming caller has correct LD_LIBRARY_PATH settings for
    # the cluster gcc. should update cmake to take care of library paths
    if hostname() == 'c2001':
        cmd.extend(['-DCMAKE_C_COMPILER:STRING=' + cluster_cc_path(),
                    '-DCMAKE_CXX_COMPILER:STRING=' + cluster_cxx_path(),
                    '-DMPI_C_COMPILER:STRING=' + cluster_mpicc_path(),
                    '-DMPI_CXX_COMPILER:STRING=' + cluster_mpicxx_path(),
                    '-DMPIEXEC:STRING=' + cluster_mpiexec_path()])
        cmd.extend(cluster_cmake_linker_flags())
    cmd.extend(['-DCMAKE_BUILD_TYPE:STRING=' + build_type,
                '-DCMAKE_INSTALL_PREFIX:STRING=' + install_dir,
                '-DCMAKE_CXX_FLAGS:STRING=-std=c++11'])
    return cmd


def build_and_install(cmake_cmd, build_dir, env=None):
    print('build dir:', build_dir)
    print(cmake_cmd)
    subprocess.check_call(cmake_cmd, cwd=build_dir, shell=False, env=env)
    subprocess.check_call(['make', '-j{}'.format(multiprocessing.cpu_count()),
                           'install'], cwd=build_dir, shell=False, env=env)


def build_gtest(build_type='Debug'):
    repo_url = 'git@github.com:google/googletest.git'
    branch = 'master'
    tag = 'release-1.8.0'
    repo_dir = os.path.join(common_dirs.full_external_repos_dir(), 'gtest')
    install_dir = os.path.join(common_dirs.third_party_dir(), 'gtest')
    shutil.rmtree(install_dir, ignore_errors=True)

    clone_or_update_git_repository(repo_dir, repo_url)
    export_git_repository(repo_dir, install_dir, branch=branch, tag=tag)
    build_dir = make_build_dir(install_dir)
    try:
        cmake_cmd = cmake_common_cmd(install_dir, build_type=build_type)
        if build_type == 'Debug':
            cmake_cmd.extend(['-Dgtest_build_samples=ON'])
        cmake_cmd.extend(['-DGTEST_HAS_PTHREAD=1', install_dir])
        build_and_install(cmake_cmd, build_dir)
    finally:
        shutil.rmtree(build_dir, ignore_errors=False)


def build_nlohmann_json(build_type='Release'):
    repo_url = 'git@github.com:nlohmann/json.git'
    branch = 'develop'
    tag = 'v3.5.0'
    repo_dir = os.path.join(common_dirs.full_external_repos_dir(), 'nlohmann')
    install_dir = os.path.join(common_dirs.third_party_dir(), 'nlohmann')
    shutil.rmtree(install_dir, ignore_errors=True)

    clone_or_update_git_repository(repo_dir, repo_url)
    export_git_repository(repo_dir, install_dir, branch=branch, tag=tag)
    build_dir = make_build_dir(install_dir)
    try:
        cmake_cmd = cmake_common_cmd(install_dir, build_type=build_type)
        if build_type == 'Release':
            cmake_cmd.extend(['-DBUILD_TESTING:BOOL=OFF'])
        cmake_cmd.extend([install_dir])
        build_and_install(cmake_cmd, build_dir)
    finally:
        shutil.rmtree(build_dir, ignore_errors=False)


def build_hdf5(mode='serial', build_type='Release'):
    external_package_dir = common_dirs.external_packages_dir()
    tarfile_name = 'hdf5-1.10.1.tar.bz2'
    package_name = 'hdf5-1.10.1'
    third_party_dir = common_dirs.third_party_dir()
    unpack_package(os.path.join(external_package_dir, tarfile_name),
                   third_party_dir)

    install_dir = os.path.join(common_dirs.third_party_dir(), 'hdf5')
    shutil.rmtree(install_dir, ignore_errors=True)
    shutil.move(os.path.join(third_party_dir, package_name), install_dir)
    build_dir = make_build_dir(install_dir)

    try:
        cmake_cmd = cmake_common_cmd(install_dir, build_type=build_type)
        build_example = 'ON' if build_type == 'Debug' else 'OFF'
        enable_debug = 'ON' if build_type == 'Debug' else 'OFF'
        # parallel hdf can not be built with c++ binding
        if mode == 'parallel':
            cmake_cmd.extend(['-DBUILD_SHARED_LIBS:BOOL=ON',
                              '-DHDF5_BUILD_CPP_LIB:BOOL=OFF',
                              '-DHDF5_BUILD_HL_LIB:BOOL=ON',
                              '-DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON',
                              '-DHDF5_ENABLE_PARALLEL:BOOL=ON',
                              '-DMPIEXEC_MAX_NUMPROCS:STRING=1024',
                              '-DHDF5_ENABLE_THREADSAFE:BOOL=OFF',
                              '-DHDF5_BUILD_EXAMPLES:BOOL={}'
                              .format(build_example),
                              '-DHDF5_ENABLE_DEBUG:BOOL={}'
                              .format(enable_debug),
                              '-DHDF5_BUILD_TOOLS:BOOL=ON',
                              install_dir])
        elif mode == 'serial':
            cmake_cmd.extend(['-DBUILD_SHARED_LIBS:BOOL=ON',
                              '-DHDF5_BUILD_CPP_LIB:BOOL=ON',
                              '-DHDF5_BUILD_HL_LIB:BOOL=ON',
                              '-DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON',
                              '-DHDF5_ENABLE_PARALLEL:BOOL=OFF',
                              '-DHDF5_ENABLE_THREADSAFE:BOOL=OFF',
                              '-DHDF5_BUILD_EXAMPLES:BOOL={}'
                              .format(build_example),
                              '-DHDF5_ENABLE_DEBUG:BOOL={}'
                              .format(enable_debug),
                              '-DHDF5_BUILD_TOOLS:BOOL=ON',
                              install_dir])
        build_and_install(cmake_cmd, build_dir)
    finally:
        shutil.rmtree(build_dir, ignore_errors=False)


def build_eigen(build_type='Release'):
    repo_url = 'git@github.com:eigenteam/eigen-git-mirror.git'
    branch = 'master'
    tag = '3.3.4'
    repo_dir = os.path.join(common_dirs.full_external_repos_dir(), 'eigen')
    install_dir = os.path.join(common_dirs.third_party_dir(), 'eigen')
    shutil.rmtree(install_dir, ignore_errors=True)

    clone_or_update_git_repository(repo_dir, repo_url)
    export_git_repository(repo_dir, install_dir, branch=branch, tag=tag)


def build_opencv(build_type='Release'):
    repo_url = 'https://github.com/opencv/opencv.git'
    branch = 'master'
    tag = '4.0.0'
    repo_dir = os.path.join(common_dirs.full_external_repos_dir(), 'opencv')
    install_dir = os.path.join(common_dirs.third_party_dir(), 'opencv')
    shutil.rmtree(install_dir, ignore_errors=True)

    clone_or_update_git_repository(repo_dir, repo_url)
    export_git_repository(repo_dir, install_dir, branch=branch, tag=tag)
    build_dir = make_build_dir(install_dir)

    try:
        cmake_cmd = cmake_common_cmd(install_dir, build_type=build_type)
        switch = 'OFF' if build_type == 'Release' else 'ON'
        cmake_cmd.extend(['-DBUILD_ZLIB:BOOL=ON', '-DBUILD_PNG:BOOL=ON',
                          '-DBUILD_JPEG:BOOL=ON', '-DBUILD_TIFF:BOOL=ON',
                          '-DBUILD_JASPER:BOOL=OFF', '-DBUILD_WEBP:BOOL=OFF',
                          '-DBUILD_OPENEXR:BOOL=OFF', '-DWITH_WEBP:BOOL=OFF',
                          '-DWITH_OPENEXR:BOOL=OFF', '-DWITH_JASPER:BOOL=OFF',
                          # if java build is not turned off explicitly, on
                          # the cluster the cmake file will hang trying to
                          # find JNI components
                          '-DWITH_LAPACK:BOOL=OFF', '-DBUILD_JAVA:BOOL=OFF',
                          '-DBUILDD_DOCS:BOOL={}'.format(switch),
                          '-DBUILD_EXAMPLES:BOOL={}'.format(switch)])
        cmake_cmd.extend([install_dir])
        build_and_install(cmake_cmd, build_dir)
    finally:
        shutil.rmtree(build_dir, ignore_errors=False)


def build_boost(build_type='Release'):
    external_package_dir = common_dirs.external_packages_dir()
    tarfile_name = 'boost_1_65_1.tar.bz2'
    package_name = 'boost_1_65_1'
    third_party_dir = common_dirs.third_party_dir()
    unpack_package(os.path.join(external_package_dir, tarfile_name),
                   third_party_dir)
    src_dir = os.path.join(third_party_dir, package_name)
    install_dir = os.path.join(third_party_dir, 'boost')
    shutil.rmtree(install_dir, ignore_errors=True)
    try:
        subprocess.check_call(['bash', 'bootstrap_wrapper.sh'],
                              cwd=common_dirs.curr_dir(), shell=False)
    finally:
        shutil.rmtree(src_dir, ignore_errors=False)


def main():
    build_libs = {
        'gtest': False,
        'nlohmann_json': False,
        'hdf5-serial': False,
        'hdf5-parallel': False,
        'eigen': False,
        'opencv': False,
        'boost': False
    }
    usage = '''python build_ext_libraries.py all | lib1 lib2 ... libn 
                      --build_type Debug | Release'''
    parser = argparse.ArgumentParser(description="build external libraies",
                                     usage=usage)
    parser.add_argument("libs", nargs='+', action="store",
                        help="build selected external libraries")
    parser.add_argument('--build_type', action='store', default='Release')
    args = parser.parse_args()
    if args.libs[0] == 'all':
        for lib in build_libs:
            build_libs[lib] = True
        build_libs['hdf5-parallel'] = False
    else:
        if 'hdf5-serial' in args.libs and 'hdf-parallel' in args.libs:
            raise ValueError(
                'can only build either hdf5-serial or hdf5-parallel')
        for lib in args.libs:
            if lib not in build_libs:
                if lib == 'hdf5':
                    raise ValueError(
                        'unknown library {}. did you mean hdf5-serial '
                        'or hdf5-parallel'.format(lib))
                raise ValueError('unknown library {}'.format(lib))
            build_libs[lib] = True
    build_type = args.build_type
    if build_type.lower() != 'release' and build_type.lower() != 'debug':
        raise ValueError('unknown build type {}'.format(build_type))
    if build_libs['gtest']:
        build_gtest(build_type=build_type)
    if build_libs['nlohmann_json']:
        build_nlohmann_json(build_type=build_type)
    if build_libs['hdf5-serial']:
        build_hdf5(mode='serial', build_type=build_type)
    if build_libs['hdf5-parallel']:
        build_hdf5(mode='parallel', build_type=build_type)
    if build_libs['eigen']:
        build_eigen(build_type=build_type)
    if build_libs['opencv']:
        build_opencv(build_type=build_type)
    if build_libs['boost']:
        build_boost(build_type=build_type)


if __name__ == '__main__':
    main()

    """
    mcp3d
        ---3rd_party
           ---gtest: install_dir
               ---build_dir, created before build and removed after
           ---nlohman_json: install_dir
              ---build_dir, created before build and removed after  
    mcp3d_external_repos
        ---gtest: repo, export to mcp3d/3rd_party/gtest  
        ---nlohman_json: repo, export to mcp3d/3rd_party/nlohman_json 
    mcp3d_external_packages
        ---hdf5-1.10.1.tar.bz2   
    """