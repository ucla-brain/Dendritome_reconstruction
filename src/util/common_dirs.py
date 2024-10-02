import os


def mcp_dir():
    d = '/ifs/loni/faculty/dong/mcp'
    assert(os.path.exists(d))
    return d


def mcp_util_dir():
    d = os.path.join(mcp_dir(), 'utils')
    assert(os.path.exists(d))
    return d


def curr_dir():
    return os.path.abspath(os.path.dirname(__file__))


def proj_dir():
    d = os.path.normpath(os.path.join(curr_dir(), '..'))
    assert(os.path.exists(d))
    return d


def proj_parent_dir():
    d = os.path.normpath(os.path.join(proj_dir(), '..'))
    assert (os.path.exists(d))
    return d

def proj_grandparent_dir():
    d = os.path.normpath(os.path.join(proj_parent_dir(), '..'))
    assert (os.path.exists(d))
    return d


def third_party_dir():
    d = os.path.normpath(os.path.join(proj_dir(), '3rd_party'))
    if not os.path.exists(d):
        os.makedirs(d)
    assert (os.path.exists(d))
    return d


def full_external_repos_dir():
    d = os.path.normpath(os.path.join(proj_grandparent_dir(),
                                      'mcp3d_external_repos'))
    if not os.path.exists(d):
        os.makedirs(d)
    assert (os.path.exists(d))
    return d


def external_packages_dir():
    d = os.path.normpath(os.path.join(proj_grandparent_dir(),
                                      'mcp3d_external_packages'))
    assert (os.path.exists(d))
    return d
