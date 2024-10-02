import os
from pathlib import Path


def file_path():
    """
    Path to this file
    """
    return Path(__file__)


def project_dir():
   """
   mcp3d project directory string
   """
   return str(file_path().parents[2])


def src_test_data_dir():
    """
    test data directory of c++ libraries
    """
    return os.path.join(project_dir(), 'test_data', 'src')


def python_test_data_dir():
    """
    test data directory of python code
    """
    return os.path.join(project_dir(), 'test_data', 'python')


def image_reader_test_data_dir():
    """
    test data directory for neural_networks.data_io.image_reader
    """
    return os.path.join(python_test_data_dir(), 'image_reader')
