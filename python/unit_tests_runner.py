from unittest import TestSuite, TestLoader, TextTestRunner
from unit_tests import test_image_reader, test_bounding_box


def main():
    # https://www.internalpointers.com/post/run-painless-test-suites-python-unittest
    loader = TestLoader()
    suite = TestSuite()
    suite.addTests(loader.loadTestsFromModule(test_image_reader))
    suite.addTests(loader.loadTestsFromModule(test_bounding_box))
    runner = TextTestRunner(verbosity=3)
    result = runner.run(suite)


if __name__ == '__main__':
    main()
