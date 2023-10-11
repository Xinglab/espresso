import unittest

import tests.alignments.test as alignments_test
import tests.cd44_example.test as cd44_example_test
import tests.high_confidence_sjs.test as high_confidence_sjs_test
import tests.isoform_assignment.test as isoform_assignment_test
import tests.read_filters.test as read_filters_test
import tests.sirv_example.test as sirv_example_test
import tests.sorted_input.test as sorted_input_test


def build_test_suite():
    loader = unittest.defaultTestLoader
    suite = unittest.TestSuite()
    suite.addTest(loader.loadTestsFromModule(alignments_test))
    suite.addTest(loader.loadTestsFromModule(cd44_example_test))
    suite.addTest(loader.loadTestsFromModule(high_confidence_sjs_test))
    suite.addTest(loader.loadTestsFromModule(isoform_assignment_test))
    suite.addTest(loader.loadTestsFromModule(read_filters_test))
    suite.addTest(loader.loadTestsFromModule(sirv_example_test))
    suite.addTest(loader.loadTestsFromModule(sorted_input_test))
    return suite


def main():
    runner = unittest.TextTestRunner(verbosity=2)
    suite = build_test_suite()
    runner.run(suite)


if __name__ == '__main__':
    main()
