library(testthat)
library(logger)
log_threshold(FATAL)

test_file('tests/test-assign_cell_types.R', reporter = 'progress') # 'summary')
