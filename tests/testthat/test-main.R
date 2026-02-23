# tests/testthat/test-main.R

test_that("Demo function runs and generates a valid dataframe", {
  # 1. Execute the demo function
  runDemoPlotMatch()

  # 2. Locate the expected output file
  out_path <- file.path(tempdir(), "OUT_Ground_Positions_trafo.txt")

  # 3. Assert the file actually exists
  expect_true(file.exists(out_path))

  # 4. Read the generated file to verify it is a valid data frame
  result_df <- read.table(out_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # 5. Assert the data frame is not empty
  expect_true(nrow(result_df) > 0)
  expect_true(ncol(result_df) > 0)

  # 6. Print the dataframe (head) to the test console
  cat("\n\n--- First few rows of the transformed plot dataframe ---\n")
  print(head(result_df))
  cat("--------------------------------------------------------\n\n")
})
