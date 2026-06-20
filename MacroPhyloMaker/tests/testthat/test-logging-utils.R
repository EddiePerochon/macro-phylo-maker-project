library(testthat)

# Access internal functions safely from package namespace
ns <- asNamespace("MacroPhyloMaker")

box_chars    <- get(".box_chars", envir = ns)
make_logger  <- get(".make_logger", envir = ns)
close_logger <- get(".close_logger", envir = ns)
log_msg      <- get("log_msg", envir = ns)
log_section  <- get("log_section", envir = ns)

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------

read_log <- function(path) {
  if (!file.exists(path)) return("")
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

capture_cat <- function(expr) {
  paste(capture.output(force(expr), type = "output"), collapse = "\n")
}

# -------------------------------------------------------------------
# .box_chars()
# -------------------------------------------------------------------

test_that(".box_chars returns UTF-8 characters by default", {
  old <- options(MacroPhyloMaker.ascii = FALSE, AntPhyloMaker.ascii = FALSE)
  on.exit(options(old), add = TRUE)

  chars <- box_chars()

  expect_equal(chars$h, "─")
  expect_equal(chars$tl, "┌")
  expect_equal(chars$tr, "┐")
  expect_equal(chars$bl, "└")
  expect_equal(chars$br, "┘")
  expect_equal(chars$v, "│")
})

test_that(".box_chars returns ASCII fallback when requested", {
  old <- options(MacroPhyloMaker.ascii = TRUE, AntPhyloMaker.ascii = FALSE)
  on.exit(options(old), add = TRUE)

  chars <- box_chars()

  expect_equal(chars$h, "-")
  expect_equal(chars$tl, "+")
  expect_equal(chars$tr, "+")
  expect_equal(chars$bl, "+")
  expect_equal(chars$br, "+")
  expect_equal(chars$v, "|")
})

test_that(".box_chars honors backward-compatible AntPhyloMaker option", {
  old <- options(MacroPhyloMaker.ascii = NULL, AntPhyloMaker.ascii = TRUE)
  on.exit(options(old), add = TRUE)

  chars <- box_chars()

  expect_equal(chars$h, "-")
  expect_equal(chars$tl, "+")
  expect_equal(chars$v, "|")
})

# -------------------------------------------------------------------
# .make_logger()
# -------------------------------------------------------------------

test_that(".make_logger creates output directory, logfile, and header when enabled", {
  dir <- tempfile("logger_dir_")

  console <- capture_cat({
    logger <- make_logger(
      outdir = dir,
      genus = "TestGenus",
      mode = "extract",
      enable = TRUE,
      console = TRUE
    )
  })
  on.exit(try(close_logger(logger), silent = TRUE), add = TRUE)

  expect_true(dir.exists(dir))
  expect_s3_class(logger, "smart_logger")
  expect_true(file.exists(logger$file))

  expect_match(basename(logger$file), "^TestGenus_extract_\\d{8}-\\d{6}\\.log$")

  txt <- read_log(logger$file)
  expect_match(txt, "=== Genus extraction log: TestGenus")
  expect_match(console, "=== Genus extraction log: TestGenus")
})

test_that(".make_logger respects file_prefix over mode-derived stem", {
  dir <- tempfile("logger_dir_")

  logger <- make_logger(
    outdir = dir,
    genus = "IgnoredGenus",
    mode = "graft",
    file_prefix = "custom_prefix",
    enable = TRUE,
    console = FALSE
  )
  on.exit(try(close_logger(logger), silent = TRUE), add = TRUE)

  expect_match(basename(logger$file), "^custom_prefix_\\d{8}-\\d{6}\\.log$")
  expect_true(file.exists(logger$file))
})

test_that(".make_logger with file logging disabled reserves path but does not create file", {
  dir <- tempfile("logger_dir_")

  console <- capture_cat({
    logger <- make_logger(
      outdir = dir,
      genus = "DisabledFile",
      mode = "generic",
      enable = FALSE,
      console = TRUE
    )
  })
  on.exit(try(close_logger(logger), silent = TRUE), add = TRUE)

  expect_false(logger$enabled)
  expect_true(is.character(logger$file))
  expect_false(file.exists(logger$file))

  # Console header still appears because console = TRUE
  expect_match(console, "=== MacroPhyloMaker log: DisabledFile")
})

test_that(".make_logger can suppress console while still writing header to file", {
  dir <- tempfile("logger_dir_")

  console <- capture_cat({
    logger <- make_logger(
      outdir = dir,
      genus = "SilentConsole",
      mode = "generic",
      enable = TRUE,
      console = FALSE
    )
  })
  on.exit(try(close_logger(logger), silent = TRUE), add = TRUE)

  expect_equal(console, "")
  expect_true(file.exists(logger$file))

  txt <- read_log(logger$file)
  expect_match(txt, "=== MacroPhyloMaker log: SilentConsole")
})

# -------------------------------------------------------------------
# log_msg()
# -------------------------------------------------------------------

test_that("log_msg writes to console and file when both are enabled", {
  dir <- tempfile("logger_dir_")
  logger <- make_logger(
    outdir = dir,
    genus = "MsgTest",
    mode = "generic",
    enable = TRUE,
    console = TRUE
  )
  on.exit(try(close_logger(logger), silent = TRUE), add = TRUE)

  console <- capture_cat({
    log_msg(logger, "hello ", "world")
  })

  expect_match(console, "hello world")

  txt <- read_log(logger$file)
  expect_match(txt, "hello world")
})

test_that("log_msg writes only to file when console mirroring is disabled", {
  dir <- tempfile("logger_dir_")
  logger <- make_logger(
    outdir = dir,
    genus = "FileOnly",
    mode = "generic",
    enable = TRUE,
    console = FALSE
  )
  on.exit(try(close_logger(logger), silent = TRUE), add = TRUE)

  before <- read_log(logger$file)

  console <- capture_cat({
    log_msg(logger, "file only message")
  })

  after <- read_log(logger$file)

  expect_equal(console, "")
  expect_false(identical(before, after))
  expect_match(after, "file only message")
})

test_that("log_msg writes only to console when file logging is disabled", {
  dir <- tempfile("logger_dir_")
  logger <- make_logger(
    outdir = dir,
    genus = "ConsoleOnly",
    mode = "generic",
    enable = FALSE,
    console = TRUE
  )
  on.exit(try(close_logger(logger), silent = TRUE), add = TRUE)

  console <- capture_cat({
    log_msg(logger, "console only message")
  })

  expect_match(console, "console only message")
  expect_false(file.exists(logger$file))
})

test_that("log_msg can be fully silenced with disabled file logging and console=FALSE", {
  dir <- tempfile("logger_dir_")
  logger <- make_logger(
    outdir = dir,
    genus = "Silent",
    mode = "generic",
    enable = FALSE,
    console = FALSE
  )
  on.exit(try(close_logger(logger), silent = TRUE), add = TRUE)

  console <- capture_cat({
    log_msg(logger, "you should not see this")
  })

  expect_equal(console, "")
  expect_false(file.exists(logger$file))
})

test_that("log_msg adds timestamp when requested", {
  dir <- tempfile("logger_dir_")
  logger <- make_logger(
    outdir = dir,
    genus = "TimeTest",
    mode = "generic",
    enable = TRUE,
    console = TRUE
  )
  on.exit(try(close_logger(logger), silent = TRUE), add = TRUE)

  console <- capture_cat({
    log_msg(logger, "timed message", .timestamp = TRUE)
  })

  expect_match(console, "^\\[[0-9]{2}:[0-9]{2}:[0-9]{2}\\] timed message")

  txt <- read_log(logger$file)
  expect_match(txt, "\\[[0-9]{2}:[0-9]{2}:[0-9]{2}\\] timed message")
})

test_that("log_msg works with non-logger input by writing to console only", {
  console <- capture_cat({
    log_msg(NULL, "console only")
  })

  expect_match(console, "console only")
})

# -------------------------------------------------------------------
# log_section()
# -------------------------------------------------------------------

test_that("log_section writes boxed UTF-8 section to console and file", {
  old <- options(MacroPhyloMaker.ascii = FALSE, AntPhyloMaker.ascii = FALSE)
  on.exit(options(old), add = TRUE)

  dir <- tempfile("logger_dir_")
  logger <- make_logger(
    outdir = dir,
    genus = "SectionTest",
    mode = "generic",
    enable = TRUE,
    console = TRUE
  )
  on.exit(try(close_logger(logger), silent = TRUE), add = TRUE)

  console <- capture_cat({
    log_section(logger, "Settings")
  })

  expect_match(console, "┌")
  expect_match(console, "│ Settings │")
  expect_match(console, "└")

  txt <- read_log(logger$file)
  expect_match(txt, "│ Settings │")
})

test_that("log_section respects ASCII fallback option", {
  old <- options(MacroPhyloMaker.ascii = TRUE, AntPhyloMaker.ascii = FALSE)
  on.exit(options(old), add = TRUE)

  dir <- tempfile("logger_dir_")
  logger <- make_logger(
    outdir = dir,
    genus = "AsciiSection",
    mode = "generic",
    enable = TRUE,
    console = TRUE
  )
  on.exit(try(close_logger(logger), silent = TRUE), add = TRUE)

  console <- capture_cat({
    log_section(logger, "Input")
  })

  expect_match(console, "\\+")
  expect_match(console, "\\| Input \\|")

  txt <- read_log(logger$file)
  expect_match(txt, "\\| Input \\|")
})

test_that("log_section respects file/console conditional logging", {
  dir <- tempfile("logger_dir_")
  logger <- make_logger(
    outdir = dir,
    genus = "FileOnlySection",
    mode = "generic",
    enable = TRUE,
    console = FALSE
  )
  on.exit(try(close_logger(logger), silent = TRUE), add = TRUE)

  console <- capture_cat({
    log_section(logger, "Output")
  })

  expect_equal(console, "")

  txt <- read_log(logger$file)
  expect_match(txt, "Output")
})

# -------------------------------------------------------------------
# .close_logger()
# -------------------------------------------------------------------

test_that(".close_logger is safe and idempotent", {
  dir <- tempfile("logger_dir_")
  logger <- make_logger(
    outdir = dir,
    genus = "CloseTest",
    mode = "generic",
    enable = TRUE,
    console = FALSE
  )

  expect_true(inherits(logger$con, "connection"))
  expect_true(isOpen(logger$con))

  expect_no_error(close_logger(logger))
  expect_no_error(close_logger(logger))
})

test_that(".close_logger is safe for loggers without open file connections", {
  dir <- tempfile("logger_dir_")
  logger <- make_logger(
    outdir = dir,
    genus = "NoFileTest",
    mode = "generic",
    enable = FALSE,
    console = FALSE
  )

  expect_silent(close_logger(logger))
})