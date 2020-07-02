library(Rfastp)
context("rgencore function")

inputbamfile <-  system.file("extdata", "ex1_sorted.bam", package="Rfastp")
outputbamfile <- "ex1_rgencore.bam"
reference <- system.file("extdata", "myreference.fa", package="Rfastp")

test_that("malformed input gives error", {
    expect_error(
        rgencore(inBam=inputbamfile, verbose=FALSE)
    )
})

test_that("correctly formated input works", {
    report <- rgencore(inBam=inputbamfile, outBam=outputbamfile, 
        refFile=reference, verbose=FALSE)
    expect_is(report, "list")
})
