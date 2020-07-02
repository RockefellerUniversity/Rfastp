library(Rfastp)
context("rfastp function")

se_read1 <- system.file("extdata","Fox3_Std_small.fq.gz",package="Rfastp")

test_that("malformed input gives error", {
    expect_error(
        rfastp(read1=se_read1, verbose=FALSE)
    )
})

test_that("correctly formatted input  works", {
    report <- rfastp(read1 = se_read1, outputFastq = './rfastp_test_se',
        thread = 4, verbose=FALSE)
    expect_is(report, "list")
})
