library(Rfastp)
context("rfastp function")

se_read1 <- system.file("extdata","Fox3_Std_small.fq.gz",package="Rfastp")
pe001_read1 <- system.file("extdata","splited_001_R1.fastq.gz",package="Rfastp")
pe002_read1 <- system.file("extdata","splited_002_R1.fastq.gz",package="Rfastp")
pe003_read1 <- system.file("extdata","splited_003_R1.fastq.gz",package="Rfastp")
pe004_read1 <- system.file("extdata","splited_004_R1.fastq.gz",package="Rfastp")
allR1 <- c(pe001_read1, pe002_read1, pe003_read1, pe004_read1)

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


test_that("catfastq: malformed input gives error", {
    expect_error(
        catfastq(output = "merged1_R1.fastq.gz",
             inputFiles = "myfiles1.fastq.gz,myfiles2.fastq.gz")
    )
})

test_that("catfastq: correctly formated input works", {
   exitcode <- catfastq(output = "merged1_R1.fastq.gz", inputFiles = allR1)
   expect_is(exitcode, "integer")
})
