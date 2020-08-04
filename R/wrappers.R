#' concatenate fastq file in R
#'
#' concatenate multiple files into a big file.
#'
#' @param output output file name [string]
#' @param inputFiles a vector of input file names [vector]
#' @param append a logical indicating append the files to a file already
#'      exists.
#' @param paired a logical indicating split the input files into two halves.
#'      the first half merged into read1, the second half merged into read2.
#' @param shuffled a logical indicating split the input file list into two 
#'      halves. The R1/R2 files are inteleaved in the inputFiles vector. 
#'     
#'
#' @return no returns.
#' @author Wei Wang
#' @export
#'
#' @examples
#'
#' pe001_read1 <- system.file("extdata","splited_001_R1.fastq.gz",
#'      package="Rfastp")
#' pe002_read1 <- system.file("extdata","splited_002_R1.fastq.gz",
#'      package="Rfastp")
#' pe003_read1 <- system.file("extdata","splited_003_R1.fastq.gz",
#'      package="Rfastp")
#' pe004_read1 <- system.file("extdata","splited_004_R1.fastq.gz",
#'      package="Rfastp")
#'
#' pe001_read2 <- system.file("extdata","splited_001_R2.fastq.gz",
#'      package="Rfastp")
#' pe002_read2 <- system.file("extdata","splited_002_R2.fastq.gz",
#'      package="Rfastp")
#' pe003_read2 <- system.file("extdata","splited_003_R2.fastq.gz",
#'      package="Rfastp")
#' pe004_read2 <- system.file("extdata","splited_004_R2.fastq.gz",
#'      package="Rfastp")
#'
#' allR1 <- c(pe001_read1, pe002_read1, pe003_read1, pe004_read1)
#' allR2 <- c(pe001_read2, pe002_read2, pe003_read2, pe004_read2)
#'
#' allreads <- c(allR1, allR2)
#' allreads_shuffled <- c(pe001_read1, pe001_read2, pe002_read1, pe002_read2,
#'                pe003_read1, pe003_read2, pe004_read1, pe004_read2)
#'
#' # a normal single-end concatenation.
#'
#' catfastq(output = "merged1_R1.fastq.gz", inputFiles = allR1)
#'
#' # a normal paired-end concatenation.
#'
#' catfastq(output = "merged_paired", inputFiles = allreads, paired=TRUE)
#'
#' # append to exists files (paired-end)
#'
#' catfastq(output="append_paired", inputFiles=allreads, append=TRUE,
#'      paired=TRUE)
#'
#' # input paired-end files are shuffled.
#'
#' catfastq(output="shuffled_paired", inputFiles=allreads_shuffled,
#'      paired=TRUE, shuffled=TRUE)

catfastq <- function(output, inputFiles, append=FALSE, paired=FALSE,
        shuffled=FALSE) {
    if (paired) {
        if (length(inputFiles) %% 2 != 0) {
            stop("the number of input files is not an even number!")
        }
        if (shuffled) {
            r1files <- inputFiles[seq(1, length(inputFiles), 2)]
            r2files <- inputFiles[seq(2, length(inputFiles), 2)]
        }
        else {
            r1files <- inputFiles[seq(1, length(inputFiles)/2)]
            r2files <- inputFiles[seq(length(inputFiles)/2+1,
                        length(inputFiles))]
        }

        if (append) {
            exitcode <- rcat(output=paste0(output, "_R1.fastq.gz"), 
                r1files, length(inputFiles)/2)
            exitcode <- rcat(output=paste0(output, "_R2.fastq.gz"),
                r2files, length(inputFiles)/2)
        }
        else {
            if (file.exists(paste0(output, "_R1.fastq.gz")) | 
                file.exists(paste0(output, "_R2.fastq.gz")) ) {
                stop("output file exists already! please change it!")
            }
            exitcode <- rcat(output=paste0(output, "_R1.fastq.gz"), 
                r1files, length(inputFiles)/2)
            exitcode <- rcat(output=paste0(output, "_R2.fastq.gz"),
                r2files, length(inputFiles)/2)
        }

    }
    else {
        if (append) {
            exitcode <- rcat(output=output, inputFiles, length(inputFiles))
        }
        else {
            if (file.exists(output)) {
                stop("output file exists already! please change it!")
            }
            exitcode <- rcat(output=output, inputFiles, length(inputFiles))
        }
    }
}


#' Summary of Quality Control
#'
#' generate a data frame of the QC summary. 
#'
#' @param json the output json of function rfastq. [json]
#'     
#' @return a data frame.
#' @author Wei Wang
#' @export
#'
#' @examples
#'
#' se_read1 <- system.file("extdata","Fox3_Std_small.fq.gz",package="Rfastp")
#' se_json_report <- rfastp(read1 = se_read1, outputFastq = './rfastp_test_se',
#'    thread = 4)
#' df_summary <- qcSummary(se_json_report)
#'

qcSummary <- function(json) {
    if (! "merged_and_filtered" %in% names(json)) {
        data.frame("Before_QC" = unlist(json$summary$before_filtering),
            "After_QC" = unlist(json$summary$after_filtering), 
            row.names = names(json$summary$before_filtering)) 
    }
    else {
        data.frame("Before_QC" = unlist(json$summary$before_filtering),
            "After_QC" = c(unlist(json$summary$after_filtering)[seq(1,7)],
            NA, unlist(json$summary$after_filtering)[8]), 
            row.names = names(json$summary$before_filtering)) 
    }
}


#' Summary of trimming
#'
#' generate a data frame of the trim summary. 
#'
#' @param json the output json of function rfastq. [json]
#'     
#' @return a data frame.
#' @author Wei Wang
#' @export
#'
#' @examples
#'
#' se_read1 <- system.file("extdata","Fox3_Std_small.fq.gz",package="Rfastp")
#' se_json_report <- rfastp(read1 = se_read1, outputFastq = './rfastp_test_se',
#'    thread = 4)
#' trim_summary <- trimSummary(se_json_report)
#'

trimSummary <- function(json) {
    if (! "read2_before_filtering" %in% names(json)) {
        data.frame("Count" = c(unlist(json$filtering_result),
                unlist(json$adapter_cutting[seq(1,3)]),
                unlist(json$adapter_cutting$read1_adapter_counts)),
            row.names = c(names(json$filtering_result),
                names(json$adapter_cutting)[seq(1,3)],
                names(json$adapter_cutting$read1_adapter_counts)) ) 
    }
    else {
        data.frame("Count" = c(unlist(json$filtering_result),
                unlist(json$adapter_cutting[seq(1,4)]),
                unlist(json$adapter_cutting$read1_adapter_counts)),
            row.names = c(names(json$filtering_result),
                names(json$adapter_cutting)[seq(1,4)],
                names(json$adapter_cutting$read1_adapter_counts)) ) 
    }
}


#' Plot of Base Quality and GC Content.
#'
#' generate a ggplot2 object of Base Quality/GC content before and after. 
#'
#' @param json the output json of function rfastq. [json]
#' @param curves plots for Base Quality("quality_curves") 
#'      or GC content("content_curves"). default is "quality_curves"
#'     
#' @return a ggplot2 object.
#' @author Wei Wang
#' @import ggplot2 reshape2
#' @export
#'
#' @examples
#'
#' se_read1 <- system.file("extdata","Fox3_Std_small.fq.gz",package="Rfastp")
#' se_json_report <- rfastp(read1 = se_read1, outputFastq = './rfastp_test_se',
#'    thread = 4)
#' # Base Quality plot is the default output:
#' p1 <- curvePlot(se_json_report)
#' p1
#'
#' p2 <- curvePlot(se_json_report, curves = "content_curves")

curvePlot <- function(json, curves = "quality_curves") {
    #globalVariables(c("position", "value", "variable"))
    if ("merged_and_filtered" %in% names(json)) {
        df1bf <- data.frame(json$read1_before_filtering[[curves]])
        df2bf <- data.frame(json$read2_before_filtering[[curves]])
        dfmerged <- data.frame(json$merged_and_filtered[[curves]])
        df1bf$position <- as.integer(rownames(df1bf))
        df2bf$position <- as.integer(rownames(df2bf))
        dfmerged$position <- as.integer(rownames(dfmerged))
        df1bf$readtype <- factor("Read1 Before QC", 
            levels=c("Read1 Before QC", "Read2 Before QC", "Merged After QC"))
        df2bf$readtype <- factor("Read2 Before QC", 
            levels=c("Read1 Before QC", "Read2 Before QC", "Merged After QC"))
        dfmerged$readtype <- factor("Merged After QC", 
            levels=c("Read1 Before QC", "Read2 Before QC", "Merged After QC"))
        tbl4plot <- rbind(melt(df1bf, c("position", "readtype")),
            melt(df2bf, c("position", "readtype")),
            melt(dfmerged, c("position", "readtype"))) 
    }
    else if (! "read2_before_filtering" %in% names(json)) {
        df1bf <- data.frame(json$read1_before_filtering[[curves]])
        df1af <- data.frame(json$read1_after_filtering[[curves]])
        df1bf$position <- as.integer(rownames(df1bf))
        df1af$position <- as.integer(rownames(df1af))
        df1bf$readtype <- factor("Read1 Before QC",
                levels=c("Read1 Before QC", "Read1 After QC"))
        df1af$readtype <- factor("Read1 After QC",
                levels=c("Read1 Before QC", "Read1 After QC"))
        tbl4plot <- rbind(melt(df1bf, c("position", "readtype")),
                melt(df1af, c("position", "readtype"))) 
    }
    else {
        df1bf <- data.frame(json$read1_before_filtering[[curves]])
        df1af <- data.frame(json$read1_after_filtering[[curves]])
        df2bf <- data.frame(json$read2_before_filtering[[curves]])
        df2af <- data.frame(json$read2_after_filtering[[curves]])
        df1bf$position <- as.integer(rownames(df1bf))
        df1af$position <- as.integer(rownames(df1af))
        df2bf$position <- as.integer(rownames(df2bf))
        df2af$position <- as.integer(rownames(df2af))
        df1bf$readtype <- factor("Read1 Before QC", levels=c("Read1 Before QC",
                    "Read1 After QC", "Read2 Before QC", "Read2 After QC"))
        df1af$readtype <- factor("Read1 After QC", levels=c("Read1 Before QC",
                    "Read1 After QC", "Read2 Before QC", "Read2 After QC"))
        df2bf$readtype <- factor("Read2 Before QC", levels=c("Read1 Before QC",
                    "Read1 After QC", "Read2 Before QC", "Read2 After QC"))
        df2af$readtype <- factor("Read2 After QC", levels=c("Read1 Before QC",
                    "Read1 After QC", "Read2 Before QC", "Read2 After QC"))
        tbl4plot <- rbind(melt(df1bf, c("position", "readtype")),
                    melt(df1af, c("position", "readtype")),
                    melt(df2bf, c("position", "readtype")),
                    melt(df2af, c("position", "readtype"))) 
    }
    ggplot(tbl4plot, aes_string(x="position", y="value", group=1)) + 
            geom_line(aes_string(color="variable")) + ylab("Base Quality") + 
            facet_wrap(~ readtype)
}


#' R wrap of fastp
#'
#' Quality control (Cut adapter, low quality trimming, UMI handling, and etc.)
#' of fastq files.
#'
#' @param read1 read1 input file name(s). [vector]
#' @param read2 read2 input file name(s). [vector]
#' @param outputFastq string of /path/prefix for output fastq [string]
#' @param unpaired for PE input, output file name for reads which the mate
#'      reads failed to pass the QC [string], default NULL, discard it. [string]
#' @param failedOut file to store reads that cannot pass the filters default
#'       NULL, discard it. [string]
#' @param merge for PE input, A logical(1) indicating whether merge each pair
#'      of reads into a single read if they are overlaped, unmerged reads will
#'      be write to `output` file. Default is FALSE. the `mergeOut` must be
#'      set if TRUE.
#' @param mergeOut under `merge` mode, file to store the merged reads. [string]
#' @param phred64 A logical indicating whether the input is using phred64
#'      scoring (it will be converted to phred33, so the output will still be
#'.      phred33)
#' @param interleaved A logical indicating whether <read1> is an interleaved
#'      FASTQ which contains both read1 and read2. Default is FALSE.
#' @param fixMGIid the MGI FASTQ ID format is not compatible with many BAM
#'      operation tools, enable this option to fix it. Default is FALSE
#' @param adapterTrimming A logical indicating whether run adapter trimming.
#'      Default is `TRUE`
#' @param adapterSequenceRead1 the adapter for read1. For SE data, if not
#'      specified, the adapter will be auto-detected. For PE data, this is used
#'      if R1/R2 are found not overlapped.
#' @param adapterSequenceRead2 the adapter for read2 (PE data only). This is
#'      used if R1/R2 are found not overlapped. If not specified, it will be the
#'      same as <adapterSequenceRead1>
#' @param adapterFasta specify a FASTA file to trim both read1 and read2 (if
#'      PE) by all the sequences in this FASTA file.
#' @param trimFrontRead1 trimming how many bases in front for read1, default 
#'      is 0.
#' @param trimTailRead1 trimming how many bases in tail for read1, default is 0'
#' @param trimFrontRead2 trimming how many bases in front for read2. If it's not
#'      specified, it will follow read1's settings
#' @param trimTailRead2 trimming how many bases in tail for read2. If it's not
#'      specified, it will follow read1's settings
#' @param maxLengthRead1 if read1 is longer than maxLengthRead1, then trim read1
#'      at its tail to make it as long as maxLengthRead1 Default 0 means no
#'      limitation.
#' @param maxLengthRead2 if read2 is longer than maxLengthRead2, then trim read2
#'      at its tail to make it as long as maxLengthRead2. Default 0 means no
#'      limitation. If it's not specified, it will follow read1's settings.
#' @param forceTrimPolyG A logical indicating force polyG tail trimming,
#'      trimming is only automatically enabled for Illumina NextSeq/NovaSeq
#'.     data.
#' @param disableTrimPolyG A logical indicating disable polyG tail trimming.
#' @param minLengthPolyG the minimum length to detect polyG in the read tail.
#'      10 by default.
#' @param trimPolyX A logical indicating force polyX tail trimming.
#' @param minLengthPolyX the minimum length to detect polyX in the read tail.
#'      10 by default.
#' @param cutLowQualFront A logical indiccating move a sliding window from
#'      front (5') to tail, drop the bases in the window if its mean quality
#'      < threshold, stop otherwise. Default is `FALSE`
#' @param cutLowQualTail A logical indiccating move a sliding window from
#'      tail (3') to front, drop the bases in the window if its mean quality
#'      < threshold, stop otherwise. Default is `FALSE`
#' @param cutSlideWindowRight A logical indicating move a sliding window from
#'      front to tail, if meet one window with mean quality < threshold, drop
#'      the bases in the window and the right part, and then stop. Default is
#'      `FALSE`
#' @param cutWindowSize the window size option shared by cutLowQualFront,
#'      cutLowQualTail, or cutSlideWindowRight. Range: 1~1000, default: 4
#' @param cutMeanQual the mean quality requirement option shared by
#'      cutLowQualFront, cutLowQualTail or cutSlideWindowRight. Range: 1~36,
#'      default: 20
#' @param cutFrontWindowSize the window size option of cutLowQualFront, default
#'      to cutWindowSize if not specified. default: 4
#' @param cutFrontMeanQual the mean quality requirement option for
#'      cutLowQualFront, default to cutMeanQual if not specified. default: 20
#' @param cutTailWindowSize the window size option of cutLowQualTail, default
#'      to cutWindowSize if not specified. default: 4
#' @param cutTailMeanQual the mean quality requirement option for
#'      cutLowQualTail, default to cutMeanQual if not specified. default: 20
#' @param cutSlideWindowSize the window size option of cutSlideWindowRight,
#'      default to cutWindowSize if not specified. default: 4
#' @param cutSlideWindowQual the mean quality requirement option for
#'      cutSlideWindowRight, default to cutMeanQual if not specified. default:
#'      20
#' @param qualityFiltering A logical indicating run quality filtering.
#'      Default is `TRUE`.
#' @param qualityFilterPhred the minimum quality value that a base is
#'      qualified. Default 15 means phred quality >=Q15 is qualified.
#' @param qualityFilterPercent Maximum percents of bases are allowed to be
#'      unqualified (0~100). Default 40 means 40\%
#' @param maxNfilter maximum number of N allowed in the sequence. read/pair is
#'      discarded if failed to pass this filter. Default is 5
#' @param averageQualFilter if one read's average quality score < 
#'      `averageQualFilter`, then this read/pair is discarded. Default 0 means
#'       no requirement.
#' @param lengthFiltering A logical indicating whether run lenght filtering.
#'      Default: TRUE
#' @param minReadLength reads shorter than minReadLength will be discarded,
#'      default is 15.
#' @param maxReadLength reads longer than maxReadLength will be discarded,
#'      default 0 means no limitation.
#' @param lowComplexityFiltering A logical indicating whethere run low
#'      complexity filter. The complexity is defined as the percentage of base
#'      that is different from its next base (base[i] != base[i+1]). Default is
#'      `FALSE`
#' @param minComplexity the threshold for low complexity filter (0~100).
#'      Default is 30, which means 30\% complexity is required. (int [=30])
#' @param index1Filter specify a file contains a list of barcodes of index1
#'      to be filtered out, one barcode per line.
#' @param index2Filter specify a file contains a list of barcodes of index2
#'      to be filtered out, one barcode per line.
#' @param maxIndexMismatch the allowed difference of index barcode for
#'      index filtering, default 0 means completely identical.
#' @param correctionOverlap A logical indicating run base correction in
#'      overlapped regions (only for PE data), default is `FALSE`
#' @param minOverlapLength the minimum length to detect overlapped region of
#'      PE reads. This will affect overlap analysis based PE merge, adapter
#'      trimming and correction. 30 by default.
#' @param maxOverlapMismatch the maximum number of mismatched bases to detect
#'      overlapped region of PE reads. This will affect overlap analysis
#'      based PE merge, adapter trimming and correction. 5 by default.
#' @param maxOverlapMismatchPercentage the maximum percentage of mismatched
#'      bases to detect overlapped region of PE reads. This will affect
#'      overlap analysis based PE merge, adapter trimming and correction.
#'      Default 20 means 20\%
#' @param umi A logical indicating whethere preprocessing unique molecular
#'      identifier (UMI). Default: `FALSE`
#' @param umiLoc specify the location of UMI, can be
#'      (index1/index2/read1/read2/per_index/per_read)
#' @param umiLength length of UMI if the UMI is in read1/read2.
#' @param umiPrefix an string indication the following string is UMI
#'      (i.e. prefix=UMI, UMI=AATTCG, final=UMIAATTCG). Only letters,
#'      numbers, and '#" allowed. No prefix by default.
#' @param umiNoConnection an logical indicating remove "_" between the UMI
#'      prefix string and the UMI string. Default is FALSE.
#' @param umiIgnoreSeqNameSpace an logical indicating ignore the space
#'      in the sequence name. Default is FALSE, the umi tag will be
#'      inserted into the sequence name before the first SPACE.
#' @param umiSkipBaseLength if the UMI is in read1/read2, skip
#'      `umiSkipBaseLength` bases following UMI, default is 0.
#' @param overrepresentationAnalysis A logical indicating overrepresentation
#'      analysis. Default is `FALSE`
#' @param overrepresentationSampling one in `overrepresentationSampling`
#'      reads will be computed for overrepresentation analysis (1~10000),
#'      smaller is slower, default is 20.
#' @param splitOutput number of files to be splitted (2~999). a sequential
#'      number prefix will be added to output name. Default is 0 (no split)
#' @param splitByLines split output by limiting lines of each file(>=1000), a
#'      sequential number prefix will be added to output name ( 0001.out.fq,
#'      0002.out.fq...), default is 0 (disabled).
#' @param thread owrker thread number, default is 2
#' @param verbose output verbose log information
#'
#' @return returns a json object of the report.
#' @author Thomas Carroll, Wei Wang
#' @importFrom rjson fromJSON
#' @export
#'
#' @examples
#'
#' # preprare for the input and output files.
#' # if the output file exists, it will be OVERWRITEN.
#'
#' se_read1 <- system.file("extdata","Fox3_Std_small.fq.gz",package="Rfastp")
#' pe_read1 <- system.file("extdata","reads1.fastq.gz",package="Rfastp")
#' pe_read2 <- system.file("extdata","reads2.fastq.gz",package="Rfastp")
#' output_fastq  <- './rfastp_test'
#'
#'
#' # a normal single-end file
#'
#' se_json_report <- rfastp(read1 = se_read1, outputFastq = './rfastp_test_se',
#'    thread = 4)
#'
#'
#' # merge paired-end data by overlap:
#'
#' pe_json_report <- rfastp(read1 = pe_read1, read2 = pe_read2, merge = TRUE,
#'    outputFastq = './rfastp_pe_test_unpaired',
#'    mergeOut = './rfastp_pe_test_merged.fastq.gz')
#'
#'
#' # a clipr example
#'
#' clipr_json_report <- rfastp(read1 = se_read1, 
#'    outputFastq = './rfastp_test_clipr',
#'    disableTrimPolyG = TRUE,
#'    cutLowQualFront = TRUE,
#'    cutFrontWindowSize = 29,
#'    cutFrontMeanQual = 20,
#'    cutLowQualTail = TRUE,
#'    cutTailWindowSize = 1,
#'    cutTailMeanQual = 5,
#'    minReadLength = 29,
#'    adapterSequenceRead1 = 'GTGTCAGTCACTTCCAGCGG'
#')


rfastp <- function(read1, read2="", outputFastq, unpaired="",
    failedOut="", merge=FALSE, mergeOut="", phred64=FALSE, interleaved=FALSE,
    fixMGIid=FALSE, adapterTrimming=TRUE, adapterSequenceRead1="auto",
    adapterSequenceRead2="auto", adapterFasta="", trimFrontRead1=0,
    trimTailRead1=0, trimFrontRead2=0, trimTailRead2=0, maxLengthRead1=0,
    maxLengthRead2=0, forceTrimPolyG=FALSE, disableTrimPolyG=FALSE,
    minLengthPolyG=10, trimPolyX=FALSE, minLengthPolyX=10, cutWindowSize=4,
    cutLowQualTail=FALSE, cutSlideWindowRight=FALSE, cutLowQualFront=FALSE,
    cutMeanQual=20, cutFrontWindowSize=4, cutFrontMeanQual=20, 
    cutTailWindowSize=4, cutTailMeanQual=20, cutSlideWindowSize=4,
    cutSlideWindowQual=20, qualityFiltering=TRUE, qualityFilterPhred=15,
    qualityFilterPercent=40, maxNfilter=5, averageQualFilter=0, 
    lengthFiltering=TRUE, minReadLength=15, maxReadLength=0,
    lowComplexityFiltering=FALSE, minComplexity=30, index1Filter="",
    index2Filter="", maxIndexMismatch=0, correctionOverlap=FALSE,
    minOverlapLength=30, maxOverlapMismatch=5, maxOverlapMismatchPercentage=20,
    umi=FALSE, umiLoc="", umiLength=0, umiPrefix="", umiSkipBaseLength=0,
    umiNoConnection=FALSE, umiIgnoreSeqNameSpace=FALSE,
    overrepresentationAnalysis=FALSE, overrepresentationSampling=20,
    splitOutput=0, splitByLines=0, thread=2, verbose=TRUE) {

    multipleInput = length(read1) > 1   
    if ( multipleInput ) {
        infilesR1 = read1
        ramstr <-  rawToChar(as.raw(sample(c(65:90,97:122), 5, replace=TRUE)))
        read1 <- paste0("catInput_", ramstr, "_R1.fastq.gz")
        exitcode <- rcat(output=read1, infilesR1, length(infilesR1))
        if (read2 != "") {
            infilesR2 = as.list(unlist(strsplit(read2,",")))
            if (length(infilesR2) != length(infilesR2)) {
                stop("the file number of Read1 and Read2 are not identical!")
            }
            read2 <- paste0("catInput_", ramstr, "_R2.fastq.gz")
            if (file.exists(read2)) {
                stop("the tmp concatenated R2 file exists already!")
            }
            exitcode <- rcat(output=read2, infilesR2, length(infilesR1))
        }
    }
    else if (length(read2) > 1) {
        stop("please double check the read2 file names, there is only one
            input file in read1.")
    }        

    if (umi & umiPrefix != "" & !umiNoConnection) {
        umiPrefix <- paste0(umiPrefix, "_")
    }

    exitcode <- runFastp(read1=read1, read2=read2, outputFastq=outputFastq,
        unpaired=unpaired, failedOut=failedOut, merge=merge, mergeOut=mergeOut,
        phred64=phred64, interleaved=interleaved, fixMGIid=fixMGIid, 
        adapterTrimming=adapterTrimming, 
        adapterSequenceRead1=adapterSequenceRead1, 
        adapterSequenceRead2=adapterSequenceRead2,
        adapterFasta=adapterFasta, trimFrontRead1=trimFrontRead1, 
        trimTailRead1=trimTailRead1, trimFrontRead2=trimFrontRead2, 
        trimTailRead2=trimTailRead2, maxLengthRead1=maxLengthRead1,
        maxLengthRead2=maxLengthRead2, forceTrimPolyG=forceTrimPolyG, 
        disableTrimPolyG=disableTrimPolyG, minLengthPolyG=minLengthPolyG, 
        trimPolyX=trimPolyX, minLengthPolyX=minLengthPolyX,
        cutLowQualFront=cutLowQualFront, cutLowQualTail=cutLowQualTail, 
        cutSlideWindowRight=cutSlideWindowRight, cutWindowSize=cutWindowSize,
        cutMeanQual=cutMeanQual, cutFrontWindowSize=cutFrontWindowSize,
        cutFrontMeanQual=cutFrontMeanQual, cutTailWindowSize=cutTailWindowSize,
        cutTailMeanQual=cutTailMeanQual, cutSlideWindowSize=cutSlideWindowSize,
        cutSlideWindowQual=cutSlideWindowQual, maxReadLength=maxReadLength,
        qualityFiltering=qualityFiltering,qualityFilterPhred=qualityFilterPhred,
        qualityFilterPercent=qualityFilterPercent, maxNfilter=maxNfilter, 
        averageQualFilter=averageQualFilter,
        lengthFiltering=lengthFiltering, minReadLength=minReadLength,
        lowComplexityFiltering=lowComplexityFiltering, 
        minComplexity=minComplexity, index1Filter=index1Filter, 
        index2Filter=index2Filter, maxIndexMismatch=maxIndexMismatch,
        correctionOverlap=correctionOverlap, minOverlapLength=minOverlapLength,
        maxOverlapMismatch=maxOverlapMismatch, 
        maxOverlapMismatchPercentage=maxOverlapMismatchPercentage, 
        umi=umi, umiLoc=umiLoc, umiLength=umiLength, umiPrefix=umiPrefix, 
        umiSkipBaseLength=umiSkipBaseLength,
        umiIgnoreSeqNameSpace=umiIgnoreSeqNameSpace, 
        overrepresentationAnalysis=overrepresentationAnalysis,
        overrepresentationSampling=overrepresentationSampling,
        splitOutput=splitOutput, splitByLines=splitByLines, thread=thread, 
        verbose=verbose)

    if (multipleInput) {
        file.remove(read1)
        if (read2 != "") {
            file.remove(read2)
        }
    }
    return(fromJSON(file = paste0(outputFastq, ".json")))
}
