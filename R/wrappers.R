
#' Call the fastp binary with additional arguments.
#'
#' @keywords internal
#'
#' @param bin The name of the fastp binary
#' @param execute Logical scalar, whether to execute the command. If FALSE,
#'   return a string with the shell command.
#' @param args A character string containing the arguments that will be passed
#'   to the binary.
#'
#' @return If \code{execute} is TRUE, returns the console output of running the
#'   fastp command. If \code{execute} is FALSE, returns the shell command.
#'

.fastpBin <- function(bin="fastp", args="", execute=TRUE) {
  if (is.null(args) || args=="") {
    stop("The fastp binary needs to be called with additional arguments")
  }
  args <- gsub("^ *| *$", "", args)
  bin <- match.arg(bin)
  call <- paste(shQuote(file.path(system.file(package="Rfastp"), bin)), args)
  if (!execute) {
    return(call)
  }
  output <- system(call, intern=TRUE)
  return(output)
}

#' Call the gencore binary with additional arguments.
#'
#' Adapted from the gencore package
#'
#' @keywords internal
#'
#' @param bin The name of the gencore binary
#' @param execute Logical scalar, whether to execute the command. If FALSE,
#'   return a string with the shell command.
#' @param args A character string containing the arguments that will be passed
#'   to the binary.
#'
#' @return If \code{execute} is TRUE, returns the console output of running the
#'   gencore command. If \code{execute} is FALSE, returns the shell command.
#'
.gencoreBin <- function(bin="gencore", args="", execute=TRUE) {
  if (is.null(args) || args=="") {
    stop("The fastp binary needs to be called with additional arguments")
  }
  args <- gsub("^ *| *$", "", args)
  bin <- match.arg(bin)
  call <- paste(shQuote(file.path(system.file(package="Rfastp"), bin)), args)
  if (!execute) {
    return(call)
  }
  output <- system(call, intern=TRUE)
  return(output)
}

#' R wrap of fastp
#'
#' Quality control (Cut adapter, low quality trimming, UMI handling, and etc.) of fastq files.
#'
#' @param read1 read1 input file name [string]
#' @param read2 read2 input file name [string]
#' @param outputFastq string of /path/prefix for output fastq [string]
#' @param unpaired for PE input, output file name for reads which the mate
#'      reads failed to pass the QC [string], default NULL, discard it. [string]
#' @param failedOut file to store reads that cannot pass the filters default
#       NULL, discard it. [string]
#' @param merge for PE input, A logical(1) indicating whether merge each pair 
#'      of reads into a single read if they are overlaped, unmerged reads will
#'      be write to `output` file. Default is FALSE. the `mergeOut` must be 
#'      set if TRUE.
#' @param mergeOut under `merge` mode, file to store the merged reads. [string]
#' @param compressLevel compression level for gzip output (1 ~ 9). 1 is fastest, 
#'      9 is smallest, default is 4. 
#' @param phred64 A logical indicating whether the input is using phred64 scoring
#'      (it will be converted to phred33, so the output will still be phred33)
#' @param interleaved A logical indicating whether <read1> is an interleaved FASTQ 
#'      which contains both read1 and read2. Default is FALSE.
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
#' @param adapterFasta specify a FASTA file to trim both read1 and read2 (if PE) 
#'      by all the sequences in this FASTA file.
#' @param trimFrontRead1 trimming how many bases in front for read1, default is 0
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
#'      trimming is only automatically enabled for Illumina NextSeq/NovaSeq data.
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
#' @param cutMeanQual the mean quality requirement option shared by cutLowQualFront, 
#'      cutLowQualTail or cutSlideWindowRight. Range: 1~36 default: 20
#' @param cutFrontWindowSize the window size option of cutLowQualFront, default to 
#'      cutWindowSize if not specified. default: 4
#' @param cutFrontMeanQual the mean quality requirement option for cutLowQualFront, 
#'      default to cutMeanQual if not specified. default: 20
#' @param cutTailWindowSize the window size option of cutLowQualTail, default to 
#'      cutWindowSize if not specified. default: 4
#' @param cutTailMeanQual the mean quality requirement option for cutLowQualTail, 
#'      default to cutMeanQual if not specified. default: 20
#' @param cutSlideWindowSize the window size option of cutSlideWindowRight, 
#'      default to cutWindowSize if not specified. default: 4
#' @param cutSlideWindowQual the mean quality requirement option for 
#'      cutSlideWindowRight, default to cutMeanQual if not specified. default: 20
#' @param qualityFiltering A logical indicating run quality filtering. 
#'      Default is `TRUE`.
#' @param qualityFilterPhred the minimum quality value that a base is qualified. 
#'      Default 15 means phred quality >=Q15 is qualified.
#' @param qualityFilterPercent Maximum percents of bases are allowed to be 
#'      unqualified (0~100). Default 40 means 40\%
#' @param maxNfilter maximum number of N allowed in the sequence. read/pair is 
#'      discarded if failed to pass this filter. Default is 5
#' @param averageQualFilter if one read's average quality score < `averageQualFilter`,
#'       then this read/pair is discarded. Default 0 means no requirement.
#' @param lengthFiltering A logical indicating whether run lenght filtering. 
#'      Default: TRUE   
#' @param minReadLength reads shorter than minReadLength will be discarded,     
#'      default is 15.
#' @param maxReadLength reads longer than maxReadLength will be discarded, 
#'      default 0 means no limitation.
#' @param lowComplexityFiltering A logical indicating whethere run low complexity 
#'      filter. The complexity is defined as the percentage of base that is 
#'      different from its next base (base[i] != base[i+1]). Default is `FALSE`
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
#' @param umiPrefix an underline will be used to connect prefix and UMI 
#'      (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default
#' @param umiSkipBaseLength if the UMI is in read1/read2, skip 
#'      `umiSkipBaseLength` bases following UMI, default is 0.
#' @param overrepresentationAnalysis A logical indicating overrepresentation
#'      analysis. Default is `FALSE`
#' @param overrepresentationSampling one in `overrepresentationSampling` 
#'      reads will be computed for overrepresentation analysis (1~10000), 
#'      smaller is slower, default is 20.
#' @param reportTitle Title of the report, default: "Rfastp report"
#' @param splitOutput number of files to be splitted (2~999). a sequential number
#'      prefix will be added to output name. Default is 0 (no split)
#' @param splitByLines split output by limiting lines of each file(>=1000), a 
#'      sequential number prefix will be added to output name ( 0001.out.fq, 
#'      0002.out.fq...), default is 0 (disabled).
#' @param splitPrefixPaddingNum the digits for the sequential number padding (1~10), 
#'      default is 4, so the filename will be padded as 0001.xxx, 0 to disable 
#'      padding
#' @param thread owrker thread number, default is 2
#' @param verbose output verbose log information
#'
#' @return returns a json object of the report.
#' @author Thomas Carroll, Wei Wang
#' @export
#' 
#' @examples
#' 
#' require(rjson)
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
#'          thread = 4)
#' 
#'
#' # merge paired-end data by overlap: 
#'
#' pe_json_report <- rfastp(read1 = pe_read1, read2 = pe_read2, merge = TRUE, 
#'          outputFastq = './rfastp_pe_test_unpaired',
#'          mergeOut = './rfastp_pe_test_merged.fastq.gz')
#'
#'
#' # a clipr example
#'
#' clipr_json_report <- rfastp(read1 = se_read1, outputFastq = './rfastp_test_clipr',
#'          disableTrimPolyG = TRUE,
#'          cutLowQualFront = TRUE,
#'          cutFrontWindowSize = 29,
#'          cutFrontMeanQual = 20,
#'          cutLowQualTail = TRUE,
#'          cutTailWindowSize = 1,
#'          cutTailMeanQual = 5,
#'          minReadLength = 29,
#'          adapterSequenceRead1 = 'GTGTCAGTCACTTCCAGCGG'
#')


rfastp <- function( read1="", read2="", 　outputFastq="", unpaired="",
    failedOut="", merge=FALSE, mergeOut="", compressLevel=4, phred64=FALSE,
    interleaved=FALSE, fixMGIid=FALSE, adapterTrimming=TRUE, 
    adapterSequenceRead1="", adapterSequenceRead2="", adapterFasta="",
    trimFrontRead1=0, trimTailRead1=0, trimFrontRead2=0, trimTailRead2=0,
    maxLengthRead1=0, maxLengthRead2=0, 
    forceTrimPolyG=FALSE, disableTrimPolyG=FALSE, minLengthPolyG=10,
    trimPolyX=FALSE, minLengthPolyX=10, 
    cutLowQualFront=FALSE, cutLowQualTail=FALSE,
    cutSlideWindowRight=FALSE, cutWindowSize=4, cutMeanQual=20, 
    cutFrontWindowSize=4, cutFrontMeanQual=20, 
    cutTailWindowSize=4, cutTailMeanQual=20,
    cutSlideWindowSize=4, cutSlideWindowQual=20,
    qualityFiltering=TRUE, qualityFilterPhred=15, qualityFilterPercent=40,
    maxNfilter=5, averageQualFilter=0, 
    lengthFiltering=TRUE, minReadLength=15, maxReadLength=0,
    lowComplexityFiltering=FALSE, minComplexity=30, 
    index1Filter="", index2Filter="", maxIndexMismatch=0,
    correctionOverlap=FALSE, minOverlapLength=30, maxOverlapMismatch=5,
    maxOverlapMismatchPercentage=20,
    umi=FALSE, umiLoc="", umiLength=0, umiPrefix="", umiSkipBaseLength=0,
    overrepresentationAnalysis=FALSE, overrepresentationSampling=20,
    splitOutput=0, splitByLines=0, splitPrefixPaddingNum=4,  
    reportTitle="Rfastp Report", thread=2, verbose=TRUE) { 

    require("rjson")

    args <- ""

    if (read1 == "" | outputFastq == "") {
        stop("Please specify the read1 file or the output file name.")
    }
    else {
        args <- paste0(args, "-i ", read1, " -o ", outputFastq, "_R1.fastq.gz ")
    }

    if ( read2 == "" ) {
        if (interleaved) {
            args <- paste0(args, " -O ", outputFastq, "_R2.fastq.gz")
        }
	else if (merge & mergeOut == "") {
            stop("Please specify the read2 file if you want to do merge, and specify the output merged file.")
        }
    }
    else {
        if (interleaved) {
            stop("You can't specify read2 when you set interleaved as TRUE")
        }

        args <- paste0(args, "-I ", read2, " -O ", outputFastq, "_R2.fastq.gz")
        if (!unpaired == "") {
            args <- paste0(args, " --unpaired1 ", unpaired, " --unpaired2 ", unpaired)
        }

        if (adapterSequenceRead1 == "" & adapterFasta == "") {
            args <- paste0(args, " --detect_adapter_for_pe")
        }

        if (merge) {
            args <- paste0(args, " --merge --merged_out ", mergeOut)
        }
    }

    if (compressLevel != 4) {
        args <- paste0(args, " -z ", compressLevel)
    }

    if (phred64) {
        args <- paste0(args, " --phred64")
    }

    if (fixMGIid) {
        args <- paste0(args, " --fix_mgi_id")
    }

    if (! adapterTrimming) {
        args <- paste0(args, " --disable_adapter_trimming")
    }

    if (adapterSequenceRead1 != "") {
        args <- paste0(args, " --adapter_sequence " , adapterSequenceRead1)
    }

    if (adapterSequenceRead2 != "") {
        args <- paste0(args, " --adapter_sequence_r2 " , adapterSequenceRead2)
    }

    if (adapterFasta != "") {
        args <- paste0(args, " --adapter_fasta")
    }

    if (trimFrontRead1 > 0) {
        args <- paste0(args," --trim_front1 ", trimFrontRead1)
    }

    if (trimTailRead1 > 0) {
        args <- paste0(args," --trim_tail1 ", trimTailRead1) 
    }
    if (trimFrontRead2 > 0) {
        args <- paste0(args," --trim_front1 ", trimFrontRead2) 
    } 
    if (trimTailRead2 >0) {
        args <- paste0(args," --trim_tail2 ", trimTailRead2) 
    }

    if (forceTrimPolyG) {
        args <- paste0(args," -g ")
        if (minLengthPolyG != 10) {
            args <- paste0(args, " --poly_g_min_len ", minLengthPolyG) 
        } 
    }
    else if (disableTrimPolyG) {
        args <- paste0(args," -G ")
    }

    if (trimPolyX) {
        args <- paste0(args," -x ")
        if (minLengthPolyX != 10) {
            args <- paste0(args, " --poly_x_min_len ", minLengthPolyX)
        }
    }

    if (cutLowQualFront) {
        args <- paste0(args, " -5")
        if (cutFrontWindowSize != 4) {
            args <- paste0(args, " --cut_front_window_size ", cutFrontWindowSize)
        }
        if (cutFrontMeanQual != 20) {
            args <- paste0(args, " --cut_front_mean_quality ", cutFrontMeanQual)
        }
    }

    if (cutLowQualTail) {
        args <- paste0(args, " -3")
        if (cutTailWindowSize != 4) {
            args <- paste0(args, " --cut_tail_window_size ", cutTailWindowSize)
        }
        if (cutTailMeanQual != 20) {
            args <- paste0(args, " --cut_tail_mean_quality ", cutTailMeanQual)
        }
    }

    if (cutSlideWindowRight) {
        args <- paste0(args, " -r")
        if (cutSlideWindowSize != 4) {
            args <- paste0(args, " --cut_right_window_size ", cutSlideWindowSize)
        }
        if (cutSlideWindowQual != 20) {
            args <- paste0(args, " --cut_right_mean_quality ", cutSlideWindowQual)
        }
    }

    if (cutWindowSize != 4) {
        args <- paste0(args, " --cut_window_size ", cutWindowSize)
    }
    if (cutMeanQual != 20) {
        args <- paste0(args, " --cut_mean_quality ", cutMeanQual)
    }

    if (! qualityFiltering) {
        args <- paste0(args, " --disable_quality_filtering")
    }

    if (qualityFilterPhred != 15) {
        args <- paste0(args, " --qualified_quality_phred ", qualityFilterPhred)
    }

    if (qualityFilterPercent != 40) {
        args <- paste0(args, " --unqualified_percent_limit ", qualityFilterPercent)
    }

    if (maxNfilter != 5) {
        args <- paste0(args, " --n_base_limit ", maxNfilter)
    }

    if (averageQualFilter != 0) {
        args <- paste0(args, " --average_qual ", averageQualFilter)
    }

    if (! lengthFiltering) {
        args <- paste0(args, " --disable_length_filtering")
    }

    if (minReadLength != 15) {
        args <- paste0(args, " --length_required ", minReadLength)
    }

    if (maxReadLength > 0) {
        args <- paste0(args, " --length_limit ", maxReadLength)
    } 

    if (lowComplexityFiltering) {
        args <- paste0(args, " --low_complexity_filter")
        if (minComplexity != 30) {
            args <- paste0(args, " --complexity_threshold ", minComplexity)
        }
    }

    if (index1Filter != "") {
        args <- paste0(args, " --filter_by_index1 ", index1Filter)
        if (maxIndexMismatch != 0) {
            args <- paste0(args, " --filter_by_index_threshold ", maxIndexMismatch)
        }
    }

    if (index2Filter != "") {
        args <- paste0(args, " --filter_by_index2 ", index2Filter)
    } 

    if (correctionOverlap) {
        args <- paste0(args, " --correction")
    }

    if (minOverlapLength != 30) {
        args <- paste0(args, " --overlap_len_require ", minOverlapLength)
    }

    if (maxOverlapMismatch != 5) {
        args <- paste0(args, " --overlap_diff_limit ", maxOverlapMismatch)
    }

    if (maxOverlapMismatchPercentage != 20) {
        args <- paste0(args, " --overlap_diff_percent_limit ", maxOverlapMismatchPercentage)
    }

    if (umi) {
        args <- paste0(args, " --umi --umi_loc ", umiLoc)
        if (umiLoc %in% c("read1", "read2")) {
            if (umiLength > 0) {
                args <- paste0(args, " --umi_len ", umiLength)
            }
            else {
                stop("You must set length of UMI if the UMI is in read1/read2.")
            }

            if (umiSkipBaseLength > 0) {
                args <- paste0(args, " --umi_skip ", umiSkipBaseLength)
            }
        }

        if (umiPrefix != "") {
            args <- paste0(args, " --umi_prefix ", umiPrefix)
        }
    }

    if (overrepresentationAnalysis) {
        args <- paste0(args, " -p")
        if (overrepresentationSampling != 20) {
            args <- paste0(args, " -P ", overrepresentationSampling)
        }
    }

    if (splitOutput > 0) {
        args <- paste0(args, " --split ", splitOutput)
    }
    else if (splitByLines >0) {
        args <- paste0(args, " --split_by_lines ", splitByLines)
    }

    if (splitPrefixPaddingNum != 4) {
        args <- paste0(args, " --split_prefix_digits ", splitPrefixPaddingNum)
    }

    if (verbose) {
        args <- paste0(args, " -V")
    }

    args <- paste0(args, ' --report_title "Rfastp Report" -w ', thread, " -h ", outputFastq, ".html -j ", outputFastq, ".json" )
    call <- paste(shQuote(file.path(system.file(package="Rfastp"), "fastp")), args)    
    system(call, intern=TRUE)
    return(fromJSON(file = paste0(outputFastq, ".json")))
}
