
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
#' @param splitOutput number of files to be splitted (2~999). a sequential number
#'      prefix will be added to output name. Default is 0 (no split)
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


rfastp <- function(read1="", read2="", outputFastq="", unpaired="", failedOut="",
    merge=FALSE, mergeOut="", phred64=FALSE, interleaved=FALSE,
    fixMGIid=FALSE, adapterTrimming=TRUE, adapterSequenceRead1="auto",
    adapterSequenceRead2="auto", adapterFasta="", trimFrontRead1=0,
    trimTailRead1=0, trimFrontRead2=0, trimTailRead2=0, maxLengthRead1=0,
    maxLengthRead2=0, forceTrimPolyG=FALSE, disableTrimPolyG=FALSE,
    minLengthPolyG=10, trimPolyX=FALSE, minLengthPolyX=10,
    cutLowQualFront=FALSE, cutLowQualTail=FALSE, cutSlideWindowRight=FALSE,
    cutWindowSize=4, cutMeanQual=20, cutFrontWindowSize=4, cutFrontMeanQual=20,
    cutTailWindowSize=4, cutTailMeanQual=20, cutSlideWindowSize=4,
    cutSlideWindowQual=20, qualityFiltering=TRUE, qualityFilterPhred=15,
    qualityFilterPercent=40, maxNfilter=5, averageQualFilter=0,
    lengthFiltering=TRUE, minReadLength=15, maxReadLength=0,
    lowComplexityFiltering=FALSE, minComplexity=30, index1Filter="",
    index2Filter="", maxIndexMismatch=0, correctionOverlap=FALSE,
    minOverlapLength=30, maxOverlapMismatch=5,
    maxOverlapMismatchPercentage=20, umi=FALSE, umiLoc="",
    umiLength=0, umiPrefix="", umiSkipBaseLength=0,
    umiNoConnection=FALSE, umiIgnoreSeqNameSpace=FALSE,
    overrepresentationAnalysis=FALSE, overrepresentationSampling=20,
    splitOutput=0, splitByLines=0, thread=2, verbose=TRUE) {

    if (umi & umiPrefix != "" & !umiNoConnection) {
        umiPrefix <- paste0(umiPrefix, "_")
    }

    exitcode <- runFastp(read1=read1, read2=read2, outputFastq=outputFastq,
        unpaired=unpaired, failedOut=failedOut, merge=merge, mergeOut=mergeOut, phred64=phred64,
        interleaved=interleaved, fixMGIid=fixMGIid, adapterTrimming=adapterTrimming,
        adapterSequenceRead1=adapterSequenceRead1, adapterSequenceRead2=adapterSequenceRead2,
        adapterFasta=adapterFasta, trimFrontRead1=trimFrontRead1, trimTailRead1=trimTailRead1,
        trimFrontRead2=trimFrontRead2, trimTailRead2=trimTailRead2, maxLengthRead1=maxLengthRead1,
        maxLengthRead2=maxLengthRead2, forceTrimPolyG=forceTrimPolyG, disableTrimPolyG=disableTrimPolyG,
        minLengthPolyG=minLengthPolyG, trimPolyX=trimPolyX, minLengthPolyX=minLengthPolyX,
        cutLowQualFront=cutLowQualFront, cutLowQualTail=cutLowQualTail, cutSlideWindowRight=cutSlideWindowRight,
        cutWindowSize=cutWindowSize, cutMeanQual=cutMeanQual, cutFrontWindowSize=cutFrontWindowSize,
        cutFrontMeanQual=cutFrontMeanQual, cutTailWindowSize=cutTailWindowSize, cutTailMeanQual=cutTailMeanQual,
        cutSlideWindowSize=cutSlideWindowSize, cutSlideWindowQual=cutSlideWindowQual,
        qualityFiltering=qualityFiltering, qualityFilterPhred=qualityFilterPhred,
        qualityFilterPercent=qualityFilterPercent, maxNfilter=maxNfilter, averageQualFilter=averageQualFilter,
        lengthFiltering=lengthFiltering, minReadLength=minReadLength, maxReadLength=maxReadLength,
        lowComplexityFiltering=lowComplexityFiltering, minComplexity=minComplexity,
        index1Filter=index1Filter, index2Filter=index2Filter, maxIndexMismatch=maxIndexMismatch,
        correctionOverlap=correctionOverlap, minOverlapLength=minOverlapLength, maxOverlapMismatch=maxOverlapMismatch,
        maxOverlapMismatchPercentage=maxOverlapMismatchPercentage, umi=umi, umiLoc=umiLoc,
        umiLength=umiLength, umiPrefix=umiPrefix, umiSkipBaseLength=umiSkipBaseLength,
        umiIgnoreSeqNameSpace=umiIgnoreSeqNameSpace,
        overrepresentationAnalysis=overrepresentationAnalysis,
        overrepresentationSampling=overrepresentationSampling,
        splitOutput=splitOutput, splitByLines=splitByLines, thread=thread, verbose=verbose)
    return(fromJSON(file = paste0(outputFastq, ".json")))
}




#' R wrap of gencore
#'
#' Quality control (Cut adapter, low quality trimming, UMI handling, and etc.) of fastq files.
#'
#' @param inBam read1 input file name of sorted bam/sam [string]
#' @param outBam file name of output bam/sam [string]
#' @param refFile Path to reference fasta file. [string]
#' @param bedFile bedfile to specify the capturing region. [string]
#' @param umiPrefix the prefix for UMI, if it has. None by default.
#' @param numSupportingReads only output consensus reads/pairs that
#'     merged by >= `supportingRead` reads/pairs. The value should be 1~10, and the
#'     default value is 1.
#' @param majorBaseRatio if the ratio of the major base in a cluster is less than
#'     `ratioMajorBase`, it will be further compared to the reference. The value should
#'     be 0.5~1.0, and the default value is 0.8
#' @param majorBaseScore if the score of the major base in a cluster is less than
#'     `scoreMajorBase`, it will be further compared to the reference. The value should
#'     be 1~20, and the default value is 6
#' @param highQual the threshold for a quality score to be considered as high qulity.
#'     Default 30 means Q30
#' @param moderateQual the threshold for a quality score to be considered as moderate
#'     qulity. Default 20 means Q20.
#' @param lowQual the threshold for a quality score to be considered as low qulity.
#'     Default 15 means Q15.
#' @param debug a logical indicating output some debug information to STDERR.
#' @param coverageSampling the sampling rate for genome scale coverage statistics.
#'     Default 10000 means 1/10000.
#' @param quitAfterContig stop when `quitAfterContig` contigs are processed.
#'     Only used for fast debugging. Default 0 means no limitation.
#' @param verbose output verbose log information
#'
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
#' inputbamfile <-  system.file("extdata", "ex1_sorted.bam", package="Rfastp")
#' outputbamfile <- "ex1_rgencore.bam"
#' reference <- system.file("extdata", "myreference.fa", package="Rfastp")
#'
#' # run rgencore
#'
#' rgencore_json_report <- rgencore(inBam=inputbamfile,
#'          outBam=outputbamfile,
#'          ref=reference
#')

rgencore <- function(inBam="", outBam="", refFile="", bedFile="", umiPrefix="",
    numSupportingReads=1, majorBaseScore=6, majorBaseRatio=0.8, quitAfterContig=0,
    highQual=30, moderateQual=20, lowQual=15, coverageSampling=10000,
    debug=FALSE, verbose=TRUE ) {
    exitcode <- runGencore(inBam=inBam, outBam=outBam, refFile=refFile, bedFile=bedFile, umiPrefix = umiPrefix,
        numSupportingReads = numSupportingReads, majorBaseScore = majorBaseScore, majorBaseRatio = majorBaseRatio,
        quitAfterContig = quitAfterContig, highQual = highQual, moderateQual = moderateQual, lowQual = lowQual,
        coverageSampling = coverageSampling, debug = debug, verbose = verbose)
    return(fromJSON(file = paste0(outBam, ".json")))
}

