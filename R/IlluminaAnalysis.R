#' @title  Analysis of Illumina Files
#'
#' @description Analyses the Illumina Files present within the provided folder
#’
#' @param input_path a path to a folder containing fastq.gz or fastq files
#'        the default is set to "./files"
#' @param trunc_forw An integer
#'        determines the length at which forward reads are truncated
#'        the default is set to 270
#' @param trunc_rev An integer
#'        determines the length at which reverse reads are truncated
#'        the default is set to 220
#' @param seq_length An integer range
#'        determines the length the aligned forward an reverse reads are allowed
#'        to have - sequences that are not within the range will be discarded
#'        default is set to 400:460
#' @param EEforward An integer value
#'        determines the amount of expected errors within each forward sequence
#'        all sequences with more errors will be discarded
#'        default is set to 2
#' @param EEreverse An integer value
#'        determines the amount of expected errors within each reverse sequence
#'        all sequences with more errors will be discarded
#'        default is set to 2
#' @param tax_trainer a PATH to the file that will be used as a trainer for the
#'        taxonomic determination of the sequences
#'        the default is set to:
#'        ".~/Gregor/Project Biofilms/NGS-analysis/files/tax/
#'        silva_nr_v138_train_set.fa.gz"
#' @param resort a logical value
#'        if true the names of the samples will be resorted
#'        (experimental feature - might break the algorithm)
#'        works only if the samplenames are of the sort
#'        LettersNumbers
#' @param multithreading a logical value
#'        is parsed to the multithread option of the underlying dada functions
#'        will enable multi-processor core usage - does not work on windows OS
#' @param skipTree a logical value
#'        if TRUE the phylogenetic TREE will not be calculated and a marker will be
#'        created to annouce this within the displayApp - used as an option to save
#'        time as processing time increases exponetially with sample number
#'
#' @import cooccur
#' @import dada2
#' @import phyloseq
#' @import phangorn
#' @import filesstrings
#' @import xlsx
#' @return
#'
#' @examples
#'
#' @export


illumina_analysis <- function(input_path = "./files",
                              # input_path is the path to a fastq file folder
                              # the default is set to "./files"
                              trunc_forw = 270,
                              # trunc_forw is integer to determine length at
                              # which forward reads are to be trunctuated
                              # default is set to 270
                              trim_left_forw = 10,
                              # trim_left_forw is an integer that determines how
                              # many bases at the beginning of the forward sequence
                              # are trimmed
                              # default is set to 10
                              trunc_rev = 220,
                              # trunc_rev is integer to determine length at
                              # which reverse reads are to be trunctuated
                              # default is set to 220
                              trim_left_rev = 10,
                              # trim_left_rev is an integer that determines how
                              # many bases at the beginning of the reverse sequence
                              # are trimmed
                              # default is set to 10
                              seq_length = 400:460,
                              # seq_length is range of the expected length of
                              # aligned forward and reverse reads
                              # default is set to range from 400 to 460
                              EEforward = 2,
                              # EEforward is an integer for the accepted
                              # expected forward errors
                              # the default is set to 2
                              EEreverse = 2,
                              # EEreverse is an integer for the accepted
                              # expected reverse errors
                              # the default is set to 2
                              tax_trainer = "~/Gregor/Project Biofilms
                              /NGS-analysis/files/tax
                              /silva_nr_v138_train_set.fa.gz",
                              # tax_trainer is the path to the training fasta
                              # file that is used for taxonomy identification
                              # compressed files can be used see Dada2
                              # documentation on how such a file should be
                              resort = TRUE,
                              # takes a logical value
                              # if true the sampleNames will be resorted
                              multithreading = TRUE,
                              # logical value
                              # if TRUE multiple processor cores will be used
                              # does not work on windows OS
                              skipTree = FALSE
                              # logical value
                              # if TRUE the calculations for the phylogenetic tree
                              # will be skipped in order to save time
){
  # list as variables - this is a necessity for use within shiny
  lsorig <- ls()
  userinput <- c(input_path = input_path,
             trim_left_forw = trim_left_forw,
             trunc_forw = trunc_forw,
             trim_left_rev = trim_left_rev,
             trunc_rev = trunc_rev,
             seq_length = seq_length,
             EEforward = EEforward,
             EEreverse = EEreverse,
             tax_trainer = tax_trainer,
             resort = resort,
             multithreading = multithreading)
  #----------------------------------------------------------------------------->Initial directory setup


  # take the time
  message(Sys.time())
  startTime <- Sys.time()

  # access the timelog.csv file if it exists
  # and assigns it to the timeLog variable
  # otherwise creates a timelog variable from scratch
  if (file.exists("./internal_files/timelog.csv")){
    timeLog <- read.csv("./internal_files/timelog.csv")
  } else {
    timeLog <- data.frame(time = integer(),
                          number_of_samples = integer(),
                          ID = integer())
  }
  ID <- nrow(timeLog) + 1


  # create subfolders within the input_path if these folders do not exist yet
  message("Creating folder structure")
  create_folder(input_path = input_path,
                folder_name = "SeqFiles/Originals",
                variable_name = "SFO")
  create_folder(input_path = input_path,
                folder_name = "SeqFiles/Filtered",
                variable_name = "SFF")
  create_folder(input_path = input_path,
                folder_name = "RData",
                variable_name = "RDF")
  create_folder(input_path = input_path,
                folder_name = "Plots",
                variable_name = "plots")
  create_folder(input_path = input_path,
                folder_name = "Tables",
                variable_name = "tabs")


  # creates a logframe dataframe
  logframe <- data.frame(Finished = character(),
                         Date = character(),
                         Time = character(),
                         stringsAsFactors = FALSE)

  # updates the logframe
  logframe <- addLog("Directory setup", logframe)

  # save all variables under BackUp.RData
  lsnow <- ls()
  lsnow <- lsnow[!(lsnow %in% lsorig)]
  save(list = lsnow,
       file = paste0(RDF,
                     "/BackUp.RData"))
  rm(lsnow)
  #----------------------------------------------------------------------------->Filtering and sorting of sequences


  # create a list of all fastq.gz files in input folder
  # full.names = T because the full path is needed for file movement
  filelist <- list.files(path = input_path,
                         pattern = ".fastq.gz",
                         full.names = TRUE)
  message(paste0("You are processing ",
                 length(filelist) / 2,
                 " sample files."))

  if (nrow(timeLog) == 0){
    message("This is your first time executing this App
            - it might take a while.")
  } else {
    message(paste0("Previously you processed and analysed ",
                   sum(timeLog$number_of_samples),
                   " samplefiles within ",
                   sum(timeLog$time),
                   " minutes."))
    message(paste0("That results in an average time of ",
                   sum(timeLog$time) / sum(timeLog$number_of_samples),
                   " minutes per samplefile."))
  }


  # moves the files of the input folder
  # into the created original Sequences sub-folder
  filesstrings::move_files(filelist, SFO)


  # read in the filepath to the fastq files as a sorted list
  # full.names = T is required because the dada2 functions
  # require the full file path
  message("Files are being read.")
  FastqForward <- gtools::mixedsort(list.files(SFO,
                                               pattern = "_R1_001.fastq.gz",
                                               full.names = TRUE))
  FastqReverse <- gtools::mixedsort(list.files(SFO,
                                               pattern = "_R2_001.fastq.gz",
                                               full.names = TRUE))

  # creates a list that contains the elements
  # of each Illumina file name split by the "_"
  name_elements <- strsplit(basename(FastqForward), "_")

  # create a list of the sample names from that by
  # taking the first character of every object in the list
  SampleN <- sapply(name_elements, `[`, 1)

  # resorts the SampleNames by splitting them between Letters and Numbers
  # counts the first occurring amount of numbers per samplename and adds to all
  # of them zeros in front until they are all of the same length
  # combines the letters and the extended numbers into the new samplename
  if (resort){
    SampleLetter <- str_extract(SampleN, "([[:alpha:]])+")
    SampleNumber <- str_extract(SampleN, "([[:digit:]])+")
    SampleNames <- NULL
    for (i in 1:length(SampleN)){
      RunningName <- paste0(
        SampleLetter[as.numeric(i)],
        str_pad(
          SampleNumber[as.numeric(i)],
          nchar(SampleNumber[length(SampleNumber)]),
          "left",
          "0"
        )
      )
      SampleNames<- append(SampleNames, RunningName)

    }
    rm(RunningName)
    rm(SampleLetter)
    rm(SampleN)

  } else {
    SampleNames <- SampleN
    rm(SampleN)
  }

  # create a file list for the filtered and trimmed files that are to be created
  FilterForward <- file.path(SFF,
                             paste0(SampleNames, "_Forw_filt.fastq.gz"))
                            # creates file names in the folder from the
                            # generated sample names an the given string
  FilterReverse <- file.path(SFF,
                             paste0(SampleNames, "_Rev_Filt.fastq.gz"))

  # use the filter and trim function of dada2 to create new fastq files
  # this will truncate the files and filter away files that are to short
  # or below a given quality score  see the dada2 documentation for
  # more information
  message("Files are being filtered.")
  FilterOutput <- dada2::filterAndTrim(FastqForward,
                                       FilterForward,
                                       FastqReverse,
                                       FilterReverse,
                                       trimLeft = c(trim_left_forw,
                                                    trim_left_rev),
                                       # trimming the left side of each
                                       # sequence at the given value - happens
                                       # after truncLen is applied
                                       truncQ = 2,
                                       # truncates the reads at the first
                                       # instance where the q score is equal or
                                       # lower to the given value - happens
                                       # before truncLen
                                       truncLen = c(trunc_forw, trunc_rev),
                                       # truncuating the forward and reverse
                                       # sequence at the given value
                                       # sequences that are shorter will be
                                       # discarded
                                       maxN = 0,
                                       maxEE = c(EEforward,
                                                 EEreverse),
                                       # read with more expected error than
                                       # expected will be discarded - this is a
                                       # better filter than truncQ - see dada2
                                       # documentation
                                       rm.phix = TRUE,
                                       # remove sequences that match to the
                                       # phiX genome that is used as an internal
                                       # control by Illumina
                                       compress = TRUE,
                                       multithread = multithreading,
                                       # is set to TRUE to use multiple
                                       # processor cores does not work on
                                       # Microsoft machines
                                       verbose = TRUE)
  # updating logframe
  logframe <- addLog("Filtering Sequences", logframe)

  # save all variables under BackUp.RData
  lsnow <- ls()
  lsnow <- lsnow[!(lsnow %in% lsorig)]
  save(list = lsnow,
       file = paste0(RDF,
                     "/BackUp.RData"))
  rm(lsnow)

  #----------------------------------------------------------------------------->Learning of Errorrates

  # use the dada2 learnErrors function to determine the error rate of the
  # sequences this is a self learning algorithm and the main advantage of the
  # dada2 package the provided ErrorScore is needed for the dada function
  message("Learning errorrate.")
  ErrorForward <- dada2::learnErrors(FilterForward,
                                     multithread = multithreading)
  ErrorReverse <- dada2::learnErrors(FilterReverse,
                                     multithread = multithreading)

  # apply dada core function to detect sample interference
  DadaForward <- dada2::dada(FilterForward,
                             ErrorForward,
                             multithread = multithreading,
                             pool = TRUE)
  DadaReverse <- dada2::dada(FilterReverse,
                             ErrorReverse,
                             multithread = multithreading,
                             pool = TRUE)

  # merge the forward and the matching reverse strands of the processed data
  MergerDada <- dada2::mergePairs(DadaForward,
                                  FilterForward,
                                  DadaReverse,
                                  FilterReverse,
                                  verbose = TRUE)

  # create a sequence table of the merged pairs
  SeqTab <- dada2::makeSequenceTable(MergerDada)
  rownames(SeqTab) <- SampleNames

  # create a new sequence table on the basis of the current one
  # including only sequences of the specified length range
  SeqTabRange <- SeqTab[, nchar(colnames(SeqTab)) %in% seq_length]


  # remove identifiable chimeras
  NonChime <- dada2::removeBimeraDenovo(SeqTabRange,
                                        method = "consensus",
                                        multithread = multithreading,
                                        verbose = TRUE)

  # updating logframe
  logframe <- addLog("Learning Errorrates", logframe)

  # save all variables under BackUp.RData
  lsnow <- ls()
  lsnow <- lsnow[!(lsnow %in% lsorig)]
  save(list = lsnow,
       file = paste0(RDF,
                     "/BackUp.RData"))
  rm(lsnow)
  #----------------------------------------------------------------------------->Assigning taxonomy

  # creates a csv files containg the sequence table with the chimeras removed
  message("Saving sequence tables.")
  write.csv(as.data.frame(t(NonChime)),
            paste0(tabs,
                   "/Sequences.csv"))

  # assigning taxonomy to the sequences
  # use the assignTaxonomy function of dada to create a table with the assigned
  # taxonomic data this uses a fa.gz trainer file provided in the
  # function input under tax_trainer
  message("Assigning Taxonomy")
  taxonomy <- dada2::assignTaxonomy(NonChime,
                                    tax_trainer,
                                    multithread = multithreading)

  message("Saving intial taxonomic data.")
  write.csv(taxonomy,
            paste0(tabs, "/Taxa.csv"))

  # updating logframe
  logframe <- addLog("Assigning taxonomy", logframe)

  # save all variables under BackUp.RData
  lsnow <- ls()
  lsnow <- lsnow[!(lsnow %in% lsorig)]
  save(list = lsnow,
       file = paste0(RDF,
                     "/BackUp.RData"))
  rm(lsnow)
  #----------------------------------------------------------------------------->Creating abundance & Co-occurrence tables

  # create an abundance table with the taxonomic ranks
  taxaAbund <- merge(taxonomy,
                     as.data.frame(t(NonChime)),
                     by = "row.names")
  row.names(taxaAbund) <- paste0("ASV",
                                 seq(nrow(taxaAbund)))
  taxaAbund <- dplyr::select(taxaAbund, -Row.names)

  # create an object containing all existing taxonomic ranks
  taxRank <- names(dplyr::select_if(taxaAbund, is.character))

  # while loop that uses the taxRank and the taxAbund object and combines all
  # elements that are equal on the level of the specific taxonomic rank
  # the loop itself can be adjusted for all tables containing a set of pure
  # character string columns and another set of pure numeric columns
  # creates a S3 object for each taxonomic rank containing the combined
  # abundance table as well a a tidy version, an Occurance table and the
  # Co_occurance data
  while(length(taxRank) > 0)
  {
    # save the amount of entries within the taxRank object as an integer within the
    # tRlength object
    tRlength <- length(taxRank)

    temps3 <- NULL

    tempdf <- plyr::ddply(taxaAbund, taxRank, plyr::numcolwise(sum))
    # create empty two column dataframe for the Names and Long_Names
    name <- data.frame(Name = as.character(),
                       Long_Name = as.character())
    # start a for loop to iterate through the rows within the tempdf dataframe object
    # it iterates starting from 1 with the variable dfrow
    for (dfrow in 1:nrow(tempdf)) {
      # if clause to detect if current row (dfrow) of the first column (1) of the
      # name dataframe object is NA
      if (is.na(name[dfrow, 1]))
      {
        # if clause to detect if the current row (dfrow) of the current column
        # (tRlength) of the tempdf dataframe object is NA
        # the check for the NA needs to be the first check as the string detection
        # that is used in the later statements would return NA in the case of an NA
        # entry and not TRUE/FALSE
        if (is.na(tempdf[dfrow, tRlength]))
        {
          # start of for loop to iterate from 0 to penultimate value of tRlength
          # using the TR variable
          for (TR in 0:(tRlength - 1))
          {
            # start while loop under the condition that the current row (dfrow) of the
            # first column (1) of the name dataframe object is NA
            if (is.na(name[dfrow, 1]))
            {
              # if clause checking the current row (dfrow) of the columns that come
              # before the current column (tRlength - TR) for their content - checks
              # if the content is neither NA nor contains "unknown" nor "incertae"
              if (!is.na(tempdf[dfrow, tRlength - TR]) &
                  !stringr::str_detect(tempdf[dfrow, tRlength - TR],
                                       stringr::fixed("unknown",
                                                      ignore_case = TRUE)) &
                  !stringr::str_detect(tempdf[dfrow, tRlength - TR],
                                       stringr::fixed("incertae",
                                                      ignore_case = TRUE)) )
              {
                # if the above stated if clauses are all TRUE the current row (dfrow)
                # of the first column (1) of the name dataframe object will be
                # assigned a string stating that the current taxonomic rank is NA
                # and what taxonomic rank is known
                name[dfrow, 1] <- paste0("NA_",
                                         taxRank[tRlength],
                                         "_(",
                                         taxRank[tRlength - TR],
                                         " = ",
                                         tempdf[dfrow, tRlength - TR],
                                         ")")
              }
              # if clause that only gets checked if at least one of the conditions of
              # the previous if clause resulted in a FALSE - this if clause has the
              # the condition is that tRlength - TR == 1 meaning which results in all
              # columns of the tempdf dataframe being checked and none containing an
              # entry that does not have the "unknown" or the "incertae" string within
              # them
              else if
              (
                tRlength - TR == 1
              )
              {
                # if the above if clause is fulfilled and thus no column of the current
                # row (dfrow) of the tempdf dataframe object contains a viable entry
                # then the first column of the current row (dfrow) of the name
                # dataframe object will be assigned a string containing that the
                # current taxonimic rank as well as the highest taxonomic rank are
                # NA
                name[dfrow, 1] <- paste0("NA_",
                                         taxRank[tRlength],
                                         "_(",
                                         taxRank[tRlength - TR],
                                         " = NA)")
              } # end of the if and else if statements checking the columns of the
              # tempdf object for their contents
            }# end of the while loop conditioned on the fact that current row (dfrow)
            # of the first column (1 - the name column) of the name dataframe object
            # is NA
          }# end of for loop iterating with TR from 0 to tRlength
        }# end of if to detect NA in current row of tempdf

        # if clause to detect if the current row (dfrow) of the current column
        # (tRlength) of the tempdf dataframe object contains the "unknown" string
        else if
        (
          stringr::str_detect(tempdf[dfrow, tRlength],
                              stringr::fixed("unknown",
                                             ignore_case = TRUE))
        )
        {
          # start of for loop to iterate from 0 to penultimate value of tRlength
          # using the TR variable
          for (TR in 0:(tRlength - 1))
          {
            # start while loop under the condition that the current row (dfrow) of the
            # first column (1) of the name dataframe object is NA
            if (is.na(name[dfrow, 1]))
            {
              # if clause checking the current row (dfrow) of the columns that come
              # before the current column (tRlength - TR) for their content - checks
              # if the content is neither NA nor contains "unknown" nor "incertae"
              if (!is.na(tempdf[dfrow, tRlength - TR]) &
                  !stringr::str_detect(tempdf[dfrow, tRlength - TR],
                                       stringr::fixed("unknown",
                                                      ignore_case = TRUE)) &
                  !stringr::str_detect(tempdf[dfrow, tRlength - TR],
                                       stringr::fixed("incertae",
                                                      ignore_case = TRUE)) )
              {
                # if the above stated if clauses are all TRUE the current row (dfrow)
                # of the first column (1) of the name dataframe object will be
                # assigned a string stating that the current taxonomic rank is unknown
                # and what taxonomic rank is known
                name[dfrow, 1] <- paste0("Unknown_",
                                         taxRank[tRlength],
                                         "_(",
                                         taxRank[tRlength - TR],
                                         " = ",
                                         tempdf[dfrow, tRlength - TR],
                                         ")")
              }
              # if clause that only gets checked if at least one of the conditions of
              # the previous if clause resulted in a FALSE - this if clause has the
              # the condition is that tRlength - TR == 1 meaning which results in all
              # columns of the tempdf dataframe being checked and none containing an
              # entry that does not have the "unknown" or the "incertae" string within
              # them
              else if
              (
                tRlength - TR == 1
              )
              {
                # if the above if clause is fulfilled and thus no column of the current
                # row (dfrow) of the tempdf dataframe object contains a viable entry
                # then the first column of the current row (dfrow) of the name
                # dataframe object will be assigned a string containing that the
                # current taxonimic rank as well as the highest taxonimic rank are
                # unknown
                name[dfrow, 1] <- paste0("Unknown_",
                                         taxRank[tRlength],
                                         "_(",
                                         taxRank[tRlength - TR],
                                         " = Unknown)")
              } # end of the if and else if statements checking the columns of the
              # tempdf object for their contents
            }# end of the while loop conditioned on the fact that current row (dfrow)
            # of the first column (1 - the name column) of the name dataframe object
            # is NA
          }# end of for loop iterating with TR from 0 to tRlength
        }# end of if to detect "unknown" string

        # if clause that gets checked only if the above if statement checking for the
        # unknown string returns a false
        # this if clause checks if the current row (dfrow) of the current column
        # (tRlenght) of the tempdf dataframe object contains the "incertae" string
        else if
        (
          stringr::str_detect(tempdf[dfrow, tRlength],
                              stringr::fixed("incertae",
                                             ignore_case = TRUE))
        )
        {
          # start of for loop to iterate from 0 to penultimate value of tRlength
          # using the TR variable
          for (TR in 0:(tRlength - 1))
          {
            # start while loop under the condition that the current row (dfrow) of the
            # first column (1) of the name dataframe object is NA
            if (is.na(name[dfrow, 1]))
            {
              # if clause checking the current row (dfrow) of the columns that come
              # before the current column (tRlength - TR) for their content - checks
              # if the content is neither NA nor contains "unknown" nor "incertae"
              if (!is.na(tempdf[dfrow, tRlength - TR]) &
                  !stringr::str_detect(tempdf[dfrow, tRlength - TR],
                                       stringr::fixed("unknown",
                                                      ignore_case = TRUE)) &
                  !stringr::str_detect(tempdf[dfrow, tRlength - TR],
                                       stringr::fixed("incertae",
                                                      ignore_case = TRUE)) )
              {
                # if the above stated if clauses are all TRUE the current row (dfrow)
                # of the first column (1) of the name dataframe object will be
                # assigned a string stating that the current taxonomic rank is
                # incertae and what taxonomic rank is known (certain)
                name[dfrow, 1] <- paste0("incertae_",
                                         taxRank[tRlength],
                                         "_(",
                                         taxRank[tRlength - TR],
                                         " = ",
                                         tempdf[dfrow, tRlength - TR],
                                         ")")
              }
              # if clause that only gets checked if at least one of the conditions of
              # the previous if clause resulted in a FALSE - this if clause has the
              # the condition is that tRlength - TR == 1 meaning which results in all
              # columns of the tempdf dataframe being checked and none containing an
              # entry that does not have the "unknown" or the "incertae" string within
              # them
              else if
              (
                tRlength - TR == 1
              )
              {
                # it the above if clause is fulfilled and thus no column of the current
                # row (dfrow) of the tempdf dataframe object contains a viable entry
                # then the first column of the current row (dfrow) of the name
                # dataframe object will be assigned a string containing that the
                # current taxonimic rank as well as the highest taxonimic rank are
                # unknown/incertae/uncertain
                name[dfrow, 1] <- paste0("incertae_",
                                         taxRank[tRlength],
                                         "_(",
                                         taxRank[tRlength - TR],
                                         " = Unknown)")
              } # end of the if and else if statements checking the columns of the
              # tempdf object for their contents
            } # end of the while loop conditioned on the fact that current row (dfrow)
            # of the first column (1 - the name column) of the name dataframe object
            # is NA -
          }# end of for loop iterating with TR from 0 to tRlength
        }# end of if to detect "incertae" string

        # else statement that gets executed if all the previous if statements of this
        # magnitude returned FALSE (meaning that current taxonomic rank (tRlength) of
        # the current row (dfrow) of the tempdf dataframe object is available and does
        # neither contains the "unknown" nor the "incertae string")
        else
        {
          # assign the contents of the current row of the current taxonomic rank of the
          # tempdf dataframe object as the contents of the current row of the first
          # column of the name dataframe object
          name[dfrow, 1] <- tempdf[dfrow, tRlength]
        } # end of if and else statements checking the content of current row of the
        # tempdf dataframe object
      }# end of if the current row in name object is NA

      # assignment of the Long_Name column within the name object:

      # iterates with the variable LongNameTaxRank from 1 to the length of the current
      # taxRank (tRlength)
      for (LongNameTaxRank in 1:tRlength) {
        # check if the current row (dfrow) has an entry in the Long_name column
        # (col 2) is NA
        if (is.na(name[dfrow, 2]))
        {
          # if statement checking the columns of current row (dfrow) of the tempdf
          # dataframe object - starting with the first column (utilizing the LongNameTaxRank
          # variable) - checking the column entries for NOT being NA
          if (!is.na(tempdf[dfrow, LongNameTaxRank]))
          {
            # assign the contents of the tempdf column and row as the contents of the
            # name row and column 2
            name[dfrow, 2] <- tempdf[dfrow, LongNameTaxRank]
          }
          # else statement - will only be executed if the above if statement returned
          # FALSE
          else
          {
            # assign "NA_(current taxonomic rank[LongNameTaxRank])" as the content of the current
            # row of the secondcolumn of the name dataframe object
            name[dfrow, 2] <- paste0("NA_(",
                                     taxRank[LongNameTaxRank],
                                     ")")
          } # end of else/if statements checking the contents of the current row
          # (dfrow) of the tempdf dataframe object

        } # end of if checking that the current row of the second column of the name
        # dataframe object is NA
        # else statement - will only be executed if the above if statement returned
        # FALSE
        else
        {
          # if statement checking the columns of current row (dfrow) of the tempdf
          # dataframe object - starting with the first column (utilizing the LongNameTaxRank
          # variable) - checking the column entries for being NA
          if (is.na(tempdf[dfrow, LongNameTaxRank]))
          {
            # append "_NA" to the entry of the current row of the second column of the
            # name dataframe object
            name[dfrow, 2] <- paste0(name[dfrow, 2], "_NA")
          }
          # else statement - will only be executed if the above if statement returned
          # FALSE
          else
          {
            # append the entry of the current row (dfrow) of the current column (LongNameTaxRank)
            # of the tempdf dataframe object to the entry of the current row (dfrow)
            # of the second column of the name dataframe object
            name[dfrow, 2] <- paste0(name[dfrow, 2], "_", tempdf[dfrow, LongNameTaxRank])
          }
        } # and of else/if statements checking if the current row (dfrow) of the
        # second column of the name dataframe object is NA

      } # end of for loop iterating from 1 to the length of the taxRank object


    } # end of for loop iterating through the amount of rows of the tempdf object
    # using the dfrow variable

    # bind the created elements within the name object onto the temporary df
    tempdf <- cbind(tempdf,
                    name)
    # rearranges the temporaray df to display all character columns before the
    # numeric columns
    tempdf <- tempdf[, c(names(dplyr::select_if(tempdf, is.character)),
                         names(dplyr::select_if(tempdf, is.numeric)))]

    # detect and assign duplicates within the Name column
    duplicate_names <- tempdf$Name[duplicated(tempdf$Name)]

    # if there are duplicates the first occasion of it will be appended an asteriks
    while(length(duplicate_names) > 0) {
      for (i in 1:nrow(tempdf)){
        if (tempdf$Name[i] %in% duplicate_names) {
          tempdf$Name[i] <- paste0(tempdf$Name[i], "*")
          duplicate_names <- tempdf$Name[duplicated(tempdf$Name)]
          break
        }

      }
    }

    # assign the element withing the Name column as the respective row name
    row.names(tempdf) <- tempdf$Name

    # create column that contains for each row the sum of all numeric elements
    tempdf$Total <- rowSums(dplyr::select_if(tempdf, is.numeric))

    # create a secondary temporary dataframe for the relative abundances
    temprela <- tempdf[, names(dplyr::select_if(tempdf, is.numeric))]
    for (i in 1:ncol(temprela)) {
      temprela[,i] <- temprela[,i] * 100 / sum(temprela[,i])
    }
    temprela <- cbind(tempdf[names(dplyr::select_if(tempdf, is.character))],
                      temprela)

    # create element contain the the rownames sorted by Abundance for highest
    # to lowest
    tNames <- row.names(tempdf[base::order(tempdf$Total, decreasing = TRUE),])

    # assign the dataframe and the sorted name element to a temporary S3 object
    temps3$Abundance <- tempdf
    temps3$relativeAbundance <- temprela
    temps3$Sorted_Names <- tNames


    # calculate cooccurance for all taxonomic ranks below the highest rank
    if (tRlength > 1)
    {
      # create a table containing just the information if the current taxonomic
      # rank exists within the samples
      tempOcc <- dplyr::select_if(dplyr::select(tempdf,
                                                -Total),is.numeric)
      tempOcc[tempOcc > 0] <- 1
      tempOcc <- cbind(dplyr::select(tempdf, where(is.character)), tempOcc)

      # assign the created dataframe as the occurance subobject within the temporary
      # S3 object (temps3)
      temps3$Occurance <- tempOcc

      message(paste0("Calculating Co-Occurance of ",
                     taxRank[tRlength]))
      # calculates Co-Occurance using the Cooccur package
      temps3$Co_Occurance <- cooccur::cooccur(dplyr::select(tempOcc,
                                                            where(is.numeric)),
                                              spp_names = TRUE)

      rm(tempOcc)
    }



    # create Tidy table with the Samples as names
    temps3$Tidy <- tidyr::pivot_longer(tempdf,
                                       c(colnames(dplyr::select_if(tempdf,
                                                                   where(
                                                                     is.numeric)
                                       )
                                       )),
                                       names_to = "Sample",
                                       values_to = "Abundance")

    # assign the temporary S3 object to the current taxonomic rank
    assign(taxRank[tRlength], temps3)

    # save an xls table with the abundances
    xlsx::write.xlsx(temps3$Abundance,
                     paste0(tabs,
                            "/Abundance.xls"),
                     sheetName = taxRank[tRlength],
                     append = TRUE)
    xlsx::write.xlsx(temps3$relativeAbundance,
                     paste0(tabs,
                            "/Abundance.xls"),
                     sheetName = paste0("relativ",
                                        taxRank[tRlength]),
                     append = TRUE)

    # move to the next higher taxonomic rank
    taxRank <- taxRank[0:(tRlength - 1)]

    # remove temporary files
    rm(tempdf)
    rm(temps3)
    rm(temprela)
  }

  # updating logframe
  logframe <- addLog("Abundance and Co-occrence tables", logframe)

  # save all variables under BackUp.RData
  lsnow <- ls()
  lsnow <- lsnow[!(lsnow %in% lsorig)]
  save(list = lsnow,
       file = paste0(RDF,
                     "/BackUp.RData"))
  rm(lsnow)
  #----------------------------------------------------------------------------->Generate distance matrix and phylogenetic tree
  message("Generating distance matrix")

  # create a dataframe with the sample names for each sample file
  SampleNames <- rownames(NonChime)
  SampleDataFrame <- data.frame(Samplename = SampleNames)
  rownames(SampleDataFrame) <- SampleNames

  # create the phylogenetic tree if the option to skip it was not chosen and create
  # a value denoting the chosen option for further use downstream
  if (!skipTree) {
    # create boolean Value to note if phylogenetic tree was created
    skipedTree <- FALSE
    # extract sequences
    message("Creating a distance matrix of the sequences -
          this process can take several minutes.")
    sequences <- dada2::getSequences(NonChime)
    # gives the sequences characters names that will be used as tip labels
    names(sequences) <- sequences
    # aligning sequences and creating distance matrix
    alignment <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(sequences),
                                     anchor = NA)
    fangoAlign <- phangorn::phyDat(as(alignment,
                                      "matrix"),
                                   type = "DNA")
    # compute pairwise distances from sequences
    fangoDist <- phangorn::dist.ml(fangoAlign)
    # neighbor joining tree estimation
    NJtree <- phangorn::NJ(fangoDist)
    # compute likelihood of phylogenetic tree
    fit <- phangorn::pml(NJtree,
                         data = fangoAlign)
    fit <- stats::update(fit,
                         k = 4,
                         inv = 0.2)
    fitGTR <- phangorn::optim.pml(fit,
                                  model = "GTR",
                                  optInv = TRUE,
                                  optGamma = TRUE,
                                  rearrangement = "stochastic",
                                  control = phangorn::pml.control(trace = 0)
    )

    # updating logframe
    logframe <- addLog("Creating distance matrix and phylogenetic tree",
                       logframe)
  } else {
    skipedTree <- TRUE
  }



  # save all variables under BackUp.RData
  lsnow <- ls()
  lsnow <- lsnow[!(lsnow %in% lsorig)]
  save(list = lsnow,
       file = paste0(RDF,
                     "/BackUp.RData"))
  rm(lsnow)
  #----------------------------------------------------------------------------->Phyloseq pipeline and Clean Up
  message("Creating phyloseq element")
  # use the phyloseq function with the created ASV table as the input OTU table
  if (!skipTree)
  {
    phylo <- phyloseq::phyloseq(phyloseq::otu_table(NonChime,
                                                    taxa_are_rows = FALSE),
                                phyloseq::sample_data(SampleDataFrame),
                                phyloseq::tax_table(taxonomy),
                                phyloseq::phy_tree(fitGTR$tree))
  } else
  {
    phylo <- phyloseq::phyloseq(phyloseq::otu_table(NonChime,
                                                    taxa_are_rows = FALSE),
                                phyloseq::sample_data(SampleDataFrame),
                                phyloseq::tax_table(taxonomy))
  }

  # store the DNA sequence information of the taxa on the seqref() portion
  # of the phyloseq element
  message("Storing sequence information.")
  dna <- Biostrings::DNAStringSet(phyloseq::taxa_names(phylo))
  names(dna) <- phyloseq::taxa_names(phylo)
  CompletePhyl <- phyloseq::merge_phyloseq(phylo,
                                           dna)
  phyloseq::taxa_names(CompletePhyl) <- paste0("ASV",
                                               seq(phyloseq::ntaxa(
                                                 CompletePhyl)))


  # create a table listing the amount of read ins after each processing step
  TrackAllTab <- cbind(FilterOutput,
                       # apply the getN function to the generated data sets
                       sapply(DadaForward,
                              getN),
                       sapply(DadaReverse,
                              getN),
                       sapply(MergerDada,
                              getN),
                       rowSums(NonChime))
  colnames(TrackAllTab) <- c("input seq.",
                             "filtered seq.",
                             "denoised forward seq.",
                             "denoised reverse seq.",
                             "merged seq.",
                             "NonChimera seq.")
  rownames(TrackAllTab) <- SampleNames
  # save the created table as a csv
  write.csv(TrackAllTab,
            paste0(tabs,
                   "/TrackSequenceLoss.csv"))

  #save the session info as an .RData file
  message("Saving session information.")
  SInfo <- sessionInfo()
  save(SInfo, file = paste0(RDF, "/SessionInfo.RData"))

  # updating logframe
  logframe <- addLog("Phyloseq and Clean Up", logframe)

  savelist <- names(dplyr::select_if(taxaAbund, is.character))
  savelist <- c(savelist[1:length(savelist)],
                "CompletePhyl",
                "userinput",
                "skipedTree")
  save(list = (savelist),
       file = paste0(RDF,
                     "/IlluminaAnalysis.RData"))

  # save all variables under BackUp.RData
  lsnow <- ls()
  lsnow <- lsnow[!(lsnow %in% lsorig)]
  save(list = lsnow,
       file = paste0(RDF,
                     "/BackUp.RData"))
  rm(lsnow)

  #save time info of session in timeLog
  endTime <- Sys.time()
  pTime <- round(as.numeric(difftime(endTime,
                                     startTime,
                                     units = "mins")
                            ),
                 digits = 1)
  newTRow <- data.frame(time = pTime,
                        number_of_samples = length(filelist) / 2 ,
                        ID = ID)
  timeLog <- rbind(timeLog, newTRow)
  write.csv(timeLog, "./internal_files/timelog.csv", row.names = FALSE)

  # end of process message
  message("Finished Illumina Basic Analysis at:")
  message(Sys.time())
  message(paste0("The process took ", pTime,
                 " minutes."))

}
