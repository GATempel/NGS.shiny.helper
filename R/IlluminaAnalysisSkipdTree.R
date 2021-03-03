#' @title  Analysis of Illumina Files
#'
#' @description Copy of the IlluminaAnalysis script with the creation of the
#'              phylogenetic tree deleted
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


illumina_analysis_treeskip <- function(input_path = "./files",
                              # input_path is the path to a fastq file folder
                              # the default is set to "./files"
                              trunc_forw = 270,
                              # trunc_forw is integer to determine length at
                              # which forward reads are to be trunctuated
                              # default is set to 270
                              trim_left_forw = 10,
                              trunc_rev = 220,
                              # trunc_rev is integer to determine length at
                              # which reverse reads are to be trunctuated
                              # default is set to 220
                              trim_left_rev = 10,
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
                              multithreading = TRUE
                              # logical value
                              # if TRUE multiple processor cores will be used
                              # does not work on windows OS
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
                         stringsAsFactors = F)

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
                         full.names = T)
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
                                               full.names = T))
  FastqReverse <- gtools::mixedsort(list.files(SFO,
                                               pattern = "_R2_001.fastq.gz",
                                               full.names = T))

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
  if (resort == T){
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
                             pool = T)
  DadaReverse <- dada2::dada(FilterReverse,
                             ErrorReverse,
                             multithread = multithreading,
                             pool = T)

  # merge the forward and the matching reverse strands of the processed data
  MergerDada <- dada2::mergePairs(DadaForward,
                                  FilterForward,
                                  DadaReverse,
                                  FilterReverse,
                                  verbose = T)

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
                                        verbose = T)

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
  while (length(taxRank) > 1) {
    # iterates through the elements in the taxRank object
    x <- length(taxRank)
    # create an empty temporary object
    temps3 <- NULL
    # create a temporary dataframe on the basis of the given dataframe
    # all rows that contain the same elements within the taxonomic Rank columns
    # will be added onto each other
    tempdf <- plyr::ddply(taxaAbund, taxRank, plyr::numcolwise(sum))
    # create empty two column dataframe for the Names
    name <- data.frame(Name = as.character(),
                       Long_Name = as.character())
    # iterate through all rows of the dataframe
    for (i in 1:nrow(tempdf)) {
      # iterate through the taxRank elements except the last
      for (j in 1:(x-1))
        # checks if a name is assigned for the current row - if not it will
        # assign the entry of the lowest hierarchical taxonomic rank that is
        # not NA if the lowest possible rank is NA an NA remark will be added
        # to the name indicating at which rank a non NA entry was found
        # a long name will be assigned in a similar mannerism by combining the
        # entries of each rank entries containing an "unknown" string will be
        # replaced by "Unknown_'taxonomic_rank'_'the entry of the taxonomic
        # rank above'"
      {
        if
        (
          # check if the current entry is neither NA not contains an "unknown"
          # string nor the "incertae sedis" string and then assigns it as the
          # name
          is.na(name[i, 1])  &
          !is.na(tempdf[i, x]) &
          !stringr::str_detect(tempdf[i, x],
                               stringr::fixed("unknown",
                                              ignore_case = T)) &
          !stringr::str_detect(tempdf[i, x],
                               stringr::fixed("incertae sedis",
                                              ignore_case = T))
        )
        {
          name[i, 1] <- tempdf[i, x]
        }
        else if
        (
          # check the rest of the entries if they are neither NA nor contain
          # the "incertae sedis" string and then assign
          # "Unknown_'taxonomic_rank'_'taxonomic rank above = its entry'
          is.na(name[i, 1]) &
          !is.na(tempdf[i, x]) &
          !stringr::str_detect(tempdf[i, x],
                               stringr::fixed("incertae sedis",
                                              ignore_case = T))
        )
        {
          name [i, 1] <- paste0("Unknown_",
                                taxRank[x],
                                "(",
                                taxRank[x - j],
                                "=",
                                tempdf[i , x - j],
                                ")")


        }
        else if
        (
          # check if the entry is not NA and then assign
          # 'taxonomic rank (Icertae)_taxonomic rank above = its entry'
          is.na(name[i, 1]) &
          !is.na(tempdf[i, x])
        )
        {
          name [i, 1] <- paste0(taxRank[x],
                                "(Incertae)_",
                                taxRank[x - j],
                                "=",
                                tempdf[i , x - j])
        }
        else if
        (
          is.na(name[i, 1]) &
          !is.na(tempdf[i, x - j]) &
          !stringr::str_detect(tempdf[i, x - j],
                               stringr::fixed("unknown",
                                              ignore_case = T)) &
          !stringr::str_detect(tempdf[i, x - j],
                               stringr::fixed("incertae",
                                              ignore_case = T))
        )
        {
          name[i, 1] <- paste0(tempdf[i, x - j],
                               "_NA(",
                               taxRank[x - (j - 1)],
                               ")" )

        }
        else if
        (
          is.na(name[i, 1]) &
          !is.na(tempdf[i, x - j]) &
          !stringr::str_detect(tempdf[i, x - j],
                               stringr::fixed("incertae",
                                              ignore_case = T))
        )
        {
          name [i, 1] <- paste0("Unknown_",
                                taxRank[x - j],
                                "(",
                                taxRank[x - (j + 1)],
                                "=",
                                tempdf[i , x - (j + 1)],
                                ")")
        } else if
        (
          is.na(name[i, 1]) &
          !is.na(tempdf[i, x - j])
        )
        {
          name [i, 1] <- paste0(taxRank[x - j],
                                "(Incertae)_",
                                taxRank[x - (j + 1)],
                                "=",
                                tempdf[i , x - (j + 1)])
        }
      }
      for (k in 1:x) {
        if (is.na(name[i, 2]))
        {
          name[i, 2] <- tempdf[i, k]
        } else {
          if (is.na(tempdf[i, k]))
          {
            name[i, 2] <- paste0(name[i, 2], "_NA")
          } else {
            name[i, 2] <- paste0(name[i, 2], "_", tempdf[i,k])
          }
        }

      }
    }

    # bind the created elements within the name object onto the temporary df
    tempdf <- cbind(tempdf,
                    name)
    # rearranges the temporaray df to display all character columns before the
    # numeric columns
    tempdf <- tempdf[, c(names(dplyr::select_if(tempdf, is.character)),
                         names(dplyr::select_if(tempdf, is.numeric)))]
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
    tNames <- row.names(tempdf[base::order(tempdf$Total, decreasing = T),])

    # assign the dataframe and the sorted name element to a temporary S3 object
    temps3$Abundance <- tempdf
    temps3$relativeAbundance <- temprela
    temps3$Sorted_Names <- tNames

    # create a table containing just the information if the current taxonomic
    # rank exists within the samples
    tempOcc <- dplyr::select_if(dplyr::select(tempdf,
                                              -Total),is.numeric)
    tempOcc[tempOcc > 0] <- 1
    tempOcc <- cbind(dplyr::select(tempdf, where(is.character)), tempOcc)

    temps3$Occurance <- tempOcc

    message(paste0("Calculating Co-Occurance of ",
                   taxRank[x]))
    # calculates Co-Occurance using the Cooccur package
    temps3$Co_Occurance <- cooccur::cooccur(dplyr::select(tempOcc,
                                                          where(is.numeric)),
                                            spp_names = TRUE)

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
    assign(taxRank[x], temps3)

    # save an xls table with the abundances
    xlsx::write.xlsx(temps3$Abundance,
                     paste0(tabs,
                            "/Abundance.xls"),
                     sheetName = taxRank[x],
                     append = T)
    xlsx::write.xlsx(temps3$relativeAbundance,
                     paste0(tabs,
                            "/Abundance.xls"),
                     sheetName = paste0("relativ",
                                        taxRank[x]),
                     append = T)

    # move to the next higher taxonomic rank
    taxRank <- taxRank[1:x-1]

    # remove temporary files
    rm(tempdf)
    rm(temps3)
    rm(tempOcc)
    rm(temprela)
  }

  # create an S3 object for the lat remaining taxRank
  temps3 <- NULL
  tempdf <- plyr::ddply(taxaAbund, taxRank, plyr::numcolwise(sum))
  name <- data.frame(Name = as.character(),
                     Long_Name = as.character())
  for (i in 1:nrow(tempdf)) {
    name[i, 1] <- tempdf[i, 1]
    name[i, 2] <- tempdf[i, 1]
  }
  tempdf <- cbind(tempdf,
                  name)
  tempdf <- tempdf[, c(names(dplyr::select_if(tempdf, is.character)),
                       names(dplyr::select_if(tempdf, is.numeric)))]
  row.names(tempdf) <- tempdf$Name
  tempdf$Total <- rowSums(dplyr::select_if(tempdf, is.numeric))
  temprela <- tempdf[, names(dplyr::select_if(tempdf, is.numeric))]
  for (i in 1:ncol(temprela)) {
    temprela[,i] <- temprela[,i] * 100 / sum(temprela[,i])
  }
  temprela <- cbind(tempdf[names(dplyr::select_if(tempdf, is.character))],
                    temprela)
  tNames <- row.names(tempdf[base::order(tempdf$Total, decreasing = T),])
  temps3$Abundance <- tempdf
  temps3$relativeAbundance <- temprela
  temps3$Sorted_Names <- tNames
  temps3$Tidy <- tidyr::pivot_longer(tempdf,
                                     c(colnames(dplyr::select_if(tempdf,
                                                                 where(
                                                                   is.numeric)
                                     )
                                     )),
                                     names_to = "Sample",
                                     values_to = "Abundance")
  assign(taxRank[1], temps3)
  xlsx::write.xlsx(temps3$Abundance,
                   paste0(tabs,
                          "/Abundance.xls"),
                   sheetName = taxRank[1],
                   append = T)
  xlsx::write.xlsx(temps3$relativeAbundance,
                   paste0(tabs,
                          "/Abundance.xls"),
                   sheetName = paste0("relativ",
                                      taxRank[1]),
                   append = T)
  rm(tempdf)
  rm(temps3)
  rm(temprela)

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

  # create boolean Value to note if phylogenetic tree was created
  skipdTree <- T

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
  phylo <- phyloseq::phyloseq(phyloseq::otu_table(NonChime,
                                                  taxa_are_rows = F),
                              phyloseq::sample_data(SampleDataFrame),
                              phyloseq::tax_table(taxonomy))

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
                "skipdTree")
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
  write.csv(timeLog, "./internal_files/timelog.csv", row.names = F)

  # end of process message
  message("Finished Illumina Basic Analysis at:")
  message(Sys.time())
  message(paste0("The process took ", pTime,
                 " minutes."))

}
