#remove_packages remove all user installed pakcges, without removing any base packages for R 
remove_packages<- function (){
  # create a list of all installed packages
  ip <- as.data.frame(installed.packages())
  head(ip)
  # if you use MRO, make sure that no packages in this library will be removed
  ip <- subset(ip, !grepl("MRO", ip$LibPath))
  # we don't want to remove base or recommended packages either\
  ip <- ip[!(ip[,"Priority"] %in% c("base", "recommended")),]
  # determine the library where the packages are installed
  path.lib <- unique(ip$LibPath)
  # create a vector with all the names of the packages you want to remove
  pkgs.to.remove <- ip[,1]
  head(pkgs.to.remove)
  # remove the packages
  sapply(pkgs.to.remove, remove.packages, lib = path.lib)
}

#install_packages: this function take a list of packages and install them
install_packages <- function(package_list){
  list_installed <- installed.packages()
  new_pkgs <- subset(package_list, !(package_list %in% list_installed[, "Package"]))
  already_install_pkgs <- subset(package_list, (package_list %in% list_installed[, "Package"]))
  print(c(" packages already installed:",already_install_pkgs))
  print(strrep("_", 70))
  if(length(new_pkgs)!=0){
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(new_pkgs, update=F)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}

#output_table1: export table1 to png and word 
output_table1 <- function (table_name) { 
  pdf_table_name<-paste0(table_name,".pdf")
  png_table_name<-paste0(table_name,".png")
  word_table_name<-paste0(table_name,".docx")
  set_flextable_defaults(background.color = "white")
  table_output<-t1flex(get(table_name))
  table_output<-fontsize(table_output,size=10)
  table_output<-line_spacing(table_output, space=0)
  table_output<-autofit(table_output)
  save_as_image(table_output, path=pdf_table_name,1, webshot = "webshot")
  save_as_image(table_output, path=png_table_name, webshot = "webshot")
  xlsxfile=paste0(table_name,".xlsx")
  exportxlsx(table_output,xlsxfile)
  ##save_as_docx(table_output, path=word_table_name) ##this line makes a word document output
}

#render.NEW: adjust table1 content- to include median, Q1 and Q3 for category variables
render.NEW <- function(x, name, data2, ...) {
  name<<-name #to be used with pvalue
  num_NA <-sum(is.na(x))
  if (!is.numeric(x)){ 
    out<-   c("", sapply(stats.apply.rounding(stats.default(x)), function(y) with(y,sprintf("%s (%s)", FREQ, PCT))))
  }
  else{  
    median <- median(x, na.rm = T)
    Q1 <- round(quantile(x, 0.25, na.rm = T), digits=2)
    Q3 <- round(quantile(x, 0.75, na.rm = T), digits=2)
    out <- c("",
             "Median [Q1, Q3]" = paste0(median," ", "[", Q1, ", ", Q3, "]" )
    )
  }
  if (num_NA==0) { ###if there is no missing, then remove the line for missing
    missing<-NULL
  }
  else { ###if there is missing, then add an extra line for missing 
    missing<-with(stats.apply.rounding(stats.default(is.na(x),),)$Yes, c(Missing=sprintf("%s (%s)", FREQ, PCT)))
  }
  out<- append(out,missing)
}

#pvalue_posthoc: calculates p value for the table1 with posthoc analysis follows
pvalue_posthoc <- function(x, ...) {
  if ("overall" %in% names(x)) { ###"overall" exist when Total was selected from overall with the table1 function. this is to exclude the overall column from being pulled to do the p value calculation
    x<-x[names(x) %in% "overall" == FALSE] 
  }
  if (!exists("table1_list_significant")){ ###if table1_list_significant doesnt exist 
    table1_list_significant <-vector(mode="list")
    significant_data_frame <- vector("list",20)
    significant_data_frame2 <- vector("list",20)
    significant_data_frame3 <- vector("list",20)
    j<<-1
  }
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if ((length(unique(g)))==2) { ###if there are only two groups (g) 
    if (is.numeric(y)) {
          ########################################################################
          ###have to check to see if one of the group is all missing- because wilcox will crash if one of the group is all missing
          ########################################################################
          # Create a data frame combining y and g
          df <- data.frame(y, g)
          # Check if any factor class has all NA values
          has_all_na <- any(sapply(unique(g), function(class) all(is.na(df$y[df$g == class]))))
          if (has_all_na) {
            p <-NA
          } else {
            # For numeric variables, perform wilcox rank sum
            output <- wilcox.test(y ~ g)
            p<-output$p.value #extract p value
          }
    } else {
      # For categorical variables, perform a chi-squared test of independence
      output <- chisq.test(table(y, g), correct=F)
      p<-output$p.value #extract p value
    }
  }
  if((length(unique(g)))>2) { ###if there are only than two groups (g) 
    if (is.numeric(y)) {
      output <- kruskal.test(y ~ g)
    } 
    else {
      output <- chisq.test(table(y, g), correct=F)
    }
    p<-output$p.value #extract p value
  }
  if (!is.numeric(y)) { ### this will have p value as NA if the groups have all have same response 
    if (length(output$observed)==(length(unique(g)))) {
      p <-NA
    }
  }  
  if (!is.numeric(y) & (length(output$observed)==(length(unique(g))))) {} ### this will have p value as NA if the groups have all have same response 
  else if (p<.05) {
    table1_list_significant<-append(table1_list_significant,name) ##adding on
    significant_data_frame[[j]]<-y ##adding on
    significant_data_frame2[[j]]<-g ##adding on
    significant_data_frame3[[j]]<-x ##adding on
    j<<-j+1
  }
  table1_list_significant<<-unlist(table1_list_significant) #turned nested list to a vector with a list of variables which were significant
  significant_data_frame<<-significant_data_frame
  significant_data_frame2<<-significant_data_frame2
  significant_data_frame3<<-significant_data_frame3
  rm(name, envir=.GlobalEnv) ##remove name after adding to table1_list_significant
  # The initial empty string places the output on the line below the variable label.
  c("", format.pval(p, digits=3, eps=0.001))
}

#pvalue: calculates p value for the table1
pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  if ("overall" %in% names(x) ) {x<-x[names(x) != "overall"]} ###"overall" exist when Total was selected from overall with the table1 function. this is to exclude the overall column from being pulled to do the p value calculation
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if ((length(unique(g)))==2) { ###if there are only two groups (g) 
    if (is.numeric(y)) {
                ########################################################################
                ###have to check to see if one of the group is all missing- because wilcox will crash if one of the group is all missing
                ########################################################################
                # Create a data frame combining y and g
                df <- data.frame(y, g)
                # Check if any factor class has all NA values
                has_all_na <- any(sapply(unique(g), function(class) all(is.na(df$y[df$g == class]))))
      if (has_all_na) {
            p <-NA
      } else {
        # For numeric variables, perform wilcox rank sum
        output <- wilcox.test(y ~ g)
        p<-output$p.value #extract p value
      }
    } else {
      # For categorical variables, perform a chi-squared test of independence
      output <- chisq.test(table(y, g), correct=F)
      p<-output$p.value #extract p value
    }
  }
  if((length(unique(g)))>2) { ###if there are only than two groups (g) 
    if (is.numeric(y)) {
      # For numeric variables, perform Kruskal test
      output <- kruskal.test(y ~ g)
    } 
    else {
      output <- chisq.test(table(y, g), correct=F)
    }
    p<-output$p.value #extract p value
  }
  if (!is.numeric(y)) { ### this will have p value as NA if the groups have all have same response 
    if (length(output$observed)==(length(unique(g)))) {
      p <-NA
    }
  }  
  # The initial empty string places the output on the line below the variable label.
  rm(name, envir=.GlobalEnv) ##remove name which was added from render step
  c("", format.pval(p, digits=3, eps=0.001))
}

###this function pulls all png files from the folder and put out a PDF file
##pdf_file_output: has to the the pdf file you want to save to with the path
##input_filepath: is the folder that you store all of your png files
merge_png_pdf <- function(pdf_file_output, input_filepath) {
  library(png)
  library(grid)
  # Save your current working directory
  old_wd <- getwd()
  setwd(input_filepath)
  #sorting the file lists based on date of creation
  details = file.info(list.files(pattern='.*[.]png'))
  details = details[with(details, order(as.POSIXct(mtime))), ]
  files = rownames(details)
  ######################################################################################
  pdf(pdf_file_output) ##you open the pdf file and start going thro the PNG files to append them
  lapply(files,function(x){
    img <- as.raster(readPNG(x))
    grid.newpage()
    grid.raster(img, interpolate = FALSE)
  })
  ###this deletes all for the png files after you put them into the PDF files 
  lapply(files,function(x){
    file.remove(x) 
  })
  dev.off()
  # Come back to the main directory
  setwd(old_wd)
}

#remove_png_in_folder remove all png files within the directory
remove_png_in_folder <- function(input_filepath) {
  old_wd <- getwd()
  setwd(input_filepath)
  details = file.info(list.files(pattern='.*[.]png'))
  files = rownames(details)
  ###this deletes all for the png files after you put them into the PDF files 
  lapply(files,function(x){
    file.remove(x) 
  })
  # Come back to the main directory
  setwd(old_wd)
}

####table1_posthoc is used immediately after table1 and output_table1 to perform post hoc tests. this function should only be used when there are more than 3 groups of comparsion
table1_posthoc <- function (input) {
  excelfile<- paste0(input,".xlsx")
  print(excelfile)
  wbook<-loadWorkbook (excelfile)
  #wbook<-createWorkbook()
  for (k in 1:length(table1_list_significant)){
    addWorksheet(wb=wbook, sheetName=table1_list_significant[[k]])
    if (!is.numeric(significant_data_frame[[k]])) {
      exceloutput<-chisq.posthoc.test(table(significant_data_frame[[k]], significant_data_frame2[[k]]), method="bonferroni")
      for (l in 1:length(significant_data_frame3[[1]])) {
        t <-l+2 ###the first two columns within excelouput is "dimension" and "value" 
        colnames(exceloutput)[t]<- names(significant_data_frame3[[1]][l])
      }
    }
    else { ##if the variable is continuous then will do post hoc for kruskal-> dunn test
      levels(significant_data_frame2[[k]])<-c(names(significant_data_frame3[[k]]))
      exceloutput_t<-dunnTest(significant_data_frame[[k]] ~ significant_data_frame2[[k]], method="bonferroni")
      exceloutput<-exceloutput_t$res
    }
    writeData(wb=wbook, sheet=table1_list_significant[[k]], x=exceloutput)
  }
  saveWorkbook(wb = wbook, file = excelfile, overwrite = TRUE)
  ### now that we ran the post hoc tests, will remove these objects which were created for the process
  rm(significant_data_frame, envir=.GlobalEnv)
  rm(significant_data_frame2, envir=.GlobalEnv)
  rm(significant_data_frame3, envir=.GlobalEnv)
  rm(table1_list_significant, envir=.GlobalEnv)
  rm(j, envir=.GlobalEnv)
}

###########section_separator is a function that make a png with a text label. when all the PNGs are combine together, the PNG generated will be act as a section separator
section_separator <- function(text_to_display) {
  if (exists("section_separator_num")) {
    section_separator_num<-section_separator_num+1
  }
  else { #if there is no section_separator_num, then will start with 1
    section_separator_num<-1
  }
  img_Name <- paste("section_separator", section_separator_num, ".", "png", sep = "")
  #open new file for output
  png(img_Name, res = 300, width=1000, height=200 , units='mm')
  text<-ggplot() +
    annotate("text", x = 100,  y = 100,
             size = 20,
             label = text_to_display) + theme_void()
  print(text)
  .GlobalEnv$section_separator_num<-section_separator_num
  #close image
  dev.off()
}

##barplot_abundance: makes a barplot for abundance by sample
barplot_abundance <- function(phyloseq_subset,outputfile) {
  temp<-phyloseq_to_df(phyloseq_subset)
  temp <- temp %>% dplyr::select(-c("OTU", "Kingdom","Class","Order","Family","Genus","Species"))
  temp2 <- reshape2::melt(temp,id.vars=c("Phylum"))
  temp3 <- temp2 %>% dplyr::group_by(variable, Phylum) %>% dplyr::summarize(Sum_abund = sum(value))
  outputfilepng=paste(outputfile,".png")
  png(file=outputfilepng, res = 300, width=750, height=500 , units='mm')
  print({ggplot(temp3, aes_string(x = "variable", y = "Sum_abund", fill = "Phylum"))+ 
      geom_bar(stat = "identity", position = "stack", 
               color = "black")+ 
      theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
      labs(Title="taxa abundance by phylum", x="Sample", y="Abundnace")
  })
  dev.off()
  temp4<-reshape2::dcast(temp3, variable ~ Phylum)
  temp4 <- temp4 %>% dplyr::rename(sample=variable)
  temp5 <- as.matrix(temp4)
  temp5 <- t(temp5)
  outputfilecsv=paste(outputfile,".csv")
  write.csv(temp5,outputfilecsv, row.names = T)
}

####post_match_it_graphs: used after matchit function to give you descriptive plot and chart to better visualize the match
post_match_it_graphs <- function (formula, outputpng, titleplot) {
  outputpngfile <- paste0(outputpng,".png")
  formula <- formula(formula) ## turn to formula
  png(file=outputpngfile, res = 300, width=400, height=500 , units='mm')
  p1<-as.grob(love.plot(m.out1,abs=T, shapes = c("triangle", "square"), 
                        colors = c("blue", "darkgreen"),sample.names=c("All", "Matched"), title=titleplot))
  p2<-as.grob(function() plot(m.out1, type = "histogram"))
  m.data2_t <- match.data(m.out1)
  p3<-as.grob(function() plot(autofit(t1flex(table1(formula ,
                                                    data=m.data2_t, overall= F, render=render.NEW, extra.col=list(`P-value`=pvalue))))))
  grid.arrange(p1, p2, p3, ncol = 1)
  dev.off()
}

####wilcox_contam_test: Wilcoxon comparsion used to run decontam objects 
wilcox_contam_test <- function(relab,neg,threshold){
  p.value <- apply(relab,1,function(x){return(wilcox.test(as.numeric(as.character(x[!neg])),
                                                          as.numeric(as.character(x[neg])),
                                                          alternative="less",exact=FALSE)$p.value)})
  fdr <- unlist(lapply(p.value,function(x){return(p.adjust(x, method="BH", n = length(p.value)))}))
  contam <- as.data.frame(cbind(p.value, fdr, ifelse(p.value<threshold, 'TRUE', 'FALSE')))
  rownames(contam) <- rownames(relab)
  colnames(contam) <- c("p.value","fdr","contaminant")
  return(contam)
}

####decontaminant_KW: identify decontaminant based on different method
decontaminant_KW<-function(input_phyloseq, 
                  SampleID.unique=NULL, #if empty, SampleID.unique would be the rowname of the sample_data of the input_phyloseq
                  sample_type_var_name,                              
                  sample_types=list(), 
                  sample_type_color=list(), 
                  sample_type_color_2nd=list(), 
                  negative_sample_type, ###the sample type that you want to be as negative control
                  compare_type=list(), 
                  stat_option=c("mean", "median"), ### the statistics to determine rank for the boxplot
                  display_contam_method=c("MannWhit","preval","freq","combin","none"), #if none is selected then there will be no red label
                  bacterial_load, #for frequency only. uses ddPCR to determine library size
                  graph_option=c("boxplot","mean_SE", "mean_SD"), 
                  test_threshold,
                  log_scale, #will transform abundance to log(abundance*100 + 1)
                  output_suffix, 
                  taxa_genus_output=c("yes","no")) {#if want the taxa name to be shortened to genus level
  ################################################################################################################################################################################
  ##############################################################setting up all the function argument and data frames############################################################## 
  ################################################################################################################################################################################
  ###ensure intput_phyloseq is a absolute count phyloseq, not relative abundance phyloseq. and make sure phyloseq is in right orientation
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  countData_check <- floor(as(otu_table(input_phyloseq), "matrix")) # round down to nearest integer. if this phsyloeq is a relative abundance, then the entire countData_check would be 0
  if (all(countData_check==0)) {stop("Please use phyloseq with absolute count. The current input physloeq contains relative abundnace")}

  ## evaluate choices
  stat_option <- match.arg(stat_option)
  display_contam_method<- match.arg(display_contam_method)
  print(stat_option)

  if (missing(output_suffix)) { output_suffix<-"OP"} #if output_suffix is missing, then will just give a suffix _OP at the end of the file
  if (missing(graph_option)) { graph_option<-"boxplot"} #if graph_option is missing, then will just have boxplot as default
  if (missing(log_scale)) { log_scale<-"no"} #if graph_option is missing, then will just have boxplot as default

  if (missing(taxa_genus_output)) { taxa_genus_output<-"no"} #default being yes for taxa_rank_name, which will just give the selected taxa name
  taxa_genus_output <- match.arg(taxa_genus_output)

  ##making sure sample_types and sample_type_color and sample_type_color_2nd has same number of elements
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  #if sample_type_color_2nd is not missing then make sure it has same number of element as sample_type_color
  if (!missing(sample_type_color_2nd)) {
    stopifnot("sample_types need to have same number of elements as sample_type_color_2nd"= length(sample_types)==length(sample_type_color_2nd))
  } #if graph_option is missing, then will just have boxplot as default
  if (missing(sample_type_color_2nd)) {
    sample_type_color_2nd<-rep("black",length(sample_types)) #make the secondary color to be black by default, unless specify
  } 
  
  ###generating a variable for each sample type 
  for (i in 1:length(sample_types)) {
    assign(paste0("sample_type_color",i), sample_type_color[[i]])
    assign(paste0("sample_type_color_2nd",i), sample_type_color_2nd[[i]])
    assign(paste0("sampletype",i), sample_types[[i]])
  }  

  ###this make sure that the phyloseq ONLY contains the sample_type_var_name variable with all options listed sample_types
  #in theory input_phyloseq should be same as input_phyloseq_2 if the user input all of the sample_types options
  keep_sample<- get_variable(input_phyloseq,sample_type_var_name) %in% sample_types
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  input_phyloseq_2 <- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 

  ####preparing compare_type list
  ####if the compare_type has "A_and_B" then will determine if taxa is contaminant for A and B when they are individually compared to the negative control
  ####if the compare_type has "A_or_B" then will determine if taxa is contaminant for A OR B when they are individually compared to the negative control
  ####if the compare_type has "A_combine_B" then will determine if taxa is contaminant when negative control is compare to both A and B together
  for (i in 1:length(compare_type)) {
    assign(paste0("compare_type",i), compare_type[[i]])
  }  
  
  ####getting relative abundance
  normalizeSample <- function(x){x/sum(x)}
  input_phyloseq_2_relab <- transformSampleCounts(input_phyloseq_2,normalizeSample) ##relative abundance
  ### Setting up ###
  # Counts from Phyloseq 
  counts.edit <- as.data.frame(otu_table(input_phyloseq_2))
  # Relative Abundance Table from Phyloseq 
  relab.edit <- as.data.frame(otu_table(input_phyloseq_2_relab))
  # reference table (not to be used)
  reference_table <- as.data.frame(sample_data(input_phyloseq_2))

  ##if SampleID.unique is missing, then will use the rowname of match to be the SampleID.unique given the rowname of sample_data of the phyloseq should be unique
  if (is.null(SampleID.unique)) {
        reference_table$SampleID_decontam<-rownames(reference_table) #this makes a variable in the DF match call SampleID.unique which would be rowname of the DF match
  } else {
        if (SampleID.unique %in% names(reference_table)){    
          reference_table$SampleID_decontam<-reference_table[[SampleID.unique]] #this makes a variable in the DF match call SampleID.unique which would be rowname of the DF match
        } else {stop(paste0(SampleID.unique," does not exist in the sample dataframe"))} #if the phyloseq sample dataframe doesn't have the SampleID.unique then stop
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ###adding bacterial load option- for frequency only 
  ###setting up the bacterial load data
  if (!missing(bacterial_load)){ #if want to use bacterial load
    if (!inherits(bacterial_load, c("data.frame"))) { #making sure that bacterial load is a data.frame
      stop_txt = paste0("The input data for bacterial_load need to be a data frame")
      stop(stop_txt, call. = FALSE)
    }
    if (ncol(bacterial_load)!=2) { #make sure the bacterial load dataframe is in the right structure
      stop_txt2 = paste0("The input data frame for bacterial_load needs to have two columns:first column being the SampleID.unique (or rowname of the sample_data of the input_phyloseq if SampleID.unique is not inputted)")
      stop(stop_txt2, call. = FALSE)
    }
    if(length(setdiff(reference_table$SampleID_decontam,bacterial_load[,1]))!=0) { ###if there are elements in the reference_table that is not in the bacterial_load first column, which should be 
      differentID<-paste(setdiff(reference_table$SampleID_decontam,bacterial_load[,1]), collapse=",  ")  
      stop_txt3= paste("The following SampleIDs do not exist in the bacterial_load dataframe:",
                        differentID)
      stop(stop_txt3, call. = FALSE)
    }
    colnames(bacterial_load)<-c("SampleID", "ddPCR_count")
  } 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #######################################################################################################################################################################
  taxa.table <- as.data.frame(tax_table(input_phyloseq_2_relab))
  # This will paste the names together all the way from Family to OTU, separted by "." 
  taxa.table$match <- paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
  taxa.table.clean <- subset(taxa.table, select=c(match))
  taxa<-taxa.table$match #make a vector with all the taxa names
  
  # This will merge by column=0 which is the rowname by taxa-name-###
  # so now count.match and relab.match both have the new condensed taxa name 
  count.match <- merge(counts.edit, taxa.table.clean, by=0) 
  relab.match <- merge(relab.edit, taxa.table.clean, by=0)
  
  ### Replace count.match, relab.match rownames with taxa-name-###
  # This step replaces rownames for the match with count.match$match
  rownames(count.match) <- count.match$match
  rownames(relab.match) <- relab.match$match
  #remove match since match is now the rownames
  #remove Row.names given that is the byproduct of merging by rowname on previous step 
  counts <- subset(count.match, select=-c(Row.names, match))
  relab <- subset(relab.match, select=-c(Row.names, match))
  
  ### Match the Sample ID to the sample type categories  
  ### We need to do this to set up the comparison objects  
  ###replacing the SampleID.unique with the sample type for the count and relative abundance table
  index_var_name<-grep(sample_type_var_name, colnames(reference_table)) #this determine the nth column which sample+type_var_name is located at in the reference_table dataframe
  index2_var_name<-grep("SampleID_decontam", colnames(reference_table))
  #this step replace all of the unique subject ID (SampleID.unique) with their corresponding group of sample type of interest
  names(counts) <- reference_table[[index_var_name]][match(names(counts), reference_table[[index2_var_name]])] 
  names(relab) <- reference_table[[index_var_name]][match(names(relab), reference_table[[index2_var_name]])]  
  
  #### this input the sample types which will be used
  comparison.list<- vector(mode = "list", length = length(sample_types))
  names(comparison.list) <- sample_types
  for (i in 1:length(sample_types)) {
    comparison.list [i] <- list(grep(sample_types[i],colnames(counts)))
  }
  #comparison.list contain multiple lists (each list provide the column number of the specific sample type of interest). If there are 3 sample types in the dataframe of phyloseq then comparison will have a list of 3
  #counts.subset is basically counts, except there is now an index on the column name 
  counts.subset <- counts[,unlist(comparison.list)] # head(counts.subset)
  relab.subset <- relab[,unlist(comparison.list)] # head(relab.subset)
  #libsize consider the total abundance by each subject 
  libsizes <- colSums(counts.subset)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #####getting a dataframe that match the sample_type with column index with the initial sampleID
  if (!missing(bacterial_load)){ #if want to use bacterial load
    basedf<-count.match %>% select(-c("Row.names","match"))
    basedf_t<-data.frame(t(basedf))
    basedf_t$SampleID<-row.names(basedf_t)
    
    base2_df_t<-data.frame(t(counts.subset))
    base2_df_t$sample_index<-row.names(base2_df_t)
    base_combine<-merge(basedf_t,base2_df_t)
    base_combine<-base_combine %>% select(c(SampleID, sample_index))    
    
    bacterial_load_indexed<-merge(bacterial_load, base_combine, by="SampleID")
    
    bacterial_load_indexed2<-bacterial_load_indexed[order(bacterial_load_indexed$sample_index),]
    rownames(bacterial_load_indexed2)<-bacterial_load_indexed2$sample_index
    bacterial_load_indexed2<-bacterial_load_indexed2 %>% select(ddPCR_count)
    
    bacterial_load_indexed3<-as.numeric(unlist(bacterial_load_indexed2))
    names(bacterial_load_indexed3)<-rownames(bacterial_load_indexed2)  
    libsizes<-bacterial_load_indexed3 #library size is replaced with the ddPCR data
    print(paste0(bacterial_load," was used for frequency contaminant determination"))
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  ### Create matrix of data for evaluation ### 
  ### This is to create an empty matrix for all top statistics ### 
  X<-nrow(relab.subset)
  top.relab.stats <- as.data.frame(matrix(ncol=length(sample_types),nrow=X))
  
  # Moving the rownames to top.relab.stats 
  # going create two columns for medians to plot 
  # This may not be necessary top.relab.stats <- NULL 
  rownames(top.relab.stats) <- rownames(relab.subset)

  #loop through the sample types and get statistics for each. if there are 3 sample_types, then this will fill up the first 3 column with the statistics
  for (i in 1:length(sample_types)) {
    if (stat_option=="median") {
      top.relab.stats[i] <-rowMedians(as.matrix(relab.subset[, grepl(sample_types[i], colnames(relab.subset))]))
    } else if (stat_option=="mean") {
      top.relab.stats[i] <-rowMeans2(as.matrix(relab.subset[, grepl(sample_types[i], colnames(relab.subset))]))
    }
    colnames(top.relab.stats)[i] <- c(paste0("stats.",sample_types[i]))
  }  
  #loop through the sample types and fill in the next several column with the rank order by that particular sample type. If there are 3 sample_types, then this will give rank order on 4th-6th column
  for (i in 1:length(sample_types)) {
    ##generating the rank order by sample types 
    j <- i+length(sample_types)
    column<-top.relab.stats[,i]
    top.relab.stats <- top.relab.stats[order(column, decreasing = TRUE), , drop = FALSE ]
    top.relab.stats[j] <- 1:nrow(top.relab.stats)
    colnames(top.relab.stats)[j] <- c(paste0("rank.order.",sample_types[i]))
  }  
  #this is to make taxa name with rank order in them
  for (i in 1:length(sample_types)) { 
    j <- i+length(sample_types)
    k <- i+2*length(sample_types)  
    column2<-top.relab.stats[,j]
    top.relab.stats[k] <- paste0(rownames(top.relab.stats)," ","(",column2,")")
    colnames(top.relab.stats)[k] <- c(paste0(paste0("rank.order.",sample_types[i]),".name"))
  }  
  #################################
  #setting up the negative control#
  #################################
  #subject-level total count for the negative control
  libsizes_negative <- libsizes[grepl(negative_sample_type, names(libsizes))]
  w<-match(negative_sample_type,sample_types)
  comparison.list_negative<- unlist(comparison.list[w]) #extract the list of negative control (this list give the column number which is consider negative control)
  comparison.list_negative<-ifelse(comparison.list_negative, T, F) #basically making this list turn into all 1 
  #for counts.subset_negative and relab.subset_negative, the columns are subject (they should only be in the negative controls)
  counts.subset_negative <- dplyr::select(counts.subset, starts_with(negative_sample_type))
  relab.subset_negative <- dplyr::select(relab.subset, starts_with(negative_sample_type))
  
  #top.relab.stats now has the summary statistics for each sample type and the corresponding rank order. each row is a specific taxa 
  for (rank_sample_type in sample_types) { #this determines the rank order 
    print(paste0("ranking:", rank_sample_type))
    l<-match(rank_sample_type,sample_types) ###look at the element in rank_sample_type to see what their index at sample_types
    r<-l+length(sample_types) 
    # Re-rank based upon the rank group decreasing from highest to lowest
    column3<-top.relab.stats[,r] ###r would be the column number which consider the rank order for the ith sample_types
    top.relab.stats <- top.relab.stats[order(column3, decreasing = F), , drop = FALSE ]
    top.2 <- top.relab.stats
    top.3 <- head(top.2, 100) #top.3 would contain the top 100 taxa by the rank order of ith sample_types 
    
    # Filtering ONLY the top 100 taxa 
    relab.subset.top <- relab.subset %>% 
      dplyr::filter(rownames(relab.subset) %in% rownames(top.3))
    
    # Re ordered and saved into a new group (previously filtered)
    #relab.subset.top.match has same order as top.3; otherwise relab.subset.top.match is same as relab.subset.top
    relab.subset.top.match <- relab.subset.top[match(rownames(top.3),rownames(relab.subset.top)),]  

    for (h in 1:length(compare_type)) {
      name_compare_type<-get(paste0("compare_type",h))
      print(paste0(negative_sample_type," compare:"))
      print(name_compare_type)
      AND<-grepl("_and_",name_compare_type, fixed=T)
      OR<-grepl("_or_",name_compare_type, fixed=T)
      COMBINE<-grepl("_combine_",name_compare_type, fixed=T)
      contam<-list()
      if ((AND==TRUE) | (OR==TRUE)) {
              if (AND==TRUE) {
                   compare_sublist <- unlist(strsplit(name_compare_type, "_and_"))
              }else if (OR==TRUE) {
                   compare_sublist <- unlist(strsplit(name_compare_type, "_or_"))
              }
              #loop through the two sample types 
              for (i in 1:length(compare_sublist)) { 
                assign(paste0("compare_sublist",i), compare_sublist[i])
              } 
      
              for (b in 1:length(compare_sublist)) {
                        # only keeping the compare list  
                        selected_sample_type<-get(paste0("compare_sublist",b))
                        p<-match(selected_sample_type,sample_types) 
                        comparison.list_select<- unlist(comparison.list[p]) 
                        comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
                        
                        # Remove unused from conc, in this case libsizes 
                        libsizes_select <- libsizes[grepl(selected_sample_type, names(libsizes))]
                        counts.subset_select <- dplyr::select(counts.subset, starts_with(selected_sample_type))
                        relab.subset_select <- dplyr::select(relab.subset, starts_with(selected_sample_type))
                        
                        counts.subset_comb<-merge(counts.subset_select,counts.subset_negative, by=0) #put back the negative control with the sample that you want to compare 
                        relab.subset_comb<-merge(relab.subset_select,relab.subset_negative, by=0)
                        #remove the new column "Row.names" created when merging by 0
                        rownames(counts.subset_comb) <- counts.subset_comb$Row.names
                        rownames(relab.subset_comb) <- relab.subset_comb$Row.names
                        counts.subset_comb <- subset(counts.subset_comb, select=-c(Row.names))
                        relab.subset_comb <- subset(relab.subset_comb, select=-c(Row.names))
                        
                        negative_comb<- append(comparison.list_select,comparison.list_negative)
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                        ###need to sort libsizes_select and libsizes_negative
                        ###this is to ensure that the correct library size is link to the corresponding samples 
                        libsizes_select_TEMP<-data.frame(libsizes_select)
                        libsizes_select_TEMP$index<- as.numeric(sub(".*\\.", "", row.names(libsizes_select_TEMP)))
                        libsizes_select_TEMP$index[is.na(libsizes_select_TEMP$index)]<-0 #replace NA with 0 
                        libsizes_select_TEMP <- libsizes_select_TEMP[order(libsizes_select_TEMP$index), ]
                        libsizes_select2<-as.numeric(libsizes_select_TEMP$libsizes_select)
                        names(libsizes_select2)<-row.names(libsizes_select_TEMP)
                        
                        libsizes_negative_TEMP<-data.frame(libsizes_negative)
                        libsizes_negative_TEMP$index<- as.numeric(sub(".*\\.", "", row.names(libsizes_negative_TEMP)))
                        libsizes_negative_TEMP$index[is.na(libsizes_negative_TEMP$index)]<-0 #replace NA with 0 
                        libsizes_negative_TEMP <- libsizes_negative_TEMP[order(libsizes_negative_TEMP$index), ]
                        libsizes_negative2<-as.numeric(libsizes_negative_TEMP$libsizes_negative)
                        names(libsizes_negative2)<-row.names(libsizes_negative_TEMP)
                        libsizes_comb <- append(libsizes_select2,libsizes_negative2 )
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#      
                        
                        contam$prev <- isContaminant(t(as.matrix(counts.subset_comb)),
                                                     method="prevalence",
                                                     neg=negative_comb,
                                                     threshold=test_threshold)
                        contam$freq <- isContaminant(t(as.matrix(counts.subset_comb)),
                                                     method="frequency",
                                                     neg=negative_comb,
                                                     conc=libsizes_comb,
                                                     threshold=test_threshold)
                        contam$combined <- isContaminant(t(as.matrix(counts.subset_comb)),
                                                         method="combined",
                                                         neg=negative_comb,
                                                         conc=libsizes_comb,
                                                         threshold=test_threshold)
                        contam$wilcox <- wilcox_contam_test(relab.subset_comb,
                                                            neg=negative_comb,
                                                            threshold=test_threshold)
                        contam_dataframe<- contam %>% data.frame #compress the lists from the different results
                        contam_dataframe$wilcox.contaminant<-as.logical(contam_dataframe$wilcox.contaminant)
          
                        contam_result<-subset(contam_dataframe, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
                        assign(paste0("contam_result",b),contam_result)
                }
            #combine the results for the two sample types
            contam_result_final<-merge(contam_result1,contam_result2, by=0)
            rownames(contam_result_final) <- contam_result_final$Row.names
            contam_result_final <- subset(contam_result_final, select=-c(Row.names))
            if (AND==TRUE) {
              contam_result_final <- contam_result_final %>% mutate(prev.contaminant = case_when((prev.contaminant.x==T & prev.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
              contam_result_final <- contam_result_final %>% mutate(freq.contaminant = case_when((freq.contaminant.x==T & freq.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
              contam_result_final <- contam_result_final %>% mutate(combined.contaminant = case_when((combined.contaminant.x==T & combined.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
              contam_result_final <- contam_result_final %>% mutate(wilcox.contaminant = case_when((wilcox.contaminant.x==T & wilcox.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
              contam_result_final <- subset(contam_result_final, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
            }else if (OR==TRUE) {
              contam_result_final <- contam_result_final %>% mutate(prev.contaminant = case_when((prev.contaminant.x==T | prev.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
              contam_result_final <- contam_result_final %>% mutate(freq.contaminant = case_when((freq.contaminant.x==T | freq.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
              contam_result_final <- contam_result_final %>% mutate(combined.contaminant = case_when((combined.contaminant.x==T | combined.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
              contam_result_final <- contam_result_final %>% mutate(wilcox.contaminant = case_when((wilcox.contaminant.x==T | wilcox.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))   
              contam_result_final <- subset(contam_result_final, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
            }
      }else{
            if (COMBINE==TRUE) { ###if compare_type is a combination set 
                    ###this works if A_combine_B... it basically combine subset A and subset B
                    combine_list <- unlist(strsplit(name_compare_type, "_combine_"))
                    f<-match(combine_list,sample_types) 
                    comparison.list_select<- unlist(comparison.list[f]) 
                    comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
                    combine_list_temp<-str_replace(name_compare_type, "_combine_", "|") 
                    #Setting up the dataset so that it include both sample types that are combine =
                    libsizes_select <- libsizes[grepl(combine_list_temp, names(libsizes))]
                    counts.subset_select <- dplyr::select(counts.subset, starts_with(combine_list))
                    relab.subset_select <- dplyr::select(relab.subset, starts_with(combine_list))
            } else { ###if compare_type does not contain _and_, _or_, _combine_.... then basically compare_type is just one of the sample types
            if (name_compare_type %in% sample_types) { #name_compare_type has to be one of the initial sample_types if it is not _and_, _or_, _combine_
                    f<-match(name_compare_type,sample_types) 
                    comparison.list_select<- unlist(comparison.list[f]) 
                    comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
                    #keeping only the sample type of interest
                    libsizes_select <- libsizes[grepl(name_compare_type, names(libsizes))]
                    counts.subset_select <- dplyr::select(counts.subset, starts_with(name_compare_type))
                    relab.subset_select <- dplyr::select(relab.subset, starts_with(name_compare_type))
              } else { stop(paste0(paste0(name_compare_type," must be in sample type list: "),sample_types)) }
          
            }
            #combining with the negative control group
            counts.subset_comb<-merge(counts.subset_select,counts.subset_negative, by=0) #put back the negative control with the sample that you want to compare 
            relab.subset_comb<-merge(relab.subset_select,relab.subset_negative, by=0)
            #remove the new column "Row.names" created when merging by 0
            rownames(counts.subset_comb) <- counts.subset_comb$Row.names
            rownames(relab.subset_comb) <- relab.subset_comb$Row.names
            counts.subset_comb <- subset(counts.subset_comb, select=-c(Row.names))
            relab.subset_comb <- subset(relab.subset_comb, select=-c(Row.names))
      
            negative_comb<- append(comparison.list_select,comparison.list_negative)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            ###need to sort libsizes_select and libsizes_negative
            ###this is to ensure that the correct library size is link to the corresponding samples 
            libsizes_select_TEMP<-data.frame(libsizes_select)
            libsizes_select_TEMP$index<- as.numeric(sub(".*\\.", "", row.names(libsizes_select_TEMP)))
            libsizes_select_TEMP$index[is.na(libsizes_select_TEMP$index)]<-0 #replace NA with 0 
            libsizes_select_TEMP <- libsizes_select_TEMP[order(libsizes_select_TEMP$index), ]
            libsizes_select2<-as.numeric(libsizes_select_TEMP$libsizes_select)
            names(libsizes_select2)<-row.names(libsizes_select_TEMP)
            
            libsizes_negative_TEMP<-data.frame(libsizes_negative)
            libsizes_negative_TEMP$index<- as.numeric(sub(".*\\.", "", row.names(libsizes_negative_TEMP)))
            libsizes_negative_TEMP$index[is.na(libsizes_negative_TEMP$index)]<-0 #replace NA with 0 
            libsizes_negative_TEMP <- libsizes_negative_TEMP[order(libsizes_negative_TEMP$index), ]
            libsizes_negative2<-as.numeric(libsizes_negative_TEMP$libsizes_negative)
            names(libsizes_negative2)<-row.names(libsizes_negative_TEMP)
            libsizes_comb <- append(libsizes_select2,libsizes_negative2 )
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#      
            
            contam$prev <- isContaminant(t(as.matrix(counts.subset_comb)),
                                         method="prevalence",
                                         neg=negative_comb,
                                         threshold=test_threshold)
            contam$freq <- isContaminant(t(as.matrix(counts.subset_comb)),
                                         method="frequency",
                                         neg=negative_comb,
                                         conc=libsizes_comb,
                                         threshold=test_threshold)
            contam$combined <- isContaminant(t(as.matrix(counts.subset_comb)),
                                             method="combined",
                                             neg=negative_comb,
                                             conc=libsizes_comb,
                                             threshold=test_threshold)
            contam$wilcox <- wilcox_contam_test(relab.subset_comb,
                                                neg=negative_comb,
                                                threshold=test_threshold)
            contam_dataframe<- contam %>% data.frame #compress the lists from the different results
            contam_dataframe$wilcox.contaminant<-as.logical(contam_dataframe$wilcox.contaminant)
            
            contam_result_final<-subset(contam_dataframe, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
      }
      lowerbound<-length(sample_types)+1
      upperbound<-2*length(sample_types)
      rankorder <- subset(top.relab.stats[c(lowerbound:upperbound)]) ##keeping the rank order only
      csv_output_contam <- merge(contam_result_final,rankorder,by="row.names")
      
      contam_result_final_t<-contam_result_final
      #adding the contaminant result back to the relative abundance table (top.3 which contains the top 100 taxa)
      if (display_contam_method=="MannWhit"){
            contam_result_final_t$contaminant <- contam_result_final_t$wilcox.contaminant  
      }
      else if (display_contam_method=="preval") {
            contam_result_final_t$contaminant <- contam_result_final_t$prev.contaminant
      }
      else if (display_contam_method=="freq") {
            contam_result_final_t$contaminant <- contam_result_final_t$freq.contaminant
      }
      else if (display_contam_method=="combin") {
            contam_result_final_t$contaminant <- contam_result_final_t$combined.contaminant
      }
      else if (display_contam_method=="none") {
            contam_result_final_t$contaminant <- "FALSE" #if display_contam_method is none then turn all the contaminant value to FALSE so all the label would be black
      }
      contam_result_final_t<-subset(contam_result_final_t, select=c(contaminant))
      top.4<-merge(top.3,contam_result_final_t, by=0, keep.x=T)
      rownames(top.4) <- top.4$Row.names
      top.4 <- subset(top.4, select=-c(Row.names))
      top.4$color <- ifelse(top.4$contaminant, "red", "black")

      for (a in sample_types) {
            assign(paste0("rank_order_",a),subset(top.4, select=c(paste0("rank.order.",a), "color", paste0(paste0("rank.order.",a),".name") )))
      }      
      
      # Transposition of the figure, maybe at this point replace the column names 
      relab.subset.top.match.transpose <- as.data.frame(t(relab.subset.top.match))
      
      # name a new column with categories 
      relab.subset.top.match.transpose$category <- rownames(relab.subset.top.match.transpose)
      
      for (a in sample_types) {
        relab.subset.top.match.transpose$category[grepl(a,relab.subset.top.match.transpose$category)] <- a
      }  

      # Reverse column order 
      relab.subset.top.match.transpose.reverse <- relab.subset.top.match.transpose[,order(ncol(relab.subset.top.match.transpose):1)]
      for (a in sample_types) {
        assign(paste0("relab.subset.top.match.transpose.reverse.",a),subset(relab.subset.top.match.transpose.reverse, category == c(a) ))
        assign(paste0("relab.subset.top.match.transpose.reverse.",a),subset(get(paste0("relab.subset.top.match.transpose.reverse.",a)), select=-c(category)))
      }  

      ###results plot
      plot.df <- as.data.frame(rbind(cbind("Mann-Whitney U test",rownames(contam_result_final)[contam_result_final$wilcox.contaminant == T]),
                                     cbind("decontam prevalence",rownames(contam_result_final)[contam_result_final$prev.contaminant == T]),
                                     cbind("decontam frequency",rownames(contam_result_final)[contam_result_final$freq.contaminant == T]),
                                     cbind("decontam combined",rownames(contam_result_final)[contam_result_final$combined.contaminant == T])))
      
      ### Plotting the results ### 
      # plot_contam_plot 
      top.4<-top.4[order(top.4[[paste0("rank.order.",rank_sample_type)]]),]
      plot.df <- plot.df[plot.df[,2] %in% rownames(top.3),]
      
      plot.df[,1] <- factor(plot.df[,1],levels=c("Mann-Whitney U test","decontam prevalence","decontam frequency","decontam combined"))
      plot.df[,2] <- factor(plot.df[,2],levels=rev(rownames(top.3)))
      plot1 <- ggplot()+
        geom_point(mapping=aes(x=!!plot.df[,1],y=!!plot.df[,2],color=!!plot.df[,1]),size=2)+
        guides(color = "none")+
        theme_bw()+  
        scale_color_discrete(drop=FALSE) +
        scale_x_discrete(drop = FALSE)+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1))
      #taxa naming
      if (tolower(taxa_genus_output)=="no"){
            plot1 <- plot1 + 
            scale_y_discrete(drop = FALSE, labels = function(y) {
              y_org<-y
              is_long <- nchar(y) > 35
              y[is_long] <- paste0(substr(sapply(strsplit(y[is_long],".",fixed=T),"[[",4),1,10),
                                   paste0("..",paste0(substr(sapply(strsplit(y[is_long],".",fixed=T),"[[",5),1,10),".."),
                                          paste0(substr(sapply(strsplit(y[is_long],".",fixed=T),"[[",6),1,10),
                                                 paste0("..",str_sub(y[is_long],-4,-1)))))
              is_short <- nchar(y) < 17
              y[is_short] <- paste0(paste0(substr(sapply(strsplit(y_org[is_short],".",fixed=T),"[[",3),1,15),".."),
                                    paste0(substr(sapply(strsplit(y_org[is_short],".",fixed=T),"[[",4),1,5),
                                           paste0("..",paste0(substr(sapply(strsplit(y_org[is_short],".",fixed=T),"[[",5),1,5),".."),
                                                  paste0(substr(sapply(strsplit(y_org[is_short],".",fixed=T),"[[",6),1,6),
                                                         paste0("..",str_sub(y_org[is_short],-4,-1))))))
              y}) 
      } else if (tolower(taxa_genus_output)=="yes") {
            plot1 <- plot1 + 
              scale_y_discrete(drop = FALSE, labels = function(y) {
                y_org<-y
                genus_no_NA<- sapply(strsplit(y,".",fixed=T),"[[",5)!="NA"
                genus_NA<- sapply(strsplit(y,".",fixed=T),"[[",4)!="NA" & sapply(strsplit(y,".",fixed=T),"[[",5)=="NA"
                family_NA<- sapply(strsplit(y,".",fixed=T),"[[",3)!="NA" & sapply(strsplit(y,".",fixed=T),"[[",4)=="NA"
                order_NA<- sapply(strsplit(y,".",fixed=T),"[[",2)!="NA" & sapply(strsplit(y,".",fixed=T),"[[",3)=="NA"
                class_NA<- sapply(strsplit(y,".",fixed=T),"[[",1)!="NA" & sapply(strsplit(y,".",fixed=T),"[[",2)=="NA"
                y[genus_no_NA] <- paste0("g_",sapply(strsplit(y[genus_no_NA],".",fixed=T),"[[",5))
                y[genus_NA] <- paste0(paste0("f_",sapply(strsplit(y[genus_NA],".",fixed=T),"[[",4)),".g_NA")
                y[family_NA] <- paste0(paste0("o_",sapply(strsplit(y[family_NA],".",fixed=T),"[[",3)),".f_NA.g_NA")
                y[order_NA] <- paste0(paste0("c_",sapply(strsplit(y[order_NA],".",fixed=T),"[[",2)),".o_NA.f_NA.g_NA")
                y[class_NA] <- paste0(paste0("p_",sapply(strsplit(y[class_NA],".",fixed=T),"[[",1)),".c_NA.o_NA.f_NA.g_NA")
                y}) 
      }
      myplots <- list()  # new empty list
      myplots[[1]]<- plot1
      
      ###looping thro each sample type 
      number<-1
      number2<-2
      
    for (a in sample_types) {
        #this gets relab.subset.top.match.transpose.reverse.____ from wide to long 
        temp <- get(paste0("relab.subset.top.match.transpose.reverse.",a)) %>% dplyr::select(everything()) %>% tidyr::gather("id", "value",1:100, factor_key = TRUE)
        temp <- merge(temp, top.4, by.x="id", by.y="row.names")
        temp$id2<- temp[[paste0(paste0("rank.order.",a),".name")]]
        #if log_scale is yes then will have log(abundance+1)
        if (tolower(log_scale)=="yes"){
          temp$value2<-log10(temp$value*100+1) 
        } else {
          temp$value2<-temp$value
        }
        print("sample_types")
        print(a)
        print("rank order")
        print(rank_sample_type)
        p <- temp %>%   mutate(id2 = fct_reorder(id2, -get(paste0("rank.order.",rank_sample_type)))) %>%
          ggplot(., aes(x=id2, y=value2)) + 
          coord_flip() +
          theme_bw() + 
          theme(axis.text.y=element_text(color=rev(top.4$color)),axis.title.y=element_blank()) +
          theme(axis.text.x=element_text(angle=90,hjust=1))+
          ggtitle(a)
        #taxa naming
        if (tolower(taxa_genus_output)=="no"){
          p <-  p + 
            scale_x_discrete(labels = function(x) {
                x_org<-x
                is_long <- nchar(x) > 35
                x[is_long] <- paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",4),1,10),
                                     paste0("..",paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",5),1,10),".."),
                                            paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",6),1,10),
                                                   paste0("..",str_sub(x[is_long],-4,-1)))))
                is_short <- nchar(x) < 17
                x[is_short] <- paste0(paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",3),1,15),".."),
                                      paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",4),1,5),
                                             paste0("..",paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",5),1,5),".."),
                                                    paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",6),1,5),
                                                           paste0("..",str_sub(x_org[is_short],-4,-1))))))
                x})
        } else if (tolower(taxa_genus_output)=="yes") {
          p <-  p + 
            scale_x_discrete(drop = FALSE, labels = function(x) {
              x_org<-x
              genus_no_NA<- sapply(strsplit(x,".",fixed=T),"[[",5)!="NA"
              genus_NA<- sapply(strsplit(x,".",fixed=T),"[[",4)!="NA" & sapply(strsplit(x,".",fixed=T),"[[",5)=="NA"
              family_NA<- sapply(strsplit(x,".",fixed=T),"[[",3)!="NA" & sapply(strsplit(x,".",fixed=T),"[[",4)=="NA"
              order_NA<- sapply(strsplit(x,".",fixed=T),"[[",2)!="NA" & sapply(strsplit(x,".",fixed=T),"[[",3)=="NA"
              class_NA<- sapply(strsplit(x,".",fixed=T),"[[",1)!="NA" & sapply(strsplit(x,".",fixed=T),"[[",2)=="NA"
              x[genus_no_NA] <- paste0("g_",sapply(strsplit(x[genus_no_NA],".",fixed=T),"[[",5))
              x[genus_NA] <- paste0(paste0("f_",sapply(strsplit(x[genus_NA],".",fixed=T),"[[",4)),".g_NA")
              x[family_NA] <- paste0(paste0("o_",sapply(strsplit(x[family_NA],".",fixed=T),"[[",3)),".f_NA.g_NA")
              x[order_NA] <- paste0(paste0("c_",sapply(strsplit(x[order_NA],".",fixed=T),"[[",2)),".o_NA.f_NA.g_NA")
              x[class_NA] <- paste0(paste0("p_",sapply(strsplit(x[class_NA],".",fixed=T),"[[",1)),".c_NA.o_NA.f_NA.g_NA")
              x}) 
        }
        if (tolower(log_scale)=="yes"){
          p<-p+ylab("Log 10 (Relative Abundance)")
        } else {
          p<-p+ylab("Relative Abundance")
        }
        if (graph_option=="boxplot") {
            p<-p+geom_boxplot(color=sample_type_color[[number]], alpha=0.2) 
        } else if(graph_option=="mean_SD") {
            p<-p+geom_jitter(color=sample_type_color[[number]],position=position_jitter(0), alpha=0.2) +
              stat_summary(fun= mean, 
                           geom="pointrange", 
                           fun.max = function (x) mean(x)+sd(x),
                           fun.min = function (x) ifelse( mean(x)-sd(x) < 0, 0, mean(x)-sd(x)),
                           color=sample_type_color_2nd[[number]], linewidth =1.0, size=0.4)
        } else if(graph_option=="mean_SE") {
          p<-p+geom_jitter(color=sample_type_color[[number]],position=position_jitter(0), alpha=0.2) +
            stat_summary(fun= mean, 
                         geom="pointrange", 
                         fun.max = function (x) mean(x)+sd(x),
                         fun.min = function (x) ifelse(mean(x)-sd(x)/sqrt(length(x)) < 0, 0, mean(x)-sd(x)/sqrt(length(x))),
                         color=sample_type_color_2nd[[number]], linewidth =1.0, size=0.4)
        }
        myplots[[number2]] <- p
        number<-number+1
        number2<-number2+1
      }
      
      ### GGARRANGE ALL PLOTS ### 
      plot.all <- ggarrange(plotlist=myplots, ncol=length(myplots), nrow=1, align="h")
    
      pdf_output=paste0(paste0(paste0(paste0(paste0("contaminant",paste0("_negC_",negative_sample_type)),paste0("_compare_",name_compare_type)),paste0("_rank_",rank_sample_type)),paste0("thres_",test_threshold)),paste0(paste0("_",output_suffix),".pdf"))
      pdf(file=pdf_output, width=20+(length(sample_types)-3)*5, height=12)
             show(plot.all)
      dev.off()        
      
      ####export csv
      csv_output=paste0(paste0(paste0(paste0("contaminant",paste0("_negC_",negative_sample_type)),paste0("_compare_",name_compare_type)),paste0("thres_",test_threshold)),paste0(paste0("_",output_suffix),".csv"))
      xlsx_output=paste0(paste0(paste0(paste0("contaminant",paste0("_negC_",negative_sample_type)),paste0("_compare_",name_compare_type)),paste0("thres_",test_threshold)),paste0(paste0("_",output_suffix),".xlsx"))
      write.csv(csv_output_contam,csv_output, row.names = FALSE)
      write.xlsx(csv_output_contam, xlsx_output, sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)
    }
  }
}


####decontaminant_testing_KW: identify decontaminant based on a particular method at a range of test threshold level 
decontaminant_testing_KW <- function (input_phyloseq, 
                                      SampleID.unique=NULL, #if empty, SampleID.unique would be the rowname of the sample_data of the input_phyloseq
                                      sample_type_var_name, 
                                      sample_types=list(), 
                                      negative_sample_type, ###the sample type that you want to be as negative control
                                      compare_type=list(), 
                                      bacterial_load, #for frequency only. uses ddPCR to determine library size
                                      method_type=c("MannWhit","preval","freq","combin"),
                                      stat_option=c("mean", "median"), ### the statistics to determine rank for the boxplot
                                      test_threshold_all_level=list(), ### test_threshold_all_level=c(0.1,0.3,0.5,0.7,0.9) would test for the level specified and create an excel spread sheet with the parameters
                                      output_excel)  { 
  
  ################################################################################################################################################################################
  ##############################################################setting up all the function argument and data frames############################################################## 
  ################################################################################################################################################################################
  ###ensure intput_phyloseq is a absolute count phyloseq, not relative abundance phyloseq. and make sure phyloseq is in right orientation
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  countData_check <- floor(as(otu_table(input_phyloseq), "matrix")) # round down to nearest integer. if this phsyloeq is a relative abundance, then the entire countData_check would be 0
  if (all(countData_check==0)) {stop("Please use phyloseq with absolute count. The current input physloeq contains relative abundnace")}
  
  ## evaluate choices
  stat_option <- match.arg(stat_option)
  method_type <- match.arg(method_type)
  print(stat_option)
  print(method_type)

  ###generating a variable for each sample type 
  for (i in 1:length(sample_types)) {
    assign(paste0("sampletype",i), sample_types[[i]])
  }
  
  ###this make sure that the phyloseq ONLY contains the sample_type_var_name variable with all options listed sample_types
  #in theory input_phyloseq should be same as input_phyloseq_2 if the user input all of the sample_types options
  keep_sample<- get_variable(input_phyloseq,sample_type_var_name) %in% sample_types
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  input_phyloseq_2 <- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
  
  ####preparing compare_type list
  ####if the compare_type has "A_and_B" then will determine if taxa is contaminant for A and B when they are individually compared to the negative control
  ####if the compare_type has "A_or_B" then will determine if taxa is contaminant for A OR B when they are individually compared to the negative control
  ####if the compare_type has "A_combine_B" then will determine if taxa is contaminant when negative control is compare to both A and B together
  for (i in 1:length(compare_type)) {
    assign(paste0("compare_type",i), compare_type[[i]])
  }  
  ####getting relative abundance
  normalizeSample <- function(x){x/sum(x)}
  input_phyloseq_2_relab <- transformSampleCounts(input_phyloseq_2,normalizeSample) ##relative abundance
  ### Setting up ###
  # Counts from Phyloseq 
  counts.edit <- as.data.frame(otu_table(input_phyloseq_2))
  # Relative Abundance Table from Phyloseq 
  relab.edit <- as.data.frame(otu_table(input_phyloseq_2_relab))
  # reference table (not to be used)
  reference_table <- as.data.frame(sample_data(input_phyloseq_2))

  ##if SampleID.unique is missing, then will use the rowname of match to be the SampleID.unique given the rowname of sample_data of the phyloseq should be unique
  if (is.null(SampleID.unique)) {
    reference_table$SampleID_decontam<-rownames(reference_table) #this makes a variable in the DF match call SampleID.unique which would be rowname of the DF match
  } else {
    if (SampleID.unique %in% names(reference_table)){    
      reference_table$SampleID_decontam<-reference_table[[SampleID.unique]] #this makes a variable in the DF match call SampleID.unique which would be rowname of the DF match
    } else {stop(paste0(SampleID.unique," does not exist in the sample dataframe"))} #if the phyloseq sample dataframe doesn't have the SampleID.unique then stop
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ###adding bacterial load option- for frequency only 
  ###setting up the bacterial load data
  if (!missing(bacterial_load)){ #if want to use bacterial load
    if (!inherits(bacterial_load, c("data.frame"))) { #making sure that bacterial load is a data.frame
      stop_txt = paste0("The input data for bacterial_load need to be a data frame")
      stop(stop_txt, call. = FALSE)
    }
    if (ncol(bacterial_load)!=2) { #make sure the bacterial load dataframe is in the right structure
      stop_txt2 = paste0("The input data frame for bacterial_load needs to have two columns:first column being the SampleID.unique (or rowname of the sample_data of the input_phyloseq if SampleID.unique is not inputted)")
      stop(stop_txt2, call. = FALSE)
    }
    if(length(setdiff(reference_table$SampleID_decontam,bacterial_load[,1]))!=0) { ###if there are elements in the reference_table that is not in the bacterial_load first column, which should be 
      differentID<-paste(setdiff(reference_table$SampleID_decontam,bacterial_load[,1]), collapse=",  ")  
      stop_txt3= paste("The following SampleIDs do not exist in the bacterial_load dataframe:",
                       differentID)
      stop(stop_txt3, call. = FALSE)
    }
    colnames(bacterial_load)<-c("SampleID", "ddPCR_count")
  } 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  #######################################################################################################################################################################
  taxa.table <- as.data.frame(tax_table(input_phyloseq_2_relab))
  # This will paste the names together all the way from Family to OTU, separted by "." 
  taxa.table$match <- paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
  taxa.table.clean <- subset(taxa.table, select=c(match))
  taxa<-taxa.table$match #make a vector with all the taxa names

  # This will merge by column=0 which is the rowname by taxa-name-###
  # so now count.match and relab.match both have the new condensed taxa name 
  count.match <- merge(counts.edit, taxa.table.clean, by=0) 
  relab.match <- merge(relab.edit, taxa.table.clean, by=0)
  
  ### Replace count.match, relab.match rownames with taxa-name-###
  # This step replaces rownames for the match with count.match$match
  rownames(count.match) <- count.match$match
  rownames(relab.match) <- relab.match$match
  #remove match since match is now the rownames
  #remove Row.names given that is the byproduct of merging by rowname on previous step 
  counts <- subset(count.match, select=-c(Row.names, match))
  relab <- subset(relab.match, select=-c(Row.names, match))
  
  ### Match the Sample ID to the sample type categories  
  ### We need to do this to set up the comparison objects  
  ###replacing the SampleID.unique with the sample type for the count and relative abundance table
  index_var_name<-grep(sample_type_var_name, colnames(reference_table)) #this determine the nth column which sample+type_var_name is located at in the reference_table dataframe
  index2_var_name<-grep("SampleID_decontam", colnames(reference_table))
  #this step replace all of the unique subject ID (SampleID.unique) with their corresponding group of sample type of interest
  names(counts) <- reference_table[[index_var_name]][match(names(counts), reference_table[[index2_var_name]])] 
  names(relab) <- reference_table[[index_var_name]][match(names(relab), reference_table[[index2_var_name]])]  
  
  #### this input the sample types which will be used
  comparison.list<- vector(mode = "list", length = length(sample_types))
  names(comparison.list) <- sample_types
  for (i in 1:length(sample_types)) {
    comparison.list [i] <- list(grep(sample_types[i],colnames(counts)))
  }
  #comparison.list contain multiple lists (each list provide the column number of the specific sample type of interest). If there are 3 sample types in the dataframe of phyloseq then comparison will have a list of 3
  #counts.subset is basically counts, except there is now an index on the column name 
  counts.subset <- counts[,unlist(comparison.list)] # head(counts.subset)
  relab.subset <- relab[,unlist(comparison.list)] # head(relab.subset)
  #libsize consider the total abundance by each subject 
  libsizes <- colSums(counts.subset)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #####getting a dataframe that match the sample_type with column index with the initial sampleID
  if (!missing(bacterial_load)){ #if want to use bacterial load
    basedf<-count.match %>% select(-c("Row.names","match"))
    basedf_t<-data.frame(t(basedf))
    basedf_t$SampleID<-row.names(basedf_t)
    
    base2_df_t<-data.frame(t(counts.subset))
    base2_df_t$sample_index<-row.names(base2_df_t)
    base_combine<-merge(basedf_t,base2_df_t)
    base_combine<-base_combine %>% select(c(SampleID, sample_index))    
    
    bacterial_load_indexed<-merge(bacterial_load, base_combine, by="SampleID")
    
    bacterial_load_indexed2<-bacterial_load_indexed[order(bacterial_load_indexed$sample_index),]
    rownames(bacterial_load_indexed2)<-bacterial_load_indexed2$sample_index
    bacterial_load_indexed2<-bacterial_load_indexed2 %>% select(ddPCR_count)
    
    bacterial_load_indexed3<-as.numeric(unlist(bacterial_load_indexed2))
    names(bacterial_load_indexed3)<-rownames(bacterial_load_indexed2)  
    libsizes<-bacterial_load_indexed3 #library size is replaced with the ddPCR data
    print(paste0(bacterial_load," was used for frequency contaminant determination"))
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  ### Create matrix of data for evaluation ### 
  ### This is to create an empty matrix for all top statistics ### 
  X<-nrow(relab.subset)
  top.relab.stats <- as.data.frame(matrix(ncol=length(sample_types),nrow=X))
  
  # Moving the rownames to top.relab.stats 
  # going create two columns for medians to plot 
  # This may not be necessary top.relab.stats <- NULL 
  rownames(top.relab.stats) <- rownames(relab.subset)
  
  #loop through the sample types and get statistics for each. if there are 3 sample_types, then this will fill up the first 3 column with the statistics
  for (i in 1:length(sample_types)) {
    if (stat_option=="median") {
      top.relab.stats[i] <-rowMedians(as.matrix(relab.subset[, grepl(sample_types[i], colnames(relab.subset))]))
    } else if (stat_option=="mean") {
      top.relab.stats[i] <-rowMeans2(as.matrix(relab.subset[, grepl(sample_types[i], colnames(relab.subset))]))
    }
    colnames(top.relab.stats)[i] <- c(paste0("stats.",sample_types[i]))
  }  
  #loop through the sample types and fill in the next several column with the rank order by that particular sample type. If there are 3 sample_types, then this will give rank order on 4th-6th column
  for (i in 1:length(sample_types)) {
    ##generating the rank order by sample types 
    j <- i+length(sample_types)
    column<-top.relab.stats[,i]
    top.relab.stats <- top.relab.stats[order(column, decreasing = TRUE), , drop = FALSE ]
    top.relab.stats[j] <- 1:nrow(top.relab.stats)
    colnames(top.relab.stats)[j] <- c(paste0("rank.order.",sample_types[i]))
  }  
  #################################
  #setting up the negative control#
  #################################
  #subject-level total count for the negative control
  libsizes_negative <- libsizes[grepl(negative_sample_type, names(libsizes))]
  w<-match(negative_sample_type,sample_types)
  comparison.list_negative<- unlist(comparison.list[w]) #extract the list of negative control (this list give the column number which is consider negative control)
  comparison.list_negative<-ifelse(comparison.list_negative, T, F) #basically making this list turn into all 1 
  #for counts.subset_negative and relab.subset_negative, the columns are subject (they should only be in the negative controls)
  counts.subset_negative <- dplyr::select(counts.subset, starts_with(negative_sample_type))
  relab.subset_negative <- dplyr::select(relab.subset, starts_with(negative_sample_type))
  ###keeping only the statistics for each sample type
  upperlimitchart<-length(sample_types)*2
  top.relab.stats_table <- top.relab.stats[c(1:upperlimitchart)]
  relab.subset_final<-relab.subset
  
  for (h in 1:length(compare_type)) {
    name_compare_type<-get(paste0("compare_type",h))
    print(paste0(negative_sample_type," compare:"))
    print(name_compare_type)
    AND<-grepl("_and_",name_compare_type, fixed=T)
    OR<-grepl("_or_",name_compare_type, fixed=T)
    COMBINE<-grepl("_combine_",name_compare_type, fixed=T)
    contam<-list()
    top.relab.stats_table_new<-top.relab.stats_table
    for (test_threshold in test_threshold_all_level) {
      if ((AND==TRUE) | (OR==TRUE)) {
        if (AND==TRUE) {
          compare_sublist <- unlist(strsplit(name_compare_type, "_and_"))
        }else if (OR==TRUE) {
          compare_sublist <- unlist(strsplit(name_compare_type, "_or_"))
        }
        #loop through the two sample types 
        for (i in 1:length(compare_sublist)) { 
          assign(paste0("compare_sublist",i), compare_sublist[i])
        } 

        for (b in 1:length(compare_sublist)) {
          # only keeping the compare list  
          selected_sample_type<-get(paste0("compare_sublist",b))
          p<-match(selected_sample_type,sample_types) 
          comparison.list_select<- unlist(comparison.list[p]) 
          comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
          
          # Remove unused from conc, in this case libsizes 
          libsizes_select <- libsizes[grepl(selected_sample_type, names(libsizes))]
          counts.subset_select <- dplyr::select(counts.subset, starts_with(selected_sample_type))
          relab.subset_select <- dplyr::select(relab.subset, starts_with(selected_sample_type))
          
          counts.subset_comb<-merge(counts.subset_select,counts.subset_negative, by=0) #put back the negative control with the sample that you want to compare 
          relab.subset_comb<-merge(relab.subset_select,relab.subset_negative, by=0)
          #remove the new column "Row.names" created when merging by 0
          rownames(counts.subset_comb) <- counts.subset_comb$Row.names
          rownames(relab.subset_comb) <- relab.subset_comb$Row.names
          counts.subset_comb <- subset(counts.subset_comb, select=-c(Row.names))
          relab.subset_comb <- subset(relab.subset_comb, select=-c(Row.names))
          
          negative_comb <- append(comparison.list_select,comparison.list_negative)
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
          ###need to sort libsizes_select and libsizes_negative
          ###this is to ensure that the correct library size is link to the corresponding samples 
          libsizes_select_TEMP<-data.frame(libsizes_select)
          libsizes_select_TEMP$index<- as.numeric(sub(".*\\.", "", row.names(libsizes_select_TEMP)))
          libsizes_select_TEMP$index[is.na(libsizes_select_TEMP$index)]<-0 #replace NA with 0 
          libsizes_select_TEMP <- libsizes_select_TEMP[order(libsizes_select_TEMP$index), ]
          libsizes_select2<-as.numeric(libsizes_select_TEMP$libsizes_select)
          names(libsizes_select2)<-row.names(libsizes_select_TEMP)
          
          libsizes_negative_TEMP<-data.frame(libsizes_negative)
          libsizes_negative_TEMP$index<- as.numeric(sub(".*\\.", "", row.names(libsizes_negative_TEMP)))
          libsizes_negative_TEMP$index[is.na(libsizes_negative_TEMP$index)]<-0 #replace NA with 0 
          libsizes_negative_TEMP <- libsizes_negative_TEMP[order(libsizes_negative_TEMP$index), ]
          libsizes_negative2<-as.numeric(libsizes_negative_TEMP$libsizes_negative)
          names(libsizes_negative2)<-row.names(libsizes_negative_TEMP)
          libsizes_comb <- append(libsizes_select2,libsizes_negative2 )
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#      
          
          contam$prev <- isContaminant(t(as.matrix(counts.subset_comb)),
                                       method="prevalence",
                                       neg=negative_comb,
                                       threshold=test_threshold)
          contam$freq <- isContaminant(t(as.matrix(counts.subset_comb)),
                                       method="frequency",
                                       neg=negative_comb,
                                       conc=libsizes_comb,
                                       threshold=test_threshold)
          contam$combined <- isContaminant(t(as.matrix(counts.subset_comb)),
                                           method="combined",
                                           neg=negative_comb,
                                           conc=libsizes_comb,
                                           threshold=test_threshold)
          contam$wilcox <- wilcox_contam_test(relab.subset_comb,
                                              neg=negative_comb,
                                              threshold=test_threshold)
          contam_dataframe<- contam %>% data.frame #compress the lists from the different results
          contam_dataframe$wilcox.contaminant<-as.logical(contam_dataframe$wilcox.contaminant)
          
          contam_result<-subset(contam_dataframe, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
          assign(paste0("contam_result",b),contam_result)
        }
        #combine the results for the two sample types
        contam_result_final<-merge(contam_result1,contam_result2, by=0)
        rownames(contam_result_final) <- contam_result_final$Row.names
        contam_result_final <- subset(contam_result_final, select=-c(Row.names))
        if (AND==TRUE) {
          contam_result_final <- contam_result_final %>% mutate(prev.contaminant = case_when((prev.contaminant.x==T & prev.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(freq.contaminant = case_when((freq.contaminant.x==T & freq.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(combined.contaminant = case_when((combined.contaminant.x==T & combined.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(wilcox.contaminant = case_when((wilcox.contaminant.x==T & wilcox.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- subset(contam_result_final, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
        }else if (OR==TRUE) {
          contam_result_final <- contam_result_final %>% mutate(prev.contaminant = case_when((prev.contaminant.x==T | prev.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(freq.contaminant = case_when((freq.contaminant.x==T | freq.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(combined.contaminant = case_when((combined.contaminant.x==T | combined.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(wilcox.contaminant = case_when((wilcox.contaminant.x==T | wilcox.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))   
          contam_result_final <- subset(contam_result_final, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
        }
      }else{
        if (COMBINE==TRUE) { ###if compare_type is a combination set 
          ###this works if A_combine_B... it basically combine subset A and subset B
          combine_list <- unlist(strsplit(name_compare_type, "_combine_"))
          f<-match(combine_list,sample_types) 
          comparison.list_select<- unlist(comparison.list[f]) 
          comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
          combine_list_temp<-str_replace(name_compare_type, "_combine_", "|") 
          #Setting up the dataset so that it include both sample types that are combine =
          libsizes_select <- libsizes[grepl(combine_list_temp, names(libsizes))]
          counts.subset_select <- dplyr::select(counts.subset, starts_with(combine_list))
          relab.subset_select <- dplyr::select(relab.subset, starts_with(combine_list))
        } else { ###if compare_type does not contain _and_, _or_, _combine_.... then basically compare_type is just one of the sample types
          
          if (name_compare_type %in% sample_types) { #name_compare_type has to be one of the initial sample_types if it is not _and_, _or_, _combine_
            f<-match(name_compare_type,sample_types) 
            comparison.list_select<- unlist(comparison.list[f]) 
            comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
            #keeping only the sample type of interest
            libsizes_select <- libsizes[grepl(name_compare_type, names(libsizes))]
            counts.subset_select <- dplyr::select(counts.subset, starts_with(name_compare_type))
            relab.subset_select <- dplyr::select(relab.subset, starts_with(name_compare_type))
          } else { stop(paste0(paste0(name_compare_type," must be in sample type list:"),sample_types)) }
          
        }
        #combining with the negative control group
        counts.subset_comb<-merge(counts.subset_select,counts.subset_negative, by=0) #put back the negative control with the sample that you want to compare 
        relab.subset_comb<-merge(relab.subset_select,relab.subset_negative, by=0)
        #remove the new column "Row.names" created when merging by 0
        rownames(counts.subset_comb) <- counts.subset_comb$Row.names
        rownames(relab.subset_comb) <- relab.subset_comb$Row.names
        counts.subset_comb <- subset(counts.subset_comb, select=-c(Row.names))
        relab.subset_comb <- subset(relab.subset_comb, select=-c(Row.names))
        
        negative_comb<- append(comparison.list_select,comparison.list_negative)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        ###need to sort libsizes_select and libsizes_negative
        ###this is to ensure that the correct library size is link to the corresponding samples 
        libsizes_select_TEMP<-data.frame(libsizes_select)
        libsizes_select_TEMP$index<- as.numeric(sub(".*\\.", "", row.names(libsizes_select_TEMP)))
        libsizes_select_TEMP$index[is.na(libsizes_select_TEMP$index)]<-0 #replace NA with 0 
        libsizes_select_TEMP <- libsizes_select_TEMP[order(libsizes_select_TEMP$index), ]
        libsizes_select2<-as.numeric(libsizes_select_TEMP$libsizes_select)
        names(libsizes_select2)<-row.names(libsizes_select_TEMP)
        
        libsizes_negative_TEMP<-data.frame(libsizes_negative)
        libsizes_negative_TEMP$index<- as.numeric(sub(".*\\.", "", row.names(libsizes_negative_TEMP)))
        libsizes_negative_TEMP$index[is.na(libsizes_negative_TEMP$index)]<-0 #replace NA with 0 
        libsizes_negative_TEMP <- libsizes_negative_TEMP[order(libsizes_negative_TEMP$index), ]
        libsizes_negative2<-as.numeric(libsizes_negative_TEMP$libsizes_negative)
        names(libsizes_negative2)<-row.names(libsizes_negative_TEMP)
        libsizes_comb <- append(libsizes_select2,libsizes_negative2 )
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#      
        
        contam$prev <- isContaminant(t(as.matrix(counts.subset_comb)),
                                     method="prevalence",
                                     neg=negative_comb,
                                     threshold=test_threshold)
        contam$freq <- isContaminant(t(as.matrix(counts.subset_comb)),
                                     method="frequency",
                                     neg=negative_comb,
                                     conc=libsizes_comb,
                                     threshold=test_threshold)
        contam$combined <- isContaminant(t(as.matrix(counts.subset_comb)),
                                         method="combined",
                                         neg=negative_comb,
                                         conc=libsizes_comb,
                                         threshold=test_threshold)
        contam$wilcox <- wilcox_contam_test(relab.subset_comb,
                                            neg=negative_comb,
                                            threshold=test_threshold)
        contam_dataframe<- contam %>% data.frame #compress the lists from the different results
        contam_dataframe$wilcox.contaminant<-as.logical(contam_dataframe$wilcox.contaminant)
        
        contam_result_final<-subset(contam_dataframe, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
      }
      ### Come back from the DECONTAM section to add this in ###    
      # Adding contam data (saved to the end to add)
      ###figure out which taxa to highlight in red based on which test you choice 
      if (method_type == "MannWhit") {
        contam_result_final$contamiant <- contam_result_final$wilcox.contaminant 
        print("MannWhit")
      }
      else if (method_type == "preval") {
        contam_result_final$contamiant <- contam_result_final$prev.contaminant 
        print("preval")
      }
      else if (method_type == "freq") {
        contam_result_final$contamiant <- contam_result_final$freq.contaminant 
        print("freq")
      }
      else if (method_type == "combin") {
        contam_result_final$contamiant <- contam_result_final$combined.contaminant
        print("combin")
      }
      contam_result_final <- subset(contam_result_final, select=c("contamiant"))
      print(test_threshold)
      colnames(contam_result_final)<- paste0("level_",test_threshold)
      top.relab.stats_table_new<-merge(top.relab.stats_table_new,contam_result_final, by="row.names")
      row.names(top.relab.stats_table_new)<-top.relab.stats_table_new$Row.names #previous merge step generate a new column call Row.names
      top.relab.stats_table_new<- subset(top.relab.stats_table_new,select=-c(Row.names))
    }
    
    ####export csv
    csv_output_name<-paste0(output_excel,paste0(paste0("contam",paste0(paste0('_negC_',negative_sample_type),paste0("_compare_",name_compare_type))),paste0(paste0("_method_",method_type),"_all_level.csv")))
    write.csv(top.relab.stats_table_new,csv_output_name)
    
    #keeping only the levels columns
    top.relab.stats_table_new2<-top.relab.stats_table_new[, grepl(paste0("^", "level_"), names(top.relab.stats_table_new))]
    
    ###export results as a dataframe
    assign(paste0(output_excel,paste0(paste0("contam",paste0(paste0('_negC_',negative_sample_type),paste0("_compare_",name_compare_type))),paste0(paste0("_method_",method_type),"_all_level"))), top.relab.stats_table_new2,.GlobalEnv)
    print("###########################")
    print("Available contaminant list all level:")
    print(paste0(output_excel,paste0(paste0("contam",paste0(paste0('_negC_',negative_sample_type),paste0("_compare_",name_compare_type))),paste0(paste0("_method_",method_type),"_all_level"))))
    print("###########################")
    
  }
}

####decontaminant_sensitivity_KW: 
decontaminant_sensitivity_KW <- function (input_phyloseq, 
                                      SampleID.unique=NULL, #if empty, SampleID.unique would be the rowname of the sample_data of the input_phyloseq
                                      sample_type_var_name, 
                                      sample_types=list(), 
                                      negative_sample_type, ###the sample type that you want to be as negative control
                                      compare_type=list(), 
                                      bacterial_load, #for frequency only. uses ddPCR to determine library size
                                      stat_option=c("mean", "median"), ### the statistics to determine rank for the boxplot
                                      test_threshold_all_level=list(), ### test_threshold_all_level=c(0.1,0.3,0.5,0.7,0.9) would test for the level specified and create an excel spread sheet with the parameters
                                      expected_contam_taxa=list(), ### list of taxa you expect to be contaminant
                                      expected_NOT_contam_taxa=list(), ### list of taxa you do NOT expect to be contaminant
                                      output_excel)  { 
  
  ################################################################################################################################################################################
  ##############################################################setting up all the function argument and data frames############################################################## 
  ################################################################################################################################################################################
  ###ensure intput_phyloseq is a absolute count phyloseq, not relative abundance phyloseq. and make sure phyloseq is in right orientation
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  countData_check <- floor(as(otu_table(input_phyloseq), "matrix")) # round down to nearest integer. if this phsyloeq is a relative abundance, then the entire countData_check would be 0
  if (all(countData_check==0)) {stop("Please use phyloseq with absolute count. The current input physloeq contains relative abundnace")}
  
  ## evaluate choices
  stat_option <- match.arg(stat_option)
  print(stat_option)

  ###generating a variable for each sample type 
  for (i in 1:length(sample_types)) {
    assign(paste0("sampletype",i), sample_types[[i]])
  }
  
  ###this make sure that the phyloseq ONLY contains the sample_type_var_name variable with all options listed sample_types
  #in theory input_phyloseq should be same as input_phyloseq_2 if the user input all of the sample_types options
  keep_sample<- get_variable(input_phyloseq,sample_type_var_name) %in% sample_types
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  input_phyloseq_2 <- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
  
  ####preparing compare_type list
  ####if the compare_type has "A_and_B" then will determine if taxa is contaminant for A and B when they are individually compared to the negative control
  ####if the compare_type has "A_or_B" then will determine if taxa is contaminant for A OR B when they are individually compared to the negative control
  ####if the compare_type has "A_combine_B" then will determine if taxa is contaminant when negative control is compare to both A and B together
  for (i in 1:length(compare_type)) {
    assign(paste0("compare_type",i), compare_type[[i]])
  }  
  ####getting relative abundance
  normalizeSample <- function(x){x/sum(x)}
  input_phyloseq_2_relab <- transformSampleCounts(input_phyloseq_2,normalizeSample) ##relative abundance
  ### Setting up ###
  # Counts from Phyloseq 
  counts.edit <- as.data.frame(otu_table(input_phyloseq_2))
  # Relative Abundance Table from Phyloseq 
  relab.edit <- as.data.frame(otu_table(input_phyloseq_2_relab))
  # reference table (not to be used)
  reference_table <- as.data.frame(sample_data(input_phyloseq_2))
  
  ##if SampleID.unique is missing, then will use the rowname of match to be the SampleID.unique given the rowname of sample_data of the phyloseq should be unique
  if (is.null(SampleID.unique)) {
    reference_table$SampleID_decontam<-rownames(reference_table) #this makes a variable in the DF match call SampleID.unique which would be rowname of the DF match
  } else {
    if (SampleID.unique %in% names(reference_table)){    
      reference_table$SampleID_decontam<-reference_table[[SampleID.unique]] #this makes a variable in the DF match call SampleID.unique which would be rowname of the DF match
    } else {stop(paste0(SampleID.unique," does not exist in the sample dataframe"))} #if the phyloseq sample dataframe doesn't have the SampleID.unique then stop
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ###adding bacterial load option- for frequency only 
  ###setting up the bacterial load data
  if (!missing(bacterial_load)){ #if want to use bacterial load
    if (!inherits(bacterial_load, c("data.frame"))) { #making sure that bacterial load is a data.frame
      stop_txt = paste0("The input data for bacterial_load need to be a data frame")
      stop(stop_txt, call. = FALSE)
    }
    if (ncol(bacterial_load)!=2) { #make sure the bacterial load dataframe is in the right structure
      stop_txt2 = paste0("The input data frame for bacterial_load needs to have two columns:first column being the SampleID.unique (or rowname of the sample_data of the input_phyloseq if SampleID.unique is not inputted)")
      stop(stop_txt2, call. = FALSE)
    }
    if(length(setdiff(reference_table$SampleID_decontam,bacterial_load[,1]))!=0) { ###if there are elements in the reference_table that is not in the bacterial_load first column, which should be 
      differentID<-paste(setdiff(reference_table$SampleID_decontam,bacterial_load[,1]), collapse=",  ")  
      stop_txt3= paste("The following SampleIDs do not exist in the bacterial_load dataframe:",
                       differentID)
      stop(stop_txt3, call. = FALSE)
    }
    colnames(bacterial_load)<-c("SampleID", "ddPCR_count")
  } 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  #######################################################################################################################################################################
  taxa.table <- as.data.frame(tax_table(input_phyloseq_2_relab))
  # This will paste the names together all the way from Family to OTU, separted by "." 
  taxa.table$match <- paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
  taxa.table.clean <- subset(taxa.table, select=c(match))
  taxa<-taxa.table$match #make a vector with all the taxa names
  
  # This will merge by column=0 which is the rowname by taxa-name-###
  # so now count.match and relab.match both have the new condensed taxa name 
  count.match <- merge(counts.edit, taxa.table.clean, by=0) 
  relab.match <- merge(relab.edit, taxa.table.clean, by=0)
  
  ### Replace count.match, relab.match rownames with taxa-name-###
  # This step replaces rownames for the match with count.match$match
  rownames(count.match) <- count.match$match
  rownames(relab.match) <- relab.match$match
  #remove match since match is now the rownames
  #remove Row.names given that is the byproduct of merging by rowname on previous step 
  counts <- subset(count.match, select=-c(Row.names, match))
  relab <- subset(relab.match, select=-c(Row.names, match))
  
  ### Match the Sample ID to the sample type categories  
  ### We need to do this to set up the comparison objects  
  ###replacing the SampleID.unique with the sample type for the count and relative abundance table
  index_var_name<-grep(sample_type_var_name, colnames(reference_table)) #this determine the nth column which sample+type_var_name is located at in the reference_table dataframe
  index2_var_name<-grep("SampleID_decontam", colnames(reference_table))
  #this step replace all of the unique subject ID (SampleID.unique) with their corresponding group of sample type of interest
  names(counts) <- reference_table[[index_var_name]][match(names(counts), reference_table[[index2_var_name]])] 
  names(relab) <- reference_table[[index_var_name]][match(names(relab), reference_table[[index2_var_name]])]  
  
  #### this input the sample types which will be used
  comparison.list<- vector(mode = "list", length = length(sample_types))
  names(comparison.list) <- sample_types
  for (i in 1:length(sample_types)) {
    comparison.list [i] <- list(grep(sample_types[i],colnames(counts)))
  }
  #comparison.list contain multiple lists (each list provide the column number of the specific sample type of interest). If there are 3 sample types in the dataframe of phyloseq then comparison will have a list of 3
  #counts.subset is basically counts, except there is now an index on the column name 
  counts.subset <- counts[,unlist(comparison.list)] # head(counts.subset)
  relab.subset <- relab[,unlist(comparison.list)] # head(relab.subset)
  #libsize consider the total abundance by each subject 
  libsizes <- colSums(counts.subset)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #####getting a dataframe that match the sample_type with column index with the initial sampleID
  if (!missing(bacterial_load)){ #if want to use bacterial load
    basedf<-count.match %>% select(-c("Row.names","match"))
    basedf_t<-data.frame(t(basedf))
    basedf_t$SampleID<-row.names(basedf_t)
    
    base2_df_t<-data.frame(t(counts.subset))
    base2_df_t$sample_index<-row.names(base2_df_t)
    base_combine<-merge(basedf_t,base2_df_t)
    base_combine<-base_combine %>% select(c(SampleID, sample_index))    
    
    bacterial_load_indexed<-merge(bacterial_load, base_combine, by="SampleID")
    
    bacterial_load_indexed2<-bacterial_load_indexed[order(bacterial_load_indexed$sample_index),]
    rownames(bacterial_load_indexed2)<-bacterial_load_indexed2$sample_index
    bacterial_load_indexed2<-bacterial_load_indexed2 %>% select(ddPCR_count)
    
    bacterial_load_indexed3<-as.numeric(unlist(bacterial_load_indexed2))
    names(bacterial_load_indexed3)<-rownames(bacterial_load_indexed2)  
    libsizes<-bacterial_load_indexed3 #library size is replaced with the ddPCR data
    print(paste0(bacterial_load," was used for frequency contaminant determination"))
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  ### Create matrix of data for evaluation ### 
  ### This is to create an empty matrix for all top statistics ### 
  X<-nrow(relab.subset)
  top.relab.stats <- as.data.frame(matrix(ncol=length(sample_types),nrow=X))
  
  # Moving the rownames to top.relab.stats 
  # going create two columns for medians to plot 
  # This may not be necessary top.relab.stats <- NULL 
  rownames(top.relab.stats) <- rownames(relab.subset)
  
  #loop through the sample types and get statistics for each. if there are 3 sample_types, then this will fill up the first 3 column with the statistics
  for (i in 1:length(sample_types)) {
    if (stat_option=="median") {
      top.relab.stats[i] <-rowMedians(as.matrix(relab.subset[, grepl(sample_types[i], colnames(relab.subset))]))
    } else if (stat_option=="mean") {
      top.relab.stats[i] <-rowMeans2(as.matrix(relab.subset[, grepl(sample_types[i], colnames(relab.subset))]))
    }
    colnames(top.relab.stats)[i] <- c(paste0("stats.",sample_types[i]))
  }  
  #loop through the sample types and fill in the next several column with the rank order by that particular sample type. If there are 3 sample_types, then this will give rank order on 4th-6th column
  for (i in 1:length(sample_types)) {
    ##generating the rank order by sample types 
    j <- i+length(sample_types)
    column<-top.relab.stats[,i]
    top.relab.stats <- top.relab.stats[order(column, decreasing = TRUE), , drop = FALSE ]
    top.relab.stats[j] <- 1:nrow(top.relab.stats)
    colnames(top.relab.stats)[j] <- c(paste0("rank.order.",sample_types[i]))
  }  
  #################################
  #setting up the negative control#
  #################################
  #subject-level total count for the negative control
  libsizes_negative <- libsizes[grepl(negative_sample_type, names(libsizes))]
  w<-match(negative_sample_type,sample_types)
  comparison.list_negative<- unlist(comparison.list[w]) #extract the list of negative control (this list give the column number which is consider negative control)
  comparison.list_negative<-ifelse(comparison.list_negative, T, F) #basically making this list turn into all 1 
  #for counts.subset_negative and relab.subset_negative, the columns are subject (they should only be in the negative controls)
  counts.subset_negative <- dplyr::select(counts.subset, starts_with(negative_sample_type))
  relab.subset_negative <- dplyr::select(relab.subset, starts_with(negative_sample_type))
  ###keeping only the statistics for each sample type
  upperlimitchart<-length(sample_types)*2
  top.relab.stats_table <- top.relab.stats[c(1:upperlimitchart)]
  relab.subset_final<-relab.subset
  
  for (h in 1:length(compare_type)) {
    name_compare_type<-get(paste0("compare_type",h))
    print(paste0(negative_sample_type," compare:"))
    print(name_compare_type)
    AND<-grepl("_and_",name_compare_type, fixed=T)
    OR<-grepl("_or_",name_compare_type, fixed=T)
    COMBINE<-grepl("_combine_",name_compare_type, fixed=T)
    contam<-list()
    top.relab.stats_table_new<-top.relab.stats_table
    test_num<-1
    for (test_threshold in test_threshold_all_level) {
      if ((AND==TRUE) | (OR==TRUE)) {
        if (AND==TRUE) {
          compare_sublist <- unlist(strsplit(name_compare_type, "_and_"))
        }else if (OR==TRUE) {
          compare_sublist <- unlist(strsplit(name_compare_type, "_or_"))
        }
        #loop through the two sample types 
        for (i in 1:length(compare_sublist)) { 
          assign(paste0("compare_sublist",i), compare_sublist[i])
        } 
        
        for (b in 1:length(compare_sublist)) {
          # only keeping the compare list  
          selected_sample_type<-get(paste0("compare_sublist",b))
          p<-match(selected_sample_type,sample_types) 
          comparison.list_select<- unlist(comparison.list[p]) 
          comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
          
          # Remove unused from conc, in this case libsizes 
          libsizes_select <- libsizes[grepl(selected_sample_type, names(libsizes))]
          counts.subset_select <- dplyr::select(counts.subset, starts_with(selected_sample_type))
          relab.subset_select <- dplyr::select(relab.subset, starts_with(selected_sample_type))
          
          counts.subset_comb<-merge(counts.subset_select,counts.subset_negative, by=0) #put back the negative control with the sample that you want to compare 
          relab.subset_comb<-merge(relab.subset_select,relab.subset_negative, by=0)
          #remove the new column "Row.names" created when merging by 0
          rownames(counts.subset_comb) <- counts.subset_comb$Row.names
          rownames(relab.subset_comb) <- relab.subset_comb$Row.names
          counts.subset_comb <- subset(counts.subset_comb, select=-c(Row.names))
          relab.subset_comb <- subset(relab.subset_comb, select=-c(Row.names))
          
          negative_comb <- append(comparison.list_select,comparison.list_negative)
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
          ###need to sort libsizes_select and libsizes_negative
          ###this is to ensure that the correct library size is link to the corresponding samples 
          libsizes_select_TEMP<-data.frame(libsizes_select)
          libsizes_select_TEMP$index<- as.numeric(sub(".*\\.", "", row.names(libsizes_select_TEMP)))
          libsizes_select_TEMP$index[is.na(libsizes_select_TEMP$index)]<-0 #replace NA with 0 
          libsizes_select_TEMP <- libsizes_select_TEMP[order(libsizes_select_TEMP$index), ]
          libsizes_select2<-as.numeric(libsizes_select_TEMP$libsizes_select)
          names(libsizes_select2)<-row.names(libsizes_select_TEMP)
          
          libsizes_negative_TEMP<-data.frame(libsizes_negative)
          libsizes_negative_TEMP$index<- as.numeric(sub(".*\\.", "", row.names(libsizes_negative_TEMP)))
          libsizes_negative_TEMP$index[is.na(libsizes_negative_TEMP$index)]<-0 #replace NA with 0 
          libsizes_negative_TEMP <- libsizes_negative_TEMP[order(libsizes_negative_TEMP$index), ]
          libsizes_negative2<-as.numeric(libsizes_negative_TEMP$libsizes_negative)
          names(libsizes_negative2)<-row.names(libsizes_negative_TEMP)
          libsizes_comb <- append(libsizes_select2,libsizes_negative2 )
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#      
          
          #generate output in global environment to check
          assign(paste0("counts.checking",b),counts.subset_comb, envir=.GlobalEnv)
          assign(paste0("relab.checking",b),relab.subset_comb,)
          assign(paste0("negative_comb.checking",b),negative_comb,envir=.GlobalEnv)
          assign(paste0("libsizes_comb.checking",b),libsizes_comb,envir=.GlobalEnv)
          
          contam$prev <- isContaminant(t(as.matrix(counts.subset_comb)),
                                       method="prevalence",
                                       neg=negative_comb,
                                       threshold=test_threshold)
          contam$freq <- isContaminant(t(as.matrix(counts.subset_comb)),
                                       method="frequency",
                                       neg=negative_comb,
                                       conc=libsizes_comb,
                                       threshold=test_threshold)
          contam$combined <- isContaminant(t(as.matrix(counts.subset_comb)),
                                           method="combined",
                                           neg=negative_comb,
                                           conc=libsizes_comb,
                                           threshold=test_threshold)
          contam$wilcox <- wilcox_contam_test(relab.subset_comb,
                                              neg=negative_comb,
                                              threshold=test_threshold)
          contam_dataframe<- contam %>% data.frame #compress the lists from the different results
          contam_dataframe$wilcox.contaminant<-as.logical(contam_dataframe$wilcox.contaminant)
          
          contam_result<-subset(contam_dataframe, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
          assign(paste0("contam_result",b),contam_result)
        }
        #combine the results for the two sample types
        contam_result_final<-merge(contam_result1,contam_result2, by=0)
        rownames(contam_result_final) <- contam_result_final$Row.names
        contam_result_final <- subset(contam_result_final, select=-c(Row.names))
        if (AND==TRUE) {
          contam_result_final <- contam_result_final %>% mutate(prev.contaminant = case_when((prev.contaminant.x==T & prev.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(freq.contaminant = case_when((freq.contaminant.x==T & freq.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(combined.contaminant = case_when((combined.contaminant.x==T & combined.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(wilcox.contaminant = case_when((wilcox.contaminant.x==T & wilcox.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- subset(contam_result_final, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
        }else if (OR==TRUE) {
          contam_result_final <- contam_result_final %>% mutate(prev.contaminant = case_when((prev.contaminant.x==T | prev.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(freq.contaminant = case_when((freq.contaminant.x==T | freq.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(combined.contaminant = case_when((combined.contaminant.x==T | combined.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
          contam_result_final <- contam_result_final %>% mutate(wilcox.contaminant = case_when((wilcox.contaminant.x==T | wilcox.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))   
          contam_result_final <- subset(contam_result_final, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
        }
      }else{
          if (COMBINE==TRUE) { ###if compare_type is a combination set 
            ###this works if A_combine_B... it basically combine subset A and subset B
            combine_list <- unlist(strsplit(name_compare_type, "_combine_"))
            f<-match(combine_list,sample_types) 
            comparison.list_select<- unlist(comparison.list[f]) 
            comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
            combine_list_temp<-str_replace(name_compare_type, "_combine_", "|") 
            #Setting up the dataset so that it include both sample types that are combine =
            libsizes_select <- libsizes[grepl(combine_list_temp, names(libsizes))]
            counts.subset_select <- dplyr::select(counts.subset, starts_with(combine_list))
            relab.subset_select <- dplyr::select(relab.subset, starts_with(combine_list))
          } else { ###if compare_type does not contain _and_, _or_, _combine_.... then basically compare_type is just one of the sample types
            
            if (name_compare_type %in% sample_types) { #name_compare_type has to be one of the initial sample_types if it is not _and_, _or_, _combine_
              f<-match(name_compare_type,sample_types) 
              comparison.list_select<- unlist(comparison.list[f]) 
              comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
              #keeping only the sample type of interest
              libsizes_select <- libsizes[grepl(name_compare_type, names(libsizes))]
              counts.subset_select <- dplyr::select(counts.subset, starts_with(name_compare_type))
              relab.subset_select <- dplyr::select(relab.subset, starts_with(name_compare_type))
            } else { stop(paste0(paste0(name_compare_type," must be in sample type list:"),sample_types)) }
            
          }
          #combining with the negative control group
          counts.subset_comb<-merge(counts.subset_select,counts.subset_negative, by=0) #put back the negative control with the sample that you want to compare 
          relab.subset_comb<-merge(relab.subset_select,relab.subset_negative, by=0)
          #remove the new column "Row.names" created when merging by 0
          rownames(counts.subset_comb) <- counts.subset_comb$Row.names
          rownames(relab.subset_comb) <- relab.subset_comb$Row.names
          counts.subset_comb <- subset(counts.subset_comb, select=-c(Row.names))
          relab.subset_comb <- subset(relab.subset_comb, select=-c(Row.names))
          
          negative_comb<- append(comparison.list_select,comparison.list_negative)
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
          ###need to sort libsizes_select and libsizes_negative
          ###this is to ensure that the correct library size is link to the corresponding samples 
          libsizes_select_TEMP<-data.frame(libsizes_select)
          libsizes_select_TEMP$index<- as.numeric(sub(".*\\.", "", row.names(libsizes_select_TEMP)))
          libsizes_select_TEMP$index[is.na(libsizes_select_TEMP$index)]<-0 #replace NA with 0 
          libsizes_select_TEMP <- libsizes_select_TEMP[order(libsizes_select_TEMP$index), ]
          libsizes_select2<-as.numeric(libsizes_select_TEMP$libsizes_select)
          names(libsizes_select2)<-row.names(libsizes_select_TEMP)
          
          libsizes_negative_TEMP<-data.frame(libsizes_negative)
          libsizes_negative_TEMP$index<- as.numeric(sub(".*\\.", "", row.names(libsizes_negative_TEMP)))
          libsizes_negative_TEMP$index[is.na(libsizes_negative_TEMP$index)]<-0 #replace NA with 0 
          libsizes_negative_TEMP <- libsizes_negative_TEMP[order(libsizes_negative_TEMP$index), ]
          libsizes_negative2<-as.numeric(libsizes_negative_TEMP$libsizes_negative)
          names(libsizes_negative2)<-row.names(libsizes_negative_TEMP)
          libsizes_comb <- append(libsizes_select2,libsizes_negative2 )
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#      
          
          #generate output in global environment to check
          counts.checking<<-counts.subset_comb
          relab.checking<<-relab.subset_comb
          negative_comb.checking<<-negative_comb
          libsizes_comb.checking<<-libsizes_comb
          
          contam$prev <- isContaminant(t(as.matrix(counts.subset_comb)),
                                       method="prevalence",
                                       neg=negative_comb,
                                       threshold=test_threshold)
          contam$freq <- isContaminant(t(as.matrix(counts.subset_comb)),
                                       method="frequency",
                                       neg=negative_comb,
                                       conc=libsizes_comb,
                                       threshold=test_threshold)
          contam$combined <- isContaminant(t(as.matrix(counts.subset_comb)),
                                           method="combined",
                                           neg=negative_comb,
                                           conc=libsizes_comb,
                                           threshold=test_threshold)
          contam$wilcox <- wilcox_contam_test(relab.subset_comb,
                                              neg=negative_comb,
                                              threshold=test_threshold)
          contam_dataframe<- contam %>% data.frame #compress the lists from the different results
          contam_dataframe$wilcox.contaminant<-as.logical(contam_dataframe$wilcox.contaminant)
          
          contam_result_final<-subset(contam_dataframe, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
      }
      
      #### add the test level as suffix to the contaminant variables
      colnames(contam_result_final)<-paste(colnames(contam_result_final),test_threshold,sep="_")
      assign(paste0("contam_result_final",test_num),contam_result_final)
      test_num<-test_num+1
    }

    #merge all of the contaminant results with different testing level
    contam_result_final_combine<-contam_result_final1
    for (u in 2:length(test_threshold_all_level)) {
      contam_result_final_combine<-merge(contam_result_final_combine,get(paste0("contam_result_final",u)),by="row.names")
      row.names(contam_result_final_combine)<-contam_result_final_combine$Row.names #previous merge step generate a new column call Row.names
      contam_result_final_combine<- subset(contam_result_final_combine,select=-c(Row.names))
    }
    top.relab.stats_table_combine<-merge(top.relab.stats_table_new,contam_result_final_combine, by="row.names")
    row.names(top.relab.stats_table_combine)<-top.relab.stats_table_combine$Row.names #previous merge step generate a new column call Row.names
    top.relab.stats_table_combine<- subset(top.relab.stats_table_combine,select=-c(Row.names))
    top.relab.stats_table_combine<<-top.relab.stats_table_combine  
    ####export csv
    csv_output_name<-paste0(output_excel,paste0(paste0("contam",paste0(paste0('_negC_',negative_sample_type),paste0("_compare_",name_compare_type))),"_all_level_SENSITIVITY.csv"))
    write.csv(top.relab.stats_table_combine,csv_output_name)
    ###################################################################################################
    ###################################################################################################
    ###################################################################################################
    ###################################################################################################
    ###########################################generate plot###########################################
    ###################################################################################################
    ###################################################################################################
    ###################################################################################################
    ###################################################################################################
    top.relab.stats_table_combine$taxa_full_name<-row.names(top.relab.stats_table_combine)
    total_column<-length(test_threshold_all_level)*4
    stats_column<-length(sample_types)*2
    
        ################################
        ######EXPECTED TO BE CONTAMINANT
        ################################
        #this step keep the taxa which are listed in expected_contam_taxa. And then it check to see if the sample is labeled as TRUE for each method (identify as contaminant)
        #It then goes through each taxa and calcuates a percentage of correctly identified contaminant
        top.relab.stats_table_combine$expected<-NULL
        for (taxa in expected_contam_taxa) {
          top.relab.stats_table_combine[,taxa]<-grepl(taxa, top.relab.stats_table_combine$taxa_full_name , fixed = TRUE)
          top.relab.stats_table_combine$expected[top.relab.stats_table_combine[,taxa]==TRUE]<-TRUE
        }
        top.relab.stats_table_combine$expected_NO<-NULL
        for (taxa2 in expected_NOT_contam_taxa) {
          top.relab.stats_table_combine[,taxa2]<-grepl(taxa2, top.relab.stats_table_combine$taxa_full_name , fixed = TRUE)
          top.relab.stats_table_combine$expected_NO[top.relab.stats_table_combine[,taxa2]==TRUE]<-TRUE
        }
        top.relab.stats_table_combine_t<-top.relab.stats_table_combine
        top.relab.stats_table_combine_expected<-top.relab.stats_table_combine_t[(!is.na(top.relab.stats_table_combine_t$expected)), ]
        top.relab.stats_table_combine_expected_org<-top.relab.stats_table_combine_expected
        for (level in test_threshold_all_level){ 
          for (method in c("prev.contaminant_", "freq.contaminant_","combined.contaminant_","wilcox.contaminant_")){
            top.relab.stats_table_combine_expected[,paste0(method,level)]<-as.integer(as.logical(top.relab.stats_table_combine_expected[,paste0(method,level)]))
          }
        }
        top.relab.stats_table_combine_expected <- top.relab.stats_table_combine_expected[, -c(1:stats_column)]# remove the relab abundance and rank columns
        top.relab.stats_table_combine_expected <- top.relab.stats_table_combine_expected[, c(1:total_column)] #keeping only the testings results 
        top.relab.stats_table_combine_total_expected<-top.relab.stats_table_combine_expected %>% dplyr::summarise_all(funs(mean))
        top.relab.stats_table_combine_total_expected$taxa_name<-"Contaminant"
        #for individual taxa of interest that is in the expected contaminant group
        for (taxa in expected_contam_taxa){
          top.relab.stats_table_combine_taxa<-top.relab.stats_table_combine_expected_org[(top.relab.stats_table_combine_expected_org[,taxa]), ]
          for (level in test_threshold_all_level){ 
            for (method in c("prev.contaminant_", "freq.contaminant_","combined.contaminant_","wilcox.contaminant_")){
              top.relab.stats_table_combine_taxa[,paste0(method,level)]<-as.integer(as.logical(top.relab.stats_table_combine_taxa[,paste0(method,level)]))
            }
          }
          top.relab.stats_table_combine_taxa <- top.relab.stats_table_combine_taxa[, -c(1:stats_column)]# remove the relab abundance and rank columns
          top.relab.stats_table_combine_taxa <- top.relab.stats_table_combine_taxa[, c(1:total_column)] #keeping only the testings results 
          top.relab.stats_table_combine_taxa <- top.relab.stats_table_combine_taxa %>% dplyr::summarise_all(funs(mean))
          top.relab.stats_table_combine_taxa$taxa_name<-taxa
          assign(paste0("top.relab.stats_table_combine_taxa",taxa),top.relab.stats_table_combine_taxa)
        }
        
        ####################################
        ######EXPECTED NOT TO BE CONTAMINANT
        ####################################
        #this step keep the taxa which are listed in expected_NOT_contam_taxa And then it check to see if the sample is labeled as FALSE for each method (identify as NOT being contaminant)
        #It then goes through each taxa and calculates a percentage of correctly identified NON-contaminant
        top.relab.stats_table_combine_expected_NO<-top.relab.stats_table_combine_t[(!is.na(top.relab.stats_table_combine_t$expected_NO) ), ]
        top.relab.stats_table_combine_expected_NO_org<-top.relab.stats_table_combine_expected_NO
        for (level in test_threshold_all_level){ 
          for (method in c("prev.contaminant_", "freq.contaminant_","combined.contaminant_","wilcox.contaminant_")){
            top.relab.stats_table_combine_expected_NO[,paste0(method,level)]<-as.integer(as.logical(top.relab.stats_table_combine_expected_NO[,paste0(method,level)])) #see how many are not label as contaminant for the ones which SHOULD NOT be contaminant
          }
        }
        top.relab.stats_table_combine_expected_NO <- top.relab.stats_table_combine_expected_NO[, -c(1:stats_column)]# remove the relab abundance and rank columns
        top.relab.stats_table_combine_expected_NO <- top.relab.stats_table_combine_expected_NO[, c(1:total_column)] #keeping only the testings results 
        top.relab.stats_table_combine_total_expected_NO<-top.relab.stats_table_combine_expected_NO %>% dplyr::summarise_all(funs(mean))
        for (level in test_threshold_all_level){ 
          for (method in c("prev.contaminant_", "freq.contaminant_","combined.contaminant_","wilcox.contaminant_")){
            top.relab.stats_table_combine_total_expected_NO[,paste0(method,level)]<-1-top.relab.stats_table_combine_total_expected_NO[,paste0(method,level)] #see how many are not label as contaminant for the ones which SHOULD NOT be contaminant
          }
        }
        top.relab.stats_table_combine_total_expected_NO$taxa_name<-"Non-contaminant"
        #for individual taxa of interest that is in the NOT expected contaminant group
        for (taxa in expected_NOT_contam_taxa){
          top.relab.stats_table_NO_taxa<-top.relab.stats_table_combine_expected_NO_org[(top.relab.stats_table_combine_expected_NO_org[,taxa]), ]
          for (level in test_threshold_all_level){ 
            for (method in c("prev.contaminant_", "freq.contaminant_","combined.contaminant_","wilcox.contaminant_")){
              top.relab.stats_table_NO_taxa[,paste0(method,level)]<-as.integer(as.logical(top.relab.stats_table_NO_taxa[,paste0(method,level)]))
            }
          }
          top.relab.stats_table_NO_taxa <- top.relab.stats_table_NO_taxa[, -c(1:stats_column)]# remove the relab abundance and rank columns
          top.relab.stats_table_NO_taxa <- top.relab.stats_table_NO_taxa[, c(1:total_column)] #keeping only the testings results 
          top.relab.stats_table_NO_taxa <- top.relab.stats_table_NO_taxa %>% dplyr::summarise_all(funs(mean))
          for (level in test_threshold_all_level){ 
            for (method in c("prev.contaminant_", "freq.contaminant_","combined.contaminant_","wilcox.contaminant_")){
              top.relab.stats_table_NO_taxa[,paste0(method,level)]<-1-top.relab.stats_table_NO_taxa[,paste0(method,level)]
            }
          }
          top.relab.stats_table_NO_taxa$taxa_name<-taxa
          assign(paste0("top.relab.stats_table_NO_taxa",taxa),top.relab.stats_table_NO_taxa)
        }
        #####combine all the files
        TOTAL_dataframe<-top.relab.stats_table_combine_total_expected
        for (taxa in expected_contam_taxa){
          TOTAL_dataframe<-rbind(TOTAL_dataframe, get(paste0("top.relab.stats_table_combine_taxa",taxa)))
        }
        TOTAL_dataframe<-rbind(TOTAL_dataframe, top.relab.stats_table_combine_total_expected_NO)
        for (taxa in expected_NOT_contam_taxa){
          TOTAL_dataframe<-rbind(TOTAL_dataframe, get(paste0("top.relab.stats_table_NO_taxa",taxa)))
        }
        #change wide to long for each method and stack them to graph 
        TOTAL_prev<-TOTAL_dataframe %>% dplyr::select(starts_with(c("prev","taxa")))
        TOTAL_prev_t<-reshape2::melt(TOTAL_prev, id.vars=c("taxa_name"))
        TOTAL_prev_t$variable<-as.numeric(str_sub(TOTAL_prev_t$variable, -3,-1))
        TOTAL_prev_t$method<-"prev"
        
        TOTAL_freq<-TOTAL_dataframe %>% dplyr::select(starts_with(c("freq","taxa")))
        TOTAL_freq_t<-reshape2::melt(TOTAL_freq, id.vars=c("taxa_name"))
        TOTAL_freq_t$variable<-as.numeric(str_sub(TOTAL_freq_t$variable, -3,-1))
        TOTAL_freq_t$method<-"freq"
        
        TOTAL_combined<-TOTAL_dataframe %>% dplyr::select(starts_with(c("combined","taxa")))
        TOTAL_combined_t<-reshape2::melt(TOTAL_combined, id.vars=c("taxa_name"))
        TOTAL_combined_t$variable<-as.numeric(str_sub(TOTAL_combined_t$variable, -3,-1))
        TOTAL_combined_t$method<-"combined"
        
        TOTAL_wilcox<-TOTAL_dataframe %>% dplyr::select(starts_with(c("wilcox","taxa")))
        TOTAL_wilcox_t<-reshape2::melt(TOTAL_wilcox, id.vars=c("taxa_name"))
        TOTAL_wilcox_t$variable<-as.numeric(str_sub(TOTAL_wilcox_t$variable, -3,-1))
        TOTAL_wilcox_t$method<-"wilcox"
        
        #########################################################################################################################################################
        #########################################################################################################################################################
        #####determine percent of contaminant identify (this step go back to the initial dataset with ALL taxa and then see how many taxa is identified as contaminant to get a percentage)
                  top.relab.stats_table_combine_TOTAL <- top.relab.stats_table_combine[, -c(1:stats_column)]# remove the relab abundance and rank columns
                  top.relab.stats_table_combine_TOTAL <- top.relab.stats_table_combine_TOTAL[, c(1:total_column)] #keeping only the testings results 
                  top.relab.stats_table_combine_TOTAL<-top.relab.stats_table_combine_TOTAL %>% dplyr::summarise_all(funs(mean))
                  top.relab.stats_table_combine_TOTAL$taxa_name2<-"TOTAL"
                  
                  TOTAL2_prev<-top.relab.stats_table_combine_TOTAL %>% dplyr::select(starts_with(c("prev","taxa")))
                  TOTAL2_prev_t<-reshape2::melt(TOTAL2_prev, id.vars=c("taxa_name2"))
                  TOTAL2_prev_t$variable<-as.numeric(str_sub(TOTAL2_prev_t$variable, -3,-1))
                  TOTAL2_prev_t$method<-"prev"
                  
                  TOTAL2_freq<-top.relab.stats_table_combine_TOTAL %>% dplyr::select(starts_with(c("freq","taxa")))
                  TOTAL2_freq_t<-reshape2::melt(TOTAL2_freq, id.vars=c("taxa_name2"))
                  TOTAL2_freq_t$variable<-as.numeric(str_sub(TOTAL2_freq_t$variable, -3,-1))
                  TOTAL2_freq_t$method<-"freq"
                  
                  TOTAL2_combined<-top.relab.stats_table_combine_TOTAL %>% dplyr::select(starts_with(c("combined","taxa")))
                  TOTAL2_combined_t<-reshape2::melt(TOTAL2_combined, id.vars=c("taxa_name2"))
                  TOTAL2_combined_t$variable<-as.numeric(str_sub(TOTAL2_combined_t$variable, -3,-1))
                  TOTAL2_combined_t$method<-"combined"
                  
                  TOTAL2_wilcox<-top.relab.stats_table_combine_TOTAL %>% dplyr::select(starts_with(c("wilcox","taxa")))
                  TOTAL2_wilcox_t<-reshape2::melt(TOTAL2_wilcox, id.vars=c("taxa_name2"))
                  TOTAL2_wilcox_t$variable<-as.numeric(str_sub(TOTAL2_wilcox_t$variable, -3,-1))
                  TOTAL2_wilcox_t$method<-"wilcox"
        #########################################################################################################################################################                  
        #########################################################################################################################################################                 
        
        for (method in c("prev", "freq","combined","wilcox")){
          g<-ggplot()+
            geom_line(data=get(paste0(paste0("TOTAL_",method),"_t")) %>% filter(taxa_name!="Contaminant" & taxa_name!="Non-contaminant" &  taxa_name!="TOTAL"), 
                      aes(x=variable, y=value, linetype =taxa_name), linewidth=0.3)+ 
            geom_line(data=get(paste0(paste0("TOTAL_",method),"_t")) %>% filter(taxa_name=="Contaminant" | taxa_name=="Non-contaminant" ), 
                      aes(x=variable, y=value, color=taxa_name), linewidth=2)+
            geom_line(data=get(paste0(paste0("TOTAL2_",method),"_t")), 
                      aes(x=variable, y=value), colour= "darkorchid1", linewidth=1, alpha=0.5)+
            labs(title=method,x="test threshold",y="percent identified")+
            scale_linetype_discrete(name = "GENUS")+ 
            scale_x_continuous(breaks = seq(0, 1, by = 0.1))+
            scale_y_continuous(breaks = seq(0, 1, by = 0.1))+
            scale_color_discrete(name = "Correctly Identified As Expected")+ 
            theme(panel.background = element_blank(),
                  panel.border=element_rect(fill=NA),
                  panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
                  panel.grid.minor = element_blank(),strip.background=element_blank(),
                  plot.title=element_text(face="bold",hjust=0.5, size = 28), 
                  plot.subtitle = element_text(hjust=0.5),
                  axis.title=element_text(face="bold", size = 18),
                  axis.text.x=element_text(face="bold",size = 16),
                  axis.text.y=element_text(face="bold",size = 16),
                  axis.ticks=element_blank(),
                  plot.margin=unit(c(1,1,1,1),"line"))
          assign(paste0("graph_",method),g)
        }

        remove_y <- theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())

        p<-list(graph_prev,
                graph_freq+remove_y,
                graph_combined+remove_y,
                graph_wilcox+remove_y)
        #put the graphs next to each other. add a title and collec the legend- place it at bottom and horizontal
        combine_p<-wrap_plots(p, nrow = 1) + plot_layout(guides = "collect")  &  theme(legend.position='bottom',legend.direction = "horizontal") &
          plot_annotation(title = paste0(paste0(negative_sample_type," compare with "),name_compare_type),
                          theme= theme(title = element_text(color = "red", face = "bold", size = 24))) 
        #output graphs in pdf format
        graph_output_title<-paste0(output_excel,paste0(paste0(paste0(negative_sample_type,"_compare_with_"),name_compare_type),".pdf"))
        pdf(graph_output_title, width=22, height=7)
        show(combine_p)
        dev.off()
    
  }
}


####decontaminant_subplot_KW: identify decontaminant based on a particular method and create boxplot for each sample type  
decontaminant_subplot_KW <- function (input_phyloseq, 
                                      SampleID.unique=NULL, #if empty, SampleID.unique would be the rowname of the sample_data of the input_phyloseq
                                      sample_type_var_name,                                       
                                      sample_types=list(), 
                                      sample_type_color=list(), 
                                      sample_type_color_2nd=list(), 
                                      negative_sample_type=list(), ###the sample type that you want to be as negative control
                                      compare_type=list(), 
                                      bacterial_load, #for frequency only. uses ddPCR to determine library size
                                      method_type=c("MannWhit","preval","freq","combin"),
                                      stat_option=c("mean", "median"), ### the statistics to determine rank for the plot
                                      graph_option=c("boxplot","mean_SE", "mean_SD"), 
                                      test_threshold, ###prevalence level for the contaminant test
                                      log_scale, #will transform abundance to log(abundance*100 + 1)
                                      output_suffix, 
                                      taxa_genus_output=c("yes","no")) {#if want the taxa name to be shortened to genus level
  ###ensure intput_phyloseq is a absolute count phyloseq, not relative abundance phyloseq. and make sure phyloseq is in right orientation
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  countData_check <- floor(as(otu_table(input_phyloseq), "matrix")) # round down to nearest integer. if this phsyloeq is a relative abundance, then the entire countData_check would be 0
  if (all(countData_check==0)) {stop("Please use phyloseq with absolute count. The current input physloeq contains relative abundnace")}
  
  ## evaluate choices
  stat_option <- match.arg(stat_option)
  method_type <- match.arg(method_type)
  print(method_type)
  print(stat_option)  
  
  if (missing(output_suffix)) { output_suffix<-"OP"} #if output_suffix is missing, then will just give a suffix _OP at the end of the file
  if (missing(graph_option)) { graph_option<-"boxplot"} #if graph_option is missing, then will just have boxplot as default
  if (missing(log_scale)) { log_scale<-"no"} #if graph_option is missing, then will just have boxplot as default
  
  if (missing(taxa_genus_output)) { taxa_genus_output<-"no"} #default being yes for taxa_rank_name, which will just give the selected taxa name
  taxa_genus_output <- match.arg(taxa_genus_output)
  
  ##making sure sample_types and sample_type_color has same number of elements
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  #if sample_type_color_2nd is not missing then make sure it has same number of element as sample_type_color
  if (!missing(sample_type_color_2nd)) {
    stopifnot("sample_types need to have same number of elements as sample_type_color_2nd"= length(sample_types)==length(sample_type_color_2nd))
  } #if graph_option is missing, then will just have boxplot as default
  if (missing(sample_type_color_2nd)) {
    sample_type_color_2nd<-rep("black",length(sample_types)) #make the secondary color to be black by default, unless specify
  } 
  
  ###generating a variable for each sample type 
  for (i in 1:length(sample_types)) {
    assign(paste0("sample_type_color",i), sample_type_color[[i]])
    assign(paste0("sample_type_color_2nd",i), sample_type_color_2nd[[i]])
    assign(paste0("sampletype",i), sample_types[[i]])
  }  
  
  ###this make sure that the phyloseq ONLY contains the sample_type_var_name variable with all options listed sample_types
  #in theory input_phyloseq should be same as input_phyloseq_2 if the user input all of the sample_types options
  keep_sample<- get_variable(input_phyloseq,sample_type_var_name) %in% sample_types
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  input_phyloseq_2 <- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
  
  ####preparing compare_type list
  ####if the compare_type has "A_and_B" then will determine if taxa is contaminant for A and B when they are individually compared to the negative control
  ####if the compare_type has "A_or_B" then will determine if taxa is contaminant for A OR B when they are individually compared to the negative control
  ####if the compare_type has "A_combine_B" then will determine if taxa is contaminant when negative control is compare to both A and B together
  for (i in 1:length(compare_type)) {
    assign(paste0("compare_type",i), compare_type[[i]])
  }  
  
  ####getting relative abundance
  normalizeSample <- function(x){x/sum(x)}
  input_phyloseq_2_relab <- transformSampleCounts(input_phyloseq_2,normalizeSample) ##relative abundance
  ### Setting up ###
  # Counts from Phyloseq 
  counts.edit <- as.data.frame(otu_table(input_phyloseq_2))
  # Relative Abundance Table from Phyloseq 
  relab.edit <- as.data.frame(otu_table(input_phyloseq_2_relab))
  # reference table (not to be used)
  reference_table <- as.data.frame(sample_data(input_phyloseq_2))
  
  ##if SampleID.unique is missing, then will use the rowname of match to be the SampleID.unique given the rowname of sample_data of the phyloseq should be unique
  if (is.null(SampleID.unique)) {
    reference_table$SampleID_decontam<-rownames(reference_table) #this makes a variable in the DF match call SampleID.unique which would be rowname of the DF match
  } else {
    if (SampleID.unique %in% names(reference_table)){    
      reference_table$SampleID_decontam<-reference_table[[SampleID.unique]] #this makes a variable in the DF match call SampleID.unique which would be rowname of the DF match
    } else {stop(paste0(SampleID.unique," does not exist in the sample dataframe"))} #if the phyloseq sample dataframe doesn't have the SampleID.unique then stop
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ###adding bacterial load option- for frequency only 
  ###setting up the bacterial load data
  if (!missing(bacterial_load)){ #if want to use bacterial load
    if (!inherits(bacterial_load, c("data.frame"))) { #making sure that bacterial load is a data.frame
      stop_txt = paste0("The input data for bacterial_load need to be a data frame")
      stop(stop_txt, call. = FALSE)
    }
    if (ncol(bacterial_load)!=2) { #make sure the bacterial load dataframe is in the right structure
      stop_txt2 = paste0("The input data frame for bacterial_load needs to have two columns:first column being the SampleID.unique (or rowname of the sample_data of the input_phyloseq if SampleID.unique is not inputted)")
      stop(stop_txt2, call. = FALSE)
    }
    if(length(setdiff(reference_table$SampleID_decontam,bacterial_load[,1]))!=0) { ###if there are elements in the reference_table that is not in the bacterial_load first column, which should be 
      differentID<-paste(setdiff(reference_table$SampleID_decontam,bacterial_load[,1]), collapse=",  ")  
      stop_txt3= paste("The following SampleIDs do not exist in the bacterial_load dataframe:",
                       differentID)
      stop(stop_txt3, call. = FALSE)
    }
    colnames(bacterial_load)<-c("SampleID", "ddPCR_count")
  } 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  #######################################################################################################################################################################
  taxa.table <- as.data.frame(tax_table(input_phyloseq_2_relab))
  # This will paste the names together all the way from Family to OTU, separted by "." 
  taxa.table$match <- paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
  taxa.table.clean <- subset(taxa.table, select=c(match))
  taxa<-taxa.table$match #make a vector with all the taxa names
  
  # This will merge by column=0 which is the rowname by taxa-name-###
  # so now count.match and relab.match both have the new condensed taxa name 
  count.match <- merge(counts.edit, taxa.table.clean, by=0) 
  relab.match <- merge(relab.edit, taxa.table.clean, by=0)
  
  ### Replace count.match, relab.match rownames with taxa-name-###
  # This step replaces rownames for the match with count.match$match
  rownames(count.match) <- count.match$match
  rownames(relab.match) <- relab.match$match
  #remove match since match is now the rownames
  #remove Row.names given that is the byproduct of merging by rowname on previous step 
  counts <- subset(count.match, select=-c(Row.names, match))
  relab <- subset(relab.match, select=-c(Row.names, match))
  
  ### Match the Sample ID to the sample type categories  
  ### We need to do this to set up the comparison objects  
  ###replacing the SampleID.unique with the sample type for the count and relative abundance table
  index_var_name<-grep(sample_type_var_name, colnames(reference_table)) #this determine the nth column which sample+type_var_name is located at in the reference_table dataframe
  index2_var_name<-grep("SampleID_decontam", colnames(reference_table))
  #this step replace all of the unique subject ID (SampleID.unique) with their corresponding group of sample type of interest
  names(counts) <- reference_table[[index_var_name]][match(names(counts), reference_table[[index2_var_name]])] 
  names(relab) <- reference_table[[index_var_name]][match(names(relab), reference_table[[index2_var_name]])]  
  
  #### this input the sample types which will be used
  comparison.list<- vector(mode = "list", length = length(sample_types))
  names(comparison.list) <- sample_types
  for (i in 1:length(sample_types)) {
    comparison.list [i] <- list(grep(sample_types[i],colnames(counts)))
  }
  #comparison.list contain multiple lists (each list provide the column number of the specific sample type of interest). If there are 3 sample types in the dataframe of phyloseq then comparison will have a list of 3
  #counts.subset is basically counts, except there is now an index on the column name 
  counts.subset <- counts[,unlist(comparison.list)] # head(counts.subset)
  relab.subset <- relab[,unlist(comparison.list)] # head(relab.subset)
  #libsize consider the total abundance by each subject 
  libsizes <- colSums(counts.subset)
  libsizes_loadconcentrate<<-libsizes
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #####getting a dataframe that match the sample_type with column index with the initial sampleID
  if (!missing(bacterial_load)){ #if want to use bacterial load
      basedf<-count.match %>% select(-c("Row.names","match"))
      basedf_t<-data.frame(t(basedf))
      basedf_t$SampleID<-row.names(basedf_t)
      
      base2_df_t<-data.frame(t(counts.subset))
      base2_df_t$sample_index<-row.names(base2_df_t)
      base_combine<-merge(basedf_t,base2_df_t)
      base_combine<-base_combine %>% select(c(SampleID, sample_index))    
      
      index_key<<-base_combine
      
      bacterial_load_indexed<-merge(bacterial_load, base_combine, by="SampleID")
      
      bacterial_load_indexed2<-bacterial_load_indexed[order(bacterial_load_indexed$sample_index),]
      rownames(bacterial_load_indexed2)<-bacterial_load_indexed2$sample_index
      bacterial_load_indexed2<-bacterial_load_indexed2 %>% select(ddPCR_count)
      
      bacterial_load_indexed3<-as.numeric(unlist(bacterial_load_indexed2))
      names(bacterial_load_indexed3)<-rownames(bacterial_load_indexed2)  
      libsizes<-bacterial_load_indexed3 #library size is replaced with the ddPCR data
      print(paste0(bacterial_load," was used for frequency contaminant determination"))
      libsizes_bacterialload<<-libsizes
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  ### Create matrix of data for evaluation ### 
  ### This is to create an empty matrix for all top statistics ### 
  X<-nrow(relab.subset)
  top.relab.stats <- as.data.frame(matrix(ncol=length(sample_types),nrow=X))
  
  # Moving the rownames to top.relab.stats 
  # going create two columns for medians to plot 
  # This may not be necessary top.relab.stats <- NULL 
  rownames(top.relab.stats) <- rownames(relab.subset)
  
  #loop through the sample types and get statistics for each. if there are 3 sample_types, then this will fill up the first 3 column with the statistics
  for (i in 1:length(sample_types)) {
    if (stat_option=="median") {
      top.relab.stats[i] <-rowMedians(as.matrix(relab.subset[, grepl(sample_types[i], colnames(relab.subset))]))
    } else if (stat_option=="mean") {
      top.relab.stats[i] <-rowMeans2(as.matrix(relab.subset[, grepl(sample_types[i], colnames(relab.subset))]))
    }
    colnames(top.relab.stats)[i] <- c(paste0("stats.",sample_types[i]))
  }  
  #loop through the sample types and fill in the next several column with the rank order by that particular sample type. If there are 3 sample_types, then this will give rank order on 4th-6th column
  for (i in 1:length(sample_types)) {
    ##generating the rank order by sample types 
    j <- i+length(sample_types)
    column<-top.relab.stats[,i]
    top.relab.stats <- top.relab.stats[order(column, decreasing = TRUE), , drop = FALSE ]
    top.relab.stats[j] <- 1:nrow(top.relab.stats)
    colnames(top.relab.stats)[j] <- c(paste0("rank.order.",sample_types[i]))
  }  
  #this is to make taxa name with rank order in them
  for (i in 1:length(sample_types)) { 
    j <- i+length(sample_types)
    k <- i+2*length(sample_types)  
    column2<-top.relab.stats[,j]
    top.relab.stats[k] <- paste0(rownames(top.relab.stats)," ","(",column2,")")
    colnames(top.relab.stats)[k] <- c(paste0(paste0("rank.order.",sample_types[i]),".name"))
  }  
  
  #################################
  #setting up the negative control#
  #################################
  #subject-level total count for the negative control
  libsizes_negative <- libsizes[grepl(negative_sample_type, names(libsizes))]
  w<-match(negative_sample_type,sample_types)
  comparison.list_negative<- unlist(comparison.list[w]) #extract the list of negative control (this list give the column number which is consider negative control)
  comparison.list_negative<-ifelse(comparison.list_negative, T, F) #basically making this list turn into all 1 
  #for counts.subset_negative and relab.subset_negative, the columns are subject (they should only be in the negative controls)
  counts.subset_negative <- dplyr::select(counts.subset, starts_with(negative_sample_type))
  relab.subset_negative <- dplyr::select(relab.subset, starts_with(negative_sample_type))
  
  ###keeping only the statistics for each sample type 
  top.relab.stats_table <- top.relab.stats[c(1:length(sample_types))]
  relab.subset_final<-relab.subset
  
  for (h in 1:length(compare_type)) {
    name_compare_type<-get(paste0("compare_type",h))
    print(paste0(negative_sample_type," compare:"))
    print(name_compare_type)
    AND<-grepl("_and_",name_compare_type, fixed=T)
    OR<-grepl("_or_",name_compare_type, fixed=T)
    COMBINE<-grepl("_combine_",name_compare_type, fixed=T)
    contam<-list()
    
    if ((AND==TRUE) | (OR==TRUE)) {
      if (AND==TRUE) {
        compare_sublist <- unlist(strsplit(name_compare_type, "_and_"))
      }else if (OR==TRUE) {
        compare_sublist <- unlist(strsplit(name_compare_type, "_or_"))
      }
      #loop through the two sample types 
      for (i in 1:length(compare_sublist)) { 
        assign(paste0("compare_sublist",i), compare_sublist[i])
      } 
      
      for (b in 1:length(compare_sublist)) {
        # only keeping the compare list  
        selected_sample_type<-get(paste0("compare_sublist",b))
        p<-match(selected_sample_type,sample_types) 
        comparison.list_select<- unlist(comparison.list[p]) 
        comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
        
        # Remove unused from conc, in this case libsizes 
        libsizes_select <- libsizes[grepl(selected_sample_type, names(libsizes))]
        counts.subset_select <- dplyr::select(counts.subset, starts_with(selected_sample_type))
        relab.subset_select <- dplyr::select(relab.subset, starts_with(selected_sample_type))
        
        counts.subset_comb<-merge(counts.subset_select,counts.subset_negative, by=0) #put back the negative control with the sample that you want to compare 
        relab.subset_comb<-merge(relab.subset_select,relab.subset_negative, by=0)
        #remove the new column "Row.names" created when merging by 0
        rownames(counts.subset_comb) <- counts.subset_comb$Row.names
        rownames(relab.subset_comb) <- relab.subset_comb$Row.names
        counts.subset_comb <- subset(counts.subset_comb, select=-c(Row.names))
        relab.subset_comb <- subset(relab.subset_comb, select=-c(Row.names))
        
        negative_comb<- append(comparison.list_select,comparison.list_negative)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        ###need to sort libsizes_select and libsizes_negative
        ###this is to ensure that the correct library size is link to the corresponding samples 
        libsizes_select_TEMP<-data.frame(libsizes_select)
        libsizes_select_TEMP$index<- as.numeric(sub(".*\\.", "", row.names(libsizes_select_TEMP)))
        libsizes_select_TEMP$index[is.na(libsizes_select_TEMP$index)]<-0 #replace NA with 0 
        libsizes_select_TEMP <- libsizes_select_TEMP[order(libsizes_select_TEMP$index), ]
        libsizes_select2<-as.numeric(libsizes_select_TEMP$libsizes_select)
        names(libsizes_select2)<-row.names(libsizes_select_TEMP)
        
        libsizes_negative_TEMP<-data.frame(libsizes_negative)
        libsizes_negative_TEMP$index<- as.numeric(sub(".*\\.", "", row.names(libsizes_negative_TEMP)))
        libsizes_negative_TEMP$index[is.na(libsizes_negative_TEMP$index)]<-0 #replace NA with 0 
        libsizes_negative_TEMP <- libsizes_negative_TEMP[order(libsizes_negative_TEMP$index), ]
        libsizes_negative2<-as.numeric(libsizes_negative_TEMP$libsizes_negative)
        names(libsizes_negative2)<-row.names(libsizes_negative_TEMP)
        libsizes_comb <- append(libsizes_select2,libsizes_negative2 )
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#      
        
        contam$prev <- isContaminant(t(as.matrix(counts.subset_comb)),
                                     method="prevalence",
                                     neg=negative_comb,
                                     threshold=test_threshold)

        contam$freq <- isContaminant(t(as.matrix(counts.subset_comb)),
                                     method="frequency",
                                     neg=negative_comb,
                                     conc=libsizes_comb,
                                     threshold=test_threshold)

        contam$combined <- isContaminant(t(as.matrix(counts.subset_comb)),
                                         method="combined",
                                         neg=negative_comb,
                                         conc=libsizes_comb,
                                         threshold=test_threshold)
        contam$wilcox <- wilcox_contam_test(relab.subset_comb,
                                            neg=negative_comb,
                                            threshold=test_threshold)
        contam_dataframe<- contam %>% data.frame #compress the lists from the different results
        contam_dataframe$wilcox.contaminant<-as.logical(contam_dataframe$wilcox.contaminant)
        
        contam_result<-subset(contam_dataframe, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
        assign(paste0("contam_result",b),contam_result)
      }
      #combine the results for the two sample types
      
      contam_result1a<<-contam_result1
      contam_result2a<<-contam_result2
      
      contam_result_final<-merge(contam_result1,contam_result2, by=0)
      rownames(contam_result_final) <- contam_result_final$Row.names
      contam_result_final <- subset(contam_result_final, select=-c(Row.names))
      if (AND==TRUE) {
        contam_result_final <- contam_result_final %>% mutate(prev.contaminant = case_when((prev.contaminant.x==T & prev.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
        contam_result_final <- contam_result_final %>% mutate(freq.contaminant = case_when((freq.contaminant.x==T & freq.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
        contam_result_final <- contam_result_final %>% mutate(combined.contaminant = case_when((combined.contaminant.x==T & combined.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
        contam_result_final <- contam_result_final %>% mutate(wilcox.contaminant = case_when((wilcox.contaminant.x==T & wilcox.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
        contam_result_final <- subset(contam_result_final, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
      }else if (OR==TRUE) {
        contam_result_final <- contam_result_final %>% mutate(prev.contaminant = case_when((prev.contaminant.x==T | prev.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
        contam_result_final <- contam_result_final %>% mutate(freq.contaminant = case_when((freq.contaminant.x==T | freq.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
        contam_result_final <- contam_result_final %>% mutate(combined.contaminant = case_when((combined.contaminant.x==T | combined.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))
        contam_result_final <- contam_result_final %>% mutate(wilcox.contaminant = case_when((wilcox.contaminant.x==T | wilcox.contaminant.y==T)  ~ TRUE, TRUE ~ FALSE))   
        contam_result_final <- subset(contam_result_final, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
      }
    }else{
      if (COMBINE==TRUE) { ###if compare_type is a combination set 
        ###this works if A_combine_B... it basically combine subset A and subset B
        combine_list <- unlist(strsplit(name_compare_type, "_combine_"))
        f<-match(combine_list,sample_types) 
        comparison.list_select<- unlist(comparison.list[f]) 
        comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
        combine_list_temp<-str_replace(name_compare_type, "_combine_", "|") 
        #Setting up the dataset so that it include both sample types that are combine =
        libsizes_select <- libsizes[grepl(combine_list_temp, names(libsizes))]
        counts.subset_select <- dplyr::select(counts.subset, starts_with(combine_list))
        relab.subset_select <- dplyr::select(relab.subset, starts_with(combine_list))
      } else { ###if compare_type does not contain _and_, _or_, _combine_.... then basically compare_type is just one of the sample types
        
        if (name_compare_type %in% sample_types) { #name_compare_type has to be one of the initial sample_types if it is not _and_, _or_, _combine_
          f<-match(name_compare_type,sample_types) 
          comparison.list_select<- unlist(comparison.list[f]) 
          comparison.list_select<-ifelse(comparison.list_select, F, T) #basically making this list turn into all 1 
          #keeping only the sample type of interest
          libsizes_select <- libsizes[grepl(name_compare_type, names(libsizes))]
          counts.subset_select <- dplyr::select(counts.subset, starts_with(name_compare_type))
          relab.subset_select <- dplyr::select(relab.subset, starts_with(name_compare_type))
        } else { stop(paste0(paste0(name_compare_type," must be in sample type list:"),sample_types)) }
        
      }
      #combining with the negative control group
      counts.subset_comb<-merge(counts.subset_select,counts.subset_negative, by=0) #put back the negative control with the sample that you want to compare 
      relab.subset_comb<-merge(relab.subset_select,relab.subset_negative, by=0)
      #remove the new column "Row.names" created when merging by 0
      rownames(counts.subset_comb) <- counts.subset_comb$Row.names
      rownames(relab.subset_comb) <- relab.subset_comb$Row.names
      counts.subset_comb <- subset(counts.subset_comb, select=-c(Row.names))
      relab.subset_comb <- subset(relab.subset_comb, select=-c(Row.names))
      
      negative_comb<- append(comparison.list_select,comparison.list_negative)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      ###need to sort libsizes_select and libsizes_negative
      ###this is to ensure that the correct library size is link to the corresponding samples 
      libsizes_select_TEMP<-data.frame(libsizes_select)
      libsizes_select_TEMP$index<- as.numeric(sub(".*\\.", "", row.names(libsizes_select_TEMP)))
      libsizes_select_TEMP$index[is.na(libsizes_select_TEMP$index)]<-0 #replace NA with 0 
      libsizes_select_TEMP <- libsizes_select_TEMP[order(libsizes_select_TEMP$index), ]
      libsizes_select2<-as.numeric(libsizes_select_TEMP$libsizes_select)
      names(libsizes_select2)<-row.names(libsizes_select_TEMP)
      
      libsizes_negative_TEMP<-data.frame(libsizes_negative)
      libsizes_negative_TEMP$index<- as.numeric(sub(".*\\.", "", row.names(libsizes_negative_TEMP)))
      libsizes_negative_TEMP$index[is.na(libsizes_negative_TEMP$index)]<-0 #replace NA with 0 
      libsizes_negative_TEMP <- libsizes_negative_TEMP[order(libsizes_negative_TEMP$index), ]
      libsizes_negative2<-as.numeric(libsizes_negative_TEMP$libsizes_negative)
      names(libsizes_negative2)<-row.names(libsizes_negative_TEMP)
      libsizes_comb <- append(libsizes_select2,libsizes_negative2 )
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#      
      
      contam$prev <- isContaminant(t(as.matrix(counts.subset_comb)),
                                   method="prevalence",
                                   neg=negative_comb,
                                   threshold=test_threshold)
      contam$freq <- isContaminant(t(as.matrix(counts.subset_comb)),
                                   method="frequency",
                                   neg=negative_comb,
                                   conc=libsizes_comb,
                                   threshold=test_threshold)

      contam$combined <- isContaminant(t(as.matrix(counts.subset_comb)),
                                       method="combined",
                                       neg=negative_comb,
                                       conc=libsizes_comb,
                                       threshold=test_threshold)
      contam$wilcox <- wilcox_contam_test(relab.subset_comb,
                                          neg=negative_comb,
                                          threshold=test_threshold)
      contam_dataframe<- contam %>% data.frame #compress the lists from the different results
      contam_dataframe$wilcox.contaminant<-as.logical(contam_dataframe$wilcox.contaminant)
      
      contam_result_final<-subset(contam_dataframe, select=c(prev.contaminant,freq.contaminant,combined.contaminant, wilcox.contaminant))
    }
    lowerbound<-length(sample_types)+1
    upperbound<-2*length(sample_types)
    rankorder <- subset(top.relab.stats[c(lowerbound:upperbound)]) ##keeping the rank order only
    contam_result_final_t<-contam_result_final
    if (method_type=="MannWhit"){
      contam_result_final_t$contaminant <- contam_result_final_t$wilcox.contaminant  
    }
    else if (method_type=="preval") {
      contam_result_final_t$contaminant <- contam_result_final_t$prev.contaminant
    }
    else if (method_type=="freq") {
      contam_result_final_t$contaminant <- contam_result_final_t$freq.contaminant
    }
    else if (method_type=="combin") {
      contam_result_final_t$contaminant <- contam_result_final_t$combined.contaminant
    }
    contam_result_final_t<-subset(contam_result_final_t,select=c(contaminant))
    csv_output_contam <- merge(contam_result_final_t,rankorder,by="row.names")
    
    ### Come back from the DECONTAM section to add this in ###    
    # Adding contam data (saved to the end to add)
    ###figure out which taxa to highlight in red based on which test you choice 
    if (method_type == "MannWhit") {
      top.relab.stats$contaminant <- contam_result_final$wilcox.contaminant[match(rownames(top.relab.stats), rownames(contam_result_final))] ##color based on prev decontam
      print("MannWhit")
    }
    else if (method_type == "preval") {
      top.relab.stats$contaminant <- contam_result_final$prev.contaminant[match(rownames(top.relab.stats), rownames(contam_result_final))] ##color based on prev decontam
      print("preval")
    }
    else if (method_type == "freq") {
      top.relab.stats$contaminant <- contam_result_final$freq.contaminant[match(rownames(top.relab.stats), rownames(contam_result_final))] ##color based on prev decontam
      print("freq")
    }
    else if (method_type == "combin") {
      top.relab.stats$contaminant <- contam_result_final$combined.contaminant[match(rownames(top.relab.stats), rownames(contam_result_final))] ##color based on prev decontam
      print("combin")
    }
    ###########################################################################################################
    ###########################################################################################################
    top.relab.stats_final<-top.relab.stats
    DecontKW_contaminant_list<-subset(top.relab.stats,select=c(contaminant))
    assign(paste0(paste0("Contam_list_NC_",negative_sample_type),paste0("_compare_",paste0(name_compare_type,paste0("_",output_suffix)))), DecontKW_contaminant_list,.GlobalEnv)
    print("###########################")
    print("Available contaminant list:")
    print(paste0(paste0("Contam_list_NC_",negative_sample_type),paste0("_compare_",paste0(name_compare_type,paste0("_",output_suffix)))))
    print("###########################")
    ###########################################################################################################
    ###########################################################################################################
    for (rank_sample_type in sample_types) {
      l<-match(rank_sample_type,sample_types) ###look at the element in rank_sample_type to see what their index at sample_types
      r<-l+length(sample_types)
      # Re-rank based upon the rank group decreasing from highest to lowest
      column3<-top.relab.stats[,r]
      top.relab.stats <- top.relab.stats[order(column3, decreasing = F), , drop = FALSE ]
      top.2 <- top.relab.stats
      top.3 <- head(top.2, 100)
      # Filtering ONLY the top 100 taxa 
      relab.subset.top <- relab.subset %>% 
        dplyr::filter(rownames(relab.subset) %in% rownames(top.3))
      # Re ordered and saved into a new group (previously filtered)
      # relab.subset.top.match now have only 100 taxa ordered by top background 
      relab.subset.top.match <- relab.subset.top[match(rownames(top.3),rownames(relab.subset.top)),]    
      top.3$color <- ifelse(top.3$contaminant, "red", "black")
      
      for (a in sample_types) {
        assign(paste0("rank_order_",a),subset(top.3, select=c(paste0("rank.order.",a), "color", paste0(paste0("rank.order.",a),".name") )))
      }      
      
      # Transposition of the figure, maybe at this point replace the column names 
      # The previous edit should have removed all character columns 
      relab.subset.top.match.transpose <- as.data.frame(t(relab.subset.top.match))
      # name a new column with categories 
      relab.subset.top.match.transpose$category <- rownames(relab.subset.top.match.transpose)
      
      for (a in sample_types) {
        relab.subset.top.match.transpose$category[grepl(a,relab.subset.top.match.transpose$category)] <- a
      }  
      
      # Reverse column order 
      relab.subset.top.match.transpose.reverse <- relab.subset.top.match.transpose[,order(ncol(relab.subset.top.match.transpose):1)]
      relab.subset.top.match.transpose.reverse_final<-relab.subset.top.match.transpose.reverse
      
      for (a in sample_types) {
        assign(paste0("relab.subset.top.match.transpose.reverse.",a),subset(relab.subset.top.match.transpose.reverse, category == c(a) ))
        assign(paste0("relab.subset.top.match.transpose.reverse.",a),subset(get(paste0("relab.subset.top.match.transpose.reverse.",a)), select=-c(category)))
      }  
      
      #get the first plot- rank_sample_type is the sample type that is being ranked.... we want the ranked plot on the left 
      temp <- get(paste0("relab.subset.top.match.transpose.reverse.",rank_sample_type)) %>% dplyr::select(everything()) %>% tidyr::gather("id", "value",1:100, factor_key = TRUE)
      temp <- merge(temp, top.3, by.x="id", by.y="row.names")
      #if log_scale is yes then will have log(abundance+1)
      if (tolower(log_scale)=="yes"){
        temp$value2<-log10(temp$value*100+1) 
      } else {
        temp$value2<-temp$value
      }
      p <- temp %>%   mutate(id = fct_reorder(id, -get(paste0("rank.order.",rank_sample_type)))) %>%
        ggplot(., aes(x=id, y=value2)) + 
        coord_flip() +
        theme_bw() + 
        theme(axis.text.y=element_text(color=rev(top.3$color)),axis.title.y=element_blank()) +
        theme(axis.text.x=element_text(angle=90,hjust=1))+
        ggtitle(rank_sample_type) 
      
      #taxa naming
      if (tolower(taxa_genus_output)=="no"){
        p <-  p + 
          scale_x_discrete(labels = function(x) {
            x_org<-x
            is_long <- nchar(x) > 35
            x[is_long] <- paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",4),1,10),
                                 paste0("..",paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",5),1,10),".."),
                                        paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",6),1,10),
                                               paste0("..",str_sub(x[is_long],-4,-1)))))
            is_short <- nchar(x) < 17
            x[is_short] <- paste0(paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",3),1,15),".."),
                                  paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",4),1,5),
                                         paste0("..",paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",5),1,5),".."),
                                                paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",6),1,5),
                                                       paste0("..",str_sub(x_org[is_short],-4,-1))))))
            x})
      } else if (tolower(taxa_genus_output)=="yes") {
        p <-  p + 
          scale_x_discrete(drop = FALSE, labels = function(x) {
            x_org<-x
            genus_no_NA<- sapply(strsplit(x,".",fixed=T),"[[",5)!="NA"
            genus_NA<- sapply(strsplit(x,".",fixed=T),"[[",4)!="NA" & sapply(strsplit(x,".",fixed=T),"[[",5)=="NA"
            family_NA<- sapply(strsplit(x,".",fixed=T),"[[",3)!="NA" & sapply(strsplit(x,".",fixed=T),"[[",4)=="NA"
            order_NA<- sapply(strsplit(x,".",fixed=T),"[[",2)!="NA" & sapply(strsplit(x,".",fixed=T),"[[",3)=="NA"
            class_NA<- sapply(strsplit(x,".",fixed=T),"[[",1)!="NA" & sapply(strsplit(x,".",fixed=T),"[[",2)=="NA"
            x[genus_no_NA] <- paste0("g_",sapply(strsplit(x[genus_no_NA],".",fixed=T),"[[",5))
            x[genus_NA] <- paste0(paste0("f_",sapply(strsplit(x[genus_NA],".",fixed=T),"[[",4)),".g_NA")
            x[family_NA] <- paste0(paste0("o_",sapply(strsplit(x[family_NA],".",fixed=T),"[[",3)),".f_NA.g_NA")
            x[order_NA] <- paste0(paste0("c_",sapply(strsplit(x[order_NA],".",fixed=T),"[[",2)),".o_NA.f_NA.g_NA")
            x[class_NA] <- paste0(paste0("p_",sapply(strsplit(x[class_NA],".",fixed=T),"[[",1)),".c_NA.o_NA.f_NA.g_NA")
            x}) 
      }
      if (tolower(log_scale)=="yes"){
        p<-p+ylab("Log 10 (Relative Abundance)")
      } else {
        p<-p+ylab("Relative Abundance")
      }
      if (graph_option=="boxplot") {
        p<-p+geom_boxplot(color=sample_type_color[[l]], alpha=0.2) 
      } else if(graph_option=="mean_SD") {
        p<-p+geom_jitter(color=sample_type_color[[l]],position=position_jitter(0), alpha=0.2) +
          stat_summary(fun= mean, 
                       geom="pointrange", 
                       fun.max = function (x) mean(x)+sd(x),
                       fun.min = function (x) ifelse( mean(x)-sd(x) < 0, 0, mean(x)-sd(x)),
                       color=sample_type_color_2nd[[l]], linewidth =1.0, size=0.4)
      } else if(graph_option=="mean_SE") {
        p<-p+geom_jitter(color=sample_type_color[[l]],position=position_jitter(0), alpha=0.2) +
          stat_summary(fun= mean, 
                       geom="pointrange", 
                       fun.max = function (x) mean(x)+sd(x)/sqrt(length(x)),
                       fun.min = function (x) ifelse( mean(x)-sd(x)/sqrt(length(x)) < 0,0, mean(x)-sd(x)/sqrt(length(x))),
                       color=sample_type_color_2nd[[l]], linewidth =1.0, size=0.4)
      }
      myplots <- list()  # new empty list
      myplots[[1]]<- p
      
      #remove the sample type that is being ranked. so new_sample_type_name and new_sample_type_color contain the sample types which were not ranked                
      new_sample_type_name<-sample_types[-l] #remove the lth element (lth element is the sample type that is being ranked by )
      new_sample_type_color<-sample_type_color[-l] #remove the lth element (lth element is the sample type that is being ranked by )
      new_sample_type_color_2nd<-sample_type_color_2nd[-l] #remove the lth element (lth element is the sample type that is being ranked by )
      for (q in 1:length(new_sample_type_name)) {
        assign(paste0("new_sample_type_color",q), new_sample_type_color[[q]])
        assign(paste0("new_sample_type_color_2nd",q), sample_type_color_2nd[[q]])
        assign(paste0("new_sample_type_name",q), new_sample_type_name[[q]])
      }  
      #remove the axis label for the sample types which were not ranked (since the ranked sample type would be on the left most)  
      remove_y<- theme(axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       axis.title.y = element_blank())
      
      ###looping thro each sample type 
      number<-1
      number2<-2
      for (a in new_sample_type_name) {
        temp <- get(paste0("relab.subset.top.match.transpose.reverse.",a)) %>% dplyr::select(everything()) %>% tidyr::gather("id", "value",1:100, factor_key = TRUE)
        temp <- merge(temp, top.3, by.x="id", by.y="row.names")
        #if log_scale is yes then will have log(abundance+1)
        if (tolower(log_scale)=="yes"){
          temp$value2<-log10(temp$value*100+1) 
        } else {
          temp$value2<-temp$value
        }
        p <- temp %>%   mutate(id = fct_reorder(id, -get(paste0("rank.order.",rank_sample_type)))) %>%
          ggplot(., aes(x=id, y=value2)) + 
          coord_flip() +
          theme_bw() + 
          theme(axis.text.y=element_text(color=rev(top.3$color)),axis.title.y=element_blank()) +
          theme(axis.text.x=element_text(angle=90,hjust=1))+
          ggtitle(a) 
        
        #taxa naming
        if (tolower(taxa_genus_output)=="no"){
          p <-  p + 
            scale_x_discrete(labels = function(x) {
              x_org<-x
              is_long <- nchar(x) > 35
              x[is_long] <- paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",4),1,10),
                                   paste0("..",paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",5),1,10),".."),
                                          paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",6),1,10),
                                                 paste0("..",str_sub(x[is_long],-4,-1)))))
              is_short <- nchar(x) < 17
              x[is_short] <- paste0(paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",3),1,15),".."),
                                    paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",4),1,5),
                                           paste0("..",paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",5),1,5),".."),
                                                  paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",6),1,5),
                                                         paste0("..",str_sub(x_org[is_short],-4,-1))))))
              x})
        } else if (tolower(taxa_genus_output)=="yes") {
          p <-  p + 
            scale_x_discrete(drop = FALSE, labels = function(x) {
              x_org<-x
              genus_no_NA<- sapply(strsplit(x,".",fixed=T),"[[",5)!="NA"
              genus_NA<- sapply(strsplit(x,".",fixed=T),"[[",4)!="NA" & sapply(strsplit(x,".",fixed=T),"[[",5)=="NA"
              family_NA<- sapply(strsplit(x,".",fixed=T),"[[",3)!="NA" & sapply(strsplit(x,".",fixed=T),"[[",4)=="NA"
              order_NA<- sapply(strsplit(x,".",fixed=T),"[[",2)!="NA" & sapply(strsplit(x,".",fixed=T),"[[",3)=="NA"
              class_NA<- sapply(strsplit(x,".",fixed=T),"[[",1)!="NA" & sapply(strsplit(x,".",fixed=T),"[[",2)=="NA"
              x[genus_no_NA] <- paste0("g_",sapply(strsplit(x[genus_no_NA],".",fixed=T),"[[",5))
              x[genus_NA] <- paste0(paste0("f_",sapply(strsplit(x[genus_NA],".",fixed=T),"[[",4)),".g_NA")
              x[family_NA] <- paste0(paste0("o_",sapply(strsplit(x[family_NA],".",fixed=T),"[[",3)),".f_NA.g_NA")
              x[order_NA] <- paste0(paste0("c_",sapply(strsplit(x[order_NA],".",fixed=T),"[[",2)),".o_NA.f_NA.g_NA")
              x[class_NA] <- paste0(paste0("p_",sapply(strsplit(x[class_NA],".",fixed=T),"[[",1)),".c_NA.o_NA.f_NA.g_NA")
              x}) 
        }
        
        #if log_scale is yes then will have log(abundance+1)
        if (tolower(log_scale)=="yes"){
          p<-p+ylab("Log 10 (Relative Abundance)")
        } else {
          p<-p+ylab("Relative Abundance")
        }
        if (graph_option=="boxplot") {
          p<-p+geom_boxplot(color=new_sample_type_color[[number]], alpha=0.2) 
        } else if(graph_option=="mean_SD") {
          p<-p+geom_jitter(color=new_sample_type_color[[number]],position=position_jitter(0), alpha=0.2) +
            stat_summary(fun= mean, 
                         geom="pointrange", 
                         fun.max = function (x) mean(x)+sd(x),
                         fun.min = function (x) ifelse( mean(x)-sd(x) < 0, 0, mean(x)-sd(x)),
                         color=new_sample_type_color_2nd[[number]], linewidth =1.0, size=0.4)
        } else if(graph_option=="mean_SE") {
          p<-p+geom_jitter(color=new_sample_type_color[[number]],position=position_jitter(0), alpha=0.2) +
            stat_summary(fun= mean, 
                         geom="pointrange", 
                         fun.max = function (x) mean(x)+sd(x)/sqrt(length(x)),
                         fun.min = function (x) ifelse( mean(x)-sd(x)/sqrt(length(x)) < 0,0, mean(x)-sd(x)/sqrt(length(x))),
                         color=new_sample_type_color_2nd[[number]], linewidth =1.0, size=0.4)
        }
        myplots[[number2]] <- p + remove_y
        number<-number+1
        number2<-number2+1
      }
      
      pdf_output=paste0(paste0(paste0(paste0(paste0(paste0("contaminant_",method_type),paste0("_negC_",negative_sample_type)),paste0("_compare_",name_compare_type)),paste0("_rank_",rank_sample_type)),paste0("thres_",test_threshold)),paste0(paste0("_",output_suffix),".pdf"))
      pdf(file=pdf_output, width=14+(length(sample_types)-3)*4, height=(14+(length(sample_types)-3)*4)*10/8.5)
             show(wrap_plots(myplots, nrow=1))
      dev.off()        
      #merging in the statistics
      top.relab.stats_table$Row.names<-row.names(top.relab.stats_table)
      csv_output_contam_final<-merge(csv_output_contam, top.relab.stats_table, by="Row.names")
      
      ####export csv
      csv_output=paste0(paste0(paste0(paste0(paste0("contaminant_",method_type),paste0("_negC_",negative_sample_type)),paste0("_compare_",name_compare_type)),paste0("thres_",test_threshold)),paste0(paste0("_",output_suffix),".csv"))
      write.csv(csv_output_contam_final,csv_output, row.names = FALSE)
      xlsx_output=paste0(paste0(paste0(paste0(paste0("contaminant_",method_type),paste0("_negC_",negative_sample_type)),paste0("_compare_",name_compare_type)),paste0("thres_",test_threshold)),paste0(paste0("_",output_suffix),".xlsx"))
      write.xlsx(csv_output_contam_final, xlsx_output, sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)
    }    
    ############################################################################################################################
    ############################################################################################################################
    ############################################################################################################################
    ##Generating boxplot rank by the corresponding sample type 
    count_table<-as.data.frame(relab.subset_final)
    for (g in 1:length(sample_types)) {
      assign(paste0("relab.subset_final_sampletype",g), count_table[,grepl(get(paste0("sampletype",g)), colnames(count_table))])
    }
    top.relab.stats_final$color <- ifelse(top.relab.stats_final$contaminant, "red", "black")
    
    end_bound<-3*length(sample_types)+2
    for (g in 1:length(sample_types)) {
      h<-g+length(sample_types)
      assign(paste0("rank_order_sample_type",g), as.data.frame(top.relab.stats_final[,c(h)]))
      temp_df<-get(paste0("rank_order_sample_type",g))
      rownames(temp_df)<-rownames(top.relab.stats_final)
      colnames(temp_df)<-"rank"
      temp_df3<-as.data.frame(top.relab.stats_final[,c(h,end_bound)])
      temp_df3<-temp_df3[order(temp_df3$rank, decreasing = T), ,drop = FALSE ]
      temp_df3<-tail(temp_df3, 100)
      assign(paste0("rank_order_sample_type",g), temp_df)
      assign(paste0("color_rank_order_sample_type",g), temp_df3) 
    }
    for (g in 1:length(sample_types)) {
      assign(paste0("relab.subset_final_sampletype_merge",g),merge(get(paste0("relab.subset_final_sampletype",g)),get(paste0("rank_order_sample_type",g)),by="row.names"))
      temp_df2<-get(paste0("relab.subset_final_sampletype_merge",g))
      rownames(temp_df2)=temp_df2$Row.name
      temp_df2=temp_df2[2:length(temp_df2)]
      temp_df2 <- temp_df2[order(temp_df2$rank, decreasing = F), , drop = FALSE ]
      temp_df2 <- head(temp_df2, 100)#keep top 100
      temp_df2 = subset(temp_df2, select = -c(rank))
      temp_df2 <- as.data.frame(t(temp_df2))
      temp_df2_b <- temp_df2[,order(ncol(temp_df2):1)]
      assign(paste0("relab.subset_final_sampletype_merge.reverse",g),temp_df2_b)
    }
    ###looping thro each sample type 
    number3<-1
    myplots_b<-list()
    for (g in sample_types) {
      temp_data<-get(paste0("relab.subset_final_sampletype_merge.reverse",number3)) %>% dplyr::select(everything()) %>% tidyr::gather("id", "value",1:100, factor_key = TRUE)
      #if log_scale is yes then will have log(abundance+1)
      if (tolower(log_scale)=="yes"){
        temp_data$value2<-log10(temp_data$value*100+1) 
      } else {
        temp_data$value2<-temp_data$value
      }
      pp <- temp_data  %>%
        ggplot(., aes(x=id, y=value2)) + 
        coord_flip() +
        theme_bw() + 
        theme(axis.text.y=element_text(color=get(paste0("color_rank_order_sample_type",number3))$color),axis.title.y=element_blank()) +
        theme(axis.text.x=element_text(angle=90,hjust=1))+
        ggtitle(g) 
      #taxa naming
      if (tolower(taxa_genus_output)=="no"){
        pp <-  pp + 
          scale_x_discrete(labels = function(x) {
            x_org<-x
            is_long <- nchar(x) > 35
            x[is_long] <- paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",4),1,10),
                                 paste0("..",paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",5),1,10),".."),
                                        paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",6),1,10),
                                               paste0("..",str_sub(x[is_long],-4,-1)))))
            is_short <- nchar(x) < 17
            x[is_short] <- paste0(paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",3),1,15),".."),
                                  paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",4),1,5),
                                         paste0("..",paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",5),1,5),".."),
                                                paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",6),1,5),
                                                       paste0("..",str_sub(x_org[is_short],-4,-1))))))
            x})
      } else if (tolower(taxa_genus_output)=="yes") {
        pp <-  pp + 
          scale_x_discrete(drop = FALSE, labels = function(x) {
            x_org<-x
            genus_no_NA<- sapply(strsplit(x,".",fixed=T),"[[",5)!="NA"
            genus_NA<- sapply(strsplit(x,".",fixed=T),"[[",4)!="NA" & sapply(strsplit(x,".",fixed=T),"[[",5)=="NA"
            family_NA<- sapply(strsplit(x,".",fixed=T),"[[",3)!="NA" & sapply(strsplit(x,".",fixed=T),"[[",4)=="NA"
            order_NA<- sapply(strsplit(x,".",fixed=T),"[[",2)!="NA" & sapply(strsplit(x,".",fixed=T),"[[",3)=="NA"
            class_NA<- sapply(strsplit(x,".",fixed=T),"[[",1)!="NA" & sapply(strsplit(x,".",fixed=T),"[[",2)=="NA"
            x[genus_no_NA] <- paste0("g_",sapply(strsplit(x[genus_no_NA],".",fixed=T),"[[",5))
            x[genus_NA] <- paste0(paste0("f_",sapply(strsplit(x[genus_NA],".",fixed=T),"[[",4)),".g_NA")
            x[family_NA] <- paste0(paste0("o_",sapply(strsplit(x[family_NA],".",fixed=T),"[[",3)),".f_NA.g_NA")
            x[order_NA] <- paste0(paste0("c_",sapply(strsplit(x[order_NA],".",fixed=T),"[[",2)),".o_NA.f_NA.g_NA")
            x[class_NA] <- paste0(paste0("p_",sapply(strsplit(x[class_NA],".",fixed=T),"[[",1)),".c_NA.o_NA.f_NA.g_NA")
            x}) 
      }
      #if log_scale is yes then will have log(abundance+1)
      if (tolower(log_scale)=="yes"){
        pp<-pp+ylab("Log 10 (Relative Abundance)")
      } else {
        pp<-pp+ylab("Relative Abundance")
      }
      if (graph_option=="boxplot") {
        pp<-pp+geom_boxplot(color=sample_type_color[[number3]], alpha=0.2) 
      } else if(graph_option=="mean_SD") {
        pp<-pp+geom_jitter(color=sample_type_color[[number3]],position=position_jitter(0), alpha=0.2) +
          stat_summary(fun= mean, 
                       geom="pointrange", 
                       fun.max = function (x) mean(x)+sd(x),
                       fun.min = function (x) ifelse( mean(x)-sd(x) < 0, 0, mean(x)-sd(x)),
                       color=sample_type_color_2nd[[number3]], linewidth =1.0, size=0.4)
      } else if(graph_option=="mean_SE") {
        pp<-pp+geom_jitter(color=sample_type_color[[number3]],position=position_jitter(0), alpha=0.2) +
          stat_summary(fun= mean, 
                       geom="pointrange", 
                       fun.max = function (x) mean(x)+sd(x)/sqrt(length(x)),
                       fun.min = function (x) ifelse( mean(x)-sd(x)/sqrt(length(x)) < 0,0, mean(x)-sd(x)/sqrt(length(x))),
                       color=sample_type_color_2nd[[number3]], linewidth =1.0, size=0.4)
      }
      myplots_b[[number3]] <- pp
      number3<-number3+1
    }
    all_plot.all <- ggarrange(plotlist=myplots_b, ncol=length(myplots), nrow=1, align="h")
    pdf_output=paste0(paste0(paste0(paste0(paste0(paste0("contaminant_",method_type),paste0("_negC_",negative_sample_type)),paste0("_compare_",name_compare_type)),"_top_"),paste0("thres_",test_threshold)),paste0(paste0("_",output_suffix),".pdf"))
    pdf(file=pdf_output, width=20, height=12)
    show(all_plot.all)
    dev.off()    
    #####################################################################################################
    #####################################################################################################
  }
}

#prune_phylo_KW: prune the phyloseq based on count cut off on a particular percentage of samples
prune_phylo_KW <- function (phyloseq, 
                            count_cut_off, #count cut off
                            percent_with_cut_off, #0.05 % of sample 
                            output_phyloseq
                            ) {                 
  
  #stopifnot(count_cut_off>0 & all.equal(count_cut_off, as.integer(count_cut_off)) ) #make sure count_cut_off is an integer that is larger than 0
  #stopifnot(percent_with_cut_off<1 & percent_with_cut_off>0) #make sure percent_with_cut_off is a percent (between 0 and 1)
  phyloseq_original<-phyloseq
  phyloseq.wh1 = genefilter_sample(phyloseq, filterfun_sample(function(x) x > count_cut_off), A = percent_with_cut_off * nsamples(phyloseq))
  output_phyloseq_T <- prune_taxa(phyloseq.wh1, phyloseq)
  print(output_phyloseq_T)
  assign(output_phyloseq,output_phyloseq_T,.GlobalEnv)
}  

##abundance_review_KW: look at the abundance of a particular ASV by sample_types 
abundance_review_KW <- function (input_phyloseq, 
                              sample_type_var_name, 
                              sample_types=list(),
                              count_type=c("absolute","relative"), 
                              graph_option=c("dotplot","violin","kernel","boxplot"),
                              ASV_review, 
                              output_name){
  #make sure graph option is within the choices
  graph_option <- match.arg(graph_option)
  print(graph_option)
  
  #make sure phyloseq only include the sample_types of interest
  keep_sample<- get_variable(input_phyloseq,sample_type_var_name) %in% sample_types
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  input_phyloseq_2 <- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
  
  for (i in 1:length(sample_types)) {
    assign(paste0("sampletype",i), sample_types[[i]])
  }
  
  if (count_type == "relative") {
          normalizeSample = function(x) { x/sum(x)}
          #getting the relative abundance
          countdata = data.frame(otu_table(transformSampleCounts(input_phyloseq_2, normalizeSample)))
  }          
  else if (count_type == "absolute") {          
          countdata <- as.data.frame(otu_table(input_phyloseq_2))
  }
  
  coldata<-as.data.frame(sample_data(input_phyloseq_2))
  taxa.table <- as.data.frame(tax_table(input_phyloseq_2))
  row.names(countdata)<-paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
  countdata= data.frame(lapply(countdata, function(x) as.numeric(as.character(x))),
                        check.names=F, row.names = rownames(countdata))
  coldata<- coldata[order(row.names(coldata)),]
  countdata <-countdata[, order(colnames(countdata))]
  #Subest the MetaData into the different sample type 
  for (i in 1:length(sample_types)) {
    assign(paste0("coldata_sampletype",i), coldata[coldata[[sample_type_var_name]] %in% c(get(paste0("sampletype",i))),])
  }
  for (i in 1:length(sample_types)) {
    wanted<-which(colnames(countdata) %in% row.names(get(paste0("coldata_sampletype",i))))
    assign(paste0("countdata",i), countdata[,wanted])
  }
  
  for (i in 1:length(sample_types)) {
    wanted<-which(rownames(get(paste0("countdata",i))) %in% ASV_review)
    assign(paste0("countdata_ASV",i), get(paste0("countdata",i))[wanted,])
    assign(paste0("countdata_ASVb",i), as.data.frame(t(as.matrix(get(paste0("countdata_ASV",i))))))
    assign(paste0("countdata_ASVb",i), cbind(get(paste0("countdata_ASVb",i)), get(paste0("sampletype",i))))
    d<-get(paste0("countdata_ASVb",i))
    colnames(d)<-c("abundance","sample_types")
    assign(paste0("countdata_ASVb",i),d)
  }
  combined_countdata<-countdata_ASVb1
  for (i in 2:length(sample_types)) {
    combined_countdata<- rbind(combined_countdata, get(paste0("countdata_ASVb",i)))
  }

  if (graph_option=="dotplot") {
    p<-ggplot(combined_countdata, aes(x=sample_types, y=abundance)) + 
      geom_dotplot(binwidth = 8, binaxis='y', dotsize = 5, stackdir='center')+
      coord_flip()
  }
  else if (graph_option=="violin") {
    p<-ggplot(combined_countdata, aes(x = factor(sample_types), y = log(abundance + 1))) +
      geom_violin(draw_quantiles = c(.25, .5, .75, .95))+
      coord_flip()
  }
  else if (graph_option=="kernel") {
    p<-ggplot(combined_countdata, 
              aes(x = log10(abundance*100 + 1), 
                  fill = sample_types)) + geom_density(alpha = 0.4) 
  }
  else if (graph_option=="boxplot"){
    for (i in 1:length(sample_types)) {
      assign(paste0("quant_ASVc",i), quantile(get(paste0("countdata_ASVb",i))$abundance)[4]+IQR(get(paste0("countdata_ASVb",i))$abundance)*1.7) #getting the 3rd quantile
    }
    quant_list <- quant_ASVc1
    for (i in 2:length(sample_types)) {
      quant_list<- append(quant_list,get(paste0("quant_ASVc",i))) #getting the 3rd quantile
    }
    upperlimit_boxplot<-max(quant_list)#determine where to set the limit of the axis so the boxplot is zoomed in
    p<-combined_countdata %>% 
      ggplot(aes(x = sample_types, y = abundance)) +
      stat_boxplot(geom ='errorbar', width = 0.02) +
      geom_boxplot(width = 0.8, outlier.shape = NA) + 
      coord_flip(ylim=c(0, upperlimit_boxplot))+
      ylab("Relative Abundance")
  }
  p<-p+ theme(plot.title=element_text(size=10,hjust=0.5, face="bold"),
              axis.title=element_text(face="bold", size=18),
              axis.text.x =element_text(size=18),
              axis.text.y =element_text(size=14))+
        xlab("Sample Type")+
        ggtitle(ASV_review)
  png_output<- paste0(output_name,".png")
  png(png_output, res = 300, width=250, height=200 , units='mm') #Plot of Loading vectors
  show(p)
  dev.off()
  csv_output<- paste0(output_name,".csv")
  write.csv(combined_countdata,csv_output)
}


###round p value 
###https://stackoverflow.com/questions/23018256/printing-p-values-with-0-001
pvalr <- function(pvals, sig.limit = .001, digits = 3) {
  roundr <- function(x, digits = 1) {
    res <- sprintf(paste0('%.', digits, 'f'), x)
    zzz <- paste0('0.', paste(rep('0', digits), collapse = ''))
    res[res == paste0('-', zzz)] <- zzz
    res
  }
  sapply(pvals, function(x, sig.limit) {
    if (x < sig.limit) return(sprintf('< %s', format(sig.limit)))
    if (x > .1) return(roundr(x, digits = 2)) else return(roundr(x, digits = digits))
  }, sig.limit = sig.limit)
}


##read_count_KW: boxplot for read count abundance by sample type
#read_count_KW (input_object=Transplant_PS_all, 
#                               sample_type_var_name="sample_types", 
#                               sample_types=c("Lower", "Upper", "BKG"), 
#                               sample_type_color=c("steelblue2", "purple2", "darkgoldenrod1"),
#                               pvalue="yes",
#                               output_name="read_count_all_nonegcontrol",
#                               width=100,  #dimension of the output png file
#                               height=150,
#                               font_size_ajust_percent=100) default is 100%, if want to decrease font size of axis then put a number less than 100
# this function check to see if the input_object is a phyloseq or a data.frame. if it is a data.frame, then you need to specify what the count_variable (contains read count) is call 
read_count_KW <- function(input_object, 
                          sample_type_var_name, 
                          sample_types=list(), 
                          sample_type_color=list(), 
                          p_value=c("yes","no"), 
                          output_name, 
                          width=100, 
                          height=150, 
                          xlabel, ylabel, 
                          xlabel_size=14, ylabel_size=14, axis_title_size=16,
                          y_log_scale=c("yes", "no"),
                          count_variable, 
                          pair_by_variable=NULL,...) {
  list2env(list(...), environment())
  #making sure the number of item in sample_types is same as sample_type_color
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  #making sure only yes and no is selected for p_value
  p_value <- match.arg(p_value)
  if (missing(y_log_scale)) {y_log_scale<-"yes"} else {y_log_scale<-tolower(y_log_scale)}
  y_log_scale<- match.arg(y_log_scale)
  #make sure that if paired comparison is wanted, there is only two groups being compared
  if (!missing(pair_by_variable)){stopifnot("Paired comparison only works with two groups comparison. Please limit sample_types to two elements"= length(sample_types)==2)}
  # Abundance of All Samples
  sample_type_color_t<-sample_type_color
  names(sample_type_color_t) <-sample_types 
  #setting up the label for the axis
  if (missing(xlabel)) { #if xlabel is missing then will have default  
    xlabel<-NULL
  } 
  if (missing(ylabel)) { #if ylabel is missing then will have default  
    if (y_log_scale=="yes"){
            ylabel<-"Total read count (log10)"  
    } else {
            ylabel<-"Total read count"
    }
  } 
  ### if the input_object is a phyloseq
  if(class(input_object)=="phyloseq"){
    input_phyloseq<-input_object
    #make sure input_phyloseq only include the sample_types of interest
    keep_sample<- get_variable(input_phyloseq,sample_type_var_name) %in% sample_types
    #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
    phyloseq_t <- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
    abundance_count<-colSums(abundances(phyloseq_t, transform = "identity")) 
    sample_df<- sample_data(phyloseq_t)
    #setting up input_dataframe, which will be used to create the statistics table
    sample_df2<-data.frame(sample_df)
    sample_df3<-sample_df2[c(sample_type_var_name,pair_by_variable)] # if pair_by_variable is not provided, then sample_df3 will only have sample_type_var_name
    input_dataframe <- merge(abundance_count, sample_df3,by="row.names")
    input_dataframe$abundance_count<-input_dataframe$x
  }else if (class(input_object)=="data.frame"){ #if input is a dataframe 
    #FYI: if input_object is a data.frame, then the data needs to be in long format. It needs to contain 2 variables (input the name in: sample_type_var_name and count_variable) within the data.frame"
    if (missing(count_variable)) {stop("need to specify: count_variable")}
    input_dataframe<-input_object
    input_dataframe$abundance_count <- input_dataframe[[count_variable]]
  }else {#if the input is neither data.frame nor physloeq
    stop("this function only works on phyloseq object or data frame")
  }
  input_dataframe$var1<-input_dataframe[[sample_type_var_name]] #standardizing the variable name. var1 would become the sample type separator
  #if paired comparison is wanted 
  if (!missing(pair_by_variable)){ 
    input_dataframe$group_var<-input_dataframe[[pair_by_variable]]
    input_dataframe<-input_dataframe[order(input_dataframe$var1, input_dataframe$group_var),]    
  }
  #turn into a factor with the same order that the function input has
  input_dataframe[["var1"]] <- factor(input_dataframe[["var1"]], sample_types)
  
  #setting up the graph
  p<-ggplot(input_dataframe, aes(x=var1, y=abundance_count, fill=var1)) + 
    stat_boxplot(geom = "errorbar",width = 0.3, color = "gray70")+
    geom_boxplot(outlier.shape = NA)
  #if paired comparison is wanted 
  if (!missing(pair_by_variable)){
     p <- p + geom_line(aes(group = group_var), color="gray70", alpha=0.5) +
       geom_jitter(width=0, size=1, shape= 1, color="azure4") #if there is paired line, will have jitter to be 0 so is not as messy
  } else { #if paired comparison is not needed then will not have lines connect betwen groups
     p <- p + geom_jitter(width=0.3, size=1, shape= 1, color="azure4") #if no paired line, then will have jitter 
  }
  #get csv output of the statistics
  stat_output_csv<- input_dataframe %>% dplyr::group_by(var1) %>% dplyr::summarise( mean = mean(abundance_count),
                                                                                    stddev = sd(abundance_count),
                                                                                    median = median(abundance_count),
                                                                                    val_min = min(abundance_count),
                                                                                    val_max = max(abundance_count),
                                                                                    per10 = as.numeric(quantile(abundance_count, prob = c(0.10))),
                                                                                    per25 = as.numeric(quantile(abundance_count, prob = c(0.25))),
                                                                                    per75 = as.numeric(quantile(abundance_count, prob = c(0.75))),
                                                                                    per90 = as.numeric(quantile(abundance_count, prob = c(0.90))))
  output_csv=paste0(output_name,"_stats.csv")
  write.csv(stat_output_csv,output_csv, row.names=FALSE)
  #continue with setting up graph 
  #adding the aesthetics for the graphs 
  p <- p + theme_classic() +
    scale_fill_manual(values=sample_type_color_t) +
    scale_x_discrete(labels=sample_types) +
    labs(title=NULL,
         x = xlabel, 
         y= ylabel) +
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          axis.title=element_text(face="bold", size=axis_title_size),
          axis.text.x =element_text(size=xlabel_size),
          axis.text.y =element_text(size=ylabel_size),
          legend.position="none")
  if (y_log_scale=="yes"){ #if want y axis in log scale
    p <- p  + scale_y_log10()  
  } 
  if (p_value=="yes") { ##if want p_value included int he diagram
    #if paired comparison is wanted 
    if (!missing(pair_by_variable)){
        stat.test <- input_dataframe %>%
              wilcox_test(abundance_count ~ var1, paired=T) %>% add_significance()
      }else { # if no paired comparison 
          stat.test <- input_dataframe %>%
              wilcox_test(abundance_count ~ var1) %>% add_significance()
      }
    if (y_log_scale=="yes"){
      stat.test<- stat.test %>% add_xy_position(x = "var1", y.trans = function(x){log10(x)+0.1})
    } else {
      stat.test<- stat.test %>% add_xy_position(x = "var1", y.trans = function(x){x*1.1})
    }
    if (!exists("step.increase")) {step.increase<-0.05} 
    #limiting p value to 3 digits after decimal
    stat.test$p<-pvalr(stat.test$p, digits = 3)
    p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.02, step.increase = step.increase, inherit.aes = F, size=3)
  }
  ###determine sample_types string length
  string_length <- max(nchar(as.character(input_dataframe$var1))) 
  if (string_length>20) { #if the sample type string is too long then will title the x axis text 
    p<-p+ theme(axis.text.x = element_text(angle = 8, hjust = 0.5, vjust=0.5))
  }
  png_output<-paste0(output_name,".png") 
  png(file=png_output, res = 300, width=width, height=height , units='mm')
        show(p)
  dev.off()
}

##read_count_v2_KW: boxplot for read count abundance by sample type
#read_count_v2_KW (input_object=Transplant_PS_all, 
#                               sample_type_var_name="sample_types", 
#                               sample_types=c("Lower", "Upper", "BKG"), 
#                               sample_type_color=c("steelblue2", "purple2", "darkgoldenrod1"),
#                               pvalue="yes",
#                               output_name="read_count_all_nonegcontrol",
#                               width=100,  #dimension of the output png file
#                               height=150,
#                               font_size_ajust_percent=100) default is 100%, if want to decrease font size of axis then put a number less than 100
# this function check to see if the input_object is a phyloseq or a data.frame. if it is a data.frame, then you need to specify what the count_variable (contains read count) is call 
read_count_v2_KW <- function(input_object, 
                          sample_type_var_name, 
                          sample_types=list(), 
                          sample_type_color=list(), 
                          sample_type2_var_name, 
                          sample_types2=list(), 
                          sample_type2_color=list(), 
                          p_value=c("yes","no"), 
                          p_value_select=c("inter","intra","both"),
                          output_name, 
                          width=100, height=150, 
                          xlabel_size=14, ylabel_size=14, axis_title_size=16,
                          xlabel, ylabel,
                          y_log_scale=c("yes", "no"),
                          count_variable, ...) {
  list2env(list(...), environment())
  #making sure the number of item in sample_types is same as sample_type_color
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  #making sure only yes and no is selected for p_value
  p_value <- match.arg(p_value)
  p_value_select <- match.arg(p_value_select)
  if (missing(y_log_scale)) {y_log_scale<-"yes"} else {y_log_scale<-tolower(y_log_scale)}
  y_log_scale<- match.arg(y_log_scale)
  if(missing(p_value_select)){p_value_select<-"both"} #default p_value_select is both, which calculates for inter and intra-group 
  # set up coloring for sample_type_var_name
  sample_type_color_t<-sample_type_color
  names(sample_type_color_t) <-sample_types 
  # set up coloring for sample_type2_var_name
  sample_type2_color_t<-sample_type2_color
  names(sample_type2_color_t) <-sample_types2 
  #setting up the label for the axis
  if (missing(xlabel)) { #if xlabel is missing then will have default  
    xlabel<-NULL
  } 
  if (missing(ylabel)) { #if ylabel is missing then will have default  
    if (y_log_scale=="yes"){
      ylabel<-"Total read count (log10)"  
    } else {
      ylabel<-"Total read count"
    }
  } 
  ### if the input_object is a phyloseq\
  ### setting up input_dataframe which would make sure the input dataframe into ggplot is consistent, regardless of using phyloseq or data.frame
  if(class(input_object)=="phyloseq"){
    input_phyloseq<-input_object
    #make sure input_phyloseq only include the sample_types of interest
    keep_sample<- get_variable(input_phyloseq,sample_type_var_name) %in% sample_types
    #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
    phyloseq_t <- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
    abundance_count<-colSums(abundances(phyloseq_t, transform = "identity")) 
    sample_df<- sample_data(phyloseq_t)
    #setting up input_dataframe, which will be used to create the statistics table
    sample_df2<-data.frame(sample_df)
    sample_df3<-sample_df2[c(sample_type_var_name,sample_type2_var_name)]
    input_dataframe <- merge(abundance_count, sample_df3,by="row.names")
    input_dataframe$abundance_count<-input_dataframe$x
  }else if (class(input_object)=="data.frame"){ #if input is a dataframe 
    #FYI: if input_object is a data.frame, then the data needs to be in long format. It needs to contain 2 variables (input the name in: sample_type_var_name and count_variable) within the data.frame"
    if (missing(count_variable)) {stop("need to specify: count_variable")}
    input_dataframe<-input_object
    input_dataframe$abundance_count <- input_dataframe[[count_variable]]
  }else {#if the input is neither data.frame nor physloeq
    stop("this function only works on phyloseq object or data frame")
  }
  ### setting up sample_type_var_name to be var1 and sample_type2_var_name to be var2 
  input_dataframe$var1<-as.factor(input_dataframe[[sample_type_var_name]])
  input_dataframe$var2<-as.factor(input_dataframe[[sample_type2_var_name]])
  input_dataframe2<-input_dataframe
  #turn into a factor with the same order that the function input has
  input_dataframe2[["var1"]] <- factor(input_dataframe2[["var1"]], sample_types)
  input_dataframe2[["var2"]] <- factor(input_dataframe2[["var2"]], sample_types2)
  
  p<-ggplot(input_dataframe2, aes(x=var1, y=abundance_count, 
                                  color=var1, fill=var2))+
            stat_boxplot(geom = "errorbar",width = 0.2, position = position_dodge(width=0.7)) +
            geom_boxplot(width=0.6, outlier.shape = NA, position= position_dodge(width=0.7)) + 
            geom_jitter(position=position_jitterdodge(), size=1, shape= 1, color="azure4") 
  #get csv output of the statistics
  stat_output_csv<- input_dataframe2 %>% dplyr::group_by(var1, var2) %>% dplyr::summarise( mean = mean(abundance_count),
                                                                                          stddev = sd(abundance_count),
                                                                                          median = median(abundance_count),
                                                                                          val_min = min(abundance_count),
                                                                                          val_max = max(abundance_count),
                                                                                          per10 = as.numeric(quantile(abundance_count, prob = c(0.10))),
                                                                                          per25 = as.numeric(quantile(abundance_count, prob = c(0.25))),
                                                                                          per75 = as.numeric(quantile(abundance_count, prob = c(0.75))),
                                                                                          per90 = as.numeric(quantile(abundance_count, prob = c(0.90))))
  output_csv=paste0(output_name,"_stats.csv")
  write.csv(stat_output_csv,output_csv, row.names=FALSE)
  #calculate p value
  if (p_value=="yes") { ##if want p_value included int he diagram
        #depending if inter, intra or both is selected for p_value_select
        #INTRA group comparson for var1
        stat.test <- input_dataframe2 %>%
                dplyr::group_by(var1) %>% 
                wilcox_test(abundance_count~ var2) %>%
                add_significance()
        if (p_value_select=="intra" | p_value_select=="both" ){
          if (y_log_scale=="yes"){
            stat.test <- stat.test %>% add_xy_position(x = "var1", y.trans = function(x){log10(x)})
          } else {
            stat.test <- stat.test %>% add_xy_position(x = "var1")
          }
        }
        #INTER grioup comparison for var1
        stat.test2 <- input_dataframe2 %>%
                wilcox_test(abundance_count ~ var1) %>% add_significance()
        if (p_value_select=="inter"){
          if (y_log_scale=="yes"){
            stat.test2<- stat.test2 %>% add_xy_position(x = "var1", y.trans = function(x){log10(x)})
          } else {
            stat.test2<- stat.test2 %>% add_xy_position(x = "var1")
          }
        } else if (p_value_select=="both"){
          if (y_log_scale=="yes"){
            stat.test2<- stat.test2 %>% add_xy_position(x = "var1", y.trans = function(x){log10(x)+0.3}) #make sure the inter and intra group p values dont overlap 
          } else {
            stat.test2<- stat.test2 %>% add_xy_position(x = "var1", y.trans = function(x){x*1.1})
          }
        }
  }
  #continue with setting up graph 
  #adding the aesthetics for the graphs 
  p <- p + theme_classic() +
    scale_color_manual(values=sample_type_color_t) +
    scale_fill_manual(values=sample_type2_color_t) +
    scale_x_discrete(labels=sample_types) +
    labs(title=NULL, x = xlabel, y= ylabel) +
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          axis.title=element_text(face="bold", size=axis_title_size),
          axis.text.x =element_text(size=xlabel_size),
          axis.text.y =element_text(size=ylabel_size),
          legend.position="none")
  if (y_log_scale=="yes"){
    p <- p  + scale_y_log10()  
  } 
  #add on the p value now 
  if (p_value=="yes") { ##if want p_value included int he diagram
    if (p_value_select=="intra" | p_value_select=="both" ){ #add on intragroup comparson 
        stat.test$p<-pvalr(stat.test$p, digits = 3)
        p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0, inherit.aes = F, size=3)
    }
    if (p_value_select=="inter" | p_value_select=="both" ){ #add on intergroup comparson
        if (!exists("step.increase")) {step.increase<-0.05} 
        stat.test2$p<-pvalr(stat.test2$p, digits = 3)
        p <- p + stat_pvalue_manual(stat.test2, label = "p", tip.length = 0.02, step.increase = step.increase, inherit.aes = F, size=3)
    }
  }
  ###determine sample_types_var_name string length
  string_length <- max(nchar(as.character(input_dataframe$var1))) 
  if (string_length>18) { #if the sample type string is too long then will title the x axis text 
         p<-p+ theme(axis.text.x = element_text(angle = 8, hjust = 0.5, vjust=0.5))
  }
  #save output as png file 
  png_output<-paste0(output_name,".png") 
  png(file=png_output, res = 300, width=width, height=height , units='mm')
        show(p)
  dev.off()
}

##alpha_diversity_KW: calculates the alpha diversity (shannon)
#alpha_diversity_KW (input_phyloseq=Transplant_PS_all, 
#                               sample_type_var_name="sample_type", 
#                               sample_types=c("Lower", "Upper", "BKG"), 
#                               sample_type_color=c("steelblue2", "purple2", "darkgoldenrod1"), 
#                               p_value="yes",
#                               output_name="alpha_diversity_all_nonegcontrol",
#                               width=100, #dimension of the output png file
#                               height=150,
#                               font_size_ajust_percent=100) default is 100%, if want to decrease font size of axis then put a number less than 100 
alpha_diversity_KW <- function(input_phyloseq, 
                               sample_type_var_name, 
                               sample_types=list(), 
                               sample_type_color=list(), 
                               p_value=c("yes","no"), 
                               output_name, 
                               width=100, height=150, 
                               xlabel_size=14, ylabel_size=14, axis_title_size=16,
                               xlabel, ylabel,                          
                               pair_by_variable=NULL, ...) {
  list2env(list(...), environment())
  #making sure the number of item in sample_types is same as sample_type_color
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  #making sure only yes and no is selected for p_value
  p_value <- match.arg(p_value)
  #make sure that if paired comparison is wanted, there is only two groups being compared
  if (!missing(pair_by_variable)){stopifnot("Paired comparison only works with two groups comparison. Please limit sample_types to two elements"= length(sample_types)==2)}
  # ALPHA DIVERSITY of All Samples
  sample_type_color_t<-sample_type_color
  names(sample_type_color_t) <-sample_types 
  #make sure input_phyloseq only include the sample_types of interest
  keep_sample<- get_variable(input_phyloseq,sample_type_var_name) %in% sample_types
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  phyloseq_t <- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
  
  Shannon_diversity <- data.frame(vegan::diversity(otu_table(phyloseq_t), index = "shannon", MARGIN = 2, base = exp(1)))
  colnames(Shannon_diversity)[1]<-"Shannon_diversity"
  
  if (missing(xlabel)) { #if xlabel is missing then will have default  
    xlabel<-NULL
  } 
  if (missing(ylabel)) { #if ylabel is missing then will have default  
    ylabel<-"Shannon diversity index"
  }
  alpha_dataframe<-phyloseq_t@sam_data
  alpha_dataframe<- merge(alpha_dataframe, Shannon_diversity, by="row.names")
  input_dataframe<- alpha_dataframe %>% dplyr::select(c(sample_type_var_name, "Shannon_diversity",pair_by_variable))# if pair_by_variable is not provided, then input_dataframe only have Shannon_diversity and sample_type_var_name
  #standardized variables 
  input_dataframe$var1<-input_dataframe[[sample_type_var_name]] #standardizing the variable name. var1 would become the sample type separator
  if (!missing(pair_by_variable)){ 
    input_dataframe$group_var<-input_dataframe[[pair_by_variable]]
    input_dataframe<-input_dataframe[order(input_dataframe$var1, input_dataframe$group_var),]    
  }
  #turn into a factor with the same order that the function input has
  input_dataframe[["var1"]] <- factor(input_dataframe[["var1"]], sample_types)
  
  #Setting up the graphs
  p<-ggplot(input_dataframe, aes(var1, Shannon_diversity, fill=var1)) +
    stat_boxplot(geom = "errorbar",width = 0.3, color = "gray70")+
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values=sample_type_color_t) +
    scale_x_discrete(labels=sample_types) +
    scale_y_continuous() + 
    theme_classic() +
    labs(title=NULL,
         x = xlabel, 
         y= ylabel) +
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          axis.title=element_text(face="bold", size=axis_title_size),
          axis.text.x =element_text(size=xlabel_size),
          axis.text.y =element_text(size=ylabel_size),
          legend.position="none")
  if (!missing(pair_by_variable)){
    p <- p + geom_line(aes(group = group_var), color="gray70", alpha=0.5) +
      geom_jitter(width=0, size=1, shape= 1, color="azure4") #if there is paired line, will have jitter to be 0 so is not as messy
  } else { #if paired comparison is not needed then will not have lines connect betwen groups
    p <- p + geom_jitter(width=0.3, size=1, shape= 1, color="azure4") #if no paired line, then will have jitter 
  }
  if (p_value=="yes") { ##if want p_value included in the diagram
    if (!missing(pair_by_variable)){
      stat.test <- input_dataframe %>%
        wilcox_test(Shannon_diversity ~ var1, paired=T) %>% add_significance()
    }else { # if no paired comparison 
      stat.test <- input_dataframe %>%
        wilcox_test(Shannon_diversity ~ var1) %>% add_significance()
    }
    stat.test<- stat.test %>% add_xy_position(x = "var1", y.trans = function(x){x*1.1})
    if (!exists("step.increase")) {step.increase<-0.05} 
    #limiting p value to 3 digits after decimal
    stat.test$p<-pvalr(stat.test$p, digits = 3)
    p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.02, step.increase = step.increase,  inherit.aes = F, size=3)
  }
  ###determine sample_types string length
  string_length <- max(nchar(as.character(input_dataframe$var1)))
  if (string_length>20) { #if the sample type string is too long then will title the x axis text 
    p<-p+ theme(axis.text.x = element_text(angle = 8, hjust = 0.5, vjust=0.5))
  }
  png_output<-paste0(output_name,".png")
  png(file=png_output, res = 300, width=width, height=height, units='mm')
      show(p)
  dev.off()
  csv_output<-paste0(output_name,".csv")
  write.csv(Shannon_diversity, csv_output)
}

alpha_diversity_v2_KW <- function(input_phyloseq, 
                                  sample_type_var_name, 
                                  sample_types=list(), 
                                  sample_type_color=list(), 
                                  sample_type2_var_name, 
                                  sample_types2=list(), 
                                  sample_type2_color=list(), 
                                  p_value=c("yes","no"), 
                                  p_value_select=c("inter","intra","both"),
                                  output_name, 
                                  width=100, height=150, 
                                  xlabel_size=14, ylabel_size=14, axis_title_size=16,
                                  xlabel, ylabel, ...) {
  list2env(list(...), environment())
  #making sure the number of item in sample_types is same as sample_type_color
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  #making sure only yes and no is selected for p_value
  p_value <- match.arg(p_value)
  p_value_select <- match.arg(p_value_select)
  if(missing(p_value_select)){p_value_select<-"both"} #default p_value_select is both, which calculates for inter and intra-group 
  # set up coloring for sample_type_var_name
  sample_type_color_t<-sample_type_color
  names(sample_type_color_t) <-sample_types 
  # set up coloring for sample_type2_var_name
  sample_type2_color_t<-sample_type2_color
  names(sample_type2_color_t) <-sample_types2 
  #setting up the label for the axis
  if (missing(xlabel)) { #if xlabel is missing then will have default  
    xlabel<-NULL
  } 
  if (missing(ylabel)) { #if ylabel is missing then will have default  
    ylabel<-"Shannon Diversity Index"  
  } 
  
  #make sure input_phyloseq only include the sample_types of interest
  keep_sample<- (get_variable(input_phyloseq,sample_type_var_name) %in% sample_types &
                   get_variable(input_phyloseq,sample_type2_var_name) %in% sample_types2)
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  phyloseq_t <- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
  #getting Shannon diversity
  Shannon_diversity <- vegan::diversity(otu_table(phyloseq_t), index = "shannon", MARGIN = 2, base = exp(1))
  Shannon_diversity<-data.frame(Shannon_diversity)
  colnames(Shannon_diversity)[1]<-"Shannon_diversity"
  alpha_dataframe<-phyloseq_t@sam_data
  alpha_dataframe<- merge(alpha_dataframe, Shannon_diversity, by="row.names")
  input_dataframe<- alpha_dataframe %>% dplyr::select(c(sample_type_var_name,sample_type2_var_name, "Shannon_diversity"))
  ### setting up sample_type_var_name to be var1 and sample_type2_var_name to be var2 
  input_dataframe$var1<-as.factor(input_dataframe[[sample_type_var_name]])
  input_dataframe$var2<-as.factor(input_dataframe[[sample_type2_var_name]])
  input_dataframe2<-input_dataframe
  #turn into a factor with the same order that the function input has
  input_dataframe2[["var2"]] <- factor(input_dataframe2[["var2"]], sample_types2)
  input_dataframe2[["var2"]] <- factor(input_dataframe2[["var2"]], sample_types2)
  
  p<-ggplot(input_dataframe2, aes(x=var1, y=Shannon_diversity, color=var1, fill=var2))+
    stat_boxplot(geom = "errorbar",width = 0.2, position = position_dodge(width=0.7)) +
    geom_boxplot(width=0.6, outlier.shape = NA, position= position_dodge(width=0.7)) + 
    geom_jitter(position=position_jitterdodge(), size=1, shape= 1, color="azure4") 
  #get csv output of the statistics
  stat_output_csv<- input_dataframe %>% dplyr::group_by(var1, var2) %>% dplyr::summarise( mean = mean(Shannon_diversity),
                                                                                          stddev = sd(Shannon_diversity),
                                                                                          median = median(Shannon_diversity),
                                                                                          val_min = min(Shannon_diversity),
                                                                                          val_max = max(Shannon_diversity),
                                                                                          per10 = as.numeric(quantile(Shannon_diversity, prob = c(0.10))),
                                                                                          per25 = as.numeric(quantile(Shannon_diversity, prob = c(0.25))),
                                                                                          per75 = as.numeric(quantile(Shannon_diversity, prob = c(0.75))),
                                                                                          per90 = as.numeric(quantile(Shannon_diversity, prob = c(0.90))))
  output_csv=paste0(output_name,"_stats.csv")
  write.csv(stat_output_csv,output_csv, row.names=FALSE)
  #calculate p value
  #depending if inter, intra or both is selected for p_value_select
  #INTRA group comparson for var1
  stat.test <- input_dataframe2 %>%
    dplyr::group_by(var1) %>% 
    wilcox_test(Shannon_diversity ~ var2) %>% add_significance()
  if (p_value_select=="intra" | p_value_select=="both" ){
    stat.test <- stat.test %>% add_xy_position(x = "var1")
  }
  #INTER grioup comparison for var1
  stat.test2 <- input_dataframe2 %>%
    wilcox_test(Shannon_diversity ~ var1) %>%  add_significance()
  if (p_value_select=="inter" | p_value_select=="both"){
    stat.test2<- stat.test2 %>% add_xy_position(x = "var1", y.trans = function(x){x*1.1})
  } 
  
  #continue with setting up graph 
  #adding the aesthetics for the graphs 
  p <- p + theme_classic() +
    scale_color_manual(values=sample_type_color_t) +
    scale_fill_manual(values=sample_type2_color_t) +
    scale_x_discrete(labels=sample_types) +
    labs(title=NULL, x = xlabel, y= ylabel) +
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          axis.title=element_text(face="bold", size=axis_title_size),
          axis.text.x =element_text(size=xlabel_size),
          axis.text.y =element_text(size=ylabel_size),
          legend.position="none")
  #add on the p value now 
  if (p_value=="yes") { ##if want p_value included int he diagram
    if (p_value_select=="intra" | p_value_select=="both" ){ #add on intragroup comparson 
      stat.test$p<-pvalr(stat.test$p, digits = 3)
      p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0, inherit.aes = F, size=3)
    }
    if (p_value_select=="inter" | p_value_select=="both" ){ #add on intergroup comparson
      if (!exists("step.increase")) {step.increase<-0.05} 
      stat.test2$p<-pvalr(stat.test2$p, digits = 3)
      p <- p + stat_pvalue_manual(stat.test2, label = "p", tip.length = 0.02, step.increase = step.increase, inherit.aes = F, size=3)
    }
  }
  ###determine sample_types_var_name string length
  string_length <- max(nchar(as.character(input_dataframe$var1))) 
  if (string_length>18) { #if the sample type string is too long then will title the x axis text 
    p<-p+ theme(axis.text.x = element_text(angle = 8, hjust = 0.5, vjust=0.5))
  }
  #save output as png file 
  png_output<-paste0(output_name,".png") 
  png(file=png_output, res = 300, width=width, height=height , units='mm')
  show(p)
  dev.off()
}

##beta_diversity_KW: calculates the Beta diversity (Adonis)
#beta_diversity_KW (input_phyloseq=Transplant_PS_all, 
#                               sample_type_var_name="sample_type", 
#                               sample_types=c("Lower", "Upper", "BKG"), 
#                               sample_type_color=c("steelblue2", "purple2", "darkgoldenrod1"), 
#                               p_value="yes",
#                               output_name="beta_diversity_all_nonegcontrol",
#                               width=200, #dimension of the output png file
#                               height=200,
#                               font_size_adjust_percent=100) default is 100%, if want to decrease font size of axis then put a number less than 100 
beta_diversity_KW <- function(input_phyloseq, 
                              sample_type_var_name, 
                              sample_types=list(), 
                              sample_type_color=list(), 
                              sample_type2_var_name,   #secondary variable to determine the color of the individual dots  
                              sample_types2=list(),
                              sample_type2_color=list(),
                              p_value=c("yes","no"),
                              output_name, 
                              DMM_cluster_assign=NULL,
                              DMM_cluster_color=list(),
                              width=200, height=200,
                              plot_title_size=14, axis_title_size=18, label_size=7,
                              x_axis_flip, 
                              beta_method=c("bray", "unifrac"),
                              p_value_location=c("TR","TL","BR","BL")) { #p_value_location: TR-top right; TL- top left; BR- bottom right; BL- bottom left
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  if(missing(x_axis_flip)) { x_axis_flip<-"no"} #default being no for x_axis_flip
  if(missing(beta_method)) { beta_method<-"bray"} #default being bray curtis
  if(missing(p_value_location) & tolower(p_value)=="yes"){
    p_value_location<-"BL" #default location for p_value will be in the bottom left corner
  } else if(!missing(p_value_location) & tolower(p_value)=="yes"){ #if p_value is wanted and p_value_location is not missing
    print("p_value_location: TR-top right; TL- top left; BR- bottom right; BL- bottom left")
    p_value_location<-match.arg(p_value_location)
    print(paste0("selected:",p_value_location))
  }
  #making sure the number of item in sample_types is same as sample_type_color
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  if (!missing(sample_type2_var_name)) {
          stopifnot("sample_types2 need to have same number of elements as sample_types2_color"= length(sample_types2)==length(sample_type2_color))
  }
  #making sure only yes and no is selected for p_value
  p_value <- match.arg(p_value)
  #Beta diversity
  #make sure input_phyloseq only include the sample_types of interest
  keep_sample<- get_variable(input_phyloseq,sample_type_var_name) %in% sample_types
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  phyloseq_t <- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
  sample_type_color_t<-sample_type_color
  names(sample_type_color_t) <-sample_types
  sample_type_sort<-sort(sample_types)
  #setting up the secondary color for the points
  if (!missing(sample_type2_var_name)) {
          sample_type2_color_t<-sample_type2_color
          names(sample_type2_color_t) <-sample_types2
  }
  if (beta_method=="bray"){
      #Create Distance Matrix with Bray
      vegdist=vegdist(t(otu_table(phyloseq_t)), method = "bray") #row is samplesID #column is taxa
  } else {
      set.seed(123)
      vegdist= phyloseq::UniFrac(phyloseq_t, weighted = TRUE) 
  }
  #Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
  CmdScale <- cmdscale(vegdist, k =10)
  #calculated Sample variance for each PC
  vars <- apply(CmdScale, 2, var)
  #Create Variable with the Percent Variance
  percentVar <- round(100 * (vars/sum(vars)))
  #Merge PC Data with MetaData
  newResults <- merge(x = CmdScale, y = sample_data(phyloseq_t), by = "row.names", all.x = TRUE)  #sample_data(phyloseq_t) row is the sampleID
  #Rename Variables for PC1 and PC2
  colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
  colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
  colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"
  formula <- as.formula(paste("cbind(PC1,PC2)",sample_type_var_name, sep=" ~ "))
  #Calculate the Centroid Value
  centroids <- aggregate(formula ,data= newResults, mean)
  #Merge the Centroid Data into the PCOA Data
  newResults <- merge(newResults,centroids,by=sample_type_var_name,suffixes=c("",".centroid"))
  # Calculate p-value with ADONIS
  x=adonis2(vegdist ~ phyloseq_t@sam_data[[sample_type_var_name]])
  pvalues<-x$`Pr(>F)`[[1]]
  subtitle_output<- paste0("Adonis, p=",pvalues)
  png_output<-paste0(output_name,".png")
  if (!missing(DMM_cluster_assign)) { #if DMM_cluster_assign option is selected, determine if the DMM_cluster_assign object exist
    if (exists(DMM_cluster_assign)=="FALSE"){
      print(paste0("DMM clustering assignment",paste0(DMM_cluster_assign," does not exist")))
      stop("please run function: decontaminant_subplot_KW prior to using this option- decontam") 
    } else {
      #setting up cluster color: 5 default color is set 
      #if DMM_cluster_color is missing then the default color will be used. otherwise, will use the DMM_cluster_color which is indicated
      cluster_color_t<-c("cornflowerblue","hotpink1","cyan4","goldenrod3",'azure4')
      names(cluster_color_t) <- c("1","2","3","4","5")
      if (!missing(DMM_cluster_color)) #if DMM_cluster_color is not missing then will replace the color 
          { for(i in 1:length(DMM_cluster_color)) {
            cluster_color_t[i]<-DMM_cluster_color[i] #replace the selected default color with selected 
            }
      } 
      #add the cluster color to the sample_type color 
      sample_type_color_t<- c(sample_type_color_t, cluster_color_t)
      #merging the DMM clustering assignment
      newResults<-merge(newResults,get(DMM_cluster_assign),by.x="name", by.y="row.names", all=T)
      p<-ggplot(newResults, aes(PC1, PC2)) + # Graph PC1 and PC2
        geom_point(size=2, aes(color=Cluster_num)) + # Set the size of the points
        xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2 
        labs(title=NULL, 
             fill="Centroids") +
        geom_segment(alpha=0.5,aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= get(sample_type_var_name)))+ 
        geom_point(data=centroids, aes(x=PC1, y=PC2, color= get(sample_type_var_name)), size=0) + 
        geom_label_repel(data = centroids, aes(x=PC1, y=PC2, color= get(sample_type_var_name), label=sample_type_sort), size=label_size) +
        scale_color_manual(name=NULL,  
                           values=sample_type_color_t) +  
        theme(panel.background = element_blank(),
              panel.border=element_rect(fill=NA),
              panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
              panel.grid.minor = element_blank(),strip.background=element_blank(),
              plot.title=element_text(face="bold",hjust=0.5, size = plot_title_size), 
              plot.subtitle = element_text(hjust=0.5),
              axis.title=element_text(face="bold", size =axis_title_size),
              axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
              axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
              axis.ticks=element_blank(),
              plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
    }
  } else { # if DMM_cluster_assign is missing
    if (!(missing(sample_type2_var_name)) & !(missing(sample_type2_var_name)) & !(missing(sample_type2_color))) { #if secondary color for points is needed
                #add the secondary dot color to sample_type color 
                sample_type_color_t<- c(sample_type_color_t, sample_type2_color_t)
                p<-ggplot(newResults, aes(PC1, PC2, color= get(sample_type_var_name))) + # Graph PC1 and PC2
                  geom_point(size=2, aes(color=get(sample_type2_var_name))) + # Set the size of the points
                  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
                  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2 
                  labs(title=NULL, 
                       fill="Centroids") +
                  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= get(sample_type_var_name)))+ 
                  geom_point(data=centroids, aes(x=PC1, y=PC2, color= get(sample_type_var_name)), size=0) + 
                  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=sample_type_sort), size=label_size) +
                  scale_color_manual(name=NULL,  
                                     values=sample_type_color_t) +
                  theme(panel.background = element_blank(),
                        panel.border=element_rect(fill=NA),
                        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
                        panel.grid.minor = element_blank(),strip.background=element_blank(),
                        plot.title=element_text(face="bold",hjust=0.5, size = plot_title_size), 
                        plot.subtitle = element_text(hjust=0.5),
                        axis.title=element_text(face="bold", size = axis_title_size),
                        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
                        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
                        axis.ticks=element_blank(),
                        plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
    } else {
                p<-ggplot(newResults, aes(PC1, PC2, color= get(sample_type_var_name))) + # Graph PC1 and PC2
                  geom_point(size=2) + # Set the size of the points
                  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
                  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2 
                  labs(title=NULL, 
                       fill="Centroids") +
                  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= get(sample_type_var_name)))+ 
                  geom_point(data=centroids, aes(x=PC1, y=PC2, color= get(sample_type_var_name)), size=0) + 
                  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=sample_type_sort), size=label_size) +
                  scale_color_manual(name=NULL,  
                                     values=sample_type_color_t) +
                  theme(panel.background = element_blank(),
                        panel.border=element_rect(fill=NA),
                        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
                        panel.grid.minor = element_blank(),strip.background=element_blank(),
                        plot.title=element_text(face="bold",hjust=0.5, size = plot_title_size), 
                        plot.subtitle = element_text(hjust=0.5),
                        axis.title=element_text(face="bold", size = axis_title_size),
                        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
                        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
                        axis.ticks=element_blank(),
                        plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
    }
  }  
  newResults<<-newResults
  if (tolower(x_axis_flip)=="yes") {
    p<-p+scale_x_reverse()
  }
  if (p_value=="yes") { ###when the axis is flip, the annotation for the p value need to placed on max(newResults$PC1), instead of min(newResults$PC2))
    if (tolower(x_axis_flip)=="no"){
            if (p_value_location=="BR"){
                  p <- p + annotate("text", x=max(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
            } else if (p_value_location=="BL"){
                  p <- p + annotate("text", x=min(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
            } else if (p_value_location=="TR"){
                  p <- p + annotate("text", x=max(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
            } else if (p_value_location=="TL"){
                  p <- p + annotate("text", x=min(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
            }
    } else if (tolower(x_axis_flip)=="yes") {
            if (p_value_location=="BR"){
              p <- p + annotate("text", x=min(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
            } else if (p_value_location=="BL"){
              p <- p + annotate("text", x=max(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
            } else if (p_value_location=="TR"){
              p <- p + annotate("text", x=min(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
            } else if (p_value_location=="TL"){
              p <- p + annotate("text", x=max(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
            }
    }
  }
  png(file=png_output, res = 300, width=width, height=height , units='mm')
      show(p)
  dev.off()
}

##PLS_DA_KW: PLS-DA component 1 loading 
#PLS_DA_KW (input_phyloseq=Transplant_PS_all, 
#                               sample_type_var_name="sample_type", 
#                               sample_types=c("Lower", "Upper", "BKG"), 
#                               sample_type_color=c("steelblue2", "purple2", "darkgoldenrod1"), 
#                               top_display_loading=25,       ### this option let you select how many taxa to include in the graph
#                               output_name="plsda_upper_lower",
#                               width=200, #dimension of the output png file
#                               height=200
#                               decontam="Contam_list_NC_BKG_compare_Lower") 
PLS_DA_KW <- function(input_phyloseq, 
                      sample_type_var_name, 
                      sample_types=list(), 
                      sample_type_color=list(), 
                      top_display_loading, 
                      output_name, 
                      width=500, height=500, 
                      decontam=NULL) {
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}

  #making sure the number of item in sample_types is same as sample_type_color
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  #make sure input_phyloseq only include the sample_types of interest
  keep_sample<- get_variable(input_phyloseq,sample_type_var_name) %in% sample_types
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  phyloseq_t <- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
  
  sample_type_color_t<-sample_type_color
  names(sample_type_color_t) <-sample_types  
  sample_type_t<-sample_types
  names(sample_type_t)<-sample_type_color
  sample_type_t<-sort(sample_type_t)
  sample_type_color_t2<-names(sample_type_t)
  data.OTU<-t(data.matrix(otu_table(phyloseq_t)))
  taxa_data.OTU<-as.data.frame(tax_table(phyloseq_t))
  taxa_data.OTU$new_name<-paste(taxa_data.OTU[[2]],taxa_data.OTU[[3]],taxa_data.OTU[[4]],taxa_data.OTU[[5]], taxa_data.OTU[[6]], taxa_data.OTU[[7]],rownames(taxa_data.OTU), sep=".")
  
  colnames(data.OTU)<-taxa_data.OTU$new_name
  sample_type_var<-sample_data(phyloseq_t)[[sample_type_var_name]] ###keeping the column with the variable of interest only
  X<-data.OTU
  Y<-sample_type_var
  col_name<-paste0(paste0(substr(colnames(X), 1, 20), ".."),str_sub(colnames(X),-7,-1)) #limiting the name to 20 characters
  MyResult.splsda <- mixOmics::splsda(X, Y, ncomp = 5) # 1 Run the method
  png_output1<- paste0(output_name,"_a.png")
  png(png_output1, res = 300, width=width, height=height , units='mm') #PLSDA with confidence ellipses
  p3<-plotIndiv(MyResult.splsda, comp =c(1,2), col.per.group = sample_type_color_t2,
                group = Y, ind.names = FALSE,  # colour points by class
                ellipse = TRUE, # include 95% confidence ellipse for each class
                legend = TRUE, title = '(a) PLSDA with confidence ellipses',
  )
  show(p3)
  dev.off()
  png_output3<- paste0(output_name,"_b.png")
  png(png_output3, res = 300, width=width*3/4, height=height , units='mm') #Plot of Loading vectors
  p4<-plotLoadings(MyResult.splsda, 
                   contrib = 'max', method = 'mean',
                   comp=1, 
                   ndisplay = top_display_loading, name.var=col_name) ##display top components

  if (!missing(decontam)) { #if decontam option is selected, determine if the df object exist
    if (exists(decontam)=="FALSE"){
      print(paste0("decontam list:",paste0(decontam," does not exist")))
      stop("please run function: decontaminant_subplot_KW prior to using this option- decontam") 
    } else{
      p4_original<-p4
      p4<-merge(p4, get(decontam),by="row.names", All.x=all)
      row.names(p4)<-p4$Row.names #previous merge step generate a new column call Row.names
      p4<- subset(p4,select=-c(Row.names))
      contaminant_color <- ifelse(p4$contaminant == "TRUE", "red", "black")
    }
  }  else{
    contaminant_color<-"black" ###if no decontam list was supplied, then will have all the label as black 
  }
  p <- ggplot(p4, aes(x = reorder(rownames(p4), -importance), y = importance))+
    geom_col(aes(fill=GroupContrib), width=0.6) + scale_fill_manual("sample type",values=sample_type_color_t) + 
    ###shortened taxa name###
    scale_x_discrete(
      labels = function(x) {
        x_org<-x
        is_long <- nchar(x) > 35
        x[is_long] <- paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",4),1,10),
                             paste0("..",paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",5),1,10),".."),
                                    paste0(substr(sapply(strsplit(x[is_long],".",fixed=T),"[[",6),1,10),
                                           paste0("..",str_sub(x[is_long],-4,-1)))))
        is_short <- nchar(x) < 17
        x[is_short] <- paste0(paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",3),1,15),".."),
                              paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",4),1,5),
                                     paste0("..",paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",5),1,5),".."),
                                            paste0(substr(sapply(strsplit(x_org[is_short],".",fixed=T),"[[",6),1,5),
                                                   paste0("..",str_sub(x_org[is_short],-4,-1))))))
        x
      }
    )+
    #########################
    coord_flip()+
    theme(panel.background = element_blank(),
          panel.border=element_rect(fill=NA),
          panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
          panel.grid.minor = element_blank(),strip.background=element_blank(),
          plot.title=element_text(face="bold",hjust=0.5, size = 14), 
          plot.subtitle = element_text(hjust=0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x=element_text(colour = "black", size = 14),
          axis.text.y=element_text(colour = contaminant_color, size = 22),
          axis.ticks=element_blank(),
          plot.margin=unit(c(1,1,1,1),"line"),
          legend.title = element_text(size=14), 
          legend.text = element_text(size=14))
  show(p)
  dev.off()
  
}

##Edger_phylo_KW:  
#Edger_phylo_KW(input_phyloseq=Transplant_PS_all, 
#               variable_to_compare="sample_types",
#               outcome_to_compare=c("Lower", "Upper"),
#               outcome_to_compare_color=c("steelblue2", "purple2"),
#               FDR_cut_off=0.2,
#               legend_onplot="yes",
#               output_name="Edger_bronch_Lower_Upper", 
#               decontam="Contam_list_NC_BKG_compare_Lower")
Edger_phylo_KW <- function (input_phyloseq, 
                            variable_to_compare,
                            outcome_to_compare=list(),
                            outcome_to_compare_color=list(),
                            FDR_cut_off=0.2,
                            number_display=10,
                            graph_option=c("volcano","lollipop"),
                            display_all_results_volcano=c("yes","no"), #display all the results in volcano plot (including results that is not significant: FDR not below the cut off)
                            abundance_size_volcano=c("yes","no"), #change size of the dot based on relative abundance
                            output_name,
                            legend_onplot=c("yes","no"),
                            decontam=NULL,
                            taxa_genus_output) {
  ###ensure intput_phyloseq is a absolute count phyloseq, not relative abundance phyloseq. and make sure phyloseq is in right orientation
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  countData_check <- floor(as(otu_table(input_phyloseq), "matrix")) # round down to nearest integer. if this phsyloeq is a relative abundance, then the entire countData_check would be 0
  if (all(countData_check==0)) {stop("Please use phyloseq with absolute count. The current input physloeq contains relative abundnace")}
  
  if (missing(taxa_genus_output)) { taxa_genus_output<-"no"} #default being yes for taxa_rank_name, which will just give the selected taxa name
  if (missing(graph_option)) { graph_option<-"lollipop"} #default being lollipop for graph_option
  if (missing(display_all_results_volcano)) { display_all_results_volcano<-"no"} #default being lollipop for graph_option
  if (missing(abundance_size_volcano)) { abundance_size_volcano<-"yes"} #default being lollipop for graph_option
  
  ## evaluate choices
  legend_onplot <- match.arg(legend_onplot)

  for (i in 1:length(outcome_to_compare)) {
    assign(paste0("outcome_to_compare",i), outcome_to_compare[[i]])
    assign(paste0("outcome_to_compare_color",i), outcome_to_compare_color[[i]])
  }  
  
  ##coldata is the meta data
  coldata<-as.data.frame(sample_data(input_phyloseq))
  #getting taxa name
  taxa.table <- as.data.frame(tax_table(input_phyloseq))
  countdata<-as.data.frame(otu_table(input_phyloseq))
  row.names(countdata)<-paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
  taxa.table$combine_name<-paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
  #Convert to numeric
  countdata= data.frame(lapply(countdata, function(x) as.numeric(as.character(x))),
                        check.names=F, row.names = rownames(countdata))
  #Set order
  coldata_SAMPLE <- data.frame(coldata[order(row.names(coldata)),],check.names=F)
  countdata_SAMPLE <- data.frame(countdata[, order(colnames(countdata))],check.names=F)
  
  #making a compare variable to standardize the input variable for comparsion- make sure the reference is always control
  coldata_SAMPLE$variable_to_compare<-coldata_SAMPLE[[variable_to_compare]]
  coldata_SAMPLE<- coldata_SAMPLE %>% dplyr::mutate(compare=case_when((variable_to_compare==outcome_to_compare1)~"control", (variable_to_compare==outcome_to_compare2)~"treatment"))
  coldata_SAMPLE<-coldata_SAMPLE[!is.na(coldata_SAMPLE$compare),] ###dropping subjects if they are not in the group of comparsion
  coldata_SAMPLE$compare<- as.factor(coldata_SAMPLE$compare) #converting variable to factor
  countdata_SAMPLE<-countdata_SAMPLE[,colnames(countdata_SAMPLE) %in% rownames(coldata_SAMPLE)]#making sure coldata has same subjects and countdata

  taxa_name<- subset(taxa.table, select=c(combine_name))
  
  #Set order
  coldata_SAMPLE <- data.frame(coldata_SAMPLE[order(row.names(coldata_SAMPLE)),])
  countdata_SAMPLE <- data.frame(countdata_SAMPLE[, order(colnames(countdata_SAMPLE))])
  
  #set up the phyloseq to make sure the control and treatment groups are correct. remove any observation that do not belong to the control and treatment group
  phyloseq_edited<-input_phyloseq
  sample_data(phyloseq_edited)$variable_to_compare <-sample_data(phyloseq_edited)[[variable_to_compare]]
  phyloseq_edited2 <- phyloseq_edited %>% ps_mutate(compare=case_when((variable_to_compare==outcome_to_compare1)~"control", (variable_to_compare==outcome_to_compare2)~"treatment"))
  keep_sample<- get_variable(phyloseq_edited2,"compare") %in% c("control", "treatment")
  phyloseq_edited3 <- prune_samples(keep_sample, phyloseq_edited2)  ###keeping only samples with the outcome variable of choice 

              #Run EDGER
              ######from phyloseq_to_edgeR
              #https://joey711.github.io/phyloseq-extensions/edgeR.html 
                    if( !taxa_are_rows(phyloseq_edited3) ){ phyloseq_edited3 <- t(phyloseq_edited3) }
                    x = as(otu_table(phyloseq_edited3), "matrix")
                    # Add one to protect against overflow, log(0) issues.
                    x = x + 1
                    # Check `compare` argument
                    if( identical(all.equal(length("compare"), 1), TRUE) & nsamples(phyloseq_edited3) > 1 ){
                      # Assume that compare was a sample variable name (must be categorical)
                      group = get_variable(phyloseq_edited3, "compare")
                    }
                    # Define gene annotations (`genes`) as tax_table
                    taxonomy = tax_table(phyloseq_edited3, errorIfNULL=FALSE)
                    if( !is.null(taxonomy) ){
                      taxonomy = data.frame(as(taxonomy, "matrix"))
                    } 
                    # Now turn into a DGEList
                    y = DGEList(counts=x, group=group, genes=taxonomy,remove.zeros = TRUE)
                    # Calculate the normalization factors
                    z = calcNormFactors(y, method="TMM")
                    # Check for division by zero inside `calcNormFactors`
                    if( !all(is.finite(z$samples$norm.factors)) ){
                      stop("Something wrong with edgeR::calcNormFactors on this data,
                           non-finite $norm.factors, consider changing `method` argument")
                    }
                    # Estimate dispersions
                    dgeFull<-estimateTagwiseDisp(estimateCommonDisp(z))
              dgeTest <- exactTest(dgeFull)
              dgeTest<<-dgeTest
              dgeFull<<-dgeFull
              
              
              #Create Table
              resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table))
              resNoFilt <- data.frame(resNoFilt)
              resNoFilt<-merge(resNoFilt,taxa_name, by="row.names")
              row.names(resNoFilt)<-resNoFilt$combine_name 
              resNoFilt<- subset(resNoFilt,select=c(logFC, logCPM, PValue, FDR))
              resNoFilt_export<<-resNoFilt
              #check to see if there is any taxa with FDR<FDR_cut_off. If not then will stop the function, because there will be error if it continues as the object will have no rows
                    checkifallzero<-resNoFilt[resNoFilt$FDR<FDR_cut_off,]
                    if (dim(checkifallzero)[1] == 0) {
                              print(paste0("no Taxa with FDR more than",FDR_cut_off))
                              return(NULL)
                    }
              resNoFilt[nrow(resNoFilt) + 1,] <- c(0, 0, 0, 0) ###adding this row of 0 so that even if there is no significant top 10 for upregulated or downregulated, the code will still run
              #Create Table of Down Regulated
              sigDownReg <- resNoFilt[resNoFilt$FDR<FDR_cut_off,] #keeping only FDR < 0.2
              #sigDownReg <- resNoFilt
              sigDownReg <- sigDownReg[sigDownReg$logFC<=0,] ##keeping only the negative logFC
              sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
              sigDownReg1 <- sigDownReg %>% dplyr::slice(1:number_display) #Subset Top (number_display)
              #Create Table of Down Regulated
              sigUpReg <- resNoFilt[resNoFilt$FDR<FDR_cut_off,] #keeping only FDR < 0.2
              #sigUpReg <- resNoFilt
              sigUpReg <- sigUpReg[sigUpReg$logFC>=0,] ##keeping only the positive logFC
              sigUpReg <- sigUpReg[order(sigUpReg$logFC, decreasing=TRUE),]
              sigUpReg1 <- sigUpReg %>% dplyr::slice(1:number_display) #Subset Top (number_display)
              #Merge Top 10s
              res <- rbind(sigDownReg1,sigUpReg1)
              res <- res[res$logFC!=0,] 
              #remove the 0 column which was used as a place holder
              resNoFilt <- resNoFilt[resNoFilt$logFC!=0,] 
              
              #Create Relative Abundance Table
              df <-
                countdata_SAMPLE %>% 
                rownames_to_column('gs') %>%
                dplyr::group_by(gs) %>% 
                dplyr::summarise_all(funs(sum)) %>%
                dplyr::mutate_if(is.numeric, funs(./sum(.))) %>%
                column_to_rownames('gs')
              df_export<<-df
              #Get the ColData for Each Comparison
              coldata.1 <- coldata_SAMPLE[coldata_SAMPLE[[variable_to_compare]]==outcome_to_compare1,]
              coldata.2 <- coldata_SAMPLE[coldata_SAMPLE[[variable_to_compare]]==outcome_to_compare2,]
              #keep Count data only for each comparison
              needed<-which(colnames(df) %in% rownames(coldata.1))    
              df.1 <- df[,needed]
              needed2<-which(colnames(df) %in% rownames(coldata.2))    
              df.2 <- df[,needed2]
              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
              if (display_all_results_volcano=="yes" & graph_option=="volcano"){ #if want to display all results, then will replace res by resNoFilt (which does not select the top results)
                res<-resNoFilt
              }
              
              #Convert Resuts table into a data.frame
              res <- as.data.frame(res)
              #decide what otu to save 
              otu.to.save <-as.character(rownames(res))
              #from relative table we should get the mean across the row of the otu table
              df.1.meanRA <- rowMeans(df.1)
              df.2.meanRA <- rowMeans(df.2)
              #need to subset AND reorder just the otus that we have 
              df.1.meanRA.save <- df.1.meanRA[otu.to.save]
              df.2.meanRA.save <- df.2.meanRA[otu.to.save]
              #add the abundnace data for the res dataframe
              res$abundance.1 <- df.1.meanRA.save
              res$abundance.2 <- df.2.meanRA.save
              #Set Names of Results Table
              res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Taxa","logFC", "lfcSE", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
              #Remove Any Data without LOGFC data
              res <- res[!is.na(res$logFC),]
              # Reorder Results based on FDR for comparison 1
              res = res[order(-res$logFC, na.last = TRUE), ]
              #Create Order
              res <- res %>% dplyr::mutate(start = 1:n())
              #Convert Important columns to Numeric
              res$adj.P.Val <-   as.numeric(as.character(res$adj.P.Val))
              res$logFC <-       as.numeric(as.character(res$logFC))
              #Replace NA
              res <- res %>% dplyr::mutate(adj.P.Val = if_else(is.na(adj.P.Val), 0.9, adj.P.Val))
              ##convert abundance to numeric
              res$abundance.1 <- as.numeric(as.character(res$abundance.1))
              res$abundance.2 <- as.numeric(as.character(res$abundance.2))
              #Create Variable for Color based on Comparison, FDR and LOGFC
              res$col <- ifelse(res$logFC>0, "B",
                                ifelse(res$logFC<0, "A","D"))
              res$abundance <- ifelse(res$col=="A", res$abundance.1, ifelse(res$col=="B", res$abundance.2,0))
              
              if (taxa_genus_output=="no") {
                res$Taxa2 <- paste0(substr(sapply(strsplit(res$Taxa,".",fixed=T),"[[",1),1,6),
                                    paste0("..",paste0(substr(sapply(strsplit(res$Taxa,".",fixed=T),"[[",4),1,10),".."),
                                           paste0(substr(sapply(strsplit(res$Taxa,".",fixed=T),"[[",5),1,10),paste0("..",str_sub(res$Taxa,-4,-1)))))  
              } else if (taxa_genus_output=="yes") {#keeping only genus level. if genus level is NA, then will go up a taxa rank
                res<- res %>% dplyr::mutate(Taxa2=case_when((sapply(strsplit(Taxa,".",fixed=T),"[[",1)=="NA")~"k_Bacteria.p_NA.c_NA.o_NA.f_NA.g_NA",
                                                            ((sapply(strsplit(Taxa,".",fixed=T),"[[",1)!="NA") & (sapply(strsplit(Taxa,".",fixed=T),"[[",2)=="NA"))~paste0(paste0("p_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",1)),".c_NA.o_NA.f_NA.g_NA"), 
                                                            ((sapply(strsplit(Taxa,".",fixed=T),"[[",2)!="NA") & (sapply(strsplit(Taxa,".",fixed=T),"[[",3)=="NA"))~paste0(paste0("c_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",2)),".o_NA.f_NA.g_NA"), 
                                                            ((sapply(strsplit(Taxa,".",fixed=T),"[[",3)!="NA") & (sapply(strsplit(Taxa,".",fixed=T),"[[",4)=="NA"))~paste0(paste0("o_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",3)),".f_NA.g_NA"), 
                                                            ((sapply(strsplit(Taxa,".",fixed=T),"[[",4)!="NA") & (sapply(strsplit(Taxa,".",fixed=T),"[[",5)=="NA"))~paste0(paste0("f_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",4)),".g_NA"),
                                                            ((sapply(strsplit(Taxa,".",fixed=T),"[[",5)!="NA"))~paste0("g_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",5)) ))
              }
              #########EXCEL EXPORT
              ####keeping all results 
              res_all <- rbind(sigDownReg,sigUpReg)
              res_all <- res_all[res_all$logFC!=0,] 
              
              #Convert Resuts table into a data.frame
              res_all <- as.data.frame(res_all)
              #decide what otu to save 
              otu.to.save_all <-as.character(rownames(res_all))
              #need to subset AND reorder just the otus that we have 
              df.1.meanRA.save <- df.1.meanRA[otu.to.save_all]
              df.2.meanRA.save <- df.2.meanRA[otu.to.save_all]
              #add the abundnace data for the res dataframe
              res_all$abundance.1 <- df.1.meanRA.save
              res_all$abundance.2 <- df.2.meanRA.save
              #Set Names of Results Table
              res_all <- setNames(cbind(rownames(res_all), res_all, row.names = NULL), c("Taxa","logFC", "lfcSE", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
              csv_output<-paste0(output_name,"_significant.csv")
              write.csv(res_all,file=csv_output, row.names = TRUE)
              
              #res_all contains all results including not significant results
              #Convert Resuts table into a data.frame
              res_all <- as.data.frame(resNoFilt)
              #decide what otu to save 
              otu.to.save_all <-as.character(rownames(res_all))
              #need to subset AND reorder just the otus that we have 
              df.1.meanRA.save <- df.1.meanRA[otu.to.save_all]
              df.2.meanRA.save <- df.2.meanRA[otu.to.save_all]
              #add the abundnace data for the res dataframe
              res_all$abundance.1 <- df.1.meanRA.save
              res_all$abundance.2 <- df.2.meanRA.save
              csv_output2<-paste0(output_name,"_all.csv")
              write.csv(resNoFilt_export,file=csv_output2, row.names = TRUE)
              
              ######################################################################################################################
              ######################################################################################################################
              ####merge in decontam list to change the contaminant taxa red
              if (!missing(decontam)) { #if decontam option is selected, determine if the df object exist
                if (exists(decontam)=="FALSE"){
                  print(paste0("decontam list:",paste0(decontam," does not exist")))
                  stop("please run function: decontaminant_subplot_KW prior to using this option- decontam") 
                } else{
                  res_original<-res
                  res<-merge(res, get(decontam),by.x="Taxa",by.y="row.names",All.x=all)
                  res<- res[order(res$logFC, na.last = TRUE), ]
                  contaminant_color <- ifelse(res$contaminant == "TRUE", "red", "black")
                }
              }  else{
                contaminant_color<-"black" ###if no decontam list was supplied, then will have all the label as black 
              }
              
              
              ######################################################################################################################
              ######################################################################################################################
              
              if (graph_option=="volcano") {
                          #=========================================================
                          #////////////////////VOLCANO PLOT///////////////////////
                          #=========================================================
                          # Compute significance, with a maximum of 350 for the p-values set to 0 due to limitation of computation precision
                          res$sig <- -log10(res$adj.P.Val)
                          res[is.infinite(res$sig),"sig"] <- 350
                          ## Volcano plot of adjusted p-values
                          cols <- densCols(res$logFC, res$sig)
                          cols[res$pvalue ==0] <- "purple"
                          #only significant values get display
                          cols[res$logFC > 0 & res$adj.P.Val < FDR_cut_off ] <- outcome_to_compare_color2
                          cols[res$logFC < 0 & res$adj.P.Val < FDR_cut_off ] <- outcome_to_compare_color1
                          cols[res$adj.P.Val > FDR_cut_off ] <- "grey"
                          
                          ####scaling the abundance so the circle size on the graph is standardized
                          res_keep<-res[!(is.na(res$adj.P.Val)),]
                          res_keep<-res_keep[res_keep$adj.P.Val<FDR_cut_off,]
                          min_value<-min(res_keep$abundance)
                          #make the abundance of the nonsignificant results to be minimal value
                          res$abundance[res$adj.P.Val > FDR_cut_off]<-min_value
                          
                          if (abundance_size_volcano=="yes"){ #if want point size to be scaled by relative abundance
                            p<-ggplot(res, aes(x = logFC, y = sig, label=Taxa2)) +
                              geom_point(aes(size =abundance),color=cols, alpha=0.5) + #Chose Colors and size for dots
                              scale_size(range = c(0, 40))+
                              geom_text_repel(aes(label=ifelse(res$logFC<(-1) & res$adj.P.Val < FDR_cut_off , as.character(res$Taxa2),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5,  colour = contaminant_color) +
                              geom_text_repel(aes(label=ifelse(res$logFC>1 & res$adj.P.Val < FDR_cut_off , as.character(res$Taxa2),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5,  colour = contaminant_color) +
                              geom_hline(yintercept=-log10(FDR_cut_off), color="red",linetype="dashed") +
                              xlab("Effect size: log2(fold-change)") +
                              ylab("-log10(adjusted p-value)") + 
                              theme(panel.background = element_blank(),
                                    panel.border=element_rect(fill=NA),
                                    legend.background = element_rect(color=NA),
                                    legend.key = element_rect(colour = "transparent", fill = "white"),
                                    plot.title=element_text(size=23, face="bold"))
                          } else {
                            p<-ggplot(res, aes(x = logFC, y = sig, label=Taxa2)) +
                              geom_point(color=cols, alpha=0.5) + #Chose Colors and size for dots
                              geom_text_repel(aes(label=ifelse(res$logFC<(-1) & res$adj.P.Val < FDR_cut_off , as.character(res$Taxa2),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5,  colour = contaminant_color) +
                              geom_text_repel(aes(label=ifelse(res$logFC>1 & res$adj.P.Val < FDR_cut_off , as.character(res$Taxa2),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5,  colour = contaminant_color) +
                              geom_hline(yintercept=-log10(FDR_cut_off), color="red",linetype="dashed") +
                              xlab("Effect size: log2(fold-change)") +
                              ylab("-log10(adjusted p-value)") + 
                              theme(panel.background = element_blank(),
                                    panel.border=element_rect(fill=NA),
                                    plot.title=element_text(size=23, face="bold"))
                          }
              }
              else if (graph_option=="lollipop") {
                        p<-ggplot(res, aes(y=reorder(Taxa2,-start), x=logFC,fill=col,size=abundance)) +
                          geom_point(color="black",alpha=0.8,shape=21)+
                          geom_segment(data=res[res$adj.P.Val<0.2,],aes(yend=reorder(Taxa2,-start)), xend=(-30), color= "black", linetype = "solid",linewidth=1)+ 
                          scale_fill_manual(values=c("B"=outcome_to_compare_color2,"A"=outcome_to_compare_color1,"D"="white"))+ 
                          scale_size_continuous(name="Relative Abundance",range=c(5, 20))+                
                          ggtitle("EdgeR")+
                          theme(panel.background = element_blank(),
                                panel.border=element_rect(fill=NA),
                                panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
                                panel.grid.minor = element_blank(),
                                strip.background=element_blank(),
                                axis.title=element_text(size=20,face="bold"),
                                axis.text.x=element_text(colour="black", size=18, face="bold"),
                                axis.text.y=element_text(colour=contaminant_color,face="bold",size=10),
                                axis.ticks=element_line(colour="gray70"),
                                legend.background = element_rect(color=NA),
                                legend.key = element_rect(colour = "transparent", fill = "white"),
                                plot.title=element_text(size=23, face="bold"))+
                          xlab("") +
                          ylab("")+
                          #xlim(-7,7)+
                          geom_vline(xintercept=0, color="red",linetype="dashed")+
                          guides(fill="none")
                        ##save the legend on a separate png file              
              }          
              if (legend_onplot=="no") { #if want to hide the legend on the plot, then will save the legend on separate png
                  leg <- get_legend(p)
                  png_output2<-paste0(output_name,"_legend.png")
                  png(png_output2, res = 300, width=200, height=200, units='mm')
                  legend<-as_ggplot(leg)
                  show(legend)
                  dev.off()
                  p<-p + theme(legend.position = "none") ###remove the legend
              }
              png_output<-paste0(output_name,".png")
              png(file=png_output, res = 300, width=200, height=200 , units='mm')
                    show(p)
              dev.off()
              ##########################################################################################################################
              ##########################################################################################################################
}

##Deseq_phylo_KW: 
#Deseq_phylo_KW(Transplant_PS_all, 
#               variable_to_compare="sample_type",
#               outcome_to_compare=c("Lower", "Upper"),
#               outcome_to_compare_color=c("steelblue2", "purple2"),
#               alpha_level=0.2,
#               graph_option="volcano",
#               legend_onplot="yes",
#               output_name="Deseq_bronch_Lower_Upper2", 
#               decontam="Contam_list_NC_BKG_compare_Lower")
Deseq_phylo_KW <- function (input_phyloseq, 
                            variable_to_compare,
                            outcome_to_compare=list(),
                            outcome_to_compare_color=list(),
                            alpha_level=0.2,
                            number_display=10,
                            graph_option=c("volcano","lollipop"),
                            legend_onplot=c("yes","no"),
                            output_name,
                            decontam=NULL) {
        ###ensure intput_phyloseq is a absolute count phyloseq, not relative abundance phyloseq. and make sure phyloseq is in right orientation
        # Enforce orientation. Samples are columns
        if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
        countData_check <- floor(as(otu_table(input_phyloseq), "matrix")) # round down to nearest integer. if this phsyloeq is a relative abundance, then the entire countData_check would be 0
        if (all(countData_check==0)) {stop("Please use phyloseq with absolute count. The current input physloeq contains relative abundnace")}
        
        alpha<-alpha_level
        ## evaluate choices
        graph_option <- match.arg(graph_option)
        print(graph_option)
        legend_onplot <- match.arg(legend_onplot)
        
        for (i in 1:length(outcome_to_compare)) {
          assign(paste0("outcome_to_compare",i), outcome_to_compare[[i]])
          assign(paste0("outcome_to_compare_color",i), outcome_to_compare_color[[i]])
        }  

        ##coldata is the meta data
        coldata<-as.data.frame(sample_data(input_phyloseq))
        #getting taxa name
        taxa.table <- as.data.frame(tax_table(input_phyloseq))
        countdata<-as.data.frame(otu_table(input_phyloseq))
        row.names(countdata)<-paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
        #Convert to numeric
        countdata= data.frame(lapply(countdata, function(x) as.numeric(as.character(x))),
                              check.names=F, row.names = rownames(countdata))
        #Set order
        coldata_SAMPLE <- data.frame(coldata[order(row.names(coldata)),],check.names=F)
        countdata_SAMPLE <- data.frame(countdata[, order(colnames(countdata))],check.names=F)
        
        #making a compare variable to standardize the input variable for comparsion- make sure the reference is always control
        coldata_SAMPLE$variable_to_compare <-coldata_SAMPLE[[variable_to_compare]]
        coldata_SAMPLE <- coldata_SAMPLE %>% dplyr::mutate(compare=case_when((variable_to_compare==outcome_to_compare1)~"control", (variable_to_compare==outcome_to_compare2)~"treatment"))
        coldata_SAMPLE <- coldata_SAMPLE[!is.na(coldata_SAMPLE$compare),] ###dropping subjects if they are not in the group of comparsion
        coldata_SAMPLE$compare<- as.factor(coldata_SAMPLE$compare) #convering variable to factor
        countdata_SAMPLE<-countdata_SAMPLE[,colnames(countdata_SAMPLE) %in% rownames(coldata_SAMPLE)] #this makes sure countdata_SAMPLE and coldata_SAMPLE both have same subjects 
        
        #Set order
        coldata_SAMPLE <- data.frame(coldata_SAMPLE[order(row.names(coldata_SAMPLE)),])
        countdata_SAMPLE <- data.frame(countdata_SAMPLE[, order(colnames(countdata_SAMPLE))])
        
        #set up the phyloseq to make sure the control and treatment groups are correct. remove any observation that do not belong to the control and treatment group
        phyloseq_edited<-input_phyloseq
        sample_data(phyloseq_edited)$variable_to_compare <-sample_data(phyloseq_edited)[[variable_to_compare]]
        phyloseq_edited2 <- phyloseq_edited %>% ps_mutate(compare=case_when((variable_to_compare==outcome_to_compare1)~"control", (variable_to_compare==outcome_to_compare2)~"treatment"))
        keep_sample<- get_variable(phyloseq_edited2,"compare") %in% c("control", "treatment")
        phyloseq_edited3 <- prune_samples(keep_sample, phyloseq_edited2)  ###keeping only samples with the outcome variable of choice 
        
        ###get taxa name list so it can be merge back to the result from deseq
        taxa.table_V2 <- as.data.frame(tax_table(phyloseq_edited3))
        taxa.table_V2$combine_name<-paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
        taxa_name<- subset(taxa.table_V2, select=c(combine_name))
        # Calculate geometric means prior to estimate size factors
        gm_mean = function(x, na.rm = TRUE){
          exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
        }
        
        #Create new DESEQ Objects        
        diagdds<- phyloseq_to_deseq2(phyloseq_edited3, ~ compare)
        geoMeans = apply(counts(diagdds), 1, gm_mean)
        diagdds <- estimateSizeFactors(diagdds, geoMeans = geoMeans)
        #making sure reference is control. resetting the level to make sure control is reference
        diagdds$compare <- droplevels(diagdds$compare)
        diagdds$compare <- relevel(diagdds$compare, ref ="control")
        
        diagdds = estimateDispersions(diagdds, fitType = "parametric")
        diagdds <- nbinomWaldTest(diagdds)
        #Output Result Table
        res <- results(diagdds, cooksCutoff = FALSE)
        res_org<<-res
        
        #Convert Resuts table into a data.frame
        res <- as.data.frame(res)
        res<-merge(res,taxa_name, by="row.names") #merging back the taxa name from the phyloseq
        row.names(res)<-res$combine_name 
        res<- subset(res,select=c(baseMean,log2FoldChange,lfcSE,stat,pvalue,padj))
        resNoFilt<- res #resNoFilt is backup object

        #Create Relative Abundance Table from the 
        df <-
          countdata_SAMPLE %>% 
          rownames_to_column('gs') %>%
          dplyr::group_by(gs) %>% 
          dplyr::summarise_all(funs(sum)) %>%
          dplyr::mutate_if(is.numeric, funs(./sum(.))) %>%
          column_to_rownames('gs')
        #Get the ColData for Each Comparison
        coldata.1 <- coldata_SAMPLE[coldata_SAMPLE[[variable_to_compare]]==outcome_to_compare1,]
        coldata.2 <- coldata_SAMPLE[coldata_SAMPLE[[variable_to_compare]]==outcome_to_compare2,]
        #keep Count data only for each comparison
        needed<-which(colnames(df) %in% rownames(coldata.1))    
        df.1 <- df[,needed]
        needed2<-which(colnames(df) %in% rownames(coldata.2))    
        df.2 <- df[,needed2]
        #decide what otu to save 
        otu.to.save <-as.character(rownames(res))
        #from relative table we should get the mean across the row of the otu table
        df.1.meanRA <- rowMeans(df.1)
        df.2.meanRA <- rowMeans(df.2)
        #need to subset AND reorder just the otus that we have 
        df.1.meanRA.save <- df.1.meanRA[otu.to.save]
        df.2.meanRA.save <- df.2.meanRA[otu.to.save]
        #add the abundnace data for the res dataframe
        res$abundance.1 <- df.1.meanRA.save
        res$abundance.2 <- df.2.meanRA.save
        #Set Names of Results Table
        res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Taxa","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
        #adding the shortened taxa name
        res$Taxa2 <- paste0(substr(sapply(strsplit(res$Taxa,".",fixed=T),"[[",1),1,6),
                            paste0("..",paste0(substr(sapply(strsplit(res$Taxa,".",fixed=T),"[[",4),1,10),".."),
                                   paste0(substr(sapply(strsplit(res$Taxa,".",fixed=T),"[[",5),1,10),paste0("..",str_sub(res$Taxa,-4,-1)))))
        #Write Tables of Differential Analysis
        csv_output<-paste0(output_name,".csv")
        write.csv(res,file=csv_output, row.names = TRUE)
        
        res_regulate<-as.data.frame(res) 
        
        res_regulate_testing<<-res_regulate
        res_regulate<-res_regulate %>% filter(!is.na(adj.P.Val))
        #check to see if there is any taxa with adjusted p value<alpha. If not then will stop the function, because there will be error if it continues as the object will have no rows
        checkifallzero<-res_regulate[res_regulate$adj.P.Val<alpha,]
        if (dim(checkifallzero)[1] == 0) {
          print(paste0("no Taxa with adjusted p value more than"),alpha)
          return(NULL)
        }
        res_regulate_positiveFC<-res_regulate %>% filter(res_regulate$logFC > 0 & res_regulate$adj.P.Val < alpha)
        if (dim(res_regulate_positiveFC)[1] == 0) {### if there is no upregulated with adjusted p value less than alpha, then 
          res_regulate_positiveFC[nrow(res_regulate_positiveFC) + 1,] <- replicate(ncol(res_regulate_positiveFC),0) ###adding this row of 0 so that even if there is no significant top 10 for upregulated or downregulated, the code will still run
        }
        res_regulate_positiveFC <- res_regulate_positiveFC[order(-res_regulate_positiveFC$logFC),] #order has the lowest value on the top, so need the - to reverse the order since slice cut from top down 
        res_regulate_positiveFC1 <- res_regulate_positiveFC %>% dplyr::slice(1:number_display)
        
        res_regulate_negativeFC<-res_regulate %>% filter(res_regulate$logFC < 0 & res_regulate$adj.P.Val < alpha)
        if (dim(res_regulate_negativeFC)[1] == 0) {
          res_regulate_negativeFC[nrow(res_regulate_negativeFC) + 1,] <- replicate(ncol(res_regulate_negativeFC),0) ###adding this row of 0 so that even if there is no significant top 10 for upregulated or downregulated, the code will still run
        }
        res_regulate_negativeFC <- res_regulate_negativeFC[order(res_regulate_negativeFC$logFC),]
        res_regulate_negativeFC1 <- res_regulate_negativeFC %>% dplyr::slice(1:number_display) #Subset Top 10
        #Merge Top 10s (number_display)
        res_regulate_final <- rbind(res_regulate_negativeFC1,res_regulate_positiveFC1)
        res_regulate_final <- res_regulate_final[res_regulate_final$logFC!=0,] 
        res_regulate_final$col <- ifelse(res_regulate_final$logFC>0, "B",
                                  ifelse(res_regulate_final$logFC<0, "A","D"))
        #Convert to dataframe
        res_regulate_final$abundance <- ifelse(res_regulate_final$col=="A", res_regulate_final$abundance.1, ifelse(res_regulate_final$col=="B", res$abundance.2,0))
        
        # Reorder Results based on FDR for comparison 1
        res_regulate_final = res_regulate_final[order(-res_regulate_final$logFC, na.last = TRUE), ]
        #Create Order
        res_regulate_final <- res_regulate_final %>% dplyr::mutate(start = 1:n())
        res_regulate_final$abundance<-ifelse(res_regulate_final$adj.P.Val<alpha & res_regulate_final$logFC>0, res_regulate_final$abundance.2, 
                                             ifelse(res_regulate_final$adj.P.Val<alpha & res_regulate_final$logFC<0, res_regulate_final$abundance.1,2))
        
        ######################################################################################################################
        ######################################################################################################################
        ####merge in decontam list to change the contaminant taxa red
        if (!missing(decontam)) { #if decontam option is selected, determine if the df object exist
          if (exists(decontam)=="FALSE"){
            print(paste0("decontam list:",paste0(decontam," does not exist")))
            stop("please run function: decontaminant_subplot_KW prior to using this option- decontam") 
          } else{
            res_original<-res
            res<-merge(res, get(decontam),by.x="Taxa",by.y="row.names",All.x=all)
            contaminant_color_volcano <<- ifelse(res$contaminant == "TRUE", "red", "black")
            
            res_regulate_final<-merge(res_regulate_final, get(decontam),by.x="Taxa",by.y="row.names",All.x=all)
            res_regulate_final<- res_regulate_final[order(res_regulate_final$logFC, na.last = TRUE), ]
            contaminant_color_ball_line <<- ifelse(res_regulate_final$contaminant == "TRUE", "red", "black")
          }
        }  else{
          contaminant_color_volcano<<-"black" ###if no decontam list was supplied, then will have all the label as black 
          contaminant_color_ball_line<<-"black"
        }

        ######################################################################################################################
        ######################################################################################################################
        
        if (graph_option=="volcano") {
              contaminant_color<-contaminant_color_volcano
              #=========================================================
              #////////////////////VOLCANO PLOT///////////////////////
              #=========================================================
              # Compute significance, with a maximum of 320 for the p-values set to 0 due to limitation of computation precision
              res$sig <- -log10(res$adj.P.Val)
              sum(is.infinite(res$sig))
              res[is.infinite(res$sig),"sig"] <- 350
              ## Volcano plot of adjusted p-values
              cols <- densCols(res$logFC, res$sig)
              cols[res$pvalue ==0] <- "purple"
              cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- outcome_to_compare_color2
              cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- outcome_to_compare_color1
              res$pch <- 19
              res$pch[res$pvalue ==0] <- 6
              
              ####scaling the abundance so the circle size on the graph is standardized
              res$abundance_scaled<-ifelse(res$adj.P.Val<alpha & res$logFC>0, res$abundance.2, ifelse(res$adj.P.Val<alpha & res$logFC<0, res$abundance.1,2))
              res_keep<-res[!(is.na(res$adj.P.Val)),]
              res_keep<-res_keep[res_keep$adj.P.Val<alpha,]
              max_value<-max(res_keep$abundance_scaled)
              value_max<-40/max_value
              res <- res %>% mutate(abundance_scaled2=case_when((adj.P.Val<alpha & !(is.na(adj.P.Val))) ~ abundance_scaled*value_max, T ~abundance_scaled ))
              
              #PLOT IT
              g<-ggplot(res, aes(x = logFC, y = sig, label=Taxa2)) +
                geom_point(size = res$abundance_scaled2,color=cols, alpha=0.5) + #Chose Colors and size for dots
                geom_text_repel(aes(label=ifelse(res$logFC<(-1) & res$adj.P.Val < alpha , as.character(res$Taxa2),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5,  colour = contaminant_color) +
                geom_text_repel(aes(label=ifelse(res$logFC>1 & res$adj.P.Val < alpha , as.character(res$Taxa2),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5, colour = contaminant_color) +
                theme(legend.position = "none") +
                geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
                xlab("Effect size: log2(fold-change)") +
                ylab("-log10(adjusted p-value)") + 
                #ylim(0,20)+
                ggtitle("Deseq")+
                theme(plot.title=element_text(size=23, face="bold"))
              png_output<-paste0(output_name,".png")
              png(file=png_output, res = 300, width=200, height=200 , units='mm')
              show(g)
              dev.off() 
        }
        else if (graph_option=="lollipop") {
            contaminant_color<-contaminant_color_ball_line
            p<-ggplot(res_regulate_final, aes(y=reorder(Taxa2,-start), x=logFC,fill=col,size=abundance)) +
              geom_point(color="black",alpha=0.8,shape=21)+
              geom_segment(data=res_regulate_final[res_regulate_final$adj.P.Val<alpha,],aes(yend=reorder(Taxa2,-start)), xend=(-30), color= "black", linetype = "solid",size=1)+ 
              scale_fill_manual(values=c("B"=outcome_to_compare_color2,"A"=outcome_to_compare_color1,"D"="white"))+ 
              scale_size_continuous(name="Relative Abundance",range=c(5, 20))+
              ggtitle("Deseq")+
              theme(panel.background = element_blank(),
                    panel.border=element_rect(fill=NA),
                    panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
                    panel.grid.minor = element_blank(),
                    strip.background=element_blank(),
                    axis.title=element_text(size=20,face="bold"),
                    axis.text.x=element_text(colour="black", size=18, face="bold"),
                    axis.text.y=element_text(colour=contaminant_color,face="bold",size=10),
                    axis.ticks=element_line(colour="black"),
                    legend.background = element_rect(color=NA),
                    legend.key = element_rect(colour = "transparent", fill = "white"),
                    plot.title=element_text(size=23, face="bold"))+
              xlab("") +
              ylab("")+
              #xlim(-7,7)+
              geom_vline(xintercept=0, color="red",linetype="dashed")+
              guides(fill="none")
            ##save the legend on a separate png file              
            if (legend_onplot=="no") { #if want to hide the legend on the plot, then will save the legend on separate png
              leg <- get_legend(p)
              png_output2<-paste0(output_name,"_legend.png")
              png(png_output2, res = 300, width=200, height=200, units='mm')
              legend<-as_ggplot(leg)
              show(legend)
              dev.off()
            }
            png_output<-paste0(output_name,".png")
            png(file=png_output, res = 300, width=200, height=200 , units='mm')
            if (legend_onplot=="no") { #if want to hide the legend
              p<-p + theme(legend.position = "none") ###remove the legend
            }
            show(p)
            dev.off()
        }
}

##Lefse_phylo_KW: 
#Lefse_phylo_KW(input_phyloseq=Transplant_PS_all, 
#               variable_to_compare="sample_type",
#               outcome_to_compare=c("Lower", "Upper"),
#               outcome_to_compare_color=c("steelblue2", "purple2"),
#               kw_cutoff_set=0.05,
#               wilcoxon_cutoff_set=0.05,
#               lda_cutoff_set=2,
#               graph_option="boxplot",
#               output_name="Lefse_bronch_Lower_Upper2", 
#               decontam="Contam_list_NC_BKG_compare_Lower")
Lefse_phylo_KW <- function(input_phyloseq, 
                           variable_to_compare, 
                           outcome_to_compare=list(), 
                           outcome_to_compare_color=list(), 
                           norm_method="CPM",
                           taxa_rank_select="none",
                           kw_cutoff_set=0.05,
                           wilcoxon_cutoff_set=0.05,
                           lda_cutoff_set=2,
                           graph_option=c("boxplot","bargraph", "none"), 
                           taxa_rank_name=c("yes","no"),
                           output_name, 
                           width=400, height=300, 
                           decontam=NULL, 
                           log_scale=c("yes","no"),
                           font_size_adjust_percent=100) {
  ###ensure intput_phyloseq is a absolute count phyloseq, not relative abundance phyloseq. and make sure phyloseq is in right orientation
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  countData_check <- floor(as(otu_table(input_phyloseq), "matrix")) # round down to nearest integer. if this phsyloeq is a relative abundance, then the entire countData_check would be 0
  if (all(countData_check==0)) {stop("Please use phyloseq with absolute count. The current input physloeq contains relative abundnace")}
  
  if (missing(taxa_rank_name)) { taxa_rank_name<-"yes"} #default being yes for taxa_rank_name, which will just give the selected taxa name
  if (missing(log_scale)) {log_scale<-"yes"} #default being yes for log scale with the relative abundance graph
  ####decontam list only look at ASV level... so would not be able to identify the contaminant when taxa level is selected
  if (!missing(decontam) & taxa_rank_select!="none") {
        stop(paste0("Decontam option is only available when taxa_rank_select is none"))
  } 
  
  ## evaluate choices
  graph_option <- match.arg(graph_option)
  taxa_rank_name <- match.arg(taxa_rank_name)
  log_scale<-match.arg(log_scale)
  print(graph_option)
  
  #taxa_rank option
  rank_type_list<-phyloseq::rank_names(input_phyloseq)
  for (d in 1:length(rank_type_list)) {
    assign(paste0("rank_type_list",d), rank_type_list[[d]])
  }
  
  rank_type_list2<-append(rank_type_list,c("all", "none"))
  #checking taxa_rank options
  if ((taxa_rank_select %in% rank_type_list2)==FALSE) {
    warning("Availble taxa_rank options:")
    warning( paste(rank_type_list2, collapse=", "))
    stop(paste0(taxa_rank_select," is not one of the taxa rank name of the phyloseq. Other options include 'all' and 'none' "))
  }

  for (i in 1:length(outcome_to_compare)) {
    assign(paste0("outcome_to_compare",i), outcome_to_compare[[i]])
    assign(paste0("outcome_to_compare_color",i), outcome_to_compare_color[[i]])
  }  

  #setting up the color
  outcome_to_compare_color_t<-outcome_to_compare_color
  names(outcome_to_compare_color_t) <-outcome_to_compare 
  
  print(variable_to_compare)
  print(outcome_to_compare)
  #need to get a vector of TRUE/FALSE to determine which subjects to keep in the phyloseq 
  keep_sample<- get_variable(input_phyloseq,variable_to_compare) %in% outcome_to_compare
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  phyloseq_update <- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
  
  #run lefse
  res <- run_lefse(phyloseq_update,
                   group=variable_to_compare,
                   subgroup = NULL,
                   taxa_rank = taxa_rank_select,
                   transform = "identity",
                   norm = norm_method,
                   kw_cutoff = kw_cutoff_set,
                   wilcoxon_cutoff = wilcoxon_cutoff_set,
                   lda_cutoff = lda_cutoff_set,
                   bootstrap_n = 30,
                   bootstrap_fraction = 2/3)
  marker_table(res)
  res<<-res
  #get the taxa name 
  if (taxa_rank_select %in% rank_type_list){ #if taxa_rank is selected then will subset the phyloseq to the level of taxa indicated
        phyloseq_update_sub<-tax_glom(phyloseq_update, taxrank=taxa_rank_select )
  }else{ 
        phyloseq_update_sub<-phyloseq_update
  }
  taxa.table_lefse <- as.data.frame(tax_table(phyloseq_update_sub))
  if (taxa_rank_select=="none"){
    taxa.table_lefse$feature<-row.names(taxa.table_lefse)  
    taxa.table_lefse$rowname_taxa<-row.names(taxa.table_lefse)
    taxa.table_lefse$Taxa<-paste(taxa.table_lefse[[2]],taxa.table_lefse[[3]],taxa.table_lefse[[4]],taxa.table_lefse[[5]], taxa.table_lefse[[6]], taxa.table_lefse[[7]],rownames(taxa.table_lefse), sep=".")
    #shorten taxa name
    taxa.table_lefse$Taxa2 <- paste0(substr(sapply(strsplit(taxa.table_lefse$Taxa,".",fixed=T),"[[",1),1,6),
                                     paste0("..",paste0(substr(sapply(strsplit(taxa.table_lefse$Taxa,".",fixed=T),"[[",4),1,10),".."),
                                            paste0(substr(sapply(strsplit(taxa.table_lefse$Taxa,".",fixed=T),"[[",5),1,10),paste0("..",str_sub(taxa.table_lefse$Taxa,-4,-1)))))
  } else if (taxa_rank_select=="all"){
    stop(paste0("taxa_rank_select:", paste0(taxa_rank_select," is currently not supported with this function")))
  } else { #if taxa_rank_select is any of the rank_names(physeq)
        number<-match(taxa_rank_select,rank_type_list)
        if (number==1) {
                  taxa.table_lefse$name_match<-paste(taxa.table_lefse[[1]],rownames(taxa.table_lefse), sep=".") #this is to merge with feature in marker_table(res)  
        } else if (number==2) {
                  taxa.table_lefse$name_match<-paste(taxa.table_lefse[[1]],taxa.table_lefse[[2]],rownames(taxa.table_lefse), sep=".") #this is to merge with feature in marker_table(res)  
        } else if (number==3) {
                  taxa.table_lefse$name_match<-paste(taxa.table_lefse[[1]],taxa.table_lefse[[2]],taxa.table_lefse[[3]],
                                                     rownames(taxa.table_lefse), sep=".") #this is to merge with feature in marker_table(res)  
        } else if (number==4) {
                  taxa.table_lefse$name_match<-paste(taxa.table_lefse[[1]],taxa.table_lefse[[2]],taxa.table_lefse[[3]],
                                                     taxa.table_lefse[[4]],rownames(taxa.table_lefse), sep=".") #this is to merge with feature in marker_table(res)  
        } else if (number==5) {
                  taxa.table_lefse$name_match<-paste(taxa.table_lefse[[1]],taxa.table_lefse[[2]],taxa.table_lefse[[3]],
                                                     taxa.table_lefse[[4]],taxa.table_lefse[[5]],rownames(taxa.table_lefse), sep=".") #this is to merge with feature in marker_table(res)  
        } else if (number==6) {
                  taxa.table_lefse$name_match<-paste(taxa.table_lefse[[1]],taxa.table_lefse[[2]],taxa.table_lefse[[3]],
                                                     taxa.table_lefse[[4]],taxa.table_lefse[[5]],taxa.table_lefse[[6]],rownames(taxa.table_lefse), sep=".") #this is to merge with feature in marker_table(res)  
        } else if (number==7) {
                  taxa.table_lefse$name_match<-paste(taxa.table_lefse[[1]],taxa.table_lefse[[2]],taxa.table_lefse[[3]],
                                                     taxa.table_lefse[[4]],taxa.table_lefse[[5]],taxa.table_lefse[[6]],taxa.table_lefse[[7]],rownames(taxa.table_lefse), sep=".") #this is to merge with feature in marker_table(res)  
        } else if (number==8) {
                  taxa.table_lefse$name_match<-paste(taxa.table_lefse[[1]],taxa.table_lefse[[2]],taxa.table_lefse[[3]],
                                                     taxa.table_lefse[[4]],taxa.table_lefse[[5]],taxa.table_lefse[[6]],taxa.table_lefse[[7]],taxa.table_lefse[[8]],rownames(taxa.table_lefse), sep=".") #this is to merge with feature in marker_table(res)  
        }
        max_name_match_length<-max(nchar(taxa.table_lefse$name_match)) #get the maximum length of the name_match. if the length is more than 30 then will shorten the name
        taxa.table_lefse$feature<-taxa.table_lefse[[taxa_rank_select]]
        if (max_name_match_length>30) {
            #shorten taxa name
            taxa.table_lefse$Taxa2 <- paste0(substr(sapply(strsplit(taxa.table_lefse$name_match,".",fixed=T),"[[",1),1,6),
                                             paste0("..",paste0(substr(sapply(strsplit(taxa.table_lefse$name_match,".",fixed=T),"[[",2),1,10),".."),
                                                    paste0(substr(sapply(strsplit(taxa.table_lefse$name_match,".",fixed=T),"[[",3),1,10),paste0("..",str_sub(taxa.table_lefse$name_match,-4,-1)))))
        } else{ #if the length of the new combine name is less than 30 then no need to condense the name
            taxa.table_lefse$Taxa2<-taxa.table_lefse$name_match
        }
        taxa.table_lefse$rowname_taxa<-row.names(taxa.table_lefse)
  }
  if (taxa_rank_name=="yes"){
    taxa.table_lefse$Taxa2<-taxa.table_lefse$feature
  }
  resdf <- merge(marker_table(res) %>% data.frame(), taxa.table_lefse, by="feature")
  
  taxa.table_lefse<<-taxa.table_lefse
  resdf<<-resdf
  #sort the result by the variable_to_compare (which is same as enrich_group) group and then descending order. the new variable indexx will be the sort index
  resdf<-resdf %>% group_by(enrich_group) %>% arrange(desc(ef_lda)) 
  resdf<-resdf[order(match(resdf$enrich_group, outcome_to_compare)),]
  resdf$indexx<-nrow(resdf):1
  #create an index dataframe which contain the order the graph will be display
  index_df<-subset(resdf,select=c(Taxa2,indexx))
  index_df$Taxa2<-sub("^","ID",index_df$Taxa2)
  
  index_df<<-index_df
  ######################################################################################################################
  ######################################################################################################################
  ####merge in decontam list to change the contaminant taxa red
  ####merge in decontam list to change the contaminant taxa red
  if (!missing(decontam)) { #if decontam option is selected, determine if the df object exist
    if (exists(decontam)=="FALSE"){
      print(paste0("decontam list:",paste0(decontam," does not exist")))
      stop("please run function: decontaminant_subplot_KW prior to using this option- decontam") 
    } else{
      resdf_original<-resdf
      resdf<-merge(resdf, get(decontam),by.x="Taxa",by.y="row.names",All.x=all)
      contaminant_color <- ifelse(resdf$contaminant == "TRUE", "red", "black")
    }
  }  else{
    contaminant_color<-"black" ###if no decontam list was supplied, then will have all the label as black 
  }

  ######################################################################################################################
  ######################################################################################################################
  if (taxa_rank_select=="none"){
    xlabeling<-"ASV"
  } else {
    xlabeling<-paste0("Taxa rank:",taxa_rank_select)
  }
  
  #get the bar graph of the lefse result
  g<- resdf %>% 
    mutate(Taxa2=fct_reorder(Taxa2,indexx)) %>%
    ggplot(aes(x = Taxa2, y = ef_lda)) +
    geom_col(aes(fill=enrich_group), width = 0.8)+
    scale_fill_manual(values=outcome_to_compare_color_t) +
    theme(legend.text=element_text(size=20*(font_size_adjust_percent/100)),legend.title=element_text(size=20*(font_size_adjust_percent/100)),
          axis.text.y=element_text(colour=contaminant_color,face="bold",size=16*(font_size_adjust_percent/100)), axis.title.y=element_text(size=18*(font_size_adjust_percent/100)),
          axis.text.x=element_text(size=17*(font_size_adjust_percent/100)), axis.title.x=element_text(size=18*(font_size_adjust_percent/100)),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    xlab(xlabeling) +
    ylab("LDA (log10)")+
    coord_flip()
  
  ###export result table
  csv_output<-paste0(output_name,".csv")
  write.csv(resdf,file=csv_output, row.names = TRUE)
  
  if (graph_option=="none") {
    png_output<-paste0(output_name,".png")
    png(file=png_output, res = 300, width=width, height=height , units='mm')
    show(g)
    dev.off() 
  }
  else { ##if graph_option is bargraph or boxplot then will get the relative abundance
          #####here is to get the other graph which contains the relative abundance result
          if (log_scale=="yes") {
            # need to have relative abundance multiple by 100 to get percentage. 
            # because some of the relative abundance is very small and less than 1 even after multiply by 100. will add 1 to it so log of the result number will always be positive
                    normalizeSample = function(x) { (x/sum(x))*100+1} 
          } else {
                    normalizeSample = function(x) { x/sum(x)}
          }
          
          #getting the relative abundance
          rel.otu.total = data.frame(otu_table(transformSampleCounts(phyloseq_update_sub, normalizeSample)))
          taxa_name<- subset(taxa.table_lefse,select=c(Taxa2))
          rel.otu.total<<-rel.otu.total
          taxa_name<<-taxa_name
          #get the new shortened taxa name for each ASV merge with the OTU table so the name is standardized
          rel.otu.total2 <- merge(rel.otu.total,taxa_name,by="row.names") #get the Taxa2 which was previously created
          #removing the new colum "Row.names which is generated by default as a result of merging by row name
          row.names(rel.otu.total2)<-rel.otu.total2$Row.names 
          rel.otu.total2<- subset(rel.otu.total2,select=-c(Row.names))
          rel.otu.total2<<-rel.otu.total2
          
          #keeping ASV which were significant in lefse
          for (i in 1:length(outcome_to_compare)) {
            resc<-resdf[resdf$enrich_group==outcome_to_compare[[i]],]
            rel.otuC.save1<-rel.otu.total2[resc$rowname_taxa,]
            row.names(rel.otuC.save1)<-sub("^","ID",rel.otuC.save1$Taxa2) #add ID in front of Taxa2 name to avoid issues with being column name (cannot be number to start with)
            rel.otuC.save2<-subset(rel.otuC.save1,select=-c(Taxa2))
            assign(paste0("rel.otuC.save_",i), data.frame(t(rel.otuC.save2)))
          }
          #cbind all of the rel.otuC.save#
          rel.otu.total3<-rel.otuC.save_1
          for (i in 2:length(outcome_to_compare)) {
            rel.otu.total3<-cbind(get(paste0("rel.otuC.save_",i)), rel.otu.total3)
          }
          rel.otu.total3<<-rel.otu.total3
          
          #OTU dataframe does not contain the variable_to_compare so need to get the sample_data from the phyloseq
          sample_data<-data.frame(sample_data(phyloseq_update))
          sample_data<-subset(sample_data,select=c(variable_to_compare))
          
          #merging the relative abudance OTU to the sample_data 
          rel.otu.total4<-merge(rel.otu.total3,sample_data,by="row.names")
          row.names(rel.otu.total4)<-rel.otu.total4$Row.names #previous merge step generate a new column call Row.names
          rel.otu.total4<- subset(rel.otu.total4,select=-c(Row.names))
          
          wanted<-which(colnames(rel.otu.total4) %in% variable_to_compare)
          rel.otu.total4$group_types<-rel.otu.total4[,wanted]
          rel.otu.total4<-rel.otu.total4[,-wanted]
          rel.otu.total4<<-rel.otu.total4
          #transforming the data from wide to long so we can get the average of relative abundance 
          # based on variable_to_compare and the features (which is the significant ASV based on lefse)
          rel.otu.total4_box <- melt(setDT(rel.otu.total4), id.vars = c("group_types"), variable.name = "features",value.name ="value")
          rel.otu.total4_box$value<- as.numeric(rel.otu.total4_box$value)
          rel.otu.total4_box$value[rel.otu.total4_box$value==Inf]<-0
          #the melt step above would remove some punctations but not all. this gsub step would remove all the special characters so the merge can happen without problem
          rel.otu.total4_box$features<-gsub("[[:punct:]]", "", rel.otu.total4_box$features)
          index_df$Taxa2<-gsub("[[:punct:]]", "", index_df$Taxa2)
          #merging rel.otu.total4_box with the index data frame 
          rel.otu.total4_box<-merge(rel.otu.total4_box,index_df, by.x="features", by.y="Taxa2")
          rel.otu.total4_box<-rel.otu.total4_box[order(match(rel.otu.total4_box$group_types, outcome_to_compare)),]
        
          rel.otu.total4_box<<-rel.otu.total4_box
          
          if (graph_option=="boxplot") {
            q<-rel.otu.total4_box %>% 
              mutate(features=fct_reorder(features,indexx)) %>%
              ggplot(aes(x = features, y = value, fill = group_types)) +
              stat_boxplot(geom ='errorbar', width = 0.02) +
              geom_boxplot(width = 0.8) + 
              scale_fill_manual(values=outcome_to_compare_color_t) +
              coord_flip() +
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                    axis.text.x=element_text(size=17*(font_size_adjust_percent/100)), axis.title.x=element_text(size=18*(font_size_adjust_percent/100)),
                    legend.text=element_text(size=20*(font_size_adjust_percent/100)),legend.title=element_text(size=20*(font_size_adjust_percent/100)))
              if (log_scale=="yes") {
                q<-q+ scale_y_log10()+
                  ylab("Log(Mean Relative Abundance+1)")
              } else {
                q<-q+ ylab("Mean Relative Abundance")
              }
        
          }
          
          mean<-rel.otu.total4 %>%
            group_by(group_types) %>%
            summarise_all(list(~mean(.)))
          rel.otu.total4_mean<-data.frame(t(mean))
          rel.otu.total4_mean<-rel.otu.total4_mean %>% row_to_names(row_number = 1)
          rel.otu.total4_mean$features<-row.names(rel.otu.total4_mean)
          rel.otu.total4_mean <- melt(setDT(rel.otu.total4_mean), id.vars = c("features"), variable.name = "group_types",value.name ="mean")
          rel.otu.total4_mean$mean<-as.numeric(rel.otu.total4_mean$mean)
          #the melt step above would remove some punctations but not all. this gsub step would remove all the special characters so the merge can happen without problem
          rel.otu.total4_mean$features<-gsub("[[:punct:]]", "", rel.otu.total4_mean$features)
          rel.otu.total4_mean<<-rel.otu.total4_mean
          
          #merging rel.otu.total4_box with the index data frame 
          rel.otu.total4_mean<-merge(rel.otu.total4_mean,index_df, by.x="features", by.y="Taxa2", all=T)
          rel.otu.total4_mean2<<-rel.otu.total4_mean
          
          sam_size<-ncol(rel.otu.total4)
          sd<-rel.otu.total4 %>%
            group_by(group_types) %>%
            summarise_all(list(~sd(.)))
          rel.otu.total4_sd<-data.frame(t(sd))
          rel.otu.total4_sd<-rel.otu.total4_sd %>% row_to_names(row_number = 1)
          rel.otu.total4_sd$features<-row.names(rel.otu.total4_sd)
          rel.otu.total4_sd <- melt(setDT(rel.otu.total4_sd), id.vars = c("features"), variable.name = "group_types",value.name ="sd")
          rel.otu.total4_sd$se<-as.numeric(rel.otu.total4_sd$sd)/sqrt(sam_size)
          #the melt step above would remove some punctations but not all. this gsub step would remove all the special characters so the merge can happen without problem
          rel.otu.total4_sd$features<-gsub("[[:punct:]]", "", rel.otu.total4_sd$features)
          rel.otu.total4_sd<<-rel.otu.total4_sd
          
          rel.otu.total4_final<-merge(rel.otu.total4_mean,rel.otu.total4_sd, by=c("features","group_types"), all=T)
          
          rel.otu.total4_final<<-rel.otu.total4_final
          if (graph_option=="bargraph") {
            q<-rel.otu.total4_final %>% 
              mutate(features=fct_reorder(features,indexx)) %>%
              ggplot(aes(x = features, y = mean, fill = group_types)) +
              geom_bar(stat = "identity", width = 0.8, position = position_dodge(0.8)) + 
              geom_errorbar( aes(x=features, ymin=mean-se, ymax=mean+se), width=0.5, colour="black", position = position_dodge(0.8),  size=0.8) +
              scale_fill_manual(values=outcome_to_compare_color_t) +
              coord_flip() +
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                    axis.text.x=element_text(size=17*(font_size_adjust_percent/100)), axis.title.x=element_text(size=18*(font_size_adjust_percent/100)),
                    legend.text=element_text(size=20*(font_size_adjust_percent/100)),legend.title=element_text(size=20*(font_size_adjust_percent/100)))
              if (log_scale=="yes") {
                q<-q+ scale_y_log10()+ ylab("Log(Mean Relative Abundance+1)")
              } else {
                q<-q+ ylab("Mean Relative Abundance")
              }
            
          }
          remove_y<- theme(  axis.text.y = element_blank(),
                             axis.ticks.y = element_blank(),
                             axis.title.y = element_blank())
          plot.all <- list(g, q  +remove_y)
          png_output<-paste0(output_name,".png")
          png(file=png_output, res = 300, width=width, height=height , units='mm')
          show(ggpubr::ggarrange(plotlist=plot.all, ncol=2,nrow=1, widths = c(0.5, 0.8), common.legend=T, legend="bottom"))
          dev.off() 
  }
}


##multi_compare_KW: 
#multi_compare_KW(input_phyloseq=Transplant_PS_all, 
#               variable_to_compare="sample_type",
#               control_group=,
#               treatment_group=c("Lower", "Upper"),
#               treatment_group_color=c("steelblue2", "purple2"),
#               adj_p_cut_off=0.2,
#               number_display=10,
#               legend_onplot="no",
#               diff_method="boxplot",
#               output_name="Lefse_bronch_Lower_Upper2", 
#               decontam="Contam_list_NC_BKG_compare_Lower",
#               width=500,
#               height=400,
#               font_size_adjust_percent=100)
multi_compare_KW <- function (input_phyloseq, 
                              variable_to_compare,
                              control_group,
                              control_group_color,
                              treatment_group=list(), 
                              treatment_group_color=list(),
                              adj_p_cut_off=0.2,
                              number_display=10,
                              legend_onplot=c("yes","no"),
                              diff_method=c("Edger","Deseq"),
                              output_name,
                              decontam=NULL,
                              width=500, height=400,
                              font_size_adjust_percent=100,
                              taxa_genus_output) {
  ###ensure intput_phyloseq is a absolute count phyloseq, not relative abundance phyloseq. and make sure phyloseq is in right orientation
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  countData_check <- floor(as(otu_table(input_phyloseq), "matrix")) # round down to nearest integer. if this phsyloeq is a relative abundance, then the entire countData_check would be 0
  if (all(countData_check==0)) {stop("Please use phyloseq with absolute count. The current input physloeq contains relative abundnace")}
  
  if (missing(taxa_genus_output)) { taxa_genus_output<-"no"} #default being yes for taxa_rank_name, which will just give the selected taxa name
  
  ## evaluate choices
  legend_onplot <- match.arg(legend_onplot)
  diff_method <- match.arg(diff_method)
  
  for (i in 1:length(treatment_group)) {
    assign(paste0("treatment_group",i), treatment_group[[i]])
    assign(paste0("treatment_group_color",i), treatment_group_color[[i]])
  }  
  
  #set color for each group
  treatment_group_color_t<-treatment_group_color
  names(treatment_group_color_t) <-treatment_group 
  
  ##coldata is the meta data
  coldata<-as.data.frame(sample_data(input_phyloseq))
  #getting taxa name
  taxa.table <- as.data.frame(tax_table(input_phyloseq))
  countdata<-as.data.frame(otu_table(input_phyloseq))
  row.names(countdata)<-paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
  taxa.table$combine_name<-paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
  #Convert to numeric
  countdata= data.frame(lapply(countdata, function(x) as.numeric(as.character(x))),
                        check.names=F, row.names = rownames(countdata))
  #Set order
  coldata_SAMPLE <- data.frame(coldata[order(row.names(coldata)),],check.names=F)
  countdata_SAMPLE <- data.frame(countdata[, order(colnames(countdata))],check.names=F)
  
  taxa_name<- subset(taxa.table, select=c(combine_name))
  coldata_SAMPLE_org<-coldata_SAMPLE
  countdata_SAMPLE_org<-countdata_SAMPLE
  
  #setting up a list to loop- make a graph for up regulated results and a graph for down regualated results
  regulation<-c("up", "down")
  
  if (diff_method=="Edger") {
    h<-1
    for(treat_group in treatment_group) { #loop thro treatment_group to get Edger with control vs each treatment choices within treatment_group
                print(treat_group)
                coldata_SAMPLE<-coldata_SAMPLE_org #starting with the original dataframes
                countdata_SAMPLE<-countdata_SAMPLE_org #starting with the original dataframes
                phyloseq_edited<-input_phyloseq #starting with the original phyloseq
                
                #making a compare variable to standardize the input variable for comparsion- make sure the reference is always control
                coldata_SAMPLE$variable_to_compare<-coldata_SAMPLE[[variable_to_compare]]
                coldata_SAMPLE<- coldata_SAMPLE %>% dplyr::mutate(compare=case_when((variable_to_compare==control_group)~"control", (variable_to_compare==treat_group)~"treatment"))
                coldata_SAMPLE<-coldata_SAMPLE[!is.na(coldata_SAMPLE$compare),] ###dropping subjects if they are not in the group of comparsion
                coldata_SAMPLE$compare<- as.factor(coldata_SAMPLE$compare) #converting variable to factor
                countdata_SAMPLE<-countdata_SAMPLE[,colnames(countdata_SAMPLE) %in% rownames(coldata_SAMPLE)]#making sure coldata has same subjects and countdata
                
                #Set order
                coldata_SAMPLE <- data.frame(coldata_SAMPLE[order(row.names(coldata_SAMPLE)),])
                countdata_SAMPLE <- data.frame(countdata_SAMPLE[, order(colnames(countdata_SAMPLE))])
                #set up the phyloseq to make sure the control and treatment groups are correct. remove any observation that do not belong to the control and treatment group
                sample_data(phyloseq_edited)$variable_to_compare <-sample_data(phyloseq_edited)[[variable_to_compare]]
                phyloseq_edited2 <- phyloseq_edited %>% ps_mutate(compare=case_when((variable_to_compare==control_group)~"control", (variable_to_compare==treat_group)~"treatment"))
                keep_sample<- get_variable(phyloseq_edited2,"compare") %in% c("control", "treatment")
                phyloseq_edited3 <- prune_samples(keep_sample, phyloseq_edited2)  ###keeping only samples with the outcome variable of choice 
                
                #Run EDGER
                ######from phyloseq_to_edgeR
                #https://joey711.github.io/phyloseq-extensions/edgeR.html 
                if( !taxa_are_rows(phyloseq_edited3) ){ phyloseq_edited3 <- t(phyloseq_edited3) }
                x = as(otu_table(phyloseq_edited3), "matrix")
                # Add one to protect against overflow, log(0) issues.
                x = x + 1
                # Check `compare` argument
                if( identical(all.equal(length("compare"), 1), TRUE) & nsamples(phyloseq_edited3) > 1 ){
                  # Assume that compare was a sample variable name (must be categorical)
                  group = get_variable(phyloseq_edited3, "compare")
                }
                # Define gene annotations (`genes`) as tax_table
                taxonomy = tax_table(phyloseq_edited3, errorIfNULL=FALSE)
                if( !is.null(taxonomy) ){
                  taxonomy = data.frame(as(taxonomy, "matrix"))
                } 
                # Now turn into a DGEList
                y = DGEList(counts=x, group=group, genes=taxonomy,remove.zeros = TRUE)
                # Calculate the normalization factors
                z = calcNormFactors(y, method="TMM")
                # Check for division by zero inside `calcNormFactors`
                if( !all(is.finite(z$samples$norm.factors)) ){
                  stop("Something wrong with edgeR::calcNormFactors on this data,
                                       non-finite $norm.factors, consider changing `method` argument")
                }
                # Estimate dispersions
                dgeFull<-estimateTagwiseDisp(estimateCommonDisp(z))
                dgeTest<-exactTest(dgeFull)

                #Create Table
                resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table))
                resNoFilt_org<<-resNoFilt
                resNoFilt<-merge(resNoFilt,taxa_name, by="row.names")
                row.names(resNoFilt)<-resNoFilt$combine_name 
                resNoFilt<- subset(resNoFilt,select=c(logFC, logCPM, PValue, FDR)) 
                resNoFilt$grouping<-treat_group

                #Create Table of Up Regulated
                sigUpReg <- resNoFilt[resNoFilt$FDR<adj_p_cut_off,] #keeping only adjust p value < cut off 
                sigUpReg <- sigUpReg[sigUpReg$logFC>=0,] ##keeping only the positive logFC
                sigUpReg <- sigUpReg[!is.na(sigUpReg$logFC),] #remove missing logFC
                sigUpReg <- sigUpReg[order(sigUpReg$logFC, decreasing=TRUE),]
                
                #######$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#######
                ###keeping only the down regulated results- direction of the control
                sigDownReg <- resNoFilt[resNoFilt$FDR<adj_p_cut_off,] #keeping only adj p < cut off 
                sigDownReg <- sigDownReg[sigDownReg$logFC<=0,] ##keeping only the positive logFC
                sigDownReg <- sigDownReg[!is.na(sigDownReg$logFC),] #remove missing logFC
                sigDownReg <- sigDownReg[order(sigDownReg$logFC, decreasing=FALSE),]
                #######$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#######
                #######$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#######
              
                ##################################################################################################################################
                ##################################################################################################################################
                #Create Relative Abundance Table
                df <-
                  countdata_SAMPLE %>% 
                  rownames_to_column('gs') %>%
                  dplyr::group_by(gs) %>% 
                  dplyr::summarise_all(funs(sum)) %>%
                  dplyr::mutate_if(is.numeric, funs(./sum(.))) %>%
                  column_to_rownames('gs')
                #Get the ColData for Each Comparison
                coldata.1 <- coldata_SAMPLE[coldata_SAMPLE[[variable_to_compare]]==control_group,]
                coldata.2 <- coldata_SAMPLE[coldata_SAMPLE[[variable_to_compare]]==treat_group,]
                #keep Count data only for each comparison
                needed<-which(colnames(df) %in% rownames(coldata.1))    
                df.1 <- df[,needed]
                needed2<-which(colnames(df) %in% rownames(coldata.2))    
                df.2 <- df[,needed2]
                
                for (reg_group in regulation) {
                        #depends on if we want to get the upregulated or downregulated results, will set res to the desire results within the loop 
                        if (reg_group=="up") {  
                                res <- sigUpReg %>% dplyr::slice(1:number_display) #Subset Top (number_display)
                        } else if (reg_group=="down"){
                                res <- sigDownReg %>% dplyr::slice(1:number_display) #Subset Top (number_display)
                        }   

                        #Convert Resuts table into a data.frame
                        res <- as.data.frame(res)
                        #decide what otu to save 
                        otu.to.save <-as.character(rownames(res))
                        #from relative table we should get the mean across the row of the otu table
                        df.1.meanRA <- rowMeans(df.1)
                        df.2.meanRA <- rowMeans(df.2)
                        #need to subset AND reorder just the otus that we have 
                        df.1.meanRA.save <- df.1.meanRA[otu.to.save]
                        df.2.meanRA.save <- df.2.meanRA[otu.to.save]
                        #add the abundnace data for the res dataframe
                        res$abundance.1 <- df.1.meanRA.save
                        res$abundance.2 <- df.2.meanRA.save
                        #Set Names of Results Table
                        res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Taxa","logFC", "lfcSE", "pvalue", "adj.P.Val","group","abundance.1","abundance.2")) 
                        #Remove Any Data without LOGFC data
                        res <- res[!is.na(res$logFC),]
                        # Reorder Results based on FDR for comparison 1
                        res = res[order(-res$logFC, na.last = TRUE), ]
                        #Create Order
                        res <- res %>% dplyr::mutate(start = 1:n())
                        #Convert Important columns to Numeric
                        res$adj.P.Val <-   as.numeric(as.character(res$adj.P.Val))
                        res$logFC <-       as.numeric(as.character(res$logFC))
                        #Replace NA
                        res <- res %>% dplyr::mutate(adj.P.Val = if_else(is.na(adj.P.Val), 0.9, adj.P.Val))
                        ##convert abundance to numeric
                        res$abundance.1 <- as.numeric(as.character(res$abundance.1))
                        res$abundance.2 <- as.numeric(as.character(res$abundance.2))
                        #Create Variable for Color based on Comparison, FDR and LOGFC
                        if (reg_group=="up") {  
                          res$abundance <-res$abundance.2
                        } else if (reg_group=="down"){
                          res$abundance <-res$abundance.1
                        } 
                        
                        ###reassign taxa name on the graph
                        if (taxa_genus_output=="no") {
                          res$Taxa2 <- paste0(substr(sapply(strsplit(res$Taxa,".",fixed=T),"[[",1),1,6),
                                              paste0("..",paste0(substr(sapply(strsplit(res$Taxa,".",fixed=T),"[[",4),1,10),".."),
                                                     paste0(substr(sapply(strsplit(res$Taxa,".",fixed=T),"[[",5),1,10),paste0("..",str_sub(res$Taxa,-4,-1)))))  
                        } else if (taxa_genus_output=="yes") { #keeping only genus level. if genus level is NA, then will go up a taxa rank
                          res<- res %>% dplyr::mutate(Taxa2=case_when((sapply(strsplit(Taxa,".",fixed=T),"[[",1)=="NA")~"k_Bacteria.p_NA.c_NA.o_NA.f_NA.g_NA",
                                                                      ((sapply(strsplit(Taxa,".",fixed=T),"[[",1)!="NA") & (sapply(strsplit(Taxa,".",fixed=T),"[[",2)=="NA"))~paste0(paste0("p_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",1)),".c_NA.o_NA.f_NA.g_NA"), 
                                                                      ((sapply(strsplit(Taxa,".",fixed=T),"[[",2)!="NA") & (sapply(strsplit(Taxa,".",fixed=T),"[[",3)=="NA"))~paste0(paste0("c_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",2)),".o_NA.f_NA.g_NA"), 
                                                                      ((sapply(strsplit(Taxa,".",fixed=T),"[[",3)!="NA") & (sapply(strsplit(Taxa,".",fixed=T),"[[",4)=="NA"))~paste0(paste0("o_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",3)),".f_NA.g_NA"), 
                                                                      ((sapply(strsplit(Taxa,".",fixed=T),"[[",4)!="NA") & (sapply(strsplit(Taxa,".",fixed=T),"[[",5)=="NA"))~paste0(paste0("f_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",4)),".g_NA"),
                                                                      ((sapply(strsplit(Taxa,".",fixed=T),"[[",5)!="NA"))~paste0("g_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",5)) ))
                        }
                      assign(paste0(paste0("res",h),reg_group),res)
                }
                #########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#########
                #########EXCEL EXPORT
                ####keeping all results 
                sig_select_all<-resNoFilt
                #Convert Resuts table into a data.frame
                sig_select_all <- as.data.frame(sig_select_all)
                #decide what otu to save 
                otu.to.save_all <-as.character(rownames(sig_select_all))
                #from relative table we should get the mean across the row of the otu table
                df.1.meanRA <- rowMeans(df.1)
                df.2.meanRA <- rowMeans(df.2)
                #need to subset AND reorder just the otus that we have 
                df.1.meanRA.save <- df.1.meanRA[otu.to.save_all]
                df.2.meanRA.save <- df.2.meanRA[otu.to.save_all]
                #add the abundnace data for the res dataframe
                sig_select_all$abundance.1 <- df.1.meanRA.save
                sig_select_all$abundance.2 <- df.2.meanRA.save
                #Set Names of Results Table
                sig_select_all <- setNames(cbind(rownames(sig_select_all), sig_select_all, row.names = NULL), c("Taxa","logFC", "lfcSE", "pvalue", "adj.P.Val","group","abundance.1","abundance.2")) 
                assign(paste0("sig_select_all",h),sig_select_all)
                #########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#########
                
                h<-h+1
    }
    
    
  } else if (diff_method=="Deseq") {
    h<-1
    for(treat_group in treatment_group) { #loop thro treatment_group to get Edger with control vs each treatment choices within treatment_group
                print(treat_group)
                coldata_SAMPLE<-coldata_SAMPLE_org #starting with the original dataframes
                countdata_SAMPLE<-countdata_SAMPLE_org #starting with the original dataframes
                phyloseq_edited<-input_phyloseq #starting with the original phyloseq
                
                #making a compare variable to standardize the input variable for comparsion- make sure the reference is always control
                coldata_SAMPLE$variable_to_compare<-coldata_SAMPLE[[variable_to_compare]]
                coldata_SAMPLE<- coldata_SAMPLE %>% dplyr::mutate(compare=case_when((variable_to_compare==control_group)~"control", (variable_to_compare==treat_group)~"treatment"))
                coldata_SAMPLE<-coldata_SAMPLE[!is.na(coldata_SAMPLE$compare),] ###dropping subjects if they are not in the group of comparsion
                coldata_SAMPLE$compare<- as.factor(coldata_SAMPLE$compare) #converting variable to factor
                countdata_SAMPLE<-countdata_SAMPLE[,colnames(countdata_SAMPLE) %in% rownames(coldata_SAMPLE)]#making sure coldata has same subjects and countdata
                
                #Set order
                coldata_SAMPLE <- data.frame(coldata_SAMPLE[order(row.names(coldata_SAMPLE)),])
                countdata_SAMPLE <- data.frame(countdata_SAMPLE[, order(colnames(countdata_SAMPLE))])
                #set up the phyloseq to make sure the control and treatment groups are correct. remove any observation that do not belong to the control and treatment group
                sample_data(phyloseq_edited)$variable_to_compare <-sample_data(phyloseq_edited)[[variable_to_compare]]
                phyloseq_edited2 <- phyloseq_edited %>% ps_mutate(compare=case_when((variable_to_compare==control_group)~"control", (variable_to_compare==treat_group)~"treatment"))
                keep_sample<- get_variable(phyloseq_edited2,"compare") %in% c("control", "treatment")
                phyloseq_edited3 <- prune_samples(keep_sample, phyloseq_edited2)  ###keeping only samples with the outcome variable of choice 
                
                ###get taxa name list so it can be merge back to the result from deseq
                taxa.table_V2 <- as.data.frame(tax_table(phyloseq_edited3))
                taxa.table_V2$combine_name<-paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
                taxa_name<- subset(taxa.table_V2, select=c(combine_name))
                # Calculate geometric means prior to estimate size factors
                gm_mean = function(x, na.rm = TRUE){
                  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
                }
                
                #Create new DESEQ Objects        
                diagdds<- phyloseq_to_deseq2(phyloseq_edited3, ~ compare)
                geoMeans = apply(counts(diagdds), 1, gm_mean)
                diagdds <- estimateSizeFactors(diagdds, geoMeans = geoMeans)
                #making sure reference is control. resetting the level to make sure control is reference
                diagdds$compare <- droplevels(diagdds$compare)
                diagdds$compare <- relevel(diagdds$compare, ref ="control")
                diagdds = estimateDispersions(diagdds, fitType = "parametric")
                diagdds <- nbinomWaldTest(diagdds)
                #Output Result Table
                res <- results(diagdds, cooksCutoff = FALSE)
                res_org<<-res
                
                #Convert Resuts table into a data.frame
                res <- as.data.frame(res)
                res<-merge(res,taxa_name, by="row.names") #merging back the taxa name from the phyloseq
                row.names(res)<-res$combine_name 
                res <- res %>% dplyr::rename(logFC=log2FoldChange)
                res<- subset(res,select=c(baseMean,logFC,lfcSE,stat,pvalue,padj))
                res$grouping<-treat_group #add treatment group
                resNoFilt<- res #resNoFilt is backup object
                #keeping only significant results for the graphs 
                sigUpReg <- resNoFilt[resNoFilt$padj<adj_p_cut_off,] #keeping only adj p < cut off 
                sigUpReg <- sigUpReg[sigUpReg$logFC>=0,] ##keeping only the positive logFC
                sigUpReg <- sigUpReg[!is.na(sigUpReg$logFC),] #remove missing logFC
                sigUpReg <- sigUpReg[order(sigUpReg$logFC, decreasing=TRUE),]

                #######$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#######
                ###keeping only the down regulated results- direction of the control
                sigDownReg <- resNoFilt[resNoFilt$padj<adj_p_cut_off,] #keeping only adj p < cut off 
                sigDownReg <- sigDownReg[sigDownReg$logFC<=0,] ##keeping only the positive logFC
                sigDownReg <- sigDownReg[!is.na(sigDownReg$logFC),] #remove missing logFC
                sigDownReg <- sigDownReg[order(sigDownReg$logFC, decreasing=FALSE),]
                #######$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#######
                #######$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#######
                
                ##################################################################################################################################
                ##################################################################################################################################
                #Create Relative Abundance Table
                df <-
                  countdata_SAMPLE %>% 
                  rownames_to_column('gs') %>%
                  dplyr::group_by(gs) %>% 
                  dplyr::summarise_all(funs(sum)) %>%
                  dplyr::mutate_if(is.numeric, funs(./sum(.))) %>%
                  column_to_rownames('gs')
                #Get the ColData for Each Comparison
                coldata.1 <- coldata_SAMPLE[coldata_SAMPLE[[variable_to_compare]]==control_group,]
                coldata.2 <- coldata_SAMPLE[coldata_SAMPLE[[variable_to_compare]]==treat_group,]
                #keep Count data only for each comparison
                needed<-which(colnames(df) %in% rownames(coldata.1))    
                df.1 <- df[,needed]
                needed2<-which(colnames(df) %in% rownames(coldata.2))    
                df.2 <- df[,needed2]
                for (reg_group in regulation) {
                        #depends on if we want to get the upregulated or downregulated results, will set res to the desire results within the loop 
                        if (reg_group=="up") {  
                          res <- sigUpReg %>% dplyr::slice(1:number_display) #Subset Top (number_display)
                        } else if (reg_group=="down"){
                          res <- sigDownReg %>% dplyr::slice(1:number_display) #Subset Top (number_display)
                        }   
                        #Convert Resuts table into a data.frame
                        res <- as.data.frame(res)
                        #decide what otu to save 
                        otu.to.save <-as.character(rownames(res))
                        #from relative table we should get the mean across the row of the otu table
                        df.1.meanRA <- rowMeans(df.1)
                        df.2.meanRA <- rowMeans(df.2)
                        #need to subset AND reorder just the otus that we have 
                        df.1.meanRA.save <- df.1.meanRA[otu.to.save]
                        df.2.meanRA.save <- df.2.meanRA[otu.to.save]
                        #add the abundnace data for the res dataframe
                        res$abundance.1 <- df.1.meanRA.save
                        res$abundance.2 <- df.2.meanRA.save
                        #Set Names of Results Table
                        res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Taxa","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","group","abundance.1","abundance.2"))  
                        #Remove Any Data without LOGFC data
                        res <- res[!is.na(res$logFC),]
                        # Reorder Results based on FDR for comparison 1
                        res = res[order(-res$logFC, na.last = TRUE), ]
                        #Create Order
                        res <- res %>% dplyr::mutate(start = 1:n())
                        #Convert Important columns to Numeric
                        res$adj.P.Val <-   as.numeric(as.character(res$adj.P.Val))
                        res$logFC <-       as.numeric(as.character(res$logFC))
                        #Replace NA
                        res <- res %>% dplyr::mutate(adj.P.Val = if_else(is.na(adj.P.Val), 0.9, adj.P.Val))
                        ##convert abundance to numeric
                        res$abundance.1 <- as.numeric(as.character(res$abundance.1))
                        res$abundance.2 <- as.numeric(as.character(res$abundance.2))
                        #Create Variable for Color based on Comparison, FDR and LOGFC
                        if (reg_group=="up") {  
                                res$abundance <-res$abundance.2
                        } else if (reg_group=="down"){
                                res$abundance <-res$abundance.1
                        } 
                        res$Taxa2 <- paste0(substr(sapply(strsplit(res$Taxa,".",fixed=T),"[[",1),1,6),
                                            paste0("..",paste0(substr(sapply(strsplit(res$Taxa,".",fixed=T),"[[",4),1,10),".."),
                                                   paste0(substr(sapply(strsplit(res$Taxa,".",fixed=T),"[[",5),1,10),paste0("..",str_sub(res$Taxa,-4,-1)))))
                        assign(paste0(paste0("res",h),reg_group),res)
                }           
                #########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#########
                #########EXCEL EXPORT
                ####keeping all results 
                sig_select_all<-resNoFilt
                #Convert Resuts table into a data.frame
                sig_select_all <- as.data.frame(sig_select_all)
                #decide what otu to save 
                otu.to.save_all <-as.character(rownames(sig_select_all))
                #from relative table we should get the mean across the row of the otu table
                df.1.meanRA <- rowMeans(df.1)
                df.2.meanRA <- rowMeans(df.2)
                #need to subset AND reorder just the otus that we have 
                df.1.meanRA.save <- df.1.meanRA[otu.to.save_all]
                df.2.meanRA.save <- df.2.meanRA[otu.to.save_all]
                #add the abundnace data for the res dataframe
                sig_select_all$abundance.1 <- df.1.meanRA.save
                sig_select_all$abundance.2 <- df.2.meanRA.save
                #Set Names of Results Table
                sig_select_all <- setNames(cbind(rownames(sig_select_all), sig_select_all, row.names = NULL), c("Taxa","baseMean", "logFC", "lfcSE", "stat", "pvalue",  "adj.P.Val","group","abundance.1","abundance.2")) 
                assign(paste0("sig_select_all",h),sig_select_all)
                #########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#########                
                h<-h+1
    }
  }
  C_resup<-res1up
  C_resdown<-res1down
  C_sig_select_all<-sig_select_all1
  for(l in 2:length(treatment_group)) { #loop thro treatment_group to get Edger with control vs each treatment choices within treatment_group
    C_resup<-rbind(C_resup, get(paste0(paste0("res",l),"up")))
    C_resdown<-rbind(C_resdown, get(paste0(paste0("res",l),"down")))
    C_sig_select_all<-rbind(C_sig_select_all, get(paste0("sig_select_all",l))) 
  }
  csv_output<-paste0(output_name,"_all.csv")
  write.csv(C_sig_select_all,file=csv_output, row.names = TRUE)
  
  #determine the min and max abundance for both upregulated and downregulated results
  #this is used to standardized the abundance circle size for upregulated and downregulated results
  min_C_resup<-min(C_resup$abundance)
  max_C_resup<-max(C_resup$abundance)
  min_C_resdown<-min(C_resdown$abundance)
  max_C_resdown<-max(C_resdown$abundance)
  
  min_result<-min(min_C_resup,min_C_resdown)
  max_result<-max(max_C_resup,max_C_resdown)
  
    ######%%%%%%%%%%%%%%%%%%%%%%%%%######
    ######%%%%%%%%%%%%%%%%%%%%%%%%%######
    ###### for upregulated results ######
    ######%%%%%%%%%%%%%%%%%%%%%%%%%######
    ######%%%%%%%%%%%%%%%%%%%%%%%%%######
              ######################################################################################################################
              #set order for C_resup then make an index variable
              C_resup<-C_resup[with(C_resup, order(group, start)),] 
              C_resup$indexx<-nrow(C_resup):1
              #make a spacing between group so easier to see 
              C_resup<-C_resup[order(C_resup$indexx, decreasing = FALSE), ]  
              C_resup<-setDT(C_resup)[ , ID := .GRP, by = group]   
              C_resup$ID<-C_resup$ID-1
              C_resup$indexx<-C_resup$indexx+C_resup$ID
              #getting group size for the grey line which divide the different group
              groups_size <- C_resup %>% dplyr::group_by(group)
              n_groups(groups_size)
              group_size_list<-group_size(groups_size)
              print(group_size_list)
              group_size_list<-rev(group_size_list)
              group_size_list<-cumsum(group_size_list)
              temp_list<-unlist(list(1:length(treatment_group)))
              group_size_list<-group_size_list+temp_list
              
              C_resup<-as.data.frame(C_resup) #so it can be merge with decontam. the previous step setDT made the dataframe into a data.table dataframe which has issues merging with decontam 
              ######################################################################################################################
              ####merge in decontam list to change the contaminant taxa red
              if (!missing(decontam)) { #if decontam option is selected, determine if the df object exist
                if (exists(decontam)=="FALSE"){
                  print(paste0("decontam list:",paste0(decontam," does not exist")))
                  stop("please run function: decontaminant_subplot_KW prior to using this option- decontam") 
                } else{
                  C_resup<-merge(C_resup, get(decontam),by.x="Taxa",by.y="row.names",All.x=all)
                  C_resup<- C_resup[order(C_resup$logFC, na.last = TRUE), ]
                  contaminant_color <- ifelse(C_resup$contaminant == "TRUE", "red", "black")
                }
              }  else{
                contaminant_color<-"black" ###if no decontam list was supplied, then will have all the label as black 
              }
              ######################################################################################################################
              ######################################################################################################################
              
              p<- C_resup %>%
                ggplot(aes(y=indexx, x=logFC, fill=group, size=abundance)) +
                geom_point(color="black",alpha=0.8,shape=21)+
                geom_segment(aes(yend=indexx), xend=(-30), color= "black", linetype = "solid",linewidth=1)+ 
                scale_fill_manual(values=treatment_group_color_t)+ 
                scale_size_continuous(name="Relative Abundance",range=c(5, 20), limit=c(min_result, max_result))+  
                scale_y_continuous(expand=c(0,0),limits=c(0,length(C_resup$indexx)+length(treatment_group)),breaks=C_resup$indexx,labels=C_resup$Taxa2)+
                theme(panel.background = element_blank(),
                      panel.border=element_rect(fill=NA),
                      panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
                      panel.grid.minor = element_blank(),
                      strip.background=element_blank(),
                      axis.title=element_text(size=20*(font_size_adjust_percent/100),face="bold"),
                      axis.text.x=element_text(colour="black", size=18*(font_size_adjust_percent/100), face="bold"),
                      axis.text.y=element_text(colour=contaminant_color,face="bold",size=10*(font_size_adjust_percent/100)),
                      axis.ticks=element_line(colour="black"),
                      legend.background = element_rect(color=NA),
                      legend.key = element_rect(colour = "transparent", fill = "white"),
                      plot.title=element_text(size=23*(font_size_adjust_percent/100), face="bold"))+
                xlab("logFC") +
                ylab("")+
                #xlim(-7,7)+
                geom_vline(xintercept=0, color="red",linetype="dashed") + 
                geom_hline(yintercept=group_size_list, color="grey30",linetype="dotted") + guides(fill = guide_legend(override.aes = list(size = 5) ) )
              ##save the legend on a separate png file              
              if (legend_onplot=="no") { #if want to hide the legend on the plot, then will save the legend on separate png
                leg <- get_legend(p)
                png_output2<-paste0(output_name,"_legend.png")
                png(png_output2, res = 300, width=width, height=height, units='mm')
                legend<-as_ggplot(leg)
                show(legend)
                dev.off()
              }
              if (diff_method=="Edger") {
                p<-p+  ggtitle("EdgeR")
              } else if (diff_method=="Deseq") {
                p<-p+  ggtitle("Deseq")
              }
              
              png_output<-paste0(output_name,".png")
              png(file=png_output, res = 300, width=width, height=height , units='mm')
              if (legend_onplot=="no") { #if want to hide the legend
                p<-p + theme(legend.position = "none") ###remove the legend
              }
              show(p)
              dev.off()
 
      ######%%%%%%%%%%%%%%%%%%%%%%%%%%%######
      ######%%%%%%%%%%%%%%%%%%%%%%%%%%%######
      ###### for downregulated results ######
      ######%%%%%%%%%%%%%%%%%%%%%%%%%%%######
      ######%%%%%%%%%%%%%%%%%%%%%%%%%%%######
              ######################################################################################################################
              #set order for C_resdown then make an index variable
              C_resdown<-C_resdown[with(C_resdown, order(group, start)),] 
              C_resdown$indexx<-1:nrow(C_resdown)
              #make a spacing between group so easier to see 
              C_resdown<-C_resdown[order(C_resdown$indexx, decreasing = FALSE), ]  
              C_resdown<-setDT(C_resdown)[ , ID := .GRP, by = group]   
              C_resdown$ID<-C_resdown$ID-1
              C_resdown$indexx<-C_resdown$indexx+C_resdown$ID
              #getting group size for the grey line which divide the different group
              groups_size <- C_resdown %>% dplyr::group_by(group)
              n_groups(groups_size)
              group_size_list<-group_size(groups_size)
              print(group_size_list)
              group_size_list<-rev(group_size_list)
              group_size_list<-cumsum(group_size_list)
              temp_list<-unlist(list(1:length(treatment_group)))
              group_size_list<-group_size_list+temp_list
              
              C_resdown<-as.data.frame(C_resdown) #so it can be merge with decontam. the previous step setDT made the dataframe into a data.table dataframe which has issues merging with decontam 
              ######################################################################################################################
              ####merge in decontam list to change the contaminant taxa red
              if (!missing(decontam)) { #if decontam option is selected, determine if the df object exist
                if (exists(decontam)=="FALSE"){
                  print(paste0("decontam list:",paste0(decontam," does not exist")))
                  stop("please run function: decontaminant_subplot_KW prior to using this option- decontam") 
                } else{
                  C_resdown<-merge(C_resdown, get(decontam),by.x="Taxa",by.y="row.names",All.x=all)
                  C_resdown<- C_resdown[order(C_resdown$logFC, na.last = TRUE), ]
                  contaminant_color <- ifelse(C_resdown$contaminant == "TRUE", "red", "black")
                }
              }  else{
                contaminant_color<-"black" ###if no decontam list was supplied, then will have all the label as black 
              }
              ######################################################################################################################
              ######################################################################################################################
              p<- C_resdown %>%
                ggplot(aes(y=indexx, x=logFC, fill=control_group_color, size=abundance)) +
                geom_point(color="black",fill=control_group_color,alpha=0.8,shape=21)+
                geom_segment(aes(yend=indexx), xend=(-30), color= "black", linetype = "solid",linewidth=1)+ 
                scale_size_continuous(name="Relative Abundance",range=c(5, 20), limit=c(min_result, max_result))+  
                scale_y_continuous(breaks=C_resdown$indexx,labels=C_resdown$Taxa2)+
                facet_grid(group~.,scales="free",space="free")+   
                theme(panel.background = element_blank(),
                      panel.border=element_rect(fill=NA),
                      panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
                      panel.grid.minor = element_blank(),
                      strip.background=element_blank(),
                      axis.title=element_text(size=20*(font_size_adjust_percent/100),face="bold"),
                      axis.text.x=element_text(colour="black", size=18*(font_size_adjust_percent/100), face="bold"),
                      axis.text.y=element_text(colour=contaminant_color,face="bold",size=10*(font_size_adjust_percent/100)),
                      axis.ticks=element_line(colour="black"),
                      legend.background = element_rect(color=NA),
                      legend.key = element_rect(colour = "transparent", fill = "white"),
                      plot.title=element_text(size=23*(font_size_adjust_percent/100), face="bold"),
                      panel.spacing.y = unit(0,"line"))+
                xlab("logFC") +
                ylab("")+
                #xlim(-7,7)+
                geom_vline(xintercept=0, color="red",linetype="dashed")  + guides(size = guide_legend(override.aes = list(fill = "white") ) )
                ##save the legend on a separate png file              
                if (legend_onplot=="no") { #if want to hide the legend on the plot, then will save the legend on separate png
                  leg <- get_legend(p)
                  png_output2<-paste0(output_name,"_control_legend.png")
                  png(png_output2, res = 300, width=width, height=height, units='mm')
                  legend<-as_ggplot(leg)
                  show(legend)
                  dev.off()
                }
              if (diff_method=="Edger") {
                p<-p+  ggtitle("EdgeR")
              } else if (diff_method=="Deseq") {
                p<-p+  ggtitle("Deseq")
              }
              
              png_output<-paste0(output_name,"_control.png")
              png(file=png_output, res = 300, width=width, height=height , units='mm')
              if (legend_onplot=="no") { #if want to hide the legend
                p<-p + theme(legend.position = "none") ###remove the legend
              }
              show(p)
              dev.off()              
              
}

#keep_taxa_KW finds out which ASV is contain in keppTaxa phyloseq and then keeping those particular ASVs on physeq
keep_taxa_KW <- function(physeq, keepTaxa){
  allTaxa = taxa_names(physeq) #get a list of the ASV that is in the original physloeq
  keepTaxa_name = taxa_names(keepTaxa) #get a list of the ASV that is in the keepTaxa phyloseq
  myTaxa <- allTaxa[(allTaxa %in% keepTaxa_name)] #determine which ASV to keep
  process_physeq<-prune_taxa(myTaxa, physeq) #prune the original physeq of interest based on the ASV that you want to keep 
  print(process_physeq) #display the result
  return(process_physeq) #return the result (post process physloeq )
}

##Edger_paired_phylo_KW: performs paired EdgeR analysis
#Edger_paired_phylo_KW(input_phyloseq=Transplant_PS_all, 
#               Subj_var="SubjID"
#               variable_to_compare="condition",
#               outcome_to_compare=c("conditionA", "conditionB"),
#               outcome_to_compare_color=c("steelblue2", "purple2"),
#               FDR_cut_off=0.2,
#               legend_onplot="yes",
#               output_name="Edger_bronch_Lower_Upper", 
#               decontam="Contam_list_NC_BKG_compare_Lower")
Edger_paired_phylo_KW <- function (input_phyloseq, 
                            Subj_var, 
                            variable_to_compare,
                            outcome_to_compare=list(),
                            outcome_to_compare_color=list(),
                            FDR_cut_off=0.2,
                            number_display=10,
                            output_name,
                            legend_onplot=c("yes","no"),
                            decontam=NULL,
                            taxa_genus_output) {
  ###ensure intput_phyloseq is a absolute count phyloseq, not relative abundance phyloseq. and make sure phyloseq is in right orientation
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  countData_check <- floor(as(otu_table(input_phyloseq), "matrix")) # round down to nearest integer. if this phsyloeq is a relative abundance, then the entire countData_check would be 0
  if (all(countData_check==0)) {stop("Please use phyloseq with absolute count. The current input physloeq contains relative abundnace")}
  
  if (missing(taxa_genus_output)) { taxa_genus_output<-"no"} #default being yes for taxa_rank_name, which will just give the selected taxa name
  
  ## evaluate choices
  legend_onplot <- match.arg(legend_onplot)
  
  for (i in 1:length(outcome_to_compare)) {
    assign(paste0("outcome_to_compare",i), outcome_to_compare[[i]])
    assign(paste0("outcome_to_compare_color",i), outcome_to_compare_color[[i]])
  }  
  
  ##coldata is the meta data
  coldata<-as.data.frame(sample_data(input_phyloseq))
  #getting taxa name
  taxa.table <- as.data.frame(tax_table(input_phyloseq))
  countdata<-as.data.frame(otu_table(input_phyloseq))
  row.names(countdata)<-paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
  taxa.table$combine_name<-paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
  #Convert to numeric
  countdata= data.frame(lapply(countdata, function(x) as.numeric(as.character(x))),
                        check.names=F, row.names = rownames(countdata))
  #Set order
  coldata_SAMPLE <- data.frame(coldata[order(row.names(coldata)),],check.names=F)
  countdata_SAMPLE <- data.frame(countdata[, order(colnames(countdata))],check.names=F)
  
  #making a compare variable to standardize the input variable for comparsion- make sure the reference is always control
  coldata_SAMPLE$variable_to_compare<-coldata_SAMPLE[[variable_to_compare]]
  coldata_SAMPLE<- coldata_SAMPLE %>% dplyr::mutate(compare=case_when((variable_to_compare==outcome_to_compare1)~"control", (variable_to_compare==outcome_to_compare2)~"treatment"))
  coldata_SAMPLE<-coldata_SAMPLE[!is.na(coldata_SAMPLE$compare),] ###dropping subjects if they are not in the group of comparsion
  coldata_SAMPLE$compare<- as.factor(coldata_SAMPLE$compare) #converting variable to factor
  countdata_SAMPLE<-countdata_SAMPLE[,colnames(countdata_SAMPLE) %in% rownames(coldata_SAMPLE)]#making sure coldata has same subjects and countdata
  
  taxa_name<- subset(taxa.table, select=c(combine_name))
  
  #Set order
  coldata_SAMPLE <- data.frame(coldata_SAMPLE[order(row.names(coldata_SAMPLE)),])
  countdata_SAMPLE <- data.frame(countdata_SAMPLE[, order(colnames(countdata_SAMPLE))])
  
  #set up the phyloseq to make sure the control and treatment groups are correct. remove any observation that do not belong to the control and treatment group
  phyloseq_edited<-input_phyloseq
  sample_data(phyloseq_edited)$variable_to_compare <-sample_data(phyloseq_edited)[[variable_to_compare]]
  phyloseq_edited2 <- phyloseq_edited %>% ps_mutate(compare=case_when((variable_to_compare==outcome_to_compare1)~"control", (variable_to_compare==outcome_to_compare2)~"treatment"))
  keep_sample<- get_variable(phyloseq_edited2,"compare") %in% c("control", "treatment")
  phyloseq_edited3 <- prune_samples(keep_sample, phyloseq_edited2)  ###keeping only samples with the outcome variable of choice 
          ####paired testing
          Condition <- get_variable(phyloseq_edited3, variable_to_compare)
          Subject <- get_variable(phyloseq_edited3, Subj_var)
          print(sum(tapply(Condition, Subject, length) == 2))
          sample_data(phyloseq_edited3)$Subj_var <-sample_data(phyloseq_edited3)[[Subj_var]]
          keepSubject <- names(which(tapply(Condition, Subject, function (x){length(unique(x))})==2))
          phyloseq_edited4 <- subset_samples(phyloseq_edited3, Subj_var %in% keepSubject) # only keeping subjects with the paired
          #reset Condition and Subject now that only the paired samples are keeped
          Condition = get_variable(phyloseq_edited4, variable_to_compare)
          Subject= get_variable(phyloseq_edited4, Subj_var)
          design = model.matrix( ~ Subject + Condition)
          # Add one to protect against overflow, log(0) issues.
          x <- as(otu_table(phyloseq_edited4), "matrix") + 1L
          taxonomy <- data.frame(as(tax_table(phyloseq_edited4), "matrix"))
          # Now turn into a DGEList
          x <- DGEList(counts=x, group=Condition, genes=taxonomy, remove.zeros=TRUE)
          # Calculate the normalization factors and estimate dispersion
          x <- calcNormFactors(x, method="RLE")
          x <- estimateGLMCommonDisp(x, design)
          x <- estimateGLMTrendedDisp(x, design)
          x <- estimateGLMTagwiseDisp(x, design)
          fit <- glmFit(x, design)
          lrt <- glmLRT(fit)
  
  #Create Table
  resNoFilt <- topTags(lrt, n=nrow(x), adjust.method="BH", sort.by="PValue")
  resNoFilt_org<<-resNoFilt
  resNoFilt<-merge(resNoFilt,taxa_name, by="row.names")
  row.names(resNoFilt)<-resNoFilt$combine_name 
  resNoFilt<- subset(resNoFilt,select=c(logFC, logCPM, PValue, FDR))
  resNoFilt_export<<-resNoFilt
  #check to see if there is any taxa with FDR<FDR_cut_off. If not then will stop the function, because there will be error if it continues as the object will have no rows
  checkifallzero<-resNoFilt[resNoFilt$FDR<FDR_cut_off,]
  if (dim(checkifallzero)[1] == 0) {
    print(paste0("no Taxa with FDR more than",FDR_cut_off))
    return(NULL)
  }
  resNoFilt[nrow(resNoFilt) + 1,] <- c(0, 0, 0, 0) ###adding this row of 0 so that even if there is no significant top 10 for upregulated or downregulated, the code will still run
  #Create Table of Down Regulated
  sigDownReg <- resNoFilt[resNoFilt$FDR<FDR_cut_off,] #keeping only FDR < 0.2
  #sigDownReg <- resNoFilt
  sigDownReg <- sigDownReg[sigDownReg$logFC<=0,] ##keeping only the negative logFC
  sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
  sigDownReg1 <- sigDownReg %>% dplyr::slice(1:number_display) #Subset Top (number_display)
  #Create Table of Down Regulated
  sigUpReg <- resNoFilt[resNoFilt$FDR<FDR_cut_off,] #keeping only FDR < 0.2
  #sigUpReg <- resNoFilt
  sigUpReg <- sigUpReg[sigUpReg$logFC>=0,] ##keeping only the positive logFC
  sigUpReg <- sigUpReg[order(sigUpReg$logFC, decreasing=TRUE),]
  sigUpReg1 <- sigUpReg %>% dplyr::slice(1:number_display) #Subset Top (number_display)
  #Merge Top 10s
  res <- rbind(sigDownReg1,sigUpReg1)
  res <- res[res$logFC!=0,] 
  
  #Create Relative Abundance Table
  df <-
    countdata_SAMPLE %>% 
    rownames_to_column('gs') %>%
    dplyr::group_by(gs) %>% 
    dplyr::summarise_all(funs(sum)) %>%
    dplyr::mutate_if(is.numeric, funs(./sum(.))) %>%
    column_to_rownames('gs')
  df_export<<-df
  #Get the ColData for Each Comparison
  coldata.1 <- coldata_SAMPLE[coldata_SAMPLE[[variable_to_compare]]==outcome_to_compare1,]
  coldata.2 <- coldata_SAMPLE[coldata_SAMPLE[[variable_to_compare]]==outcome_to_compare2,]
  #keep Count data only for each comparison
  needed<-which(colnames(df) %in% rownames(coldata.1))    
  df.1 <- df[,needed]
  needed2<-which(colnames(df) %in% rownames(coldata.2))    
  df.2 <- df[,needed2]
  #Convert Resuts table into a data.frame
  res <- as.data.frame(res)
  #decide what otu to save 
  otu.to.save <-as.character(rownames(res))
  #from relative table we should get the mean across the row of the otu table
  df.1.meanRA <- rowMeans(df.1)
  df.2.meanRA <- rowMeans(df.2)
  #need to subset AND reorder just the otus that we have 
  df.1.meanRA.save <- df.1.meanRA[otu.to.save]
  df.2.meanRA.save <- df.2.meanRA[otu.to.save]
  #add the abundnace data for the res dataframe
  res$abundance.1 <- df.1.meanRA.save
  res$abundance.2 <- df.2.meanRA.save
  #Set Names of Results Table
  res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Taxa","logFC", "lfcSE", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
  #Remove Any Data without LOGFC data
  res <- res[!is.na(res$logFC),]
  # Reorder Results based on FDR for comparison 1
  res = res[order(-res$logFC, na.last = TRUE), ]
  #Create Order
  res <- res %>% dplyr::mutate(start = 1:n())
  #Convert Important columns to Numeric
  res$adj.P.Val <-   as.numeric(as.character(res$adj.P.Val))
  res$logFC <-       as.numeric(as.character(res$logFC))
  #Replace NA
  res <- res %>% dplyr::mutate(adj.P.Val = if_else(is.na(adj.P.Val), 0.9, adj.P.Val))
  ##convert abundance to numeric
  res$abundance.1 <- as.numeric(as.character(res$abundance.1))
  res$abundance.2 <- as.numeric(as.character(res$abundance.2))
  #Create Variable for Color based on Comparison, FDR and LOGFC
  res$col <- ifelse(res$logFC>0, "B",
                    ifelse(res$logFC<0, "A","D"))
  res$abundance <- ifelse(res$col=="A", res$abundance.1, ifelse(res$col=="B", res$abundance.2,0))
  if (taxa_genus_output=="no") {
    res$Taxa2 <- paste0(substr(sapply(strsplit(res$Taxa,".",fixed=T),"[[",1),1,6),
                        paste0("..",paste0(substr(sapply(strsplit(res$Taxa,".",fixed=T),"[[",4),1,10),".."),
                               paste0(substr(sapply(strsplit(res$Taxa,".",fixed=T),"[[",5),1,10),paste0("..",str_sub(res$Taxa,-4,-1)))))  
  } else if (taxa_genus_output=="yes") { #keeping only genus level. if genus level is NA, then will go up a taxa rank
    res<- res %>% dplyr::mutate(Taxa2=case_when((sapply(strsplit(Taxa,".",fixed=T),"[[",1)=="NA")~"k_Bacteria.p_NA.c_NA.o_NA.f_NA.g_NA",
                                                ((sapply(strsplit(Taxa,".",fixed=T),"[[",1)!="NA") & (sapply(strsplit(Taxa,".",fixed=T),"[[",2)=="NA"))~paste0(paste0("p_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",1)),".c_NA.o_NA.f_NA.g_NA"), 
                                                ((sapply(strsplit(Taxa,".",fixed=T),"[[",2)!="NA") & (sapply(strsplit(Taxa,".",fixed=T),"[[",3)=="NA"))~paste0(paste0("c_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",2)),".o_NA.f_NA.g_NA"), 
                                                ((sapply(strsplit(Taxa,".",fixed=T),"[[",3)!="NA") & (sapply(strsplit(Taxa,".",fixed=T),"[[",4)=="NA"))~paste0(paste0("o_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",3)),".f_NA.g_NA"), 
                                                ((sapply(strsplit(Taxa,".",fixed=T),"[[",4)!="NA") & (sapply(strsplit(Taxa,".",fixed=T),"[[",5)=="NA"))~paste0(paste0("f_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",4)),".g_NA"),
                                                ((sapply(strsplit(Taxa,".",fixed=T),"[[",5)!="NA"))~paste0("g_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",5)) ))
  }
  #########EXCEL EXPORT
  ####keeping all results 
  res_all <- rbind(sigDownReg,sigUpReg)
  res_all <- res_all[res_all$logFC!=0,] 
  
  #Convert Resuts table into a data.frame
  res_all <- as.data.frame(res_all)
  #decide what otu to save 
  otu.to.save_all <-as.character(rownames(res_all))
  #from relative table we should get the mean across the row of the otu table
  df.1.meanRA <- rowMeans(df.1)
  df.2.meanRA <- rowMeans(df.2)
  #need to subset AND reorder just the otus that we have 
  df.1.meanRA.save <- df.1.meanRA[otu.to.save_all]
  df.2.meanRA.save <- df.2.meanRA[otu.to.save_all]
  #add the abundnace data for the res dataframe
  res_all$abundance.1 <- df.1.meanRA.save
  res_all$abundance.2 <- df.2.meanRA.save
  #Set Names of Results Table
  res_all <- setNames(cbind(rownames(res_all), res_all, row.names = NULL), c("Taxa","logFC", "lfcSE", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
  csv_output<-paste0(output_name,"_significant.csv")
  write.csv(res_all,file=csv_output, row.names = TRUE)
  csv_output2<-paste0(output_name,"_all.csv")
  write.csv(resNoFilt_export,file=csv_output2, row.names = TRUE)
  
  ######################################################################################################################
  ######################################################################################################################
  ####merge in decontam list to change the contaminant taxa red
  if (!missing(decontam)) { #if decontam option is selected, determine if the df object exist
    if (exists(decontam)=="FALSE"){
      print(paste0("decontam list:",paste0(decontam," does not exist")))
      stop("please run function: decontaminant_subplot_KW prior to using this option- decontam") 
    } else{
      res_original<-res
      res<-merge(res, get(decontam),by.x="Taxa",by.y="row.names",All.x=all)
      res<- res[order(res$logFC, na.last = TRUE), ]
      contaminant_color <- ifelse(res$contaminant == "TRUE", "red", "black")
    }
  }  else{
    contaminant_color<-"black" ###if no decontam list was supplied, then will have all the label as black 
  }
  ######################################################################################################################
  ######################################################################################################################
  p<-ggplot(res, aes(y=reorder(Taxa2,-start), x=logFC,fill=col,size=abundance)) +
    geom_point(color="black",alpha=0.8,shape=21)+
    geom_segment(data=res[res$adj.P.Val<0.2,],aes(yend=reorder(Taxa2,-start)), xend=(-30), color= "black", linetype = "solid",linewidth=1)+ 
    scale_fill_manual(values=c("B"=outcome_to_compare_color2,"A"=outcome_to_compare_color1,"D"="white"))+ 
    scale_size_continuous(name="Relative Abundance",range=c(5, 20))+                
    ggtitle("EdgeR")+
    theme(panel.background = element_blank(),
          panel.border=element_rect(fill=NA),
          panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
          panel.grid.minor = element_blank(),
          strip.background=element_blank(),
          axis.title=element_text(size=20,face="bold"),
          axis.text.x=element_text(colour="black", size=18, face="bold"),
          axis.text.y=element_text(colour=contaminant_color,face="bold",size=10),
          axis.ticks=element_line(colour="black"),
          legend.background = element_rect(color=NA),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          plot.title=element_text(size=23, face="bold"))+
    xlab("") +
    ylab("")+
    #xlim(-7,7)+
    geom_vline(xintercept=0, color="red",linetype="dashed")+
    guides(fill="none")
  ##save the legend on a separate png file              
  if (legend_onplot=="no") { #if want to hide the legend on the plot, then will save the legend on separate png
    leg <- get_legend(p)
    png_output2<-paste0(output_name,"_legend.png")
    png(png_output2, res = 300, width=200, height=200, units='mm')
    legend<-as_ggplot(leg)
    show(legend)
    dev.off()
  }
  png_output<-paste0(output_name,".png")
  png(file=png_output, res = 300, width=200, height=200 , units='mm')
  if (legend_onplot=="no") { #if want to hide the legend
    p<-p + theme(legend.position = "none") ###remove the legend
  }
  show(p)
  dev.off()
  
  ##########################################################################################################################
  ##########################################################################################################################
}

###DMM_cluster_KW: perform DMM clustering
DMM_cluster_KW <- function (input_phyloseq, 
                                   variable_to_compare,
                                   max_cluster=5,
                                   phyloseq_uniqueID, 
                                   output_name,
                                   DMM_cluster_color=list(), 
                                   compare_stats=c("chisquare","fisher")){

  ###ensure intput_phyloseq is a absolute count phyloseq, not relative abundance phyloseq. and make sure phyloseq is in right orientation
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  countData_check <- floor(as(otu_table(input_phyloseq), "matrix")) # round down to nearest integer. if this phsyloeq is a relative abundance, then the entire countData_check would be 0
  if (all(countData_check==0)) {stop("Please use phyloseq with absolute count. The current input physloeq contains relative abundnace")}
  
  #set random seed
  set.seed(1234)
  tse <- makeTreeSummarizedExperimentFromPhyloseq(input_phyloseq)
  ###max_cluster determine maximum cluster to test (default is 5)
  max_cluster_display<-paste0("max cluster: ", max_cluster)
  print(max_cluster_display)
  tse_dmn <- mia::runDMN(tse, name="DMN", k = 1:max_cluster)
  #display the stats for each iterations
  getDMN(tse_dmn)
  
  plotDMNfit_output<- paste0(paste0("DMM_model_fit_graph",output_name),".pdf")
  plotDMNFit<-plotDMNFit(tse_dmn, type = "laplace")
  pdf(file=plotDMNfit_output, width = 16, height = 18)
            show(plotDMNFit)     
  dev.off()
  #displaying the best fit DMM clustering
  print("Best fit:")
  print(getBestDMNFit(tse_dmn, type = "laplace"))
  DMN_best_fit<-bestDMNFit(tse_dmn, type = "laplace") #number of clusters which best fit model
  
  #create a data frame laplace which contains the goodness of fit for the DMM models with the corresponding number of Dirichlet components
  laplace <- c()
  for(i in 1:(length(getDMN(tse_dmn)))){
    laplace[i] <- getDMN(tse_dmn)[[i]]@goodnessOfFit[[3]]
  }
  laplace %>% data.frame() %>%
    dplyr::rename("Model fit" = ".") %>%
    dplyr::mutate("Number of Dirichlet components" = rownames(.))
  # Group samples and return DMNGroup object that contains a summary. pgd_highest_48_72_v2 is used for grouping.
  # this provides summary for the group samples with DMM clusters
  dmn_group <- calculateDMNgroup(tse_dmn, variable = variable_to_compare,  exprs_values = "counts",
                                 k =DMN_best_fit, seed=1234)

  # Mixture weights (rough measure of the cluster size)
  ## pi gives you the component weight
  ## theta gives you component variability (small values of theta correspond to highly variable components)
  DMM_comp_weight<-DirichletMultinomial::mixturewt(getBestDMNFit(tse_dmn))
  
  # Samples-cluster assignment 
  DMM_cluster_assign <- data.frame(DirichletMultinomial::mixture(getBestDMNFit(tse_dmn), assign=T))
  colnames(DMM_cluster_assign)<-c("Cluster_num")
  DMM_cluster_assign$Cluster_num<-as.factor(DMM_cluster_assign$Cluster_num)
  DMM_cluster_assign[phyloseq_uniqueID]<-row.names(DMM_cluster_assign) #add the unique identifier to the DMM_cluster_assign dataframe
  
  assign(paste0("DMM_cluster_KW_",output_name),DMM_cluster_assign, envir=.GlobalEnv)
  print("###########################")
  print("###########################")
  print("Available DMM assignment:")
  print(paste0("DMM_cluster_KW_",output_name))
  print("###########################")
  print("###########################")
  
  #write DMM csv file
  csv.output.DMM<-paste0(paste0("DMM_cluster_KW_",output_name),".csv")
  write.csv(DMM_cluster_assign,csv.output.DMM)
  
  #DMM_cluster_assign tells you which cluster each samples belong to
  DMM_sample_var_interest<- data.frame(sample_data(input_phyloseq))%>% dplyr::select(variable_to_compare)
  
  DMM_cluster_assign_v2 <- merge(DMM_cluster_assign, DMM_sample_var_interest, by="row.names")
  DMM_cluster_assign_v3 <- DMM_cluster_assign_v2 %>% dplyr::select(-c("Row.names"))
  DMM_cluster_assign_v4 <- DMM_cluster_assign_v3 %>% dplyr::group_by(DMM_cluster_assign_v3[[variable_to_compare]],Cluster_num, .drop=FALSE) %>% tally()
  colnames(DMM_cluster_assign_v4)<-c(variable_to_compare,"Cluster_num","n")
  DMM_cluster_assign_v5 <- data.frame(DMM_cluster_assign_v4)
  DMM_cluster_assign_v5$compare<-DMM_cluster_assign_v5[[variable_to_compare]]
  DMM_cluster_assign_v6 <- DMM_cluster_assign_v5 %>% 
    dplyr::group_by(compare) %>%
    dplyr::mutate(perc =n/sum(n)) %>%
    dplyr::arrange(perc) %>%
    dplyr::mutate(labels=scales::percent(perc)) #get percentage for the piechart
  #Graphs on the DMM clusters by groups
  #choose color for the cluster group
  cluster_color_t<-c("cornflowerblue","hotpink1","cyan4","goldenrod3",'azure4') #default color if color is not specified 
  names(cluster_color_t) <- c("1","2","3","4","5")
  if (!missing(DMM_cluster_color)){ #if DMM_cluster_color is not missing then will replace the color 
        for(i in 1:length(DMM_cluster_color)) {
          cluster_color_t[i]<-DMM_cluster_color[i] #replace the selected default color with selected 
        }
  } 
  graph_DMM_all<- c() #make an empty list
  tt<-1
  for(i in sort(unique(unlist(DMM_cluster_assign_v6[variable_to_compare])))){
    DMM_cluster_assign_v6<-DMM_cluster_assign_v6 %>% dplyr::arrange(compare, Cluster_num)
    DMM_cluster_assign_v7<-DMM_cluster_assign_v6[DMM_cluster_assign_v6[[variable_to_compare]] %in% i,]
    graph_DMM<-DMM_cluster_assign_v7 %>%
      ggplot(aes(x="", y=n, fill=Cluster_num)) +
      geom_bar(stat="identity", width=1, color="white") +
      geom_label(aes(label = labels, size=14, fontface ="bold"), color = "white", 
                 position = position_stack(vjust = 0.5), 
                 show.legend = FALSE) +
      coord_polar("y", start=0) +
      theme_void() +     
      scale_fill_manual(values=cluster_color_t, drop=FALSE) +
      theme(legend.position="none",plot.title = element_text(size = 25, face = "bold"))+
      ggtitle(i)
    assign(paste0("graph_DMM_",tt),graph_DMM)
    graph_DMM_all[[tt]] <- graph_DMM
    tt<-tt+1
  }
  #combine_DMM_plot will give you pie charts  
  combine_DMM_plot<-ggarrange(plotlist=graph_DMM_all, nrow=1, common.legend = TRUE, legend="bottom")
  combine_DMM_plot_output<- paste0(paste0("DMM_model_cluster_group",output_name),".pdf")
  pdf(file=combine_DMM_plot_output, width = 5*(tt-1), height = 5)
             show(combine_DMM_plot)              
  dev.off()
  
  compare_stats <- match.arg(compare_stats)
  print(paste0("Compare statistics:",compare_stats))
  #set compare.test to either fisher test or chi-square
  if (compare_stats=="chisquare"){
        #chi square to compare the cluster groups 
        compare.test <<- function(a, b) {
          return(chisq.test(cbind(a, b)))
        }
  } else if (compare_stats=="fisher") {
        compare.test <<- function(a, b) {
          print("a")
          print(a)
          print("b")
          print(b)
        return(fisher.test(cbind(a, b)))
    }
  }
  
  #different combination of comparison
  myComparisons<-as.list(as.data.frame(combn(unique(unlist(DMM_cluster_assign_v2[variable_to_compare])), 2))) 
  combine_DMM_plot_v2<-ggplot(DMM_cluster_assign_v6, aes(x=get(variable_to_compare), y=n))+
    geom_bar(aes(fill=Cluster_num), stat = "identity", position = "dodge")+
    geom_signif(comparisons = myComparisons,
                test = "compare.test", map_signif_level=F, step_increase=0.1)+
    scale_fill_manual(values=cluster_color_t) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.title=element_text(face="bold", size=18),
          axis.text.x =element_text(size=18),
          axis.text.y =element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position="none")
  
  combine_DMM_plot_output2<- paste0(paste0("DMM_model_cluster_group",output_name),"_V2.pdf")
  pdf(file=combine_DMM_plot_output2, width = 5*(tt-1), height = 5)
              show(combine_DMM_plot_v2)
  dev.off()
  
  # Contribution of each taxa to each component
  DMM_taxa_component<-data.frame(DirichletMultinomial::fitted(getBestDMNFit(tse_dmn)))
  
  #adding on the cluster results to the phyloseq sample data
  input_phyloseq2<-ps_join(input_phyloseq, DMM_cluster_assign, by=phyloseq_uniqueID)
 
  name_input_phyloseq<-deparse(substitute(input_phyloseq))
  output_phyloseq<- paste0(name_input_phyloseq,"_DMM") #deparse and substitute turn the phyloseq object name into a character string
  assign(output_phyloseq,input_phyloseq2,.GlobalEnv)
  print("###########################")
  print("###########################")
  print("Available phyloseq with DMM assignment:")
  print(paste0(output_phyloseq))
  print("###########################")
  print("###########################")
}

###distance_lower_to_centroid_KW: calculates distance from lower airway to other sample type centroids. 
#variable_to_compare is optional- if you want to see subgroup difference to upper airway sample centroid. 
#distance_lower_to_centroid_KW(input_phyloseq=Transplant_PS_PGD_GENUS_UL_relab_prune, 
#                           sample_type_var_name="sample_type",
#                           Lower_sample_type="Lower",
#                           Lower_sample_color="steelblue2",
#                           Other_sample_type=c("Upper", "BKG"),
#                           Other_sample_color=c("purple2","darkgoldenrod1"),
#                           output_suffix="PGD",
#                           variable_to_compare="pgd_highest_48_72_v3", 
#                           outcome_to_compare=c("PGD-0or1","PGD-2","PGD-3"), 
#                           outcome_to_compare_color=c("darkgreen","darkorange1","red3"), 
#                           p_value="yes")
distance_lower_to_centroid_KW <- function(input_phyloseq, 
                                       sample_type_var_name, 
                                       Lower_sample_type, 
                                       Lower_sample_color,
                                       Other_sample_type=list(), 
                                       Other_sample_color=list(),
                                       output_suffix,
                                       variable_to_compare,
                                       outcome_to_compare=list(),
                                       outcome_to_compare_color=list(),
                                       beta_method=c("bray", "unifrac"),...) {
  list2env(list(...), environment())
  #default value if not specified
  if(!exists("plot_title_size")) {plot_title_size<-14}#default pvalue is no
  if(!exists("label_size")) {label_size<-7}#default pvalue is no
  if(!exists("axis_title_size")) {axis_title_size<-18}#default pvalue is no
  if(missing(beta_method)) { beta_method<-"bray"} #default being bray curtis

  #display the beta method used 
  print(paste0("beta method: ",beta_method))  
  
  #need the function read_count_KW
  if(exists("read_count_KW") != T) stop("Error: need read_count_KW to use this function")
  #if(!exists("y_log_scale")) {y_log_scale<-"no"}#default y_log_scale is no unless specified
  if(!exists("p_value")) {p_value<-"no"}#default pvalue is no
  #ensure the coloring is correct: 
  stopifnot("Other_sample_type need to have same number of elements as Other_sample_color"= length(Other_sample_type)==length(Other_sample_color))

  # set up coloring for everything
  Lower_Other_color_t<-c(Lower_sample_color,Other_sample_color)
  names(Lower_Other_color_t) <- c(Lower_sample_type,Other_sample_type)
  outcome_to_compare_color_t<-outcome_to_compare_color
  names(outcome_to_compare_color_t) <- outcome_to_compare
  #combine the color for sample type and the outcome to compare variable
  sample_type_color_t<- c(Lower_Other_color_t, outcome_to_compare_color_t)
  
  #ensure that Lower_sample_type and Other_sample_type are all within the sample_type_var_name of the input_phyloseq:
  stopifnot("Input for Lower_sample_type and Other_sample_type have to be elements within sample_type_var_name"= 
              c(Lower_sample_type, Other_sample_type) %in% unique(sample_data(input_phyloseq)[[sample_type_var_name]]))
  
  ### Code with all samples: 
  #Create Distance Matrix with Bray (or wUniFrac depending what you are using)
  #Create Distance Matrix with Bray
  if (beta_method=="bray"){
    #Create Distance Matrix with Bray
    vegdist=vegdist(t(otu_table(input_phyloseq)), method = "bray") #row is samplesID #column is taxa
  } else {
    vegdist= phyloseq::UniFrac(input_phyloseq, weighted = TRUE) 
  }
  #Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
  CmdScale <- cmdscale(vegdist, k =10)
  #calculated Sample variance for each PC
  vars <- apply(CmdScale, 2, var)
  #Create Variable with the Percent Variance
  percentVar <- round(100 * (vars/sum(vars)))
  #Merge PC Data with MetaData
  newResults <- merge(x = CmdScale, y = sample_data(input_phyloseq), by = "row.names", all.x = TRUE)
  newResults$sampletype_KW<-newResults[[sample_type_var_name]]
  #Rename Variables for PC1 and PC2
  colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
  colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
  colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"
  
  # Calculate the Centroid Value
  centroids<-aggregate(cbind(PC1,PC2)~sampletype_KW,data=newResults, mean) #Here you would use your grouping variable (e.g., days)
  print("Centroids:")
  print(centroids)
  
  #PCA for Lower sample only 
  Lower_PCA <- newResults %>% dplyr::filter(sampletype_KW==Lower_sample_type)
  row.names(Lower_PCA)<-Lower_PCA$name
  Lower_PCA_matrix<- Lower_PCA %>% dplyr::select(c("PC1","PC2"))
  Lower_PCA_matrix<-as.matrix(Lower_PCA_matrix)
  
  # Merge the Centroid Data into the PCOA Data
  newResults<-merge(newResults, centroids, by="sampletype_KW", suffixes=c("",".centroid"))
  
  #keeping only lower and other sample (this step is probably not needed if there are no sample types which were not included in the other_sample_type)
  centroids_lower_Other <- centroids %>% filter(sampletype_KW %in% c(Lower_sample_type, Other_sample_type))
  #this gives you the PCA of ONLY lower airway
  newResults_lower<-newResults %>% filter(sampletype_KW==Lower_sample_type)
  
  p<-ggplot(newResults_lower, aes(PC1, PC2)) # Graph PC1 and PC2
  if(!missing(variable_to_compare)){ #if there is variable_to_compare then will need read_count_KW
    p<-p+geom_point(size=2, aes(color=get(variable_to_compare)))  # Set the size of the points
  } else {
    p<-p+geom_point(size=2, aes(color=get(sample_type_var_name)))
  }
  p<- p+ xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2 
    labs(title=NULL, fill="Centroids")+
    geom_segment(alpha=0.2,aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=sampletype_KW ))+ 
    geom_point(data=centroids_lower_Other, aes(x=PC1, y=PC2, color=sampletype_KW), size=8, alpha=5)+
    geom_label_repel(data = centroids_lower_Other, aes(x=PC1, y=PC2, color=sampletype_KW, label=sampletype_KW), size=label_size) +
    scale_color_manual(name=NULL,  values=sample_type_color_t)+
    theme(panel.background = element_blank(),
          panel.border=element_rect(fill=NA),
          panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
          panel.grid.minor = element_blank(),strip.background=element_blank(),
          plot.title=element_text(face="bold",hjust=0.5, size = plot_title_size), 
          plot.subtitle = element_text(hjust=0.5),
          axis.title=element_text(face="bold", size = axis_title_size),
          axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
          axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
          axis.ticks=element_blank(),
          plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
  PCA_output<-paste0(paste0("PCA_lower_other_centroids",output_suffix),".pdf")
  pdf(PCA_output, height = 10, width = 10)
        show(p)
  dev.off()
  
  #define function to calculate distance between distance vectors and the centroid. This function is taken from betadisper (obtained from Fares code)
  Resids <- function(x, c) {
    if(is.matrix(c))
      d <- x - c
    else
      d <- sweep(x, 2, c)
    rowSums(d^2)
  }
  
  i<-1
  #Figure out distance to other centroids 
  for (sample in Other_sample_type) {
          Other_Lower_PCA <- newResults %>% dplyr::filter(sampletype_KW==Lower_sample_type | sampletype_KW==sample)
          row.names(Other_Lower_PCA)<-Other_Lower_PCA$name
          Other_Lower_PCA_matrix<- Other_Lower_PCA %>% dplyr::select(c("PC1","PC2"))
          Other_Lower_PCA_matrix<-as.matrix(Other_Lower_PCA_matrix)
          
          centroids_Other <- centroids %>% filter(sampletype_KW ==sample)
          centroids_Other_matrix<- centroids_Other %>% dplyr::select(c("PC1","PC2"))
          centroids_Other_matrix<- as.matrix(centroids_Other_matrix)
          
          #calculate the distance to Upper centroid for Upper and Lower samples 
          dist_to_centroids_test<-Resids(Other_Lower_PCA_matrix, as.numeric(centroids_Other_matrix))
          dist_to_centroids_test<-as.data.frame(dist_to_centroids_test)
          dist_to_centroids_test$subset_sample<-Other_Lower_PCA[[sample_type_var_name]]
          dist_to_centroids_test<-dist_to_centroids_test %>% dplyr::rename("distance"="dist_to_centroids_test")
          read_count_output_test<-paste0(paste0(paste0("distance_lower_to_",sample),"_"),output_suffix)
          dist_to_centroids_test<-dist_to_centroids_test
          ylabel_t<-paste0(paste0("Distance to ",sample)," airway")
          read_count_KW(input_object=dist_to_centroids_test, 
                        sample_type_var_name="subset_sample", 
                        sample_types=c(Lower_sample_type,sample) , sample_type_color=c(Lower_sample_color, Other_sample_color[i]), y_log_scale = "no", 
                        p_value=p_value, read_count_output_test, ylabel=ylabel_t, count_variable="distance")
          output_dataframe_name<- paste0(paste0(sample,"_to_lower_dist"),output_suffix)
          #export the distance to global environment for other use
          assign(output_dataframe_name,dist_to_centroids_test, envir = .GlobalEnv )
          print("###########################")
          print(paste0(paste0("Available dataframe :distance from lower to ",sample),":"))
          print(output_dataframe_name)
          print("###########################")
          if(!missing(variable_to_compare)){ #if there is variable_to_compare then will need read_count_KW
            #calculate the distance to Upper centroid for compare_variable (stored in dist_to_centroids)
            dist_to_centroids<-Resids(Lower_PCA_matrix, as.numeric(centroids_Other_matrix))
            dist_to_centroids<-as.data.frame(dist_to_centroids)
            dist_to_centroids$subset_sample<-Lower_PCA[[variable_to_compare]]
            dist_to_centroids<-dist_to_centroids %>% dplyr::rename("distance"="dist_to_centroids")
            read_count_output<-paste0(paste0(paste0("distance_lower_to_",sample),"_comparevar_"),output_suffix)
            dist_to_centroids<<-dist_to_centroids
            ylabel_t<-paste0(paste0("Distance to ",sample)," airway")
            read_count_KW(input_object=dist_to_centroids, 
                          sample_type_var_name="subset_sample", 
                          sample_types=outcome_to_compare, sample_type_color=outcome_to_compare_color, y_log_scale = "no",
                          p_value=p_value, output_name=read_count_output, ylabel=ylabel_t, count_variable="distance")
          } 
    i<-i+1
  }
}




Ancombc_phylo_KW <- function (input_phyloseq, 
                            variable_to_compare,
                            outcome_to_compare=list(),
                            outcome_to_compare_color=list(),
                            FDR_cut_off=0.2,
                            number_display=10,
                            graph_option=c("volcano","lollipop"),
                            display_all_results_volcano=c("yes","no"), #display all the results in volcano plot (including results that is not significant: FDR not below the cut off)
                            abundance_size_volcano=c("yes","no"), #change size of the dot based on relative abundance
                            output_name,
                            legend_onplot=c("yes","no"),
                            decontam=NULL,
                            taxa_genus_output, ...) {
  ###ensure intput_phyloseq is a absolute count phyloseq, not relative abundance phyloseq. and make sure phyloseq is in right orientation
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  countData_check <- floor(as(otu_table(input_phyloseq), "matrix")) # round down to nearest integer. if this phsyloeq is a relative abundance, then the entire countData_check would be 0
  if (all(countData_check==0)) {stop("Please use phyloseq with absolute count. The current input physloeq contains relative abundnace")}
  
  if (missing(taxa_genus_output)) { taxa_genus_output<-"no"} #default being yes for taxa_rank_name, which will just give the selected taxa name
  if (missing(graph_option)) { graph_option<-"lollipop"} #default being lollipop for graph_option
  if (missing(display_all_results_volcano)) { display_all_results_volcano<-"no"} #default being lollipop for graph_option
  if (missing(abundance_size_volcano)) { abundance_size_volcano<-"yes"} #default being lollipop for graph_option
  
  ## evaluate choices
  legend_onplot <- match.arg(legend_onplot)
  
  for (i in 1:length(outcome_to_compare)) {
    assign(paste0("outcome_to_compare",i), outcome_to_compare[[i]])
    assign(paste0("outcome_to_compare_color",i), outcome_to_compare_color[[i]])
  }  
  
  ##coldata is the meta data
  coldata<-as.data.frame(sample_data(input_phyloseq))
  #getting taxa name
  taxa.table <- as.data.frame(tax_table(input_phyloseq))
  countdata<-as.data.frame(otu_table(input_phyloseq))
  row.names(countdata)<-paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
  taxa.table$combine_name<-paste(taxa.table[[2]],taxa.table[[3]],taxa.table[[4]],taxa.table[[5]], taxa.table[[6]], taxa.table[[7]],rownames(taxa.table), sep=".")
  #Convert to numeric
  countdata= data.frame(lapply(countdata, function(x) as.numeric(as.character(x))),
                        check.names=F, row.names = rownames(countdata))
  #Set order
  coldata_SAMPLE <- data.frame(coldata[order(row.names(coldata)),],check.names=F)
  countdata_SAMPLE <- data.frame(countdata[, order(colnames(countdata))],check.names=F)
  
  #making a compare variable to standardize the input variable for comparsion- make sure the reference is always control
  coldata_SAMPLE$variable_to_compare<-coldata_SAMPLE[[variable_to_compare]]
  coldata_SAMPLE<- coldata_SAMPLE %>% dplyr::mutate(compare=case_when((variable_to_compare==outcome_to_compare1)~"control", (variable_to_compare==outcome_to_compare2)~"treatment"))
  coldata_SAMPLE<-coldata_SAMPLE[!is.na(coldata_SAMPLE$compare),] ###dropping subjects if they are not in the group of comparsion
  coldata_SAMPLE$compare<- as.factor(coldata_SAMPLE$compare) #converting variable to factor
  countdata_SAMPLE<-countdata_SAMPLE[,colnames(countdata_SAMPLE) %in% rownames(coldata_SAMPLE)]#making sure coldata has same subjects and countdata
  
  taxa_name<- subset(taxa.table, select=c(combine_name))
  
  #Set order
  coldata_SAMPLE <- data.frame(coldata_SAMPLE[order(row.names(coldata_SAMPLE)),])
  countdata_SAMPLE <- data.frame(countdata_SAMPLE[, order(colnames(countdata_SAMPLE))])
  
  #set up the phyloseq to make sure the control and treatment groups are correct. remove any observation that do not belong to the control and treatment group
  phyloseq_edited<-input_phyloseq
  sample_data(phyloseq_edited)$variable_to_compare <-sample_data(phyloseq_edited)[[variable_to_compare]]
  phyloseq_edited2 <- phyloseq_edited %>% ps_mutate(compare=case_when((variable_to_compare==outcome_to_compare1)~"control", (variable_to_compare==outcome_to_compare2)~"treatment"))
  keep_sample<- get_variable(phyloseq_edited2,"compare") %in% c("control", "treatment")
  phyloseq_edited3 <- prune_samples(keep_sample, phyloseq_edited2)  ###keeping only samples with the outcome variable of choice 

  phyloseq_edited4<<-phyloseq_edited3
  
  output=ancombc2(data=phyloseq_edited3, fix_formula="compare", 
                p_adj_method = "holm",
                group="compare", # for the purpose performing structural zeros testing
                prv_cut=0, lib_cut=0, #avoid pruning in this function 
                alpha=FDR_cut_off
           )
  
  raw_final_result<-output$res
  
  raw_final_result_updated<-merge(raw_final_result,taxa.table, by.x="taxon", by.y="row.names")
  row.names(raw_final_result_updated)<-raw_final_result_updated$combine_name
  
  resNoFilt<- subset(raw_final_result_updated,select=c(lfc_comparetreatment, se_comparetreatment, p_comparetreatment, q_comparetreatment))
  resNoFilt_export<<-resNoFilt
  
  #check to see if there is any taxa with FDR<FDR_cut_off. If not then will stop the function, because there will be error if it continues as the object will have no rows
  checkifallzero<-resNoFilt[resNoFilt$q_comparetreatment<FDR_cut_off,]
  if (dim(checkifallzero)[1] == 0) {
    print(paste0("no Taxa with FDR less than",FDR_cut_off))
    return(NULL)
  }
  resNoFilt[nrow(resNoFilt) + 1,] <- c(0, 0, 0, 0) ###adding this row of 0 so that even if there is no significant top 10 for upregulated or downregulated, the code will still run
  #Create Table of Down Regulated
  sigDownReg <- resNoFilt[resNoFilt$q_comparetreatment<FDR_cut_off,] #keeping only FDR < 0.2
  sigDownReg <- sigDownReg[sigDownReg$lfc_comparetreatment<=0,] ##keeping only the negative logFC
  sigDownReg <- sigDownReg[order(sigDownReg$lfc_comparetreatment),]
  sigDownReg1 <- sigDownReg %>% dplyr::slice(1:number_display) #Subset Top (number_display)
  #Create Table of Up Regulated
  sigUpReg <- resNoFilt[resNoFilt$q_comparetreatment<FDR_cut_off,] #keeping only FDR < 0.2
  sigUpReg <- sigUpReg[sigUpReg$lfc_comparetreatment>=0,] ##keeping only the positive logFC
  sigUpReg <- sigUpReg[order(sigUpReg$lfc_comparetreatment, decreasing=TRUE),]
  sigUpReg1 <- sigUpReg %>% dplyr::slice(1:number_display) #Subset Top (number_display)
  #Merge Top results
  res <- rbind(sigDownReg1,sigUpReg1)
  res <- res[res$lfc_comparetreatment!=0,] 

  #Create Relative Abundance Table
  df <-
    countdata_SAMPLE %>% 
    rownames_to_column('gs') %>%
    dplyr::group_by(gs) %>% 
    dplyr::summarise_all(funs(sum)) %>%
    dplyr::mutate_if(is.numeric, funs(./sum(.))) %>%
    column_to_rownames('gs')
  df_export<<-df
  #Get the ColData for Each Comparison
  coldata.1 <- coldata_SAMPLE[coldata_SAMPLE[[variable_to_compare]]==outcome_to_compare1,]
  coldata.2 <- coldata_SAMPLE[coldata_SAMPLE[[variable_to_compare]]==outcome_to_compare2,]
  #keep Count data only for each comparison
  needed<-which(colnames(df) %in% rownames(coldata.1))    
  df.1 <- df[,needed]
  needed2<-which(colnames(df) %in% rownames(coldata.2))    
  df.2 <- df[,needed2]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
  if (display_all_results_volcano=="yes" & graph_option=="volcano"){ #if want to display all results, then will replace res by resNoFilt (which does not select the top results)
    res<-resNoFilt  #replace with resNoFilt which has all taxa regardless if they are significant or not
  }
  
  #Convert Resuts table into a data.frame
  res <- as.data.frame(res)
  #decide what otu to save 
  otu.to.save <-as.character(rownames(res))
  #from relative table we should get the mean across the row of the otu table
  df.1.meanRA <- rowMeans(df.1)
  df.2.meanRA <- rowMeans(df.2)
  #need to subset AND reorder just the otus that we have 
  df.1.meanRA.save <- df.1.meanRA[otu.to.save]
  df.2.meanRA.save <- df.2.meanRA[otu.to.save]
  #add the abundnace data for the res dataframe
  res$abundance.1 <- df.1.meanRA.save
  res$abundance.2 <- df.2.meanRA.save
  #Set Names of Results Table
  res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Taxa","logFC", "lfcSE", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
  #Remove Any Data without LOGFC data
  res <- res[!is.na(res$logFC),]
  # Reorder Results based on FDR for comparison 1
  res = res[order(-res$logFC, na.last = TRUE), ]
  #Create Order
  res <- res %>% dplyr::mutate(start = 1:n())
  #Convert Important columns to Numeric
  res$adj.P.Val <-   as.numeric(as.character(res$adj.P.Val))
  res$logFC <-       as.numeric(as.character(res$logFC))
  #Replace NA
  res <- res %>% dplyr::mutate(adj.P.Val = if_else(is.na(adj.P.Val), 0.9, adj.P.Val))
  ##convert abundance to numeric
  res$abundance.1 <- as.numeric(as.character(res$abundance.1))
  res$abundance.2 <- as.numeric(as.character(res$abundance.2))
  #Create Variable for Color based on Comparison, FDR and LOGFC
  res$col <- ifelse(res$logFC>0, "B",
                    ifelse(res$logFC<0, "A","D"))
  res$abundance <- ifelse(res$col=="A", res$abundance.1, ifelse(res$col=="B", res$abundance.2,0))
  
  if (taxa_genus_output=="no") {
    res$Taxa2 <- paste0(substr(sapply(strsplit(res$Taxa,".",fixed=T),"[[",1),1,6),
                        paste0("..",paste0(substr(sapply(strsplit(res$Taxa,".",fixed=T),"[[",4),1,10),".."),
                               paste0(substr(sapply(strsplit(res$Taxa,".",fixed=T),"[[",5),1,10),paste0("..",str_sub(res$Taxa,-4,-1)))))  
  } else if (taxa_genus_output=="yes") {#keeping only genus level. if genus level is NA, then will go up a taxa rank
    res<- res %>% dplyr::mutate(Taxa2=case_when((sapply(strsplit(Taxa,".",fixed=T),"[[",1)=="NA")~"k_Bacteria.p_NA.c_NA.o_NA.f_NA.g_NA",
                                                ((sapply(strsplit(Taxa,".",fixed=T),"[[",1)!="NA") & (sapply(strsplit(Taxa,".",fixed=T),"[[",2)=="NA"))~paste0(paste0("p_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",1)),".c_NA.o_NA.f_NA.g_NA"), 
                                                ((sapply(strsplit(Taxa,".",fixed=T),"[[",2)!="NA") & (sapply(strsplit(Taxa,".",fixed=T),"[[",3)=="NA"))~paste0(paste0("c_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",2)),".o_NA.f_NA.g_NA"), 
                                                ((sapply(strsplit(Taxa,".",fixed=T),"[[",3)!="NA") & (sapply(strsplit(Taxa,".",fixed=T),"[[",4)=="NA"))~paste0(paste0("o_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",3)),".f_NA.g_NA"), 
                                                ((sapply(strsplit(Taxa,".",fixed=T),"[[",4)!="NA") & (sapply(strsplit(Taxa,".",fixed=T),"[[",5)=="NA"))~paste0(paste0("f_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",4)),".g_NA"),
                                                ((sapply(strsplit(Taxa,".",fixed=T),"[[",5)!="NA"))~paste0("g_",sapply(strsplit(res$Taxa,".",fixed=T),"[[",5)) ))
  }
  #########EXCEL EXPORT
  ####keeping all results 
  res_all <- rbind(sigDownReg,sigUpReg)
  res_all <- res_all[res_all$logFC!=0,] 

  #Convert Resuts table into a data.frame
  res_all <- as.data.frame(res_all)
  #decide what otu to save 
  otu.to.save_all <-as.character(rownames(res_all))
  #need to subset AND reorder just the otus that we have 
  df.1.meanRA.save <- df.1.meanRA[otu.to.save_all]
  df.2.meanRA.save <- df.2.meanRA[otu.to.save_all]
  #add the abundnace data for the res dataframe
  res_all$abundance.1 <- df.1.meanRA.save
  res_all$abundance.2 <- df.2.meanRA.save
  #Set Names of Results Table
  res_all <- setNames(cbind(rownames(res_all), res_all, row.names = NULL), c("Taxa","logFC", "lfcSE", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
  csv_output<-paste0(output_name,"_significant.csv")
  write.csv(res_all,file=csv_output, row.names = TRUE)
  
  #res_all contains all results including not significant results
  #Convert Resuts table into a data.frame
  res_all <- as.data.frame(resNoFilt)
  #decide what otu to save 
  otu.to.save_all <-as.character(rownames(res_all))
  #need to subset AND reorder just the otus that we have 
  df.1.meanRA.save <- df.1.meanRA[otu.to.save_all]
  df.2.meanRA.save <- df.2.meanRA[otu.to.save_all]
  #add the abundnace data for the res dataframe
  res_all$abundance.1 <- df.1.meanRA.save
  res_all$abundance.2 <- df.2.meanRA.save
  csv_output2<-paste0(output_name,"_all.csv")
  write.csv(resNoFilt_export,file=csv_output2, row.names = TRUE)
  
  ######################################################################################################################
  ######################################################################################################################
  ####merge in decontam list to change the contaminant taxa red
  if (!missing(decontam)) { #if decontam option is selected, determine if the df object exist
    if (exists(decontam)=="FALSE"){
      print(paste0("decontam list:",paste0(decontam," does not exist")))
      stop("please run function: decontaminant_subplot_KW prior to using this option- decontam") 
    } else{
      res_original<-res
      res<-merge(res, get(decontam),by.x="Taxa",by.y="row.names",All.x=all)
      res<- res[order(res$logFC, na.last = TRUE), ]
      contaminant_color <- ifelse(res$contaminant == "TRUE", "red", "black")
    }
  }  else{
    contaminant_color<-"black" ###if no decontam list was supplied, then will have all the label as black 
  }
  ######################################################################################################################
  ######################################################################################################################
  
  if (graph_option=="volcano") {
    #=========================================================
    #////////////////////VOLCANO PLOT///////////////////////
    #=========================================================
    # Compute significance, with a maximum of 350 for the p-values set to 0 due to limitation of computation precision
    res$sig <- -log10(res$adj.P.Val)
    res[is.infinite(res$sig),"sig"] <- 350
    ## Volcano plot of adjusted p-values
    cols <- densCols(res$logFC, res$sig)
    cols[res$pvalue ==0] <- "purple"
    #only significant values get display
    cols[res$logFC > 0 & res$adj.P.Val < FDR_cut_off ] <- outcome_to_compare_color2
    cols[res$logFC < 0 & res$adj.P.Val < FDR_cut_off ] <- outcome_to_compare_color1
    cols[res$adj.P.Val > FDR_cut_off ] <- "grey"
    
    ####scaling the abundance so the circle size on the graph is standardized
    res_keep<-res[!(is.na(res$adj.P.Val)),]
    res_keep<-res_keep[res_keep$adj.P.Val<FDR_cut_off,]
    min_value<-min(res_keep$abundance)
    #make the abundance of the nonsignificant results to be minimal value
    res$abundance[res$adj.P.Val > FDR_cut_off]<-min_value
    
    if (abundance_size_volcano=="yes"){ #if want point size to be scaled by relative abundance
      p<-ggplot(res, aes(x = logFC, y = sig, label=Taxa2)) +
        geom_point(aes(size =abundance),color=cols, alpha=0.5) + #Chose Colors and size for dots
        scale_size(range = c(0, 40))+
        geom_text_repel(aes(label=ifelse(res$logFC<(-1) & res$adj.P.Val < FDR_cut_off , as.character(res$Taxa2),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5,  colour = contaminant_color) +
        geom_text_repel(aes(label=ifelse(res$logFC>1 & res$adj.P.Val < FDR_cut_off , as.character(res$Taxa2),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5,  colour = contaminant_color) +
        geom_hline(yintercept=-log10(FDR_cut_off), color="red",linetype="dashed") +
        xlab("Effect size: log2(fold-change)") +
        ylab("-log10(adjusted p-value)") + 
        theme(panel.background = element_blank(),
              panel.border=element_rect(fill=NA),
              legend.background = element_rect(color=NA),
              legend.key = element_rect(colour = "transparent", fill = "white"),
              plot.title=element_text(size=23, face="bold"))
    } else {
      p<-ggplot(res, aes(x = logFC, y = sig, label=Taxa2)) +
        geom_point(color=cols, alpha=0.5) + #Chose Colors and size for dots
        geom_text_repel(aes(label=ifelse(res$logFC<(-1) & res$adj.P.Val < FDR_cut_off , as.character(res$Taxa2),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5,  colour = contaminant_color) +
        geom_text_repel(aes(label=ifelse(res$logFC>1 & res$adj.P.Val < FDR_cut_off , as.character(res$Taxa2),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5,  colour = contaminant_color) +
        geom_hline(yintercept=-log10(FDR_cut_off), color="red",linetype="dashed") +
        xlab("Effect size: log2(fold-change)") +
        ylab("-log10(adjusted p-value)") + 
        theme(panel.background = element_blank(),
              panel.border=element_rect(fill=NA),
              plot.title=element_text(size=23, face="bold"))
    }
  }
  else if (graph_option=="lollipop") {
    p<-ggplot(res, aes(y=reorder(Taxa2,-start), x=logFC,fill=col,size=abundance)) +
      geom_point(color="black",alpha=0.8,shape=21)+
      geom_segment(data=res[res$adj.P.Val<0.2,],aes(yend=reorder(Taxa2,-start)), xend=(-30), color= "black", linetype = "solid",linewidth=1)+ 
      scale_fill_manual(values=c("B"=outcome_to_compare_color2,"A"=outcome_to_compare_color1,"D"="white"))+ 
      scale_size_continuous(name="Relative Abundance",range=c(5, 20))+                
      ggtitle("ANCOMBC")+
      theme(panel.background = element_blank(),
            panel.border=element_rect(fill=NA),
            panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
            panel.grid.minor = element_blank(),
            strip.background=element_blank(),
            axis.title=element_text(size=20,face="bold"),
            axis.text.x=element_text(colour="black", size=18, face="bold"),
            axis.text.y=element_text(colour=contaminant_color,face="bold",size=10),
            axis.ticks=element_line(colour="gray70"),
            legend.background = element_rect(color=NA),
            legend.key = element_rect(colour = "transparent", fill = "white"),
            plot.title=element_text(size=23, face="bold"))+
      xlab("") +
      ylab("")+
      #xlim(-7,7)+
      geom_vline(xintercept=0, color="red",linetype="dashed")+
      guides(fill="none")
    ##save the legend on a separate png file              
  }          
  if (legend_onplot=="no") { #if want to hide the legend on the plot, then will save the legend on separate png
    leg <- get_legend(p)
    png_output2<-paste0(output_name,"_legend.png")
    png(png_output2, res = 300, width=200, height=200, units='mm')
    legend<-as_ggplot(leg)
    show(legend)
    dev.off()
    p<-p + theme(legend.position = "none") ###remove the legend
  }
  png_output<-paste0(output_name,".png")
  png(file=png_output, res = 300, width=200, height=200 , units='mm')
  show(p)
  dev.off()
  ##########################################################################################################################
  ##########################################################################################################################
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

###Non phyloseq functions
boxplot_nonPS_KW <- function(input_table, 
                             sample_type_var_name, 
                             sample_types=list(), 
                             sample_type_color=list(), 
                             sample_type_color_2nd=list(),
                             compare_stat_option=c("mean","median"),
                             log_scale=c("yes","no"),#will transform abundance to log(abundance*100 + 1)
                             graph_option=c("boxplot","mean_SE", "mean_SD"), 
                             output_prefix) {
  
  #print instruction for the function
  print("This function requires your input_table to be a dataframe with subject ID as column name")
  print("First row of the input_table is for classifying the subjects into different sample types")
  print("First column of the input_table is the names of the metabolites/transcripts")
  
  #making sure the number of item in sample_types is same as sample_type_color
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  
  #if sample_type_color_2nd is not missing then make sure it has same number of element as sample_type_color
  if (!missing(sample_type_color_2nd)) {
    stopifnot("sample_types need to have same number of elements as sample_type_color_2nd"= length(sample_types)==length(sample_type_color_2nd))
  } #if graph_option is missing, then will just have boxplot as default
  if (missing(sample_type_color_2nd)) {
    sample_type_color_2nd<-rep("black",length(sample_types)) #make the secondary color to be black by default, unless specify
  } 
  
  #adding default setting
  if (missing(graph_option)) { graph_option<-"boxplot"} #if graph_option is missing, then will just have boxplot as default
  if (missing(log_scale)) { log_scale<-"no"} #if graph_option is missing, then will just have boxplot as default
  
  #making sure the input is within allowed 
  compare_stat_option <- match.arg(compare_stat_option)
  log_scale <- match.arg(log_scale)
  graph_option <- match.arg(graph_option)
  
  #ensure the first row contains the sample type information in the dataframe
  sample_type_row<- head(input_table,1)
  sample_type<- as.character(sample_type_row[1,1])
  stopifnot("The first row needs to be your sample type information. first cell on the row needs to be the input sample type variable (sample_type_var_name)"= sample_type_var_name==sample_type)
  #checking what elements are in the first row 
  sample_type_elements<- unlist(unique(as.list(sample_type_row[,-1])))
  print("the following sample types are contain in the table:")
  print(sample_type_elements)
  input_table2<-input_table
  
  #keeping only sample types we care about 
  input_table3<-input_table2[,(input_table2[1,]) %in%  c(sample_type_var_name, sample_types)]
  
  input_table4<-input_table3[,-1]
  rownames(input_table4)<-input_table3[,1]
  
  input_table5<-input_table4[-1,]
  #changing all the columns to numeric
  input_table5 <- input_table5 %>% mutate_if(is.character, as.numeric)
  
  #sample_type_table only contains subjectID and sample type
  sample_type_table<-data.frame(t(input_table3[1,-1]))
  colnames(sample_type_table)<-"category"
  
  stopifnot("there is no sample in the data table with the sample types indicated on first row"= ncol(input_table3)!=1)
  top.stats <- as.data.frame(matrix(ncol=length(sample_types),nrow=nrow(input_table5)))
  rownames(top.stats) <- rownames(input_table5)

  for (i in 1:length(sample_types)) {
    if (compare_stat_option=="median") {
      top.stats[i] <-rowMedians(as.matrix(input_table5[, (input_table4[1,]) %in% sample_types[i]]))
    } else if (compare_stat_option=="mean") {
      top.stats[i] <-rowMeans2(as.matrix(input_table5[,  (input_table4[1,]) %in% sample_types[i]]))
    }
    colnames(top.stats)[i] <- c(paste0("stats.",sample_types[i]))
  }  
  
  #loop through the sample types and fill in the next several column with the rank order by that particular sample type. If there are 3 sample_types, then this will give rank order on 4th-6th column
  for (i in 1:length(sample_types)) {
    ##generating the rank order by sample types 
    j <- i+length(sample_types)
    column<-top.stats[,i]
    top.stats <- top.stats[order(column, decreasing = TRUE), , drop = FALSE ]
    top.stats[j] <- 1:nrow(top.stats)
    colnames(top.stats)[j] <- c(paste0("rank.order.",sample_types[i]))
  }  
  top.stats_output<-top.stats
  
  #this is to make taxa name with rank order in them
  for (i in 1:length(sample_types)) { 
    j <- i+length(sample_types)
    k <- i+2*length(sample_types)  
    column2<-top.stats[,j]
    top.stats[k] <- paste0(rownames(top.stats)," ","(",column2,")")
    colnames(top.stats)[k] <- c(paste0(paste0("rank.order.",sample_types[i]),".name"))
  }  

  for (rank_sample_type in sample_types) {
    l<-match(rank_sample_type,sample_types) ###look at the element in rank_sample_type to see what their index at sample_types
    r<-l+length(sample_types)
    # Re-rank based upon the rank group decreasing from highest to lowest
    column3<-top.stats[,r]
    top.stats <- top.stats[order(column3, decreasing = F), , drop = FALSE ]
    top.2 <- top.stats
    top.3 <- head(top.2, 100)
    # Filtering ONLY the top 100 taxa 
    subset.top <- input_table5 %>% 
      dplyr::filter(rownames(input_table5) %in% rownames(top.3))
    # Re ordered and saved into a new group (previously filtered)
    # subset.top.match now have only 100 taxa ordered by top background 
    subset.top.match <- subset.top[match(rownames(top.3),rownames(subset.top)),]    
    
    for (a in sample_types) {
      assign(paste0("rank_order_",a),subset(top.3, select=c(paste0("rank.order.",a), paste0(paste0("rank.order.",a),".name") )))
    }      
    
    # Transposition of the figure, maybe at this point replace the column names 
    # The previous edit should have removed all character columns 
    subset.top.match.transpose <- as.data.frame(t(subset.top.match))

    # name a new column with categories 
    subset.top.match.transpose<-merge(subset.top.match.transpose,sample_type_table, by="row.names")
    row.names(subset.top.match.transpose)<-subset.top.match.transpose$Row.names

    # Reverse column order 
    subset.top.match.transpose.reverse <- subset.top.match.transpose[,order(ncol(subset.top.match.transpose):1)]
    subset.top.match.transpose.reverse_final<<-subset.top.match.transpose.reverse
    
    for (a in sample_types) {
      assign(paste0("subset.top.match.transpose.reverse.",a),subset(subset.top.match.transpose.reverse, category == c(a) ))
      assign(paste0("subset.top.match.transpose.reverse.",a),subset(get(paste0("subset.top.match.transpose.reverse.",a)), select=-c(category)))
    }  
    
    #get the first plot- rank_sample_type is the sample type that is being ranked.... we want the ranked plot on the left 
    temp <- get(paste0("subset.top.match.transpose.reverse.",rank_sample_type)) %>% dplyr::select(everything()) %>% tidyr::gather("id", "value",1:100, factor_key = TRUE)
    temp <- merge(temp, top.3, by.x="id", by.y="row.names")
    #if log_scale is yes then will have log(abundance+1)
    if (tolower(log_scale)=="yes"){
      temp$value2<-log10(temp$value*100+1) 
    } else {
      temp$value2<-temp$value
    }
    p <- temp %>%   mutate(id = fct_reorder(id, -get(paste0("rank.order.",rank_sample_type)))) %>%
      ggplot(., aes(x=id, y=value2)) + 
      coord_flip() +
      theme_bw() + 
      theme(axis.text.y=element_text(),axis.title.y=element_blank()) +
      theme(axis.text.x=element_text(angle=90,hjust=1))+
      ggtitle(rank_sample_type) 
    
    if (tolower(log_scale)=="yes"){
      p<-p+ylab("Log 10 (Intensity)")
    } else {
      p<-p+ylab("Intensity")
    }
    if (graph_option=="boxplot") {
      p<-p+geom_boxplot(color=sample_type_color[[l]], alpha=0.2) 
    } else if(graph_option=="mean_SD") {
      p<-p+geom_jitter(color=sample_type_color[[l]],position=position_jitter(0), alpha=0.2) +
        stat_summary(fun= mean, 
                     geom="pointrange", 
                     fun.max = function (x) mean(x)+sd(x),
                     fun.min = function (x) ifelse( mean(x)-sd(x) < 0, 0, mean(x)-sd(x)),
                     color=sample_type_color_2nd[[l]], linewidth =1.0, size=0.4)
    } else if(graph_option=="mean_SE") {
      p<-p+geom_jitter(color=sample_type_color[[l]],position=position_jitter(0), alpha=0.2) +
        stat_summary(fun= mean, 
                     geom="pointrange", 
                     fun.max = function (x) mean(x)+sd(x)/sqrt(length(x)),
                     fun.min = function (x) ifelse( mean(x)-sd(x)/sqrt(length(x)) < 0,0, mean(x)-sd(x)/sqrt(length(x))),
                     color=sample_type_color_2nd[[l]], linewidth =1.0, size=0.4)
    }
    myplots <- list()  # new empty list
    myplots[[1]]<- p
    
    #remove the sample type that is being ranked. so new_sample_type_name and new_sample_type_color contain the sample types which were not ranked                
    new_sample_type_name<-sample_types[-l] #remove the lth element (lth element is the sample type that is being ranked by )
    new_sample_type_color<-sample_type_color[-l] #remove the lth element (lth element is the sample type that is being ranked by )
    new_sample_type_color_2nd<-sample_type_color_2nd[-l] #remove the lth element (lth element is the sample type that is being ranked by )
    
    #remove the axis label for the sample types which were not ranked (since the ranked sample type would be on the left most)  
    remove_y<- theme(axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.y = element_blank())
    
    ###looping thro each sample type 
    number<-1
    number2<-2
    for (a in new_sample_type_name) {
      temp <- get(paste0("subset.top.match.transpose.reverse.",a)) %>% dplyr::select(everything()) %>% tidyr::gather("id", "value",1:100, factor_key = TRUE)
      temp <- merge(temp, top.3, by.x="id", by.y="row.names")
      #if log_scale is yes then will have log(abundance+1)
      if (tolower(log_scale)=="yes"){
        temp$value2<-log10(temp$value*100+1) 
      } else {
        temp$value2<-temp$value
      }
      p <- temp %>%   mutate(id = fct_reorder(id, -get(paste0("rank.order.",rank_sample_type)))) %>%
        ggplot(., aes(x=id, y=value2)) + 
        coord_flip() +
        theme_bw() + 
        theme(axis.text.y=element_text(color=rev(top.3$color)),axis.title.y=element_blank()) +
        theme(axis.text.x=element_text(angle=90,hjust=1))+
        ggtitle(a) 
      
      #if log_scale is yes then will have log(abundance+1)
      if (tolower(log_scale)=="yes"){
        p<-p+ylab("Log 10 (Intensity)")
      } else {
        p<-p+ylab("Intensity")
      }
      if (graph_option=="boxplot") {
        p<-p+geom_boxplot(color=new_sample_type_color[[number]], alpha=0.2) 
      } else if(graph_option=="mean_SD") {
        p<-p+geom_jitter(color=new_sample_type_color[[number]],position=position_jitter(0), alpha=0.2) +
          stat_summary(fun= mean, 
                       geom="pointrange", 
                       fun.max = function (x) mean(x)+sd(x),
                       fun.min = function (x) ifelse( mean(x)-sd(x) < 0, 0, mean(x)-sd(x)),
                       color=new_sample_type_color_2nd[[number]], linewidth =1.0, size=0.4)
      } else if(graph_option=="mean_SE") {
        p<-p+geom_jitter(color=new_sample_type_color[[number]],position=position_jitter(0), alpha=0.2) +
          stat_summary(fun= mean, 
                       geom="pointrange", 
                       fun.max = function (x) mean(x)+sd(x)/sqrt(length(x)),
                       fun.min = function (x) ifelse( mean(x)-sd(x)/sqrt(length(x)) < 0,0, mean(x)-sd(x)/sqrt(length(x))),
                       color=new_sample_type_color_2nd[[number]], linewidth =1.0, size=0.4)
      }
      myplots[[number2]] <- p + remove_y
      number<-number+1
      number2<-number2+1
    }
    
    ####export csv
    csv_output=paste0(paste0(paste0(paste0(output_prefix,"boxplot_rank_by_"),rank_sample_type),paste0("_",compare_stat_option)),".csv")
    write.csv(top.stats_output,csv_output)
    
    pdf_output_boxplot<-paste0(paste0(paste0(paste0(output_prefix,"boxplot_rank_by_"),rank_sample_type),paste0("_",compare_stat_option)),".pdf")
    pdf(file=pdf_output_boxplot, width=14+(length(sample_types)-3)*4, height=(14+(length(sample_types)-3)*4)*10/8.5)
           show(patchwork::wrap_plots(myplots, nrow=1))
    dev.off()        
  }   
  
}

alpha_diversity_nonPS_KW <- function(input_table, 
                                     sample_type_var_name, 
                                     sample_types=list(), 
                                     sample_type_color=list(), 
                                     p_value=c("yes","no"), 
                                     output_name, 
                                     width=100, height=150, 
                                     xlabel_size=14, ylabel_size=14, axis_title_size=16,
                                     xlabel, ylabel,                          
                                     pair_by_variable=NULL, ...) {
  list2env(list(...), environment())
  #print instruction for the function
  print("This function requires your input_table to be a dataframe with subject ID as column name")
  print("First row of the input_table is for classifying the subjects into different sample types")
  print("First column of the input_table is the names of the metabolites/transcripts")
  
  #making sure the number of item in sample_types is same as sample_type_color
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  #making sure only yes and no is selected for p_value
  p_value <- match.arg(p_value)
  
  #making sure the coloring and order of the input for the sample_types are going to align with the output graph
  sample_type_color_t<-sample_type_color
  names(sample_type_color_t) <-sample_types
  sample_type2<-sort(sample_types)
  
  #ensure the first row contains the sample type information in the dataframe
  sample_type_row<- head(input_table,1)
  sample_type<- as.character(sample_type_row[1,1])
  stopifnot("The first row needs to be your sample type information. first cell on the row needs to be the input sample type variable (sample_type_var_name)"= sample_type_var_name==sample_type)
  #checking what elements are in the first row 
  sample_type_elements<- unlist(unique(as.list(sample_type_row[,-1])))
  print("the following sample types are contain in the table:")
  print(sample_type_elements)
  input_table2<-input_table
  
  #keeping only sample types we care about 
  input_table3<-input_table2[,(input_table2[1,]) %in%  c(sample_type_var_name, sample_types)]
  
  if (!missing(pair_by_variable)){ 
    stopifnot("you can only have two groups for paired analysis"= length(sample_types)==2)
    stopifnot("Please make sure the second row of the dataframe is the variable used for pairing"= input_table[2,1]==pair_by_variable)
  }
  
  #get the first column to be row names
  input_table4<-input_table3[,-1]
  rownames(input_table4)<-input_table3[,1]
  
  if (!missing(pair_by_variable)){ 
    #sample_type_table only contains subjectID, sample type, pair_by_varible
    sample_type_table<-data.frame(t(input_table3[c(1,2),-1])) 
    colnames(sample_type_table)<-c(sample_type_var_name,pair_by_variable) 
    input_table4<-input_table4[-c(1,2),]
  } else{
    #sample_type_table only contains subjectID and sample type
    sample_type_table<-data.frame(t(input_table3[1,-1]))
    colnames(sample_type_table)<-sample_type_var_name
    input_table4<-input_table4[-1,]
  }
  stopifnot("there is no sample in the data table with the sample types indicated on first row"= ncol(input_table3)!=1)
  #removing the first column which contains the name of the metabolite/transcripts and first row which contains sample type information
  #changing all the columns to numeric
  input_table4 <- input_table4 %>% mutate_if(is.character, as.numeric)
  Shannon_diversity <- data.frame(vegan::diversity(input_table4, index = "shannon", MARGIN = 2, base = exp(1)))
  colnames(Shannon_diversity)[1]<-"Shannon_diversity"
  
  if (missing(xlabel)) { #if xlabel is missing then will have default  
    xlabel<-NULL
  } 
  if (missing(ylabel)) { #if ylabel is missing then will have default  
    ylabel<-"Shannon diversity index"
  }
  alpha_dataframe<-sample_type_table
  alpha_dataframe<- merge(alpha_dataframe, Shannon_diversity, by="row.names")
  
  input_dataframe<- alpha_dataframe %>% dplyr::select(c(sample_type_var_name, "Shannon_diversity",pair_by_variable))# if pair_by_variable is not provided, then input_dataframe only have Shannon_diversity and sample_type_var_name
  
  #standardized variables 
  input_dataframe$var1<-input_dataframe[[sample_type_var_name]] #standardizing the variable name. var1 would become the sample type separator
  if (!missing(pair_by_variable)){ 
    input_dataframe$group_var<-input_dataframe[[pair_by_variable]]
    input_dataframe<-input_dataframe[order(input_dataframe$var1, input_dataframe$group_var),]    
  }
  #turn into a factor with the same order that the function input has
  input_dataframe[["var1"]] <- factor(input_dataframe[["var1"]], sample_types)
  
  #Setting up the graphs
  p<-ggplot(input_dataframe, aes(var1, Shannon_diversity, fill=var1)) +
    stat_boxplot(geom = "errorbar",width = 0.3, color = "gray70")+
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values=sample_type_color_t) +
    scale_x_discrete(labels=sample_types) +
    scale_y_continuous() + 
    theme_classic() +
    labs(title=NULL,
         x = xlabel, 
         y= ylabel) +
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          axis.title=element_text(face="bold", size=axis_title_size),
          axis.text.x =element_text(size=xlabel_size),
          axis.text.y =element_text(size=ylabel_size),
          legend.position="none")
  if (!missing(pair_by_variable)){
    p <- p + geom_line(aes(group = group_var), color="gray70", alpha=0.5) +
      geom_jitter(width=0, size=1, shape= 1, color="azure4") #if there is paired line, will have jitter to be 0 so is not as messy
  } else { #if paired comparison is not needed then will not have lines connect betwen groups
    p <- p + geom_jitter(width=0.3, size=1, shape= 1, color="azure4") #if no paired line, then will have jitter 
  }
  if (p_value=="yes") { ##if want p_value included in the diagram
    if (!missing(pair_by_variable)){
      stat.test <- input_dataframe %>%
        wilcox_test(Shannon_diversity ~ var1, paired=T) %>% add_significance()
    }else { # if no paired comparison 
      stat.test <- input_dataframe %>%
        wilcox_test(Shannon_diversity ~ var1) %>% add_significance()
    }
    stat.test<- stat.test %>% add_xy_position(x = "var1", y.trans = function(x){x*1.1})
    if (!exists("step.increase")) {step.increase<-0.05} 
    #limiting p value to 3 digits after decimal
    stat.test$p<-pvalr(stat.test$p, digits = 3)
    p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.02, step.increase = step.increase,  inherit.aes = F, size=3)
  }
  ###determine sample_types string length
  string_length <- max(nchar(as.character(input_dataframe$var1)))
  if (string_length>20) { #if the sample type string is too long then will title the x axis text 
    p<-p+ theme(axis.text.x = element_text(angle = 8, hjust = 0.5, vjust=0.5))
  }
  png_output<-paste0(output_name,".png")
  png(file=png_output, res = 300, width=width, height=height, units='mm')
  show(p)
  dev.off()
  csv_output<-paste0(output_name,".csv")
  write.csv(Shannon_diversity, csv_output)
}

beta_diversity_nonPS_KW <- function(input_table, 
                                    sample_type_var_name, 
                                    sample_types=list(), 
                                    sample_type_color=list(), 
                                    p_value=c("yes","no"),
                                    output_name, 
                                    width=200, height=200,
                                    plot_title_size=14, axis_title_size=18, label_size=7,
                                    x_axis_flip, 
                                    p_value_location=c("TR","TL","BR","BL")) { #p_value_location: TR-top right; TL- top left; BR- bottom right; BL- bottom left
  # Enforce orientation. Samples are columns
  if(missing(x_axis_flip)) { x_axis_flip<-"no"} #default being no for x_axis_flip
  if(missing(p_value_location) & tolower(p_value)=="yes"){
    p_value_location<-"BL" #default location for p_value will be in the bottom left corner
  } else if(!missing(p_value_location) & tolower(p_value)=="yes"){ #if p_value is wanted and p_value_location is not missing
    print("p_value_location: TR-top right; TL- top left; BR- bottom right; BL- bottom left")
    p_value_location<-match.arg(p_value_location)
    print(paste0("selected:",p_value_location))
  }
  #print instruction for the function
  print("This function requires your input_table to be a dataframe with subject ID as column name")
  print("First row of the input_table is for classifying the subjects into different sample types")
  print("First column of the input_table is the names of the metabolites/transcripts")
  
  #making sure the number of item in sample_types is same as sample_type_color
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  #making sure only yes and no is selected for p_value
  p_value <- match.arg(p_value)
  
  #making sure the coloring and order of the input for the sample_types are going to align with the output graph
  sample_type_color_t<-sample_type_color
  names(sample_type_color_t) <-sample_types
  sample_type2<-sort(sample_types)
  
  #ensure the first row contains the sample type information in the dataframe
  sample_type_row<- head(input_table,1)
  sample_type<- as.character(sample_type_row[1,1])
  stopifnot("The first row needs to be your sample type information. first cell on the row needs to be the input sample type variable (sample_type_var_name)"= sample_type_var_name==sample_type)
  #checking what elements are in the first row 
  sample_type_elements<- unlist(unique(as.list(sample_type_row[,-1])))
  print("the following sample types are contain in the table:")
  print(sample_type_elements)
  input_table2<-input_table
  
  #keeping only sample types we care about 
  input_table3<-input_table2[,(input_table2[1,]) %in%  c(sample_type_var_name, sample_types)]
  
  #sample_type_table only contains subjectID and sample type
  sample_type_table<-data.frame(t(input_table3[1,-1]))
  colnames(sample_type_table)<-sample_type_var_name
  
  stopifnot("there is no sample in the data table with the sample types indicated on first row"= ncol(input_table3)!=1)
  #removing the first column which contains the name of the metabolite/transcripts and first row which contains sample type information
  input_table4<-input_table3[,-1]
  rownames(input_table4)<-input_table3[,1]
  input_table4<-input_table4[-1,]
  
  #changing all the columns to numeric
  input_table4 <- input_table4 %>% mutate_if(is.character, as.numeric)
  
  #Beta diversity
  #Create Distance Matrix with Bray
  vegdist=vegdist(t(input_table4), method = "bray") #row is samplesID #column is now metabolite/transcript
  #Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
  CmdScale <- cmdscale(vegdist, k =10)
  #calculated Sample variance for each PC
  vars <- apply(CmdScale, 2, var)
  #Create Variable with the Percent Variance
  percentVar <- round(100 * (vars/sum(vars)))
  
  #Merge PC Data with MetaData
  newResults <- merge(x = CmdScale, y = sample_type_table, by = "row.names", all.x = TRUE)  #sample_data(phyloseq_t) row is the sampleID
  
  #Rename Variables for PC1 and PC2
  colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
  colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
  colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"
  newResults_backup<-newResults
  
  formula <- as.formula(paste("cbind(PC1,PC2)",sample_type_var_name, sep=" ~ "))
  #Calculate the Centroid Value
  centroids <- aggregate(formula ,data= newResults, mean)
  #Merge the Centroid Data into the PCOA Data
  newResults <- merge(newResults,centroids,by=sample_type_var_name,suffixes=c("",".centroid"))
  
  # Calculate p-value with ADONIS
  x=adonis2(vegdist ~ sample_type_table[[sample_type_var_name]])
  pvalues<-x$`Pr(>F)`[[1]]
  subtitle_output<- paste0("Adonis, p=",pvalues)
  
  p<-ggplot(newResults, aes(PC1, PC2, color= get(sample_type_var_name))) + # Graph PC1 and PC2
    geom_point(size=2) + # Set the size of the points
    xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2 
    labs(title=NULL, 
         fill="Centroids") +
    geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= get(sample_type_var_name)))+ 
    geom_point(data=centroids, aes(x=PC1, y=PC2, color= get(sample_type_var_name)), size=0) + 
    geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=sample_type2), size=label_size) +
    scale_color_manual(name=NULL,  
                       values=sample_type_color_t) +
    theme(panel.background = element_blank(),
          panel.border=element_rect(fill=NA),
          panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
          panel.grid.minor = element_blank(),strip.background=element_blank(),
          plot.title=element_text(face="bold",hjust=0.5, size = plot_title_size), 
          plot.subtitle = element_text(hjust=0.5),
          axis.title=element_text(face="bold", size = axis_title_size),
          axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
          axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
          axis.ticks=element_blank(),
          plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
  p<-p
  newResults<<-newResults
  if (tolower(x_axis_flip)=="yes") {
    p<-p+scale_x_reverse()
  }
  if (p_value=="yes") { ###when the axis is flip, the annotation for the p value need to placed on max(newResults$PC1), instead of min(newResults$PC2))
    if (tolower(x_axis_flip)=="no"){
      if (p_value_location=="BR"){
        p <- p + annotate("text", x=max(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
      } else if (p_value_location=="BL"){
        p <- p + annotate("text", x=min(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
      } else if (p_value_location=="TR"){
        p <- p + annotate("text", x=max(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
      } else if (p_value_location=="TL"){
        p <- p + annotate("text", x=min(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
      }
    } else if (tolower(x_axis_flip)=="yes") {
      if (p_value_location=="BR"){
        p <- p + annotate("text", x=min(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
      } else if (p_value_location=="BL"){
        p <- p + annotate("text", x=max(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
      } else if (p_value_location=="TR"){
        p <- p + annotate("text", x=min(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
      } else if (p_value_location=="TL"){
        p <- p + annotate("text", x=max(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
      }
    }
  }
  
  png_output<-paste0(output_name,".png")
  png(file=png_output, res = 300, width=width, height=height , units='mm')
  show(p)
  dev.off()
  
}

PLS_DA_nonPS_KW <- function(input_table, 
                      sample_type_var_name, 
                      sample_types=list(), 
                      sample_type_color=list(), 
                      top_display_loading=20, 
                      output_name, 
                      width=500, height=500) {
  
  #print instruction for the function
  print("This function requires your input_table to be a dataframe with subject ID as column name")
  print("First row of the input_table is for classifying the subjects into different sample types")
  print("First column of the input_table is the names of the metabolites/transcripts")
  
  #making sure the number of item in sample_types is same as sample_type_color
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))

  #setting up color
  sample_type_color_t<-sample_type_color
  names(sample_type_color_t) <-sample_types  
  sample_type_t<-sample_types
  names(sample_type_t)<-sample_type_color
  sample_type_t<-sort(sample_type_t)
  sample_type_color_t2<-names(sample_type_t)

  #ensure the first row contains the sample type information in the dataframe
  sample_type_row<- head(input_table,1)
  sample_type<- as.character(sample_type_row[1,1])
  stopifnot("The first row needs to be your sample type information. first cell on the row needs to be the input sample type variable (sample_type_var_name)"= sample_type_var_name==sample_type)
  #checking what elements are in the first row 
  sample_type_elements<- unlist(unique(as.list(sample_type_row[,-1])))
  print("the following sample types are contain in the table:")
  print(sample_type_elements)
  input_table2<-input_table
  
  #keeping only sample types we care about 
  input_table3<-input_table2[,(input_table2[1,]) %in%  c(sample_type_var_name, sample_types)]
  
  input_table4<-input_table3[,-1]
  rownames(input_table4)<-input_table3[,1]
  
  input_table5<-input_table4[-1,]
  #changing all the columns to numeric
  input_table5 <- input_table5 %>% mutate_if(is.character, as.numeric)
  input_table6<-t(input_table5)
  
  #sample_type_table only contains subjectID and sample type
  sample_type_table<-data.frame(t(input_table3[1,-1]))
  colnames(sample_type_table)<-sample_type_var_name
  
  sample_type_table[[1]]<-as.factor(sample_type_table[[1]])

  X<-input_table6
  Y<-sample_type_table[[1]]
  col_name<-paste0(paste0(substr(colnames(X), 1, 10), ".."),str_sub(colnames(X),-7,-1)) #limiting the name to 20 characters
  
  MyResult.splsda <- mixOmics::splsda(X, Y, ncomp = 5) # 1 Run the method
  p3<-plotIndiv(MyResult.splsda, comp =c(1,2), col.per.group = sample_type_color_t2,
                group = Y, ind.names = FALSE,  # colour points by class
                ellipse = TRUE, # include 95% confidence ellipse for each class
                legend = TRUE, title = '(a) PLSDA with confidence ellipses',
  )
  png_output1<- paste0(output_name,"_a.png")
  png(png_output1, res = 300, width=width, height=height , units='mm') #PLSDA with confidence ellipses
          show(p3)
  dev.off()

  p4<-plotLoadings(MyResult.splsda, 
                   contrib = 'max', method = 'mean',
                   comp=1, 
                   ndisplay = top_display_loading, name.var=col_name) ##display top components
  p <- ggplot(p4, aes(x = reorder(rownames(p4), -importance), y = importance))+
    geom_col(aes(fill=GroupContrib), width=0.6) + 
    scale_fill_manual("sample type",values=sample_type_color_t) +
    #########################
  coord_flip()+
    theme(panel.background = element_blank(),
          panel.border=element_rect(fill=NA),
          panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
          panel.grid.minor = element_blank(),strip.background=element_blank(),
          plot.title=element_text(face="bold",hjust=0.5, size = 14), 
          plot.subtitle = element_text(hjust=0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x=element_text(size = 14),
          axis.text.y=element_text(size = 22),
          axis.ticks=element_blank(),
          plot.margin=unit(c(1,1,1,1),"line"),
          legend.title = element_text(size=14), 
          legend.text = element_text(size=14))
  png_output3<- paste0(output_name,"_b.png")
  png(png_output3, res = 300, width=width*3/4, height=height , units='mm') #Plot of Loading vectors
        show(p)
  dev.off()
}

Edger_nonPS_KW <- function (input_table, 
                            sample_type_var_name, 
                            sample_types=list(), 
                            sample_type_color=list(),                             
                            FDR_cut_off=0.2,
                            number_display=10,
                            graph_option=c("volcano","lollipop"),
                            display_all_results_volcano=c("yes","no"), #display all the results in volcano plot (including results that is not significant: FDR not below the cut off)
                            abundance_size_volcano=c("yes","no"), #change size of the dot based on relative abundance
                            output_name,
                            legend_onplot=c("yes","no")) {
  
  #print instruction for the function
  print("This function requires your input_table to be a dataframe with subject ID as column name")
  print("First row of the input_table is for classifying the subjects into different sample types")
  print("First column of the input_table is the names of the metabolites/transcripts")
  
  #making sure the number of item in sample_types is same as sample_type_color
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  
  if (missing(graph_option)) { graph_option<-"lollipop"} #default being lollipop for graph_option
  if (missing(display_all_results_volcano)) { display_all_results_volcano<-"no"} #default being lollipop for graph_option
  if (missing(abundance_size_volcano)) { abundance_size_volcano<-"yes"} #default being lollipop for graph_option
  
  #setting up color
  sample_type_color_t<-sample_type_color
  names(sample_type_color_t) <-sample_types  
  sample_type_t<-sample_types
  names(sample_type_t)<-sample_type_color
  sample_type_t<-sort(sample_type_t)
  sample_type_color_t2<-names(sample_type_t)
  
  for (i in 1:length(sample_types)) {
    assign(paste0("sample_types",i), sample_types[[i]])
    assign(paste0("sample_type_color",i), sample_type_color[[i]])
  }  
  
  #ensure the first row contains the sample type information in the dataframe
  sample_type_row<- head(input_table,1)
  sample_type<- as.character(sample_type_row[1,1])
  stopifnot("The first row needs to be your sample type information. first cell on the row needs to be the input sample type variable (sample_type_var_name)"= sample_type_var_name==sample_type)
  #checking what elements are in the first row 
  sample_type_elements<- unlist(unique(as.list(sample_type_row[,-1])))
  print("the following sample types are contain in the table:")
  print(sample_type_elements)
  input_table2<-input_table
  
  #keeping only sample types we care about 
  input_table3<-input_table2[,(input_table2[1,]) %in%  c(sample_type_var_name, sample_types)]
  
  input_table4<-input_table3[,-1]
  rownames(input_table4)<-input_table3[,1]
  
  input_table5<-input_table4[-1,]
  #changing all the columns to numeric
  input_table5 <- input_table5 %>% mutate_if(is.character, as.numeric)
  input_table6<-t(input_table5)
  
  #sample_type_table only contains subjectID and sample type
  sample_type_table<-data.frame(t(input_table3[1,-1]))
  colnames(sample_type_table)<-sample_type_var_name
  
  ## making sure function input is appropiate
  legend_onplot <- match.arg(legend_onplot)
  if (missing(legend_onplot)) { legend_onplot<-"yes"} #default being yes for legend_onplot, which will just display legend on the graph
  
  #Set order
  coldata_SAMPLE <- data.frame(sample_type_table[order(row.names(sample_type_table)),,drop=FALSE],check.names=F)
  countdata_SAMPLE <- data.frame(input_table5[, order(colnames(input_table5)),drop=FALSE],check.names=F)
  
  #making a compare variable to standardize the input variable for comparsion- make sure the reference is always control
  coldata_SAMPLE$sample_type_var_name<-coldata_SAMPLE[[sample_type_var_name]]
  coldata_SAMPLE<- coldata_SAMPLE %>% dplyr::mutate(compare=case_when((sample_type_var_name==sample_types1)~"control", (sample_type_var_name==sample_types2)~"treatment"))
  coldata_SAMPLE<-coldata_SAMPLE[!is.na(coldata_SAMPLE$compare),] ###dropping subjects if they are not in the group of comparsion
  coldata_SAMPLE$compare<- as.factor(coldata_SAMPLE$compare) #converting variable to factor
  countdata_SAMPLE<-countdata_SAMPLE[,colnames(countdata_SAMPLE) %in% rownames(coldata_SAMPLE)]#making sure coldata has same subjects and countdata
  
  #Run EDGER
  dgeFull <- DGEList(countdata_SAMPLE, group=coldata_SAMPLE$compare)
  dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
                     group=dgeFull$samples$group)
  dgeFull <- calcNormFactors(dgeFull, method="TMM")
  dgeFull <- estimateCommonDisp(dgeFull)
  dgeFull <- estimateTagwiseDisp(dgeFull)
  dgeTest <- exactTest(dgeFull)
  
  #Create Table
  resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")
  resNoFilt_export<<-resNoFilt
  resNoFilt <- data.frame(resNoFilt)
  #check to see if there is any taxa with FDR<FDR_cut_off. If not then will stop the function, because there will be error if it continues as the object will have no rows
  checkifallzero<-resNoFilt[resNoFilt$FDR<FDR_cut_off,]
  #export all the results without filtering the nonsignificant results
  csv_output2<-paste0(output_name,"_all.csv")
  write.csv(resNoFilt_export,file=csv_output2, row.names = TRUE)
  if (dim(checkifallzero)[1] == 0) {
    print(paste0("no Taxa with FDR more than",FDR_cut_off))
    return(NULL)
  }
  resNoFilt[nrow(resNoFilt) + 1,] <- c(0, 0, 0, 0) ###adding this row of 0 so that even if there is no significant top 10 for upregulated or downregulated, the code will still run
  #Create Table of Down Regulated
  sigDownReg <- resNoFilt[resNoFilt$FDR<FDR_cut_off,] #keeping only FDR < 0.2
  #sigDownReg <- resNoFilt
  sigDownReg <- sigDownReg[sigDownReg$logFC<=0,] ##keeping only the negative logFC
  sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
  sigDownReg1 <- sigDownReg %>% dplyr::slice(1:number_display) #Subset Top (number_display)
  #Create Table of Down Regulated
  sigUpReg <- resNoFilt[resNoFilt$FDR<FDR_cut_off,] #keeping only FDR < 0.2
  #sigUpReg <- resNoFilt
  sigUpReg <- sigUpReg[sigUpReg$logFC>=0,] ##keeping only the positive logFC
  sigUpReg <- sigUpReg[order(sigUpReg$logFC, decreasing=TRUE),]
  sigUpReg1 <- sigUpReg %>% dplyr::slice(1:number_display) #Subset Top (number_display)
  #Merge Top N results
  #res contains the top N results we care about
  res <- rbind(sigDownReg1,sigUpReg1)
  res <- res[res$logFC!=0,] 
  #remove the 0 column which was used as a place holder
  resNoFilt <- resNoFilt[resNoFilt$logFC!=0,] 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #Create Relative Abundance Table: which will be merge back to the results
  df <-
    countdata_SAMPLE %>% 
    rownames_to_column('gs') %>%
    dplyr::group_by(gs) %>% 
    dplyr::summarise_all(funs(sum)) %>%
    dplyr::mutate_if(is.numeric, funs(./sum(.))) %>%
    column_to_rownames('gs')
  df_export<<-df
  #Get the ColData for Each Comparison
  coldata.1 <- coldata_SAMPLE[coldata_SAMPLE[[sample_type_var_name]]==sample_types1,]
  coldata.2 <- coldata_SAMPLE[coldata_SAMPLE[[sample_type_var_name]]==sample_types2,]
  #keep Count data only for each comparison
  needed<-which(colnames(df) %in% rownames(coldata.1))    
  df.1 <- df[,needed]
  needed2<-which(colnames(df) %in% rownames(coldata.2))    
  df.2 <- df[,needed2]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
  if (display_all_results_volcano=="yes" & graph_option=="volcano"){ #if want to display all results, then will replace res by resNoFilt (which does not select the top results)
          res<-resNoFilt
  }
  
  #Convert Resuts table into a data.frame
  res <- as.data.frame(res)
  #decide what otu to save 
  otu.to.save <-as.character(rownames(res))
  #from relative table we should get the mean across the row of the otu table
  df.1.meanRA <- rowMeans(df.1)
  df.2.meanRA <- rowMeans(df.2)
  #need to subset AND reorder just the otus that we have 
  df.1.meanRA.save <- df.1.meanRA[otu.to.save]
  df.2.meanRA.save <- df.2.meanRA[otu.to.save]
  #add the abundnace data for the res dataframe
  res$abundance.1 <- df.1.meanRA.save
  res$abundance.2 <- df.2.meanRA.save
  #Set Names of Results Table
  res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Name","logFC", "lfcSE", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
  #Remove Any Data without LOGFC data
  res <- res[!is.na(res$logFC),]
  # Reorder Results based on FDR for comparison 1
  res = res[order(-res$logFC, na.last = TRUE), ]
  #Create Order
  res <- res %>% dplyr::mutate(start = 1:n())
  #Convert Important columns to Numeric
  res$adj.P.Val <-   as.numeric(as.character(res$adj.P.Val))
  res$logFC <-       as.numeric(as.character(res$logFC))
  #Replace NA
  res <- res %>% dplyr::mutate(adj.P.Val = if_else(is.na(adj.P.Val), 0.9, adj.P.Val))
  ##convert abundance to numeric
  res$abundance.1 <- as.numeric(as.character(res$abundance.1))
  res$abundance.2 <- as.numeric(as.character(res$abundance.2))
  #Create Variable for Color based on Comparison, FDR and LOGFC
  res$col <- ifelse(res$logFC>0, "B",
                    ifelse(res$logFC<0, "A","D"))
  res$abundance <- ifelse(res$col=="A", res$abundance.1, ifelse(res$col=="B", res$abundance.2,0))
  
  #########EXCEL EXPORT
  #res_sign contains all the significant results
  res_sign <- rbind(sigDownReg,sigUpReg)
  res_sign <- res_sign[res_sign$logFC!=0,] 
  #Convert Resuts table into a data.frame
  res_sign <- as.data.frame(res_sign)
  #decide what otu to save 
  otu.to.save_sign <-as.character(rownames(res_sign))
  #need to subset AND reorder just the otus that we have 
  df.1.meanRA.save <- df.1.meanRA[otu.to.save_sign]
  df.2.meanRA.save <- df.2.meanRA[otu.to.save_sign]
  #add the abundnace data for the res dataframe
  res_sign$abundance.1 <- df.1.meanRA.save
  res_sign$abundance.2 <- df.2.meanRA.save
  #Export result table
  csv_output<-paste0(output_name,"_significant.csv")
  write.csv(res_sign,file=csv_output, row.names = TRUE)

  #res_all contains all results including not significant results
  #Convert Resuts table into a data.frame
  res_all <- as.data.frame(resNoFilt)
  #decide what otu to save 
  otu.to.save_all <-as.character(rownames(res_all))
  #need to subset AND reorder just the otus that we have 
  df.1.meanRA.save <- df.1.meanRA[otu.to.save_all]
  df.2.meanRA.save <- df.2.meanRA[otu.to.save_all]
  #add the abundnace data for the res dataframe
  res_all$abundance.1 <- df.1.meanRA.save
  res_all$abundance.2 <- df.2.meanRA.save
  #Export result table
  csv_output2<-paste0(output_name,"_all.csv")
  write.csv(res_all,file=csv_output2, row.names = TRUE)
  
  if (graph_option=="volcano") {
    #=========================================================
    #////////////////////VOLCANO PLOT///////////////////////
    #=========================================================
    # Compute significance, with a maximum of 350 for the p-values set to 0 due to limitation of computation precision
    res$sig <- -log10(res$adj.P.Val)
    res[is.infinite(res$sig),"sig"] <- 350
    ## Volcano plot of adjusted p-values
    cols <- densCols(res$logFC, res$sig)
    cols[res$pvalue ==0] <- "purple"
    #only significant values get display
    cols[res$logFC > 0 & res$adj.P.Val < FDR_cut_off ] <- sample_type_color2
    cols[res$logFC < 0 & res$adj.P.Val < FDR_cut_off ] <- sample_type_color1
    cols[res$adj.P.Val > FDR_cut_off ] <- "grey"
    
    ####scaling the abundance so the circle size on the graph is standardized
    res_keep<-res[!(is.na(res$adj.P.Val)),]
    res_keep<-res_keep[res_keep$adj.P.Val<FDR_cut_off,]
    min_value<-min(res_keep$abundance)
    #make the abundance of the nonsignificant results to be minimal value
    res$abundance[res$adj.P.Val > FDR_cut_off]<-min_value
    if (abundance_size_volcano=="yes"){ #if want point size to be scaled by relative abundance
        p<-ggplot(res, aes(x = logFC, y = sig, label=Name)) +
          geom_point(aes(size =abundance),color=cols, alpha=0.5) + #Chose Colors and size for dots
          scale_size(range = c(0, 40))+
          geom_text_repel(aes(label=ifelse(res$logFC<(-1) & res$adj.P.Val < FDR_cut_off , as.character(res$Name),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5) +
          geom_text_repel(aes(label=ifelse(res$logFC>1 & res$adj.P.Val < FDR_cut_off , as.character(res$Name),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5) +
          geom_hline(yintercept=-log10(FDR_cut_off), color="red",linetype="dashed") +
          xlab("Effect size: log2(fold-change)") +
          ylab("-log10(adjusted p-value)") + 
          theme(panel.background = element_blank(),
                panel.border=element_rect(fill=NA),
                legend.background = element_rect(color=NA),
                legend.key = element_rect(colour = "transparent", fill = "white"),
                plot.title=element_text(size=23, face="bold"))
    } else {
      p<-ggplot(res, aes(x = logFC, y = sig, label=Name)) +
        geom_point(color=cols, alpha=0.5) + #Chose Colors and size for dots
        geom_text_repel(aes(label=ifelse(res$logFC<(-1) & res$adj.P.Val < FDR_cut_off , as.character(res$Name),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5) +
        geom_text_repel(aes(label=ifelse(res$logFC>1 & res$adj.P.Val < FDR_cut_off , as.character(res$Name),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5) +
        geom_hline(yintercept=-log10(FDR_cut_off), color="red",linetype="dashed") +
        xlab("Effect size: log2(fold-change)") +
        ylab("-log10(adjusted p-value)") + 
        theme(panel.background = element_blank(),
              panel.border=element_rect(fill=NA),
              plot.title=element_text(size=23, face="bold"))
    }
  }
  else if (graph_option=="lollipop") {
    #=========================================================
    #////////////////////Lollipop PLOT///////////////////////
    #=========================================================
      p<-ggplot(res, aes(y=reorder(Name,-start), x=logFC,fill=col,size=abundance)) +
                geom_point(color="black",alpha=0.8,shape=21)+
                geom_segment(data=res[res$adj.P.Val<0.2,],aes(yend=reorder(Name,-start)), xend=(-30), color= "black", linetype = "solid",linewidth=1)+ 
                scale_fill_manual(values=c("B"=sample_type_color2,"A"=sample_type_color1,"D"="white"))+ 
                scale_size_continuous(name="Relative Abundance",range=c(5, 20))+                
                ggtitle("EdgeR")+
                theme(panel.background = element_blank(),
                      panel.border=element_rect(fill=NA),
                      panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
                      panel.grid.minor = element_blank(),
                      strip.background=element_blank(),
                      axis.title=element_text(size=20,face="bold"),
                      axis.text.x=element_text(size=18, face="bold"),
                      axis.text.y=element_text(face="bold",size=10),
                      axis.ticks=element_line(colour="gray70"),
                      legend.background = element_rect(color=NA),
                      legend.key = element_rect(colour = "transparent", fill = "white"),
                      plot.title=element_text(size=23, face="bold"))+
                xlab("") +
                ylab("")+
                geom_vline(xintercept=0, color="red",linetype="dashed")+
                guides(fill="none")
  }              
  ##save the legend on a separate png file              
  if (legend_onplot=="no") { #if want to hide the legend on the plot, then will save the legend on separate png
    leg <- get_legend(p)
    legend<-as_ggplot(leg)
    png_output2<-paste0(output_name,"_legend.png")
    png(png_output2, res = 300, width=200, height=200, units='mm')
        show(legend)
    dev.off()
    p<-p + theme(legend.position = "none") ###remove the legend
  }
  png_output<-paste0(output_name,".png")
  png(file=png_output, res = 300, width=200, height=200 , units='mm')
      show(p)
  dev.off()
  
}

Edger_paired_nonPS_KW <- function (input_table, 
                                   pair_by_variable, 
                                   sample_type_var_name, 
                                   sample_types=list(), 
                                   sample_type_color=list(),                             
                                   FDR_cut_off=0.2,
                                   number_display=10,
                                   output_name,
                                   legend_onplot=c("yes","no")) {
  
  #print instruction for the function
  print("This function requires your input_table to be a dataframe with subject ID as column name")
  print("First row of the input_table is for classifying the subjects into different sample types")
  print("First column of the input_table is the names of the metabolites/transcripts")
  
  #making sure the number of item in sample_types is same as sample_type_color
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  
  #making sure pair_by_variable is not missing
  stopifnot("There must be a pair_by_variable to identify the group for the paired analysis"= !missing(pair_by_variable))
  
  #making sure there is only two groups for paired analysis
  if (!missing(pair_by_variable)){ 
    stopifnot("you can only have two groups for paired analysis"= length(sample_types)==2)
    stopifnot("Please make sure the second row of the dataframe is the variable used for pairing"= input_table[2,1]==pair_by_variable)
  }
  
  #setting up color
  sample_type_color_t<-sample_type_color
  names(sample_type_color_t) <-sample_types  
  sample_type_t<-sample_types
  names(sample_type_t)<-sample_type_color
  sample_type_t<-sort(sample_type_t)
  sample_type_color_t2<-names(sample_type_t)
  
  for (i in 1:length(sample_types)) {
    assign(paste0("sample_types",i), sample_types[[i]])
    assign(paste0("sample_type_color",i), sample_type_color[[i]])
  }  
  
  #ensure the first row contains the sample type information in the dataframe
  sample_type_row<- head(input_table,1)
  sample_type<- as.character(sample_type_row[1,1])
  stopifnot("The first row needs to be your sample type information. first cell on the row needs to be the input sample type variable (sample_type_var_name)"= sample_type_var_name==sample_type)
  #checking what elements are in the first row 
  sample_type_elements<- unlist(unique(as.list(sample_type_row[,-1])))
  print("the following sample types are contain in the table:")
  print(sample_type_elements)
  input_table2<-input_table
  
  #keeping only sample types we care about 
  input_table3<-input_table2[,(input_table2[1,]) %in%  c(sample_type_var_name, sample_types)]
  
  input_table4<-input_table3[,-1]
  rownames(input_table4)<-input_table3[,1]
  
  input_table5<-input_table4[-c(1,2),]
  #changing all the columns to numeric
  input_table5 <- input_table5 %>% mutate_if(is.character, as.numeric)

  #sample_type_table only contains subjectID and sample type
  sample_type_table<-data.frame(t(input_table3[c(1,2),-1])) 
  colnames(sample_type_table)<-c(sample_type_var_name,pair_by_variable)   

  ## making sure function input is appropiate
  legend_onplot <- match.arg(legend_onplot)
  if (missing(legend_onplot)) { legend_onplot<-"yes"} #default being yes for legend_onplot, which will just display legend on the graph
  
  #Set order
  coldata_SAMPLE <- data.frame(sample_type_table[order(row.names(sample_type_table)),,drop=FALSE],check.names=F)
  countdata_SAMPLE <- data.frame(input_table5[, order(colnames(input_table5)),drop=FALSE],check.names=F)
  
  #making a compare variable to standardize the input variable for comparsion- make sure the reference is always control
  coldata_SAMPLE$sample_type_var_name<-coldata_SAMPLE[[sample_type_var_name]]
  coldata_SAMPLE<- coldata_SAMPLE %>% dplyr::mutate(compare=case_when((sample_type_var_name==sample_types1)~"control", (sample_type_var_name==sample_types2)~"treatment"))
  coldata_SAMPLE<-coldata_SAMPLE[!is.na(coldata_SAMPLE$compare),] ###dropping subjects if they are not in the group of comparsion
  coldata_SAMPLE$compare<- as.factor(coldata_SAMPLE$compare) #converting variable to factor
  countdata_SAMPLE<-countdata_SAMPLE[,colnames(countdata_SAMPLE) %in% rownames(coldata_SAMPLE)]#making sure coldata has same subjects and countdata
  
  coldata_SAMPLE<-coldata_SAMPLE
  countdata_SAMPLE<-countdata_SAMPLE
  
  transcript<-row.names(countdata_SAMPLE)
  Patient<- factor(coldata_SAMPLE$match_group)
  study_sample_type<-factor(coldata_SAMPLE$compare)
  design<-model.matrix(~Patient + study_sample_type)
  rownames(design) <- colnames(countdata_SAMPLE)
  
  #https://gtpb.github.io/ADER19F/pages/tutorial2
  y <- DGEList(counts=as.matrix(countdata_SAMPLE), genes=transcript)
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust=TRUE)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit)

  resNoFilt <- topTags(lrt, n=nrow(y), adjust.method="BH", sort.by="PValue")
  resNoFilt<- data.frame(resNoFilt)
  
  resNoFilt<- subset(resNoFilt,select=c(logFC, logCPM, PValue, FDR))
  resNoFilt_export<<-resNoFilt
  #check to see if there is any taxa with FDR<FDR_cut_off. If not then will stop the function, because there will be error if it continues as the object will have no rows
  checkifallzero<-resNoFilt[resNoFilt$FDR<FDR_cut_off,]
  if (dim(checkifallzero)[1] == 0) {
    csv_output2<-paste0(output_name,"_all.csv")
    write.csv(resNoFilt_export,file=csv_output2, row.names = TRUE)
    print(paste0("no Taxa with FDR more than",FDR_cut_off))
    return(NULL)
  }
  resNoFilt[nrow(resNoFilt) + 1,] <- c(0, 0, 0, 0) ###adding this row of 0 so that even if there is no significant top 10 for upregulated or downregulated, the code will still run
  #Create Table of Down Regulated
  sigDownReg <- resNoFilt[resNoFilt$FDR<FDR_cut_off,] #keeping only FDR < 0.2
  #sigDownReg <- resNoFilt
  sigDownReg <- sigDownReg[sigDownReg$logFC<=0,] ##keeping only the negative logFC
  sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
  sigDownReg1 <- sigDownReg %>% dplyr::slice(1:number_display) #Subset Top (number_display)
  #Create Table of Down Regulated
  sigUpReg <- resNoFilt[resNoFilt$FDR<FDR_cut_off,] #keeping only FDR < 0.2
  #sigUpReg <- resNoFilt
  sigUpReg <- sigUpReg[sigUpReg$logFC>=0,] ##keeping only the positive logFC
  sigUpReg <- sigUpReg[order(sigUpReg$logFC, decreasing=TRUE),]
  sigUpReg1 <- sigUpReg %>% dplyr::slice(1:number_display) #Subset Top (number_display)
  #Merge Top 10s
  res <- rbind(sigDownReg1,sigUpReg1)
  res <- res[res$logFC!=0,] 
  
  #Create Relative Abundance Table
  df <-
    countdata_SAMPLE %>% 
    rownames_to_column('gs') %>%
    dplyr::group_by(gs) %>% 
    dplyr::summarise_all(funs(sum)) %>%
    dplyr::mutate_if(is.numeric, funs(./sum(.))) %>%
    column_to_rownames('gs')
  df_export<<-df
  #Get the ColData for Each Comparison
  coldata.1 <- coldata_SAMPLE[coldata_SAMPLE[[sample_type_var_name]]==sample_types1,]
  coldata.2 <- coldata_SAMPLE[coldata_SAMPLE[[sample_type_var_name]]==sample_types2,]
  #keep Count data only for each comparison
  needed<-which(colnames(df) %in% rownames(coldata.1))    
  df.1 <- df[,needed]
  needed2<-which(colnames(df) %in% rownames(coldata.2))    
  df.2 <- df[,needed2]
  #Convert Resuts table into a data.frame
  res <- as.data.frame(res)
  #decide what otu to save 
  otu.to.save <-as.character(rownames(res))
  #from relative table we should get the mean across the row of the otu table
  df.1.meanRA <- rowMeans(df.1)
  df.2.meanRA <- rowMeans(df.2)
  #need to subset AND reorder just the otus that we have 
  df.1.meanRA.save <- df.1.meanRA[otu.to.save]
  df.2.meanRA.save <- df.2.meanRA[otu.to.save]
  #add the abundnace data for the res dataframe
  res$abundance.1 <- df.1.meanRA.save
  res$abundance.2 <- df.2.meanRA.save
  #Set Names of Results Table
  res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Name","logFC", "lfcSE", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
  #Remove Any Data without LOGFC data
  res <- res[!is.na(res$logFC),]
  # Reorder Results based on FDR for comparison 1
  res = res[order(-res$logFC, na.last = TRUE), ]
  #Create Order
  res <- res %>% dplyr::mutate(start = 1:n())
  #Convert Important columns to Numeric
  res$adj.P.Val <-   as.numeric(as.character(res$adj.P.Val))
  res$logFC <-       as.numeric(as.character(res$logFC))
  #Replace NA
  res <- res %>% dplyr::mutate(adj.P.Val = if_else(is.na(adj.P.Val), 0.9, adj.P.Val))
  ##convert abundance to numeric
  res$abundance.1 <- as.numeric(as.character(res$abundance.1))
  res$abundance.2 <- as.numeric(as.character(res$abundance.2))
  #Create Variable for Color based on Comparison, FDR and LOGFC
  res$col <- ifelse(res$logFC>0, "B",
                    ifelse(res$logFC<0, "A","D"))
  res$abundance <- ifelse(res$col=="A", res$abundance.1, ifelse(res$col=="B", res$abundance.2,0))
  
  #########EXCEL EXPORT
  ####keeping all results 
  res_all <- rbind(sigDownReg,sigUpReg)
  res_all <- res_all[res_all$logFC!=0,] 
  
  #Convert Resuts table into a data.frame
  res_all <- as.data.frame(res_all)
  #decide what otu to save 
  otu.to.save_all <-as.character(rownames(res_all))
  #from relative table we should get the mean across the row of the otu table
  df.1.meanRA <- rowMeans(df.1)
  df.2.meanRA <- rowMeans(df.2)
  #need to subset AND reorder just the otus that we have 
  df.1.meanRA.save <- df.1.meanRA[otu.to.save_all]
  df.2.meanRA.save <- df.2.meanRA[otu.to.save_all]
  #add the abundnace data for the res dataframe
  res_all$abundance.1 <- df.1.meanRA.save
  res_all$abundance.2 <- df.2.meanRA.save
  
  #Export result table
  csv_output<-paste0(output_name,"_significant.csv")
  write.csv(res_all,file=csv_output, row.names = TRUE)
  
  csv_output2<-paste0(output_name,"_all.csv")
  write.csv(resNoFilt_export,file=csv_output2, row.names = TRUE)
  
  p<-ggplot(res, aes(y=reorder(Name,-start), x=logFC,fill=col,size=abundance)) +
    geom_point(color="black",alpha=0.8,shape=21)+
    geom_segment(data=res[res$adj.P.Val<0.2,],aes(yend=reorder(Name,-start)), xend=(-30), color= "black", linetype = "solid",linewidth=1)+ 
    scale_fill_manual(values=c("B"=sample_type_color2,"A"=sample_type_color1,"D"="white"))+ 
    scale_size_continuous(name="Relative Abundance",range=c(5, 20))+                
    ggtitle("EdgeR (Paired)")+
    theme(panel.background = element_blank(),
          panel.border=element_rect(fill=NA),
          panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
          panel.grid.minor = element_blank(),
          strip.background=element_blank(),
          axis.title=element_text(size=20,face="bold"),
          axis.text.x=element_text(size=18, face="bold"),
          axis.text.y=element_text(face="bold",size=10),
          axis.ticks=element_line(colour="gray70"),
          legend.background = element_rect(color=NA),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          plot.title=element_text(size=23, face="bold"))+
    xlab("") +
    ylab("")+
    geom_vline(xintercept=0, color="red",linetype="dashed")+
    guides(fill="none")
  ##save the legend on a separate png file              
  if (legend_onplot=="no") { #if want to hide the legend on the plot, then will save the legend on separate png
    leg <- get_legend(p)
    png_output2<-paste0(output_name,"_legend.png")
    png(png_output2, res = 300, width=200, height=200, units='mm')
    legend<-as_ggplot(leg)
    show(legend)
    dev.off()
    p<-p + theme(legend.position = "none") ###remove the legend
  }
  png_output<-paste0(output_name,".png")
  png(file=png_output, res = 300, width=200, height=200 , units='mm')
  show(p)
  dev.off()
  
  ##########################################################################################################################
  ##########################################################################################################################
}


Deseq_nonPS_KW <- function(input_table, 
                           sample_type_var_name, 
                           sample_types=list(), 
                           sample_type_color=list(),                             
                           alpha_level=0.2,
                           number_display=10,
                           graph_option=c("volcano","lollipop"),
                           legend_onplot=c("yes","no"),
                           output_name) {
  #print instruction for the function
  print("This function requires your input_table to be a dataframe with subject ID as column name")
  print("First row of the input_table is for classifying the subjects into different sample types")
  print("First column of the input_table is the names of the metabolites/transcripts")
  
  #initially function options
  alpha<-alpha_level
  ## evaluate choices
  graph_option <- match.arg(graph_option)
  legend_onplot <- match.arg(legend_onplot)
  print(graph_option)
  
  #making sure the number of item in sample_types is same as sample_type_color
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  
  #setting up color
  sample_type_color_t<-sample_type_color
  names(sample_type_color_t) <-sample_types  
  sample_type_t<-sample_types
  names(sample_type_t)<-sample_type_color
  sample_type_t<-sort(sample_type_t)
  
  for (i in 1:length(sample_types)) {
    assign(paste0("sample_types",i), sample_types[[i]])
    assign(paste0("sample_type_color",i), sample_type_color[[i]])
  }  
  
  #ensure the first row contains the sample type information in the dataframe
  sample_type_row<- head(input_table,1)
  sample_type<- as.character(sample_type_row[1,1])
  stopifnot("The first row needs to be your sample type information. first cell on the row needs to be the input sample type variable (sample_type_var_name)"= sample_type_var_name==sample_type)
  #checking what elements are in the first row 
  sample_type_elements<- unlist(unique(as.list(sample_type_row[,-1])))
  print("the following sample types are contain in the table:")
  print(sample_type_elements)
  input_table2<-input_table
  
  #keeping only sample types we care about 
  input_table3<-input_table2[,(input_table2[1,]) %in%  c(sample_type_var_name, sample_types)]
  
  input_table4<-input_table3[,-1]
  rownames(input_table4)<-input_table3[,1]
  
  input_table5<-input_table4[-1,]
  #changing all the columns to numeric
  input_table5 <- input_table5 %>% mutate_if(is.character, as.numeric)

  #sample_type_table only contains subjectID and sample type
  sample_type_table<-data.frame(t(input_table3[1,-1]))
  colnames(sample_type_table)<-sample_type_var_name
  
  ## making sure function input is appropiate
  legend_onplot <- match.arg(legend_onplot)
  if (missing(legend_onplot)) { legend_onplot<-"yes"} #default being yes for legend_onplot, which will just display legend on the graph
  
  #Set order
  coldata_SAMPLE <- data.frame(sample_type_table[order(row.names(sample_type_table)),,drop=FALSE],check.names=F)
  countdata_SAMPLE <- data.frame(input_table5[, order(colnames(input_table5)),drop=FALSE],check.names=F)
  
  #making a compare variable to standardize the input variable for comparsion- make sure the reference is always control
  coldata_SAMPLE$sample_type_var_name<-coldata_SAMPLE[[sample_type_var_name]]
  coldata_SAMPLE<- coldata_SAMPLE %>% dplyr::mutate(compare=case_when((sample_type_var_name==sample_types1)~"control", (sample_type_var_name==sample_types2)~"treatment"))
  coldata_SAMPLE$sample_type_var_name<-as.factor(coldata_SAMPLE$sample_type_var_name)
  coldata_SAMPLE<-coldata_SAMPLE[!is.na(coldata_SAMPLE$compare),] ###dropping subjects if they are not in the group of comparsion
  coldata_SAMPLE$compare<- as.factor(coldata_SAMPLE$compare) #converting variable to factor
  countdata_SAMPLE<-countdata_SAMPLE[,colnames(countdata_SAMPLE) %in% rownames(coldata_SAMPLE)]#making sure coldata has same subjects and countdata
  
  #Deseq analysis
  diagdds <- DESeqDataSetFromMatrix(countData = countdata_SAMPLE,
                                    colData = coldata_SAMPLE,
                                    design= ~sample_type_var_name )
  
  ##Calculate geometric means prior to estimate size factor
  gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
  geoMeans = apply(counts(diagdds), 1, gm_mean)        
  #Estimate Factors of DESeq Object
  diagdds <- estimateSizeFactors(diagdds, geoMeans = geoMeans)
  #making sure reference is control. resetting the level to make sure control is reference
  diagdds$compare <- droplevels(diagdds$compare)
  diagdds$compare <- relevel(diagdds$compare, ref ="control")

  diagdds <- DESeq(diagdds)
  #Output Result Table
  res <- DESeq2::results(diagdds, cooksCutoff = FALSE)
  res_org<<-res
  
  #Convert Resuts table into a data.frame
  res <- as.data.frame(res)
  
  res<- subset(res,select=c(baseMean,log2FoldChange,lfcSE,stat,pvalue,padj))
  resNoFilt<- res #resNoFilt is backup object
  
  #Create Relative Abundance Table from the 
  df <-
    countdata_SAMPLE %>% 
    rownames_to_column('gs') %>%
    dplyr::group_by(gs) %>% 
    dplyr::summarise_all(funs(sum)) %>%
    dplyr::mutate_if(is.numeric, funs(./sum(.))) %>%
    column_to_rownames('gs')
  #Get the ColData for Each Comparison
  coldata.1 <- coldata_SAMPLE[coldata_SAMPLE[[sample_type_var_name]]==sample_types1,]
  coldata.2 <- coldata_SAMPLE[coldata_SAMPLE[[sample_type_var_name]]==sample_types2,]
  #keep Count data only for each comparison
  needed<-which(colnames(df) %in% rownames(coldata.1))    
  df.1 <- df[,needed]
  needed2<-which(colnames(df) %in% rownames(coldata.2))    
  df.2 <- df[,needed2]
  #decide what otu to save 
  otu.to.save <-as.character(rownames(res))
  #from relative table we should get the mean across the row of the otu table
  df.1.meanRA <- rowMeans(df.1)
  df.2.meanRA <- rowMeans(df.2)
  #need to subset AND reorder just the otus that we have 
  df.1.meanRA.save <- df.1.meanRA[otu.to.save]
  df.2.meanRA.save <- df.2.meanRA[otu.to.save]
  #add the abundnace data for the res dataframe
  res$abundance.1 <- df.1.meanRA.save
  res$abundance.2 <- df.2.meanRA.save
  #Set Names of Results Table
  res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Name","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
  #Write Tables of Differential Analysis
  csv_output<-paste0(output_name,".csv")
  write.csv(res,file=csv_output, row.names = TRUE)
  
  res_regulate<-as.data.frame(res) 
  
  res_regulate<-res_regulate %>% filter(!is.na(adj.P.Val))
  #check to see if there is any taxa with adjusted p value<alpha. If not then will stop the function, because there will be error if it continues as the object will have no rows
  checkifallzero<-res_regulate[res_regulate$adj.P.Val<alpha,]
  if (dim(checkifallzero)[1] == 0) {
    print(paste0("no Taxa with adjusted p value more than"),alpha)
    return(NULL)
  }
  res_regulate_positiveFC<-res_regulate %>% filter(res_regulate$logFC > 0 & res_regulate$adj.P.Val < alpha)
  if (dim(res_regulate_positiveFC)[1] == 0) {### if there is no upregulated with adjusted p value less than alpha, then 
    res_regulate_positiveFC[nrow(res_regulate_positiveFC) + 1,] <- replicate(ncol(res_regulate_positiveFC),0) ###adding this row of 0 so that even if there is no significant top 10 for upregulated or downregulated, the code will still run
  }
  res_regulate_positiveFC <- res_regulate_positiveFC[order(-res_regulate_positiveFC$logFC),] #order has the lowest value on the top, so need the - to reverse the order since slice cut from top down 
  res_regulate_positiveFC1 <- res_regulate_positiveFC %>% dplyr::slice(1:number_display)
  
  res_regulate_negativeFC<-res_regulate %>% filter(res_regulate$logFC < 0 & res_regulate$adj.P.Val < alpha)
  if (dim(res_regulate_negativeFC)[1] == 0) {
    res_regulate_negativeFC[nrow(res_regulate_negativeFC) + 1,] <- replicate(ncol(res_regulate_negativeFC),0) ###adding this row of 0 so that even if there is no significant top 10 for upregulated or downregulated, the code will still run
  }
  res_regulate_negativeFC <- res_regulate_negativeFC[order(res_regulate_negativeFC$logFC),]
  res_regulate_negativeFC1 <- res_regulate_negativeFC %>% dplyr::slice(1:number_display) #Subset Top 10
  #Merge Top 10s (number_display)
  res_regulate_final <- rbind(res_regulate_negativeFC1,res_regulate_positiveFC1)
  res_regulate_final <- res_regulate_final[res_regulate_final$logFC!=0,] 
  res_regulate_final$col <- ifelse(res_regulate_final$logFC>0, "B",
                                   ifelse(res_regulate_final$logFC<0, "A","D"))
  #Convert to dataframe
  res_regulate_final$abundance <- ifelse(res_regulate_final$col=="A", res_regulate_final$abundance.1, ifelse(res_regulate_final$col=="B", res$abundance.2,0))
  
  # Reorder Results based on FDR for comparison 1
  res_regulate_final = res_regulate_final[order(-res_regulate_final$logFC, na.last = TRUE), ]
  #Create Order
  res_regulate_final <- res_regulate_final %>% dplyr::mutate(start = 1:n())
  res_regulate_final$abundance<-ifelse(res_regulate_final$adj.P.Val<alpha & res_regulate_final$logFC>0, res_regulate_final$abundance.2, 
                                       ifelse(res_regulate_final$adj.P.Val<alpha & res_regulate_final$logFC<0, res_regulate_final$abundance.1,2))
  

  ######################################################################################################################
  ######################################################################################################################
  if (graph_option=="volcano") {
    #=========================================================
    #////////////////////VOLCANO PLOT///////////////////////
    #=========================================================
    # Compute significance, with a maximum of 320 for the p-values set to 0 due to limitation of computation precision
    res$sig <- -log10(res$adj.P.Val)
    sum(is.infinite(res$sig))
    res[is.infinite(res$sig),"sig"] <- 350
    ## Volcano plot of adjusted p-values
    cols <- densCols(res$logFC, res$sig)
    cols[res$pvalue ==0] <- "purple"
    cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- sample_type_color2
    cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- sample_type_color1
    res$pch <- 19
    res$pch[res$pvalue ==0] <- 6
    
    ####scaling the abundance so the circle size on the graph is standardized
    res$abundance_scaled<-ifelse(res$adj.P.Val<alpha & res$logFC>0, res$abundance.2, ifelse(res$adj.P.Val<alpha & res$logFC<0, res$abundance.1,2))
    res_keep<-res[!(is.na(res$adj.P.Val)),]
    res_keep<-res_keep[res_keep$adj.P.Val<alpha,]
    max_value<-max(res_keep$abundance_scaled)
    value_max<-40/max_value
    res <- res %>% mutate(abundance_scaled2=case_when((adj.P.Val<alpha & !(is.na(adj.P.Val))) ~ abundance_scaled*value_max, T ~abundance_scaled ))
    
    #PLOT IT
    g<-ggplot(res, aes(x = logFC, y = sig, label=Name)) +
      geom_point(size = res$abundance_scaled2,color=cols, alpha=0.5) + #Chose Colors and size for dots
      geom_text_repel(aes(label=ifelse(res$logFC<(-1) & res$adj.P.Val < alpha , as.character(res$Name),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5) +
      geom_text_repel(aes(label=ifelse(res$logFC>1 & res$adj.P.Val < alpha , as.character(res$Name),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5) +
      theme(legend.position = "none") +
      geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
      xlab("Effect size: log2(fold-change)") +
      ylab("-log10(adjusted p-value)") + 
      #ylim(0,20)+
      ggtitle("Deseq")+
      theme(plot.title=element_text(size=23, face="bold"))
    png_output<-paste0(output_name,".png")
    png(file=png_output, res = 300, width=200, height=200 , units='mm')
    show(g)
    dev.off() 
  }
  else if (graph_option=="lollipop") {
    p<-ggplot(res_regulate_final, aes(y=reorder(Name,-start), x=logFC,fill=col,size=abundance)) +
      geom_point(color="black",alpha=0.8,shape=21)+
      geom_segment(data=res_regulate_final[res_regulate_final$adj.P.Val<alpha,],aes(yend=reorder(Name,-start)), xend=(-30), color= "black", linetype = "solid",size=1)+ 
      scale_fill_manual(values=c("B"=sample_type_color2,"A"=sample_type_color1,"D"="white"))+ 
      scale_size_continuous(name="Relative Abundance",range=c(5, 20))+
      ggtitle("Deseq")+
      theme(panel.background = element_blank(),
            panel.border=element_rect(fill=NA),
            panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
            panel.grid.minor = element_blank(),
            strip.background=element_blank(),
            axis.title=element_text(size=20,face="bold"),
            axis.text.x=element_text(colour="black", size=18, face="bold"),
            axis.text.y=element_text(face="bold",size=10),
            axis.ticks=element_line(colour="black"),
            legend.background = element_rect(color=NA),
            legend.key = element_rect(colour = "transparent", fill = "white"),
            plot.title=element_text(size=23, face="bold"))+
      xlab("") +
      ylab("")+
      #xlim(-7,7)+
      geom_vline(xintercept=0, color="red",linetype="dashed")+
      guides(fill="none")
    ##save the legend on a separate png file              
    if (legend_onplot=="no") { #if want to hide the legend on the plot, then will save the legend on separate png
      leg <- get_legend(p)
      png_output2<-paste0(output_name,"_legend.png")
      png(png_output2, res = 300, width=200, height=200, units='mm')
      legend<-as_ggplot(leg)
      show(legend)
      dev.off()
    }
    png_output<-paste0(output_name,".png")
    png(file=png_output, res = 300, width=200, height=200 , units='mm')
    if (legend_onplot=="no") { #if want to hide the legend
      p<-p + theme(legend.position = "none") ###remove the legend
    }
    show(p)
    dev.off()
  }
  
}

Prune_w_Deseq_nonPS_KW<- function(input_table, 
                                  sample_type_var_name, 
                                  sample_types=list(), 
                                  normalized_count=100) {
  #print instruction for the function
  print("This function requires your input_table to be a dataframe with subject ID as column name")
  print("First row of the input_table is for classifying the subjects into different sample types")
  print("First column of the input_table is the names of the metabolites/transcripts")
 
  for (i in 1:length(sample_types)) {
    assign(paste0("sample_types",i), sample_types[[i]])
  }  
  
  #ensure the first row contains the sample type information in the dataframe
  sample_type_row<- head(input_table,1)
  sample_type<- as.character(sample_type_row[1,1])
  stopifnot("The first row needs to be your sample type information. first cell on the row needs to be the input sample type variable (sample_type_var_name)"= sample_type_var_name==sample_type)
  #checking what elements are in the first row 
  sample_type_elements<- unlist(unique(as.list(sample_type_row[,-1])))
  print("the following sample types are contain in the table:")
  print(sample_type_elements)
  input_table2<-input_table
  
  #keeping only sample types we care about 
  input_table3<-input_table2[,(input_table2[1,]) %in%  c(sample_type_var_name, sample_types)]
  
  input_table4<-input_table3[,-1]
  rownames(input_table4)<-input_table3[,1]
  
  input_table5<-input_table4[-1,]
  #changing all the columns to numeric
  input_table5 <- input_table5 %>% mutate_if(is.character, as.numeric)
  
  #sample_type_table only contains subjectID and sample type
  sample_type_table<-data.frame(t(input_table3[1,-1]))
  colnames(sample_type_table)<-sample_type_var_name
  
  #Set order
  coldata_SAMPLE <- data.frame(sample_type_table[order(row.names(sample_type_table)),,drop=FALSE],check.names=F)
  countdata_SAMPLE <- data.frame(input_table5[, order(colnames(input_table5)),drop=FALSE],check.names=F)
  
  #making a compare variable to standardize the input variable for comparsion- make sure the reference is always control
  coldata_SAMPLE$sample_type_var_name<-coldata_SAMPLE[[sample_type_var_name]]
  coldata_SAMPLE<- coldata_SAMPLE %>% dplyr::mutate(compare=case_when((sample_type_var_name==sample_types1)~"control", (sample_type_var_name==sample_types2)~"treatment"))
  coldata_SAMPLE$sample_type_var_name<-as.factor(coldata_SAMPLE$sample_type_var_name)
  coldata_SAMPLE<-coldata_SAMPLE[!is.na(coldata_SAMPLE$compare),] ###dropping subjects if they are not in the group of comparsion
  coldata_SAMPLE$compare<- as.factor(coldata_SAMPLE$compare) #converting variable to factor
  countdata_SAMPLE<-countdata_SAMPLE[,colnames(countdata_SAMPLE) %in% rownames(coldata_SAMPLE)]#making sure coldata has same subjects and countdata
  
  #Deseq analysis
  diagdds <- DESeqDataSetFromMatrix(countData = countdata_SAMPLE,
                                    colData = coldata_SAMPLE,
                                    design= ~sample_type_var_name )
  #Normalization Step 
  diagdds <- estimateSizeFactors(diagdds)
  #filter genes where there are less than 3 samples with normalized counts greater than or equal to "normalized_count".
  idx <- rowSums( counts(diagdds, normalized=TRUE) >= normalized_count ) >= 3
  print(count(idx))
  idx<-idx[idx==T] #keeping only filtered genes
  idx<-data.frame(idx)
  idx$gene<-rownames(idx)
  idx<-idx %>% select(-1)
  return(idx)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
beta_div_removal_testing_KW <- function(contaminant_level_df,
                                        input_phyloseq, 
                                        sample_type_var_name, 
                                        sample_types=list(), 
                                        sample_type_color=list(), 
                                        prune_sample_types=list(), #this indicates which sample types you want to prune. if sample_types is the same as prune_sample_types then all the taxa which are consider contaminant would be removed. 
                                        p_value=c("yes","no"),
                                        output_name, 
                                        width=200, height=200,
                                        plot_title_size=14, axis_title_size=18, label_size=7,
                                        x_axis_flip, 
                                        prune_method=c("remove","zero"), #prune method- if remove then the taxa will be remove from the taxa_table. if zero
                                        beta_method=c("bray", "unifrac"),
                                        p_value_location=c("TR","TL","BR","BL")) { #p_value_location: TR-top right; TL- top left; BR- bottom right; BL- bottom left
  ###### this function is exclusively run after decontaminant_testing_KW
  ###### this function display 
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  if(missing(prune_method)) { prune_method<-"remove"}#default being remove taxa that is in the prune
  if(missing(x_axis_flip)) { x_axis_flip<-"no"} #default being no for x_axis_flip
  if(missing(beta_method)) { beta_method<-"bray"} #default being bray curtis
  if(missing(p_value_location) & tolower(p_value)=="yes"){
    p_value_location<-"BL" #default location for p_value will be in the bottom left corner
  } else if(!missing(p_value_location) & tolower(p_value)=="yes"){ #if p_value is wanted and p_value_location is not missing
    print("p_value_location: TR-top right; TL- top left; BR- bottom right; BL- bottom left")
    p_value_location<-match.arg(p_value_location)
    print(paste0("selected:",p_value_location))
  }
  #making sure the number of item in sample_types is same as sample_type_color
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  #making sure only yes and no is selected for p_value
  p_value <- match.arg(p_value)
  #setting up the coloring of the graph 
  sample_type_color_t<-sample_type_color
  names(sample_type_color_t) <-sample_types
  sample_type2<-sort(sample_types)
  
  #Beta diversity
  #make sure input_phyloseq only include the sample_types of interest
  keep_sample<- get_variable(input_phyloseq,sample_type_var_name) %in% sample_types
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  phyloseq_t <- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
  
  # Function to extract last element with delimiter
  extractLastElement <- function(text, delimiter) {
    elements <- unlist(strsplit(text, delimiter))
    lastElement <- tail(elements, 1)
    return(lastElement)
  }
  
  #number_column should be equal to the number of levels which was tested 
  number_column<-ncol(contaminant_level_df)
  
  output_level_graphs <- data.frame(level = numeric(number_column),
                                    num_contaminant = numeric(number_column),
                                    p_value = numeric(number_column))
  num<-1
  combin_p<-list()
  for (i in colnames(contaminant_level_df)){
    #the columns from the file should be level_###
    level<-as.numeric(sapply(i, extractLastElement, delimiter = "_"))
    print(level) # display level
    #input value: level for contaminant testing
    output_level_graphs[num,1]<-level
    #level will be the exact level when the contaminant testing was performed at 
    #keeping only one column at a time 
    taxa_contaminant_list<-contaminant_level_df[i]
    colnames(taxa_contaminant_list)<-"contaminant"
    taxa_contaminant_list<<-taxa_contaminant_list
    # Apply the function to the column
    taxa_contaminant_list$ASV <- sapply(row.names(taxa_contaminant_list), extractLastElement, delimiter = "\\.")
    #keeping only the noncontaminant
    taxa_contaminant_list_noncont<-taxa_contaminant_list %>% filter(contaminant==FALSE)
    NOT_contaminant_list<-taxa_contaminant_list_noncont$ASV
    noncontam_taxa_num<-length(NOT_contaminant_list)
    print(paste0("number of remaining taxa: ", noncontam_taxa_num))
    #input value: number of noncontaminant
    output_level_graphs[num,2]<-noncontam_taxa_num
    
    #number of contaminant revealed
    taxa_contaminant_list_cont<-taxa_contaminant_list %>% filter(contaminant==TRUE)
    contaminant_list<-taxa_contaminant_list_cont$ASV
    contam_taxa_num<-length(contaminant_list)
    
    #input value: number of contaminant
    output_level_graphs[num,3]<-contam_taxa_num
    
    print(paste0("prune method: ", prune_method)) #display prune method
    
    #the pruning method matters because with unifrac, the phylogenic tree is considered when calculating the distance so if a taxa is prune, it would not be entered into the calculation. So
    #if you want to prune only a certain sample types, but not all, then you would have to replace the taxa with 0, rather than simply removing them from the phyloseq. 
    if (prune_method=="remove") { #if prune_method is remove then will just use the prune_taxa function which just remove the taxa from the phyloseq 
            #keeping only taxa which is NOT consider contaminant    
            phyloseq_t2 <- prune_taxa(NOT_contaminant_list, phyloseq_t)
            print("After pruning: ")
            print(phyloseq_t2)
    } else {
      #this section select on the row number in which the contaminant taxa is located and then replacing them with 0
      #it select the columns that has the sample type of interest
      taxa_indices <- which(rownames(otu_table(phyloseq_t)) %in% contaminant_list)
      #sample_data_df is the sample data of the phyloseq. row is each unique sample
      sample_data_df <- as.data.frame(sample_data(phyloseq_t))
      #figure out which sample ID has the sample type of interest
      rows_with_target_samples_types<-rownames(sample_data_df[sample_data_df[[sample_type_var_name]] %in% prune_sample_types, ])
      #use the sample IDs identified from previous step to see which column from the otu table to select 
      prune_sample_types_indices <- which(colnames(otu_table(phyloseq_t)) %in% rows_with_target_samples_types)
      phyloseq_t2<-phyloseq_t
      #now you change the value on the otu table from the phyloeq to be 0 
      for (colnum in prune_sample_types_indices) {
              otu_table(phyloseq_t2)[taxa_indices, colnum] <- 0 #in a particular sample, replace the selected taxa to 0
      }
    }

    #recalculate the relative abundance  
    normalizeSample <- function(x){x/sum(x)}
    phyloseq_t3 = transformSampleCounts(phyloseq_t2,normalizeSample) ##relative abundance

    #Create Distance Matrix with Bray
    if (beta_method=="bray"){
      #Create Distance Matrix with Bray
      vegdist=vegdist(t(otu_table(phyloseq_t3)), method = "bray") #row is samplesID #column is taxa
    } else {
      set.seed(123)
      vegdist= phyloseq::UniFrac(phyloseq_t3, weighted = TRUE) 
    }
    #Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
    CmdScale <- cmdscale(vegdist, k =10)
    #calculated Sample variance for each PC
    vars <- apply(CmdScale, 2, var)
    #Create Variable with the Percent Variance
    percentVar <- round(100 * (vars/sum(vars)))
    #Merge PC Data with MetaData
    newResults <- merge(x = CmdScale, y = sample_data(phyloseq_t3), by = "row.names", all.x = TRUE)  #sample_data(phyloseq_t3) row is the sampleID
    #Rename Variables for PC1 and PC2
    colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
    colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
    colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"
    formula <- as.formula(paste("cbind(PC1,PC2)",sample_type_var_name, sep=" ~ "))
    #Calculate the Centroid Value
    centroids <- aggregate(formula ,data= newResults, mean)
    #Merge the Centroid Data into the PCOA Data
    newResults <- merge(newResults,centroids,by=sample_type_var_name,suffixes=c("",".centroid"))
    # Calculate p-value with ADONIS
    x=adonis2(vegdist ~ phyloseq_t3@sam_data[[sample_type_var_name]])
    pvalues<-x$`Pr(>F)`[[1]]
    
    #input value: p value
    output_level_graphs[num,4]<-pvalues
    
    subtitle_output<- paste0("Adonis, p=",pvalues)
    p<-ggplot(newResults, aes(PC1, PC2, color= get(sample_type_var_name))) + # Graph PC1 and PC2
      geom_point(size=2) + # Set the size of the points
      xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2 
      labs(title=NULL, 
           fill="Centroids") +
      ggtitle(paste0(level,paste0(paste0(" ( N=",noncontam_taxa_num),")")))+
      geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= get(sample_type_var_name)))+ 
      geom_point(data=centroids, aes(x=PC1, y=PC2, color= get(sample_type_var_name)), size=0) + 
      geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=sample_type2), size=label_size) +
      scale_color_manual(name=NULL,  
                         values=sample_type_color_t) +
      theme(panel.background = element_blank(),
            panel.border=element_rect(fill=NA),
            panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
            panel.grid.minor = element_blank(),strip.background=element_blank(),
            plot.title=element_text(face="bold",hjust=0.5, size = plot_title_size), 
            plot.subtitle = element_text(hjust=0.5),
            axis.title=element_text(face="bold", size = axis_title_size),
            axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
            axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
            axis.ticks=element_blank(),
            plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
    p<-p
    newResults<<-newResults
    if (tolower(x_axis_flip)=="yes") {
      p<-p+scale_x_reverse()
    }
    if (p_value=="yes") { ###when the axis is flip, the annotation for the p value need to placed on max(newResults$PC1), instead of min(newResults$PC2))
      if (tolower(x_axis_flip)=="no"){
        if (p_value_location=="BR"){
          p <- p + annotate("text", x=max(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
        } else if (p_value_location=="BL"){
          p <- p + annotate("text", x=min(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
        } else if (p_value_location=="TR"){
          p <- p + annotate("text", x=max(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
        } else if (p_value_location=="TL"){
          p <- p + annotate("text", x=min(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
        }
      } else if (tolower(x_axis_flip)=="yes") {
        if (p_value_location=="BR"){
          p <- p + annotate("text", x=min(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
        } else if (p_value_location=="BL"){
          p <- p + annotate("text", x=max(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
        } else if (p_value_location=="TR"){
          p <- p + annotate("text", x=min(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
        } else if (p_value_location=="TL"){
          p <- p + annotate("text", x=max(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
        }
      }
    }
    
    assign(paste0("p",num),p)
    #png_output<-paste0(paste0(paste0(output_name,"_level"),level),".png")
    #png(file=png_output, res = 300, width=width, height=height , units='mm')
    #      show(p)
    #dev.off()
    combin_p[[num]]<-p
    num<-num+1
  }
  
  combin_p2<<-combin_p
  #put the graphs next to each other. add a title and collec the legend- place it at bottom and horizontal
  combine_plot<-wrap_plots(combin_p2, nrow = ceiling(number_column/5)) + plot_layout(guides = "collect")  &  
    plot_annotation(title = output_name,
                    theme= theme(title = element_text(color = "red", face = "bold", size = 24))) 
  
  #output pdf with multiple PCA plots based on level of contaminant testing  
  pdf_output<-paste0(paste0(output_name,"_level"),".pdf")
  pdf(pdf_output, width=22, height=7*ceiling(number_column/5))
  show(combine_plot)
  dev.off()
  
  colnames(output_level_graphs)<-c("level", "non_contaminant", "contaminant", "p_value")
  
  #output csv which contains how many taxa were identified as contaminant/noncontaminant, along with the p value for the PCA plots
  csv_output= paste0(output_name,"_all_level.csv")
  write.csv(output_level_graphs,csv_output, row.names = FALSE)
}

beta_div_removal_topedgeR_KW <- function(edger_results, #CSV file for EDGER results
                                         remove_top_N,
                                         input_phyloseq, 
                                         sample_type_var_name, 
                                         sample_types=list(), 
                                         sample_type_color=list(), 
                                         p_value=c("yes","no"),
                                         output_name, 
                                         width=200, height=200,
                                         plot_title_size=14, axis_title_size=18, label_size=7,
                                         x_axis_flip, 
                                         p_value_location=c("TR","TL","BR","BL")) { #p_value_location: TR-top right; TL- top left; BR- bottom right; BL- bottom left
  
  ###### this function is exclusively run after decontaminant_testing_KW
  ###### this function display 
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  if(missing(x_axis_flip)) { x_axis_flip<-"no"} #default being no for x_axis_flip
  if(missing(p_value_location) & tolower(p_value)=="yes"){
    p_value_location<-"BL" #default location for p_value will be in the bottom left corner
  } else if(!missing(p_value_location) & tolower(p_value)=="yes"){ #if p_value is wanted and p_value_location is not missing
    print("p_value_location: TR-top right; TL- top left; BR- bottom right; BL- bottom left")
    p_value_location<-match.arg(p_value_location)
    print(paste0("selected:",p_value_location))
  }
  #making sure the number of item in sample_types is same as sample_type_color
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  #making sure only yes and no is selected for p_value
  p_value <- match.arg(p_value)
  #setting up the coloring of the graph 
  sample_type_color_t<-sample_type_color
  names(sample_type_color_t) <-sample_types
  sample_type2<-sort(sample_types)
  
  #import edger results
  results<-read.csv(edger_results)
  
  # Function to extract last element with delimiter
  extractLastElement <- function(text, delimiter) {
    elements <- unlist(strsplit(text, delimiter))
    lastElement <- tail(elements, 1)
    return(lastElement)
  }
  
  results$ASV <- sapply(results$X, extractLastElement, delimiter = "\\.")
  results<-results[order(results$logFC, decreasing = FALSE),]
  new_results <- results[-c(1:remove_top_N), ] #remove the top N row
  
  new_results<-new_results[order(new_results$logFC, decreasing = TRUE),]
  new_results2 <- new_results[-c(1:remove_top_N), ] #remove the bottom N row
  
  remaining_ASV<-as.character(new_results2$ASV)
  
  remaining_ASVsdfsdf<<-remaining_ASV
  #number of remaining taxa
  remaining_ASV_num<-length(remaining_ASV)
  print(remaining_ASV_num)
  
  #Beta diversity
  #make sure input_phyloseq only include the sample_types of interest
  keep_sample<- get_variable(input_phyloseq,sample_type_var_name) %in% sample_types
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  phyloseq_t <- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
  
  phyloseq_tsdfsdfdsfsd<<-phyloseq_t
  
  phyloseq_t2 <- prune_taxa(remaining_ASV, phyloseq_t)
  
  phyloseq_t2_TESTING<<-phyloseq_t2
  
  #recalculate the relative abundance  
  normalizeSample <- function(x){x/sum(x)}
  phyloseq_t3 = transformSampleCounts(phyloseq_t2,normalizeSample) ##relative abundance
  
  #Create Distance Matrix with Bray
  vegdist=vegdist(t(otu_table(phyloseq_t3)), method = "bray") #row is samplesID #column is taxa
  #Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
  CmdScale <- cmdscale(vegdist, k =10)
  #calculated Sample variance for each PC
  vars <- apply(CmdScale, 2, var)
  #Create Variable with the Percent Variance
  percentVar <- round(100 * (vars/sum(vars)))
  #Merge PC Data with MetaData
  newResults <- merge(x = CmdScale, y = sample_data(phyloseq_t3), by = "row.names", all.x = TRUE)  #sample_data(phyloseq_t3) row is the sampleID
  #Rename Variables for PC1 and PC2
  colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
  colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
  colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"
  formula <- as.formula(paste("cbind(PC1,PC2)",sample_type_var_name, sep=" ~ "))
  #Calculate the Centroid Value
  centroids <- aggregate(formula ,data= newResults, mean)
  #Merge the Centroid Data into the PCOA Data
  newResults <- merge(newResults,centroids,by=sample_type_var_name,suffixes=c("",".centroid"))
  # Calculate p-value with ADONIS
  x=adonis2(vegdist ~ phyloseq_t3@sam_data[[sample_type_var_name]])
  pvalues<-x$`Pr(>F)`[[1]]
  
  subtitle_output<- paste0("Adonis, p=",pvalues)
  p<-ggplot(newResults, aes(PC1, PC2, color= get(sample_type_var_name))) + # Graph PC1 and PC2
    geom_point(size=2) + # Set the size of the points
    xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2 
    labs(title=NULL, 
         fill="Centroids") +
    ggtitle(paste0("removed highest/lowest LFC taxa:",remove_top_N))+
    geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= get(sample_type_var_name)))+ 
    geom_point(data=centroids, aes(x=PC1, y=PC2, color= get(sample_type_var_name)), size=0) + 
    geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=sample_type2), size=label_size) +
    scale_color_manual(name=NULL,  
                       values=sample_type_color_t) +
    theme(panel.background = element_blank(),
          panel.border=element_rect(fill=NA),
          panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
          panel.grid.minor = element_blank(),strip.background=element_blank(),
          plot.title=element_text(face="bold",hjust=0.5, size = plot_title_size), 
          plot.subtitle = element_text(hjust=0.5),
          axis.title=element_text(face="bold", size = axis_title_size),
          axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
          axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
          axis.ticks=element_blank(),
          plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
  p<-p
  newResults<<-newResults
  if (tolower(x_axis_flip)=="yes") {
    p<-p+scale_x_reverse()
  }
  if (p_value=="yes") { ###when the axis is flip, the annotation for the p value need to placed on max(newResults$PC1), instead of min(newResults$PC2))
    if (tolower(x_axis_flip)=="no"){
      if (p_value_location=="BR"){
        p <- p + annotate("text", x=max(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
      } else if (p_value_location=="BL"){
        p <- p + annotate("text", x=min(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
      } else if (p_value_location=="TR"){
        p <- p + annotate("text", x=max(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
      } else if (p_value_location=="TL"){
        p <- p + annotate("text", x=min(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
      }
    } else if (tolower(x_axis_flip)=="yes") {
      if (p_value_location=="BR"){
        p <- p + annotate("text", x=min(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
      } else if (p_value_location=="BL"){
        p <- p + annotate("text", x=max(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
      } else if (p_value_location=="TR"){
        p <- p + annotate("text", x=min(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
      } else if (p_value_location=="TL"){
        p <- p + annotate("text", x=max(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
      }
    }
  }
  png_output<-paste0(output_name,".png")
  png(file=png_output, res = 300, width=width, height=height , units='mm')
  show(p)
  dev.off()
  
}

DMM_removal_testing_KW <- function(contaminant_level_df,
                                        input_phyloseq, 
                                        sample_type_var_name, 
                                        sample_types=list(), 
                                        sample_type_color=list(), 
                                        p_value=c("yes","no"),
                                        output_name,
                                        uniqueID, 
                                        max_cluster_DMM=5, ...) { #p_value_location: TR-top right; TL- top left; BR- bottom right; BL- bottom left
  list2env(list(...), environment())
  ###### this function is exclusively run after decontaminant_testing_KW
  ###### this function display 
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  #making sure the number of item in sample_types is same as sample_type_color
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  #setting up the coloring of the graph 
  sample_type_color_t<-sample_type_color
  names(sample_type_color_t) <-sample_types
  sample_type2<-sort(sample_types)
  
  #Beta diversity
  #make sure input_phyloseq only include the sample_types of interest
  keep_sample<- get_variable(input_phyloseq,sample_type_var_name) %in% sample_types
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  phyloseq_t <- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
  
  # Function to extract last element with delimiter
  extractLastElement <- function(text, delimiter) {
    elements <- unlist(strsplit(text, delimiter))
    lastElement <- tail(elements, 1)
    return(lastElement)
  }
  
  #number_column should be equal to the number of levels which was tested 
  number_column<-ncol(contaminant_level_df)
  
  output_level_graphs <- data.frame(level = numeric(number_column),
                                    num_contaminant = numeric(number_column),
                                    p_value = numeric(number_column))
  num<-1
  combin_p<-list()
  level_setting<-colnames(contaminant_level_df)
  
  print(level_setting)
  for (i in level_setting){
    #the columns from the file should be level_###
    level<-as.numeric(sapply(i, extractLastElement, delimiter = "_"))
    print(level) # display level
    #input value: level for contaminant testing
    output_level_graphs[num,1]<-level
    #level will be the exact level when the contaminant testing was performed at 
    #keeping only one column at a time 
    taxa_contaminant_list<-contaminant_level_df[i]
    colnames(taxa_contaminant_list)<-"contaminant"
    taxa_contaminant_list<<-taxa_contaminant_list
    # Apply the function to the column
    taxa_contaminant_list$ASV <- sapply(row.names(taxa_contaminant_list), extractLastElement, delimiter = "\\.")
    #keeping only the noncontaminant
    taxa_contaminant_list_noncont<-taxa_contaminant_list %>% filter(contaminant==FALSE)
    NOT_contaminant_list<-taxa_contaminant_list_noncont$ASV
    noncontam_taxa_num<-length(NOT_contaminant_list)
    print(paste0("number of remaining taxa: ", noncontam_taxa_num))
    
    #input value: number of noncontaminant
    output_level_graphs[num,2]<-noncontam_taxa_num
    
    #number of contaminant revealed
    taxa_contaminant_list_cont<-taxa_contaminant_list %>% filter(contaminant==TRUE)
    contaminant_list<-taxa_contaminant_list_cont$ASV
    contam_taxa_num<-length(contaminant_list)
    
    #input value: number of contaminant
    output_level_graphs[num,3]<-contam_taxa_num
    
    #keeping only taxa which is NOT consider contaminant    
    phyloseq_t2 <- prune_taxa(NOT_contaminant_list, phyloseq_t)

    output_name2<-paste0(paste0(output_name,"_"),level)
    print(output_name2)
    DMM_cluster_KW(input_phyloseq=phyloseq_t2,
                   variable_to_compare=sample_type_var_name,
                   max_cluster=max_cluster_DMM,
                   phyloseq_uniqueID=uniqueID,
                   output_name=output_name2,                
                   compare_stats="fisher")
    num<-num+1
  }
  #########################
  #results without pruning#
  #########################
  output_noprune<-paste0(output_name,"_0")
  
  DMM_cluster_KW(input_phyloseq=phyloseq_t,
                 variable_to_compare=sample_type_var_name,
                 max_cluster=max_cluster_DMM,
                 phyloseq_uniqueID=uniqueID,
                 output_name=output_noprune,                
                 compare_stats="fisher")
  
  print(level_setting)
  #combining the cluster by level_num and turning it into a grid ggplot 
    num<-1
    for (i in level_setting ){
      level_num<-as.numeric(sapply(i, extractLastElement, delimiter = "_"))
      print(level_num)
      output_name3<-paste0(paste0(output_name,"_"),level_num)
      print(output_name3)
      object_num<-paste0("object",num)
      object<-get(paste0("DMM_cluster_KW_",output_name3)) #getting the DMM cluster assignment dataframe which is exported by the DMM_Cluster_KW function      
      colnames(object)<- c(i, "sample_ID")
      assign(object_num, object)
      num<-num+1
    }
    object_combin<-object1
    num2<-2 
    level_setting2<-level_setting[-1] #remove the first element with merging, since they are merging into a dataframe which is initially made with first element
    for (i in level_setting2 ){
      object_num2<-paste0("object",num2)
      object_combin<-merge(object_combin, get(object_num2), by="sample_ID")
      num2<-num2+1
    }
    
    object_noprune<-get(paste0("DMM_cluster_KW_",output_noprune)) #getting the DMM cluster assignment dataframe which is exported by the DMM_Cluster_KW function      
    colnames(object_noprune)<- c("level_0", "sample_ID")
    
    ####merging in the results without pruning
    object_combin<- merge(object_combin,object_noprune, by="sample_ID")

    #change to long format for ggplot
    object_combin_long <- reshape2::melt(object_combin, id.vars=c("sample_ID"))
    
    #remove the prefix level_
    object_combin_long$variable <- gsub("level_", "", object_combin_long$variable)
    object_combin_long<- object_combin_long %>% dplyr::rename(level=variable)
    
    #setting up the axis labeling
    list_of_levels<-unique(object_combin_long$level)
    list_of_levels<-replace(sort(list_of_levels),1, "unpruned") #the first element should be level 0 which is unpruned
    
    #make the cluster outcome a categorical variable
    object_combin_long$value<-as.factor(object_combin_long$value)
    object_combin_long<- object_combin_long %>% dplyr::rename(Cluster=value)
    
    #setting up clusters color
    cluster_color_t<-c("cornflowerblue","hotpink1","cyan4","goldenrod3",'azure4') #default color if color is not specified 
    names(cluster_color_t) <- c("1","2","3","4","5")

    combine_color_t<-c(cluster_color_t, sample_type_color_t)
    
    #export graphs- separate based on sample types
                IDclass<-data.frame(input_phyloseq@sam_data) #sample data from phyloseq
                IDclass <- IDclass %>% select(c(uniqueID,sample_type_var_name)) #keeping only unique sample ID and the sample type
                colnames(IDclass)<-c("sample_ID","org_sample_type") #renaming it so that it matches with object_combin_long
                
                object_combin_long2<- merge(object_combin_long, IDclass, by=("sample_ID"), all.x=T)
                #making sure the order of the org_sample_type is as assigned
                object_combin_long2$org_sample_type <- factor(object_combin_long2$org_sample_type, levels = sample_types)
                
                q<- object_combin_long2 %>%
                  mutate(sample_ID = fct_reorder(sample_ID, desc(org_sample_type))) %>%
                  ggplot( aes(x = level, y = sample_ID, fill=Cluster)) +
                  geom_tile(aes(color = org_sample_type,  width=1, height=0.9), size=1) +
                  scale_fill_manual(values = cluster_color_t)+
                  scale_color_manual(values = sample_type_color_t)+
                  coord_fixed(ratio = .3) +
                  scale_x_discrete(labels=list_of_levels)+
                  theme(panel.grid.minor=element_blank(),
                        panel.grid.major=element_blank(),
                        axis.ticks = element_blank(),
                        axis.title.y = element_blank(),
                        axis.text.y = element_text(size = 7,face="bold"),        
                        axis.text.x = element_text(size = 11))
                
    #export the graph- facet by sample types
    pdf_combine_output<-paste0(output_name,"_all_level_cluster.pdf")
    pdf(file=pdf_combine_output, width=12, height=20)
          show(q)
    dev.off()
                  q2<- object_combin_long2 %>%
                    mutate(sample_ID = fct_reorder(sample_ID, desc(org_sample_type))) %>%
                    ggplot( aes(x = level, y = sample_ID, fill=Cluster)) +
                    facet_grid(org_sample_type ~ .,scales = "free") +
                    geom_tile() +
                    scale_fill_manual(values = cluster_color_t)+
                    scale_x_discrete(labels=list_of_levels)+
                    theme(panel.grid.minor=element_blank(),
                          panel.grid.major=element_blank(),
                          axis.ticks = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y = element_text(size = 7,face="bold"),        
                          axis.text.x = element_text(size = 11))+
                    theme(strip.background =element_rect(colour="black"))+
                    theme(strip.text = element_text(size = 11, face="bold"))
                  
      #changing the color of the strips            
          # Find strips glob
          gt<-ggplot_gtable(ggplot_build(q2))
          strips <- which(startsWith(gt$layout$name,'strip'))
          
          # Change the fill color of each strip
          for (s in seq_along(strips)) {
            gt$grobs[[strips[s]]]$grobs[[1]]$children[[1]]$gp$fill <- sample_type_color_t[s]
          }
          
    pdf_combine_output2<-paste0(output_name,"_all_level_cluster_versionB.pdf")
    pdf(file=pdf_combine_output2, width=12, height=20)
          show(plot(gt))
    dev.off()
    
    #export the csv
    csv_combine_output<-paste0(output_name,"_all_level_cluster.csv")
    write.csv(object_combin,csv_combine_output, row.names = FALSE) 
    
}

dist2upper_removal_testing_KW <- function(contaminant_level_df,
                                        input_phyloseq, 
                                        sample_type_org, 
                                        sample_types=list(), 
                                        sample_type_color=list(), 
                                        sample_type_var_name,  
                                        Upper_sample_type, 
                                        Upper_sample_color,
                                        p_value=c("yes","no"),
                                        output_name, 
                                        uniqueID,
                                        width=200, height=200,
                                        plot_title_size=14, axis_title_size=18, label_size=7,
                                        x_axis_flip, 
                                        correlate_2_upper,
                                        beta_method=c("bray", "unifrac"),
                                        p_value_location=c("TR","TL","BR","BL")) { #p_value_location: TR-top right; TL- top left; BR- bottom right; BL- bottom left
  ###### this function is exclusively run after decontaminant_testing_KW
  ###### this function display 
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  if(missing(x_axis_flip)) { x_axis_flip<-"no"} #default being no for x_axis_flip
  if(missing(beta_method)) { beta_method<-"bray"} #default being bray curtis
  if(missing(p_value_location) & tolower(p_value)=="yes"){
    p_value_location<-"BL" #default location for p_value will be in the bottom left corner
  } else if(!missing(p_value_location) & tolower(p_value)=="yes"){ #if p_value is wanted and p_value_location is not missing
    print("p_value_location: TR-top right; TL- top left; BR- bottom right; BL- bottom left")
    p_value_location<-match.arg(p_value_location)
    print(paste0("selected:",p_value_location))
  }
  #making sure uniqueID is not missing. need to know how the sample is identify in the phyloseq
  stopifnot("uniqueID cannot be missing. Please provide uniqueID that identify each sample within the phyloseq"= !missing(uniqueID))
  
  #making sure the number of item in sample_types is same as sample_type_color
  stopifnot("sample_types need to have same number of elements as sample_type_color"= length(sample_types)==length(sample_type_color))
  #making sure only yes and no is selected for p_value
  p_value <- match.arg(p_value)
  #setting up the coloring of the graph 
  sample_type_color_t<-sample_type_color
  names(sample_type_color_t) <-sample_types
  sample_type2<-sort(sample_types)
  
  #Beta diversity
  #make sure input_phyloseq only include the sample_types of interest
  keep_sample<- get_variable(input_phyloseq,sample_type_org) %in% sample_types
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  phyloseq_t <- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
  
  # Function to extract last element with delimiter
  extractLastElement <- function(text, delimiter) {
    elements <- unlist(strsplit(text, delimiter))
    lastElement <- tail(elements, 1)
    return(lastElement)
  }
  
  contaminant_level_df_edit<-contaminant_level_df
  contaminant_level_df_edit$level_0<-FALSE
  contaminant_level_df_edit <- contaminant_level_df_edit[, c("level_0", names(contaminant_level_df_edit)[-which(names(contaminant_level_df_edit) == "level_0")])]
  
  #number_column should be equal to the number of levels which was tested 
  number_column<-ncol(contaminant_level_df_edit)
      
  
  num<-1
  for (i in colnames(contaminant_level_df_edit)){
    #the columns from the file should be level_###
    level<-as.numeric(sapply(i, extractLastElement, delimiter = "_"))
    print(level) # display level
    #level will be the exact level when the contaminant testing was performed at 
    #keeping only one column at a time 
    taxa_contaminant_list<-contaminant_level_df_edit[i]
    colnames(taxa_contaminant_list)<-"contaminant"
    taxa_contaminant_list<<-taxa_contaminant_list
    # Apply the function to the column
    taxa_contaminant_list$ASV <- sapply(row.names(taxa_contaminant_list), extractLastElement, delimiter = "\\.")
    #keeping only the noncontaminant
    taxa_contaminant_list_noncont<-taxa_contaminant_list %>% filter(contaminant==FALSE)
    NOT_contaminant_list<-taxa_contaminant_list_noncont$ASV
    noncontam_taxa_num<-length(NOT_contaminant_list)
    print(paste0("number of remaining taxa: ", noncontam_taxa_num))
    
    #number of contaminant revealed
    taxa_contaminant_list_cont<-taxa_contaminant_list %>% filter(contaminant==TRUE)
    contaminant_list<-taxa_contaminant_list_cont$ASV
    contam_taxa_num<-length(contaminant_list)
    
    #keeping only taxa which is NOT consider contaminant    
    phyloseq_t2 <- prune_taxa(NOT_contaminant_list, phyloseq_t)
    #recalculate the relative abundance  
    normalizeSample <- function(x){x/sum(x)}
    phyloseq_t3 = transformSampleCounts(phyloseq_t2,normalizeSample) ##relative abundance

    #distance to upper
    output_modified<-paste0(paste0(output_name,"_"),level)
    print(output_modified)
    
    ### Code with all samples: 
    #Create Distance Matrix with Bray (or wUniFrac depending what you are using)
    #Create Distance Matrix with Bray
    if (beta_method=="bray"){
      #Create Distance Matrix with Bray
      vegdist=vegdist(t(otu_table(phyloseq_t3)), method = "bray") #row is samplesID #column is taxa
    } else {
      set.seed(123)
      vegdist= phyloseq::UniFrac(phyloseq_t3, weighted = TRUE) 
    }

    #Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
    CmdScale <- cmdscale(vegdist, k =10)
    #calculated Sample variance for each PC
    vars <- apply(CmdScale, 2, var)
    #Create Variable with the Percent Variance
    percentVar <- round(100 * (vars/sum(vars)))
    #Merge PC Data with MetaData
    newResults <- merge(x = CmdScale, y = sample_data(phyloseq_t3), by = "row.names", all.x = TRUE)
    newResults$sampletype_KW<-newResults[[sample_type_var_name]]
    #Rename Variables for PC1 and PC2
    colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
    colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
    colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"
    
    # Calculate the Centroid Value
    centroids<-aggregate(cbind(PC1,PC2)~sampletype_KW,data=newResults, mean) #Here you would use your grouping variable (e.g., days)
    print("Centroids:")
    print(centroids)
    
    # Merge the Centroid Data into the PCOA Data
    newResults<-merge(newResults, centroids, by="sampletype_KW", suffixes=c("",".centroid"))
    
    #define function to calculate distance between distance vectors and the centroid. This function is taken from betadisper (obtained from Fares code)
    Resids <- function(x, c) {
      if(is.matrix(c))
        d <- x - c
      else
        d <- sweep(x, 2, c)
      rowSums(d^2)
    }
    i<-1
    #Figure out distance to other centroids 
    All_other_PCA <- newResults 
    row.names(All_other_PCA)<-All_other_PCA$name
    All_other_PCA_matrix<- All_other_PCA %>% dplyr::select(c("PC1","PC2"))
    All_other_PCA_matrix<-as.matrix(All_other_PCA_matrix)
    
    centroids_upper <- centroids %>% filter(sampletype_KW ==Upper_sample_type)
    centroids_upper_matrix<- centroids_upper %>% dplyr::select(c("PC1","PC2"))
    centroids_upper_matrix<- as.matrix(centroids_upper_matrix)
    
    #calculate the distance to Upper centroid for Upper and Lower samples 
    dist_to_centroids_test<-Resids(All_other_PCA_matrix, as.numeric(centroids_upper_matrix))
    dist_to_centroids_test<-as.data.frame(dist_to_centroids_test)
    dist_to_centroids_test$subset_sample<-All_other_PCA[[sample_type_org]]
    dist_to_centroids_test<-dist_to_centroids_test %>% dplyr::rename("distance"="dist_to_centroids_test")

    output_modified<-paste0(paste0(output_name,"_"),level)
    
    output_dataframe_name<- paste0("Dist_to_upper",output_modified)
    #export the distance to global environment for other use
    assign(output_dataframe_name,dist_to_centroids_test, envir = .GlobalEnv )
    print("###########################")
    print(paste0("Available dataframe :distance to upper",":"))
    print(output_dataframe_name)
    print("###########################")

    num<-num+1
  }
  #combining the different level 
  num<-1
  for (i in colnames(contaminant_level_df_edit)) {
    level<-as.numeric(sapply(i, extractLastElement, delimiter = "_"))
    output_modified<-paste0(paste0(output_name,"_"),level)
    if (num==1) { #set the first level to the final output and then merge in the subsequent levels
      combine_final_output<- get(paste0("Dist_to_upper",output_modified))
      combine_final_output = data.frame(subset(combine_final_output, select = c(distance)))
      colnames(combine_final_output)<-level
    } else{
      tempdf<-get(paste0("Dist_to_upper",output_modified))
      tempdf<-data.frame(subset(tempdf, select = c(distance)))
      colnames(tempdf)<-level
      combine_final_output<-cbind(combine_final_output,tempdf)
    }
    num<-num+1
  }
  
  combine_final_output2<-combine_final_output 
  colnames(combine_final_output2) <- paste("level" ,colnames(combine_final_output),sep="_")
  #exporting the combine distance from lower
  combine_final_output_name<-paste0("Dist_to_upper_combined",output_name)
  assign(combine_final_output_name, combine_final_output2,envir=.GlobalEnv)
  print("###########################")
  print(paste0("Available dataframe- combined :distance to upper",":"))
  print(combine_final_output_name)
  print("###########################")
  
  #############heatmap with distance to upper##############
  combine_final_output$sample_ID<-rownames(combine_final_output)
  combine_final_output_long <- reshape2::melt(combine_final_output, id.vars=c("sample_ID"))
  colnames(combine_final_output_long)<-c('sample_ID', "level", "distance_to_upper")
  
  IDclass<-data.frame(input_phyloseq@sam_data) #sample data from phyloseq
  IDclass <- IDclass %>% select(c(uniqueID,sample_type_org)) #keeping only unique sample ID and the sample type
  colnames(IDclass)<-c("sample_ID","org_sample_type") #renaming it so that it matches with object_combin_long
  
  combine_final_output_long2<- merge(combine_final_output_long, IDclass, by=("sample_ID"), all.x=T)
  
  list_of_levels<-unique(combine_final_output_long$level)
  list_of_levels<-as.character(list_of_levels)
  list_of_levels<-replace(sort(list_of_levels),1, "unpruned") #the first element should be level 0 which is unpruned
  
  combine_final_output_long2$distance_to_upper<-round(combine_final_output_long2$distance_to_upper, digits = 2)
  my_palette <- colorRampPalette(c("white", "darkgreen"))(n = length(unique(combine_final_output_long2$distance_to_upper)))
  
  qq<-combine_final_output_long2 %>%
    mutate(sample_ID = fct_reorder(sample_ID, desc(org_sample_type))) %>%
    ggplot( aes(x = level, y = sample_ID, fill=factor(distance_to_upper))) +
    facet_grid(org_sample_type ~ .,scales = "free") +
    geom_tile() +
    scale_fill_manual(values = my_palette)+
    scale_x_discrete(labels=list_of_levels)+
    labs(fill = "Distance to Upper")+
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 7,face="bold"),        
          axis.text.x = element_text(size = 11))+
    theme(strip.background =element_rect(colour="black"))+
    theme(strip.text = element_text(size = 11, face="bold"))
  
  #changing the color of the strips            
  # Find strips glob
  gt<-ggplot_gtable(ggplot_build(qq))
  strips <- which(startsWith(gt$layout$name,'strip'))
  
  # Change the fill color of each strip
  for (s in seq_along(strips)) {
    gt$grobs[[strips[s]]]$grobs[[1]]$children[[1]]$gp$fill <- sample_type_color_t[s]
  }

  #output pdf with multiple PCA plots based on level of contaminant testing  
  pdf_output<-paste0(paste0(output_name,"_distance2upper_alllevel"),".pdf")
  pdf(file=pdf_output, width=12, height=20)
        show(qq)
  dev.off()
  
  ###if want to correlate distance to upper with inflammatory markers or cell counts 
  if (!missing(correlate_2_upper)){
    #need to make sure correlate_2_upper is within 
    stopifnot("Please make sure correlate_2_upper exist in the phyloseq"= correlate_2_upper %in% colnames(sample_data(input_phyloseq)))
    
    correlate_lung.physeq<-sample_data(input_phyloseq)[,c(uniqueID,correlate_2_upper, sample_type_org)]
    correlate_lung.physeq[[correlate_2_upper]]<-as.numeric(correlate_lung.physeq[[correlate_2_upper]])
    correlate_lung.physeq<-correlate_lung.physeq[!is.na(correlate_lung.physeq[[correlate_2_upper]]),] #remove the missing lymph
    
    correlate_lung.physeq_comb<-merge(correlate_lung.physeq,combine_final_output2, by=0, all.x=T)
    rownames(correlate_lung.physeq_comb)<-correlate_lung.physeq_comb$Row.names
    correlate_lung.physeq_comb<-subset(correlate_lung.physeq_comb,select=-c(Row.names))

    ########################
    ###correlation graphs###
    ########################
    num<-1
    combin_p<-list()
    for (i in colnames(contaminant_level_df_edit)){
          #the columns from the file should be level_###
          level<-as.numeric(sapply(i, extractLastElement, delimiter = "_"))
          print(level) # display level
          
          #level will be the exact level when the contaminant testing was performed at 
          #keeping only one column at a time 
          taxa_contaminant_list<-contaminant_level_df_edit[i]
          colnames(taxa_contaminant_list)<-"contaminant"
          taxa_contaminant_list<<-taxa_contaminant_list
          # Apply the function to the column
          taxa_contaminant_list$ASV <- sapply(row.names(taxa_contaminant_list), extractLastElement, delimiter = "\\.")
          #keeping only the noncontaminant
          taxa_contaminant_list_noncont<-taxa_contaminant_list %>% filter(contaminant==FALSE)
          NOT_contaminant_list<-taxa_contaminant_list_noncont$ASV
          noncontam_taxa_num<-length(NOT_contaminant_list)

          correlate_lung.physeq_comb2<-correlate_lung.physeq_comb %>% select(i,correlate_2_upper, uniqueID, sample_type_org)
          correlate_lung.physeq_comb2<-correlate_lung.physeq_comb2[,c(i,correlate_2_upper, uniqueID, sample_type_org)]
          colnames(correlate_lung.physeq_comb2)<-c("distance_to_upper","variable","sample_ID","sample_type")
          
          pp<-ggplot(correlate_lung.physeq_comb2, aes(x=distance_to_upper, y=variable)) +
            geom_point(size=2, aes(colour=sample_type)) +  
            geom_smooth(method=lm)+ stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) +
            xlab("distance to upper") + 
            ylab(correlate_2_upper) + 
            scale_color_manual(values=sample_type_color_t)+
            ggtitle(paste0(level,paste0(paste0(" ( N=",noncontam_taxa_num),")")))+
            theme(panel.background = element_blank(),
                  panel.border=element_rect(fill=NA),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),strip.background=element_blank(),
                  plot.title=element_text(face="bold",hjust=0.5, size = plot_title_size), 
                  plot.subtitle = element_text(hjust=0.5),
                  axis.title=element_text(face="bold", size = axis_title_size),
                  axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
                  axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
                  axis.ticks=element_blank(),
                  plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
          combin_p[[num]]<-pp
          num<-num+1
    }
    combin_p2<<-combin_p
    #put the graphs next to each other. add a title and collec the legend- place it at bottom and horizontal
    combine_plot<-wrap_plots(combin_p2, nrow = ceiling(number_column/5)) + plot_layout(guides = "collect")  &  
      plot_annotation(title = paste0(paste0(output_name,": distance to upper correlates with "), correlate_2_upper),
                      theme= theme(title = element_text(color = "red", face = "bold", size = 24))) 
    
    #output pdf with multiple PCA plots based on level of contaminant testing  
    pdf_output<-paste0(paste0(paste0(output_name,"_correlates_upper_2_"),correlate_2_upper),".pdf")
    pdf(pdf_output, width=22, height=7*ceiling(number_column/5))
    show(combine_plot)
    dev.off()
    
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Dickson_graph6A_KW <- function(contaminant_level_df,
                                      input_phyloseq, 
                                      sample_type_var_name, 
                                      sample_types=list(), 
                                      variable_to_compare,
                                      prune_sample_types=list(), #this indicates which sample types you want to prune. if sample_types is the same as prune_sample_types then all the taxa which are consider contaminant would be removed. 
                                      output_name, 
                                      prune_method=c("remove","zero")){ #prune method- if remove then the taxa will be remove from the taxa_table. if zero
  ###### this function is exclusively run after decontaminant_testing_KW
  ###### this function display 
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  if(missing(prune_method)) { prune_method<-"remove"}#default being remove taxa that is in the prune
  #check to see if the variable_to_compare exist in the phyloseq for comparison 
  stopifnot("variable_to_compare does not exist in the input_phyloseq"= variable_to_compare %in% colnames(sample_data(input_phyloseq)))
  
  #Beta diversity
  #make sure input_phyloseq only include the sample_types of interest
  keep_sample<- get_variable(input_phyloseq,sample_type_var_name) %in% sample_types
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  phyloseq_t <- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
  
  # Function to extract last element with delimiter
  extractLastElement <- function(text, delimiter) {
    elements <- unlist(strsplit(text, delimiter))
    lastElement <- tail(elements, 1)
    return(lastElement)
  }
  
  #number_column should be equal to the number of levels which was tested 
  number_column<-ncol(contaminant_level_df)
  
  output_level_graphs <- data.frame(level = numeric(number_column),
                                    num_contaminant = numeric(number_column),
                                    num_noncontaminant = numeric(number_column),
                                    p_value = numeric(number_column),
                                    r_squared= numeric(number_column))
  num<-1
  combin_p<-list()
  for (i in colnames(contaminant_level_df)){
    #the columns from the file should be level_###
    level<-as.numeric(sapply(i, extractLastElement, delimiter = "_"))
    print(level) # display level
    #input value: level for contaminant testing
    output_level_graphs[num,1]<-level
    #level will be the exact level when the contaminant testing was performed at 
    #keeping only one column at a time 
    taxa_contaminant_list<-contaminant_level_df[i]
    colnames(taxa_contaminant_list)<-"contaminant"
    taxa_contaminant_list<<-taxa_contaminant_list
    # Apply the function to the column
    taxa_contaminant_list$ASV <- sapply(row.names(taxa_contaminant_list), extractLastElement, delimiter = "\\.")
    #keeping only the noncontaminant
    taxa_contaminant_list_noncont<-taxa_contaminant_list %>% filter(contaminant==FALSE)
    NOT_contaminant_list<-taxa_contaminant_list_noncont$ASV
    noncontam_taxa_num<-length(NOT_contaminant_list)
    print(paste0("number of remaining taxa: ", noncontam_taxa_num))
    #input value: number of noncontaminant
    output_level_graphs[num,2]<-noncontam_taxa_num
    
    #number of contaminant revealed
    taxa_contaminant_list_cont<-taxa_contaminant_list %>% filter(contaminant==TRUE)
    contaminant_list<-taxa_contaminant_list_cont$ASV
    contam_taxa_num<-length(contaminant_list)
    
    #input value: number of contaminant
    output_level_graphs[num,3]<-contam_taxa_num
    
    print(paste0("prune method: ", prune_method)) #display prune method
    
    #the pruning method matters because with unifrac, the phylogenic tree is considered when calculating the distance so if a taxa is prune, it would not be entered into the calculation. So
    #if you want to prune only a certain sample types, but not all, then you would have to replace the taxa with 0, rather than simply removing them from the phyloseq. 
    if (prune_method=="remove") { #if prune_method is remove then will just use the prune_taxa function which just remove the taxa from the phyloseq 
      #keeping only taxa which is NOT consider contaminant    
      phyloseq_t2 <- prune_taxa(NOT_contaminant_list, phyloseq_t)
      print("After pruning: ")
      print(phyloseq_t2)
    } else {
      #this section select on the row number in which the contaminant taxa is located and then replacing them with 0
      #it select the columns that has the sample type of interest
      taxa_indices <- which(rownames(otu_table(phyloseq_t)) %in% contaminant_list)
      #sample_data_df is the sample data of the phyloseq. row is each unique sample
      sample_data_df <- as.data.frame(sample_data(phyloseq_t))
      #figure out which sample ID has the sample type of interest
      rows_with_target_samples_types<-rownames(sample_data_df[sample_data_df[[sample_type_var_name]] %in% prune_sample_types, ])
      #use the sample IDs identified from previous step to see which column from the otu table to select 
      prune_sample_types_indices <- which(colnames(otu_table(phyloseq_t)) %in% rows_with_target_samples_types)
      phyloseq_t2<-phyloseq_t
      #now you change the value on the otu table from the phyloeq to be 0 
      for (colnum in prune_sample_types_indices) {
        otu_table(phyloseq_t2)[taxa_indices, colnum] <- 0 #in a particular sample, replace the selected taxa to 0
      }
    }
    #recalculate the relative abundance  
    normalizeSample <- function(x){x/sum(x)}
    phyloseq_t3 = transformSampleCounts(phyloseq_t2,normalizeSample) ##relative abundance
                              
                              #calculate shannon diversity
                              Shannon_diversity <- data.frame(vegan::diversity(otu_table(phyloseq_t2), index = "shannon", MARGIN = 2, base = exp(1)))
                              
                              sample_data_Dickson<-data.frame(sample_data(phyloseq_t2))
                              sample_data_Dickson<-sample_data_Dickson %>% select(c(variable_to_compare))
                              
                              Shannon_div_Dickson<-merge(Shannon_diversity, sample_data_Dickson, by=0)
                              colnames(Shannon_div_Dickson)<-c("SampleID", "Shannon.Diversity", "compare_var")
                              
                              # Calculate the p-value and R-squared
                              lm_model <- lm(compare_var ~ Shannon.Diversity, data = Shannon_div_Dickson)
                              p_value <- sprintf("%.4f", summary(lm_model)$coefficients[2, 4])  # P-value for the slope coefficient
                              r_squared <- sprintf("%.4f", summary(lm_model)$r.squared) 
                              
                              #input value: p value
                              output_level_graphs[num,4]<-p_value
                              #input value: r_squared
                              output_level_graphs[num,5]<-r_squared
                              Shannon_div_Dickson2<<-Shannon_div_Dickson
                              p<-ggplot(Shannon_div_Dickson, aes(x = Shannon.Diversity, y = compare_var)) +
                                geom_point(shape = 21, fill = "dodgerblue4", color = "black", size = 3) +
                                geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add best-fit line
                                labs(
                                  x = "Lung community diversity",
                                  y = "Lung concentration of IL-1α"
                                ) +
                                ggtitle(paste0(level,paste0(paste0(" ( Removed=",contam_taxa_num),")")))+
                                theme(panel.background = element_blank(),
                                  panel.grid=element_blank(),
                                  panel.border=element_blank(),
                                  axis.line = element_line()) +
                                # Add the p-value and R-squared to the plot
                                annotate("text",
                                  x = max(Shannon_div_Dickson$Shannon.Diversity),  # Adjust the x and y coordinates as needed
                                  y = max(Shannon_div_Dickson$compare_var, na.rm = TRUE),
                                  label = paste("P-value =", p_value),
                                  hjust = 1,vjust = 1,
                                  color = "black") +
                                annotate("text",
                                  x = max(Shannon_div_Dickson$Shannon.Diversity),  # Adjust the x and y coordinates as needed
                                  y = max(Shannon_div_Dickson$compare_var, na.rm = TRUE)-max(Shannon_div_Dickson$compare_var, na.rm = TRUE)/20,  # Adjust the vertical position
                                  label = paste("R-squared =", r_squared),
                                  hjust = 1,vjust = 1,
                                  color = "black")
    assign(paste0("p",num),p)
    combin_p[[num]]<-p
    num<-num+1
  }
  
  combin_p2<<-combin_p
  #put the graphs next to each other. add a title and collec the legend- place it at bottom and horizontal
  combine_plot<-wrap_plots(combin_p2, nrow = ceiling(number_column/5)) + plot_layout(guides = "collect")  &  
    plot_annotation(title = output_name,
                    theme= theme(title = element_text(color = "red", face = "bold", size = 24))) 
  
  #output pdf with multiple PCA plots based on level of contaminant testing  
  pdf_output<-paste0(paste0(output_name,"_level"),".pdf")
  pdf(pdf_output, width=22, height=7*ceiling(number_column/5))
  show(combine_plot)
  dev.off()
  
  colnames(output_level_graphs)<-c("level", "non_contaminant", "contaminant", "p_value", "r_squared")
  
  #output csv which contains how many taxa were identified as contaminant/noncontaminant, along with the p value for the PCA plots
  csv_output= paste0(output_name,"_all_level.csv")
  write.csv(output_level_graphs,csv_output, row.names = FALSE)
}


Dickson_graph6B_KW <- function(contaminant_level_df,
                               input_phyloseq, 
                               sample_type_var_name, 
                               sample_types=list(), 
                               variable_to_compare,
                               prune_sample_types=list(), #this indicates which sample types you want to prune. if sample_types is the same as prune_sample_types then all the taxa which are consider contaminant would be removed. 
                               prune_method=c("remove","zero"),
                               p_value=c("yes","no"),
                               output_name, 
                               width=200, height=200,
                               plot_title_size=14, axis_title_size=18, label_size=7,
                               x_axis_flip=c("yes","no"), 
                               beta_method=c("bray", "unifrac"),
                               p_value_location=c("TR","TL","BR","BL")){ #prune method- if remove then the taxa will be remove from the taxa_table. if zero
  ###### this function is exclusively run after decontaminant_testing_KW
  ###### this function display 
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(input_phyloseq) ){ input_phyloseq <- t(input_phyloseq)}
  if(missing(prune_method)) { prune_method<-"remove"}#default being remove taxa that is in the prune
  #check to see if the variable_to_compare exist in the phyloseq for comparison 
  stopifnot("variable_to_compare does not exist in the input_phyloseq"= variable_to_compare %in% colnames(sample_data(input_phyloseq)))
  if(missing(x_axis_flip)) { x_axis_flip<-"no"} #default being no for x_axis_flip
  x_axis_flip <- match.arg(x_axis_flip)
  if(missing(beta_method)) { beta_method<-"bray"} #default being bray curtis
  #making sure only yes and no is selected for p_value
  if(missing(p_value)) { p_value<-"yes"} #default being yes
  p_value <- match.arg(p_value)
  if(missing(p_value_location) & tolower(p_value)=="yes"){
    p_value_location<-"BL" #default location for p_value will be in the bottom left corner
  } else if(!missing(p_value_location) & tolower(p_value)=="yes"){ #if p_value is wanted and p_value_location is not missing
    print("p_value_location: TR-top right; TL- top left; BR- bottom right; BL- bottom left")
    p_value_location<-match.arg(p_value_location)
    print(paste0("selected:",p_value_location))
  }

  #Beta diversity
  #make sure input_phyloseq only include the sample_types of interest
  keep_sample<- get_variable(input_phyloseq,sample_type_var_name) %in% sample_types
  #prune_samples is used rather than subset_samples because subset_samples does not work well within a function 
  phyloseq_t <- prune_samples(keep_sample, input_phyloseq)  ###keeping only samples with the outcome variable of choice 
  
  # Function to extract last element with delimiter
  extractLastElement <- function(text, delimiter) {
    elements <- unlist(strsplit(text, delimiter))
    lastElement <- tail(elements, 1)
    return(lastElement)
  }
  
  #number_column should be equal to the number of levels which was tested 
  number_column<-ncol(contaminant_level_df)
  
  output_level_graphs <- data.frame(level = numeric(number_column),
                                    num_contaminant = numeric(number_column),
                                    num_noncontaminant = numeric(number_column),
                                    p_value = numeric(number_column),
                                    r_squared= numeric(number_column))
  num<-1
  combin_p<-list()
  for (i in colnames(contaminant_level_df)){
    #the columns from the file should be level_###
    level<-as.numeric(sapply(i, extractLastElement, delimiter = "_"))
    print(level) # display level
    #input value: level for contaminant testing
    output_level_graphs[num,1]<-level
    #level will be the exact level when the contaminant testing was performed at 
    #keeping only one column at a time 
    taxa_contaminant_list<-contaminant_level_df[i]
    colnames(taxa_contaminant_list)<-"contaminant"
    taxa_contaminant_list<<-taxa_contaminant_list
    # Apply the function to the column
    taxa_contaminant_list$ASV <- sapply(row.names(taxa_contaminant_list), extractLastElement, delimiter = "\\.")
    #keeping only the noncontaminant
    taxa_contaminant_list_noncont<-taxa_contaminant_list %>% filter(contaminant==FALSE)
    NOT_contaminant_list<-taxa_contaminant_list_noncont$ASV
    noncontam_taxa_num<-length(NOT_contaminant_list)
    print(paste0("number of remaining taxa: ", noncontam_taxa_num))
    #input value: number of noncontaminant
    output_level_graphs[num,2]<-noncontam_taxa_num
    
    #number of contaminant revealed
    taxa_contaminant_list_cont<-taxa_contaminant_list %>% filter(contaminant==TRUE)
    contaminant_list<-taxa_contaminant_list_cont$ASV
    contam_taxa_num<-length(contaminant_list)
    
    #input value: number of contaminant
    output_level_graphs[num,3]<-contam_taxa_num
    
    print(paste0("prune method: ", prune_method)) #display prune method
    
    #the pruning method matters because with unifrac, the phylogenic tree is considered when calculating the distance so if a taxa is prune, it would not be entered into the calculation. So
    #if you want to prune only a certain sample types, but not all, then you would have to replace the taxa with 0, rather than simply removing them from the phyloseq. 
    if (prune_method=="remove") { #if prune_method is remove then will just use the prune_taxa function which just remove the taxa from the phyloseq 
      #keeping only taxa which is NOT consider contaminant    
      phyloseq_t2 <- prune_taxa(NOT_contaminant_list, phyloseq_t)
      print("After pruning: ")
      print(phyloseq_t2)
    } else {
      #this section select on the row number in which the contaminant taxa is located and then replacing them with 0
      #it select the columns that has the sample type of interest
      taxa_indices <- which(rownames(otu_table(phyloseq_t)) %in% contaminant_list)
      #sample_data_df is the sample data of the phyloseq. row is each unique sample
      sample_data_df <- as.data.frame(sample_data(phyloseq_t))
      #figure out which sample ID has the sample type of interest
      rows_with_target_samples_types<-rownames(sample_data_df[sample_data_df[[sample_type_var_name]] %in% prune_sample_types, ])
      #use the sample IDs identified from previous step to see which column from the otu table to select 
      prune_sample_types_indices <- which(colnames(otu_table(phyloseq_t)) %in% rows_with_target_samples_types)
      phyloseq_t2<-phyloseq_t
      #now you change the value on the otu table from the phyloeq to be 0 
      for (colnum in prune_sample_types_indices) {
        otu_table(phyloseq_t2)[taxa_indices, colnum] <- 0 #in a particular sample, replace the selected taxa to 0
      }
    }
    
    #recalculate the relative abundance  
    normalizeSample <- function(x){x/sum(x)}
    phyloseq_t3 = transformSampleCounts(phyloseq_t2,normalizeSample) ##relative abundance

    #Create Distance Matrix with Bray
    if (beta_method=="bray"){
      #Create Distance Matrix with Bray
      vegdist=vegdist(t(otu_table(phyloseq_t3)), method = "bray", na.rm=TRUE) #row is samplesID #column is taxa
    } else {
      set.seed(123)
      vegdist= phyloseq::UniFrac(phyloseq_t3, weighted = TRUE) 
    }

    #Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
    CmdScale <- cmdscale(vegdist, k =10)
    #calculated Sample variance for each PC
    vars <- apply(CmdScale, 2, var)
    #Create Variable with the Percent Variance
    percentVar <- round(100 * (vars/sum(vars)))
    #Merge PC Data with MetaData
    newResults <- merge(x = CmdScale, y = sample_data(phyloseq_t3), by = "row.names", all.x = TRUE)  #sample_data(phyloseq_t3) row is the sampleID
    #Rename Variables for PC1 and PC2
    colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
    colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
    colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"
    
    phyloseq_t344<<-phyloseq_t3
    
    # Calculate p-value with ADONIS with the variable_to_compare 
    x=adonis2(vegdist ~ phyloseq_t3@sam_data[[variable_to_compare]],na.action = na.exclude)
    pvalues<-x$`Pr(>F)`[[1]]
    
    #input value: p value
    output_level_graphs[num,4]<-pvalues
    
    subtitle_output<- paste0("Adonis, p=",pvalues)
    p<-ggplot(newResults, aes(PC1, PC2, size = get(variable_to_compare) )) + # Graph PC1 and PC2
      geom_point(color= "dodgerblue4") + # Set the size of the points
      scale_size(range = c(0.1, 20))+
      xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2 
      labs(title=NULL, 
           fill="Centroids") +
      
      ggtitle(paste0(level,paste0(paste0(" ( Removed=",contam_taxa_num),")")))+
      theme(panel.background = element_blank(),
            panel.border=element_rect(fill=NA),
            panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
            panel.grid.minor = element_blank(),strip.background=element_blank(),
            plot.title=element_text(face="bold",hjust=0.5, size = plot_title_size), 
            plot.subtitle = element_text(hjust=0.5),
            axis.title=element_text(face="bold", size = axis_title_size),
            axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
            axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
            axis.ticks=element_blank(),
            plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
    p<-p
    newResults<<-newResults
    if (tolower(x_axis_flip)=="yes") {
      p<-p+scale_x_reverse()
    }
    if (p_value=="yes") { ###when the axis is flip, the annotation for the p value need to placed on max(newResults$PC1), instead of min(newResults$PC2))
      if (tolower(x_axis_flip)=="no"){
        if (p_value_location=="BR"){
          p <- p + annotate("text", x=max(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
        } else if (p_value_location=="BL"){
          p <- p + annotate("text", x=min(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
        } else if (p_value_location=="TR"){
          p <- p + annotate("text", x=max(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
        } else if (p_value_location=="TL"){
          p <- p + annotate("text", x=min(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
        }
      } else if (tolower(x_axis_flip)=="yes") {
        if (p_value_location=="BR"){
          p <- p + annotate("text", x=min(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
        } else if (p_value_location=="BL"){
          p <- p + annotate("text", x=max(newResults$PC1),y=min(newResults$PC2), vjust=0, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
        } else if (p_value_location=="TR"){
          p <- p + annotate("text", x=min(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=1, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
        } else if (p_value_location=="TL"){
          p <- p + annotate("text", x=max(newResults$PC1),y=max(newResults$PC2), vjust=1, hjust=0, label = subtitle_output, fontface = "italic", col = "black", size = label_size) ###adding the p value
        }
      }
    }
    
    assign(paste0("p",num),p)
    combin_p[[num]]<-p
    num<-num+1
  }
  
  combin_p2<<-combin_p
  #put the graphs next to each other. add a title and collec the legend- place it at bottom and horizontal
  combine_plot<-wrap_plots(combin_p2, nrow = ceiling(number_column/5)) + plot_layout(guides = "collect")  &  
    plot_annotation(title = output_name,
                    theme= theme(title = element_text(color = "red", face = "bold", size = 24))) 
  
  #output pdf with multiple PCA plots based on level of contaminant testing  
  pdf_output<-paste0(paste0(output_name,"_level"),".pdf")
  pdf(pdf_output, width=22, height=7*ceiling(number_column/5))
  show(combine_plot)
  dev.off()
  
  colnames(output_level_graphs)<-c("level", "non_contaminant", "contaminant", "p_value")
  
  #output csv which contains how many taxa were identified as contaminant/noncontaminant, along with the p value for the PCA plots
  csv_output= paste0(output_name,"_all_level.csv")
  write.csv(output_level_graphs,csv_output, row.names = FALSE)
}


    
    
    
    
    
    



