###################### Alternative assignment .TXT to .FASTQ ####################### 
# Author: Anne Manders                                                             #
# Date: 06/07/2022 - 07/07/2022                                                    #
# Description: Converts DNA sequence of the allele HLA-A stored in .txt format     #
#   into two fastq format files (R1 and R2).                                       #
# Time spent: 6/7(20:30-22:30), 7/7(8:00-11:30 (functional by 11:00), 20:00-21:30  #
#                                       (finishing touches and upload to GitHub))  #
####################################################################################

## Install and load packages
BiocManager::install("readr")
library("readr") # v 2.1.2
BiocManager::install("devtools")
library("devtools") # v 2.4.3
BiocManager::install("ShortRead")
library("ShortRead") # v 1.52.0
# installed ShortRead to test if the created files were readable by read-write 
#   fastq packages.


#'Creates two files (R1 and R2) in fastq format using generated R1 and R2 reads.
#'1) Run the entire script to load the packages and functions.
#'2) In the console, execute create_fastq() with the chosen values for the parameters 
#'   overlap and read_length.
#'
#'@usage create_fastq(overlap,read_length)
#'@param overlap The numbers of bp the reads should overlap.
#'@param read_length The length  the reads should be (in bp).
#'@examples 
#'create_fastq(overlap = 30,read_length = 151).Takes less than a second.
create_fastq = function(overlap, read_length){
  # Load the genomic file in txt format (delim = new line).
  # Allele name is column name and genomic sequence is value 
  input = read_delim("./Input/HLA-ASequence.txt", delim = "\n", show_col_types = FALSE)
  
  # Get allele name
  allele = names(input)
  
  # Get genomic sequence
  sequence = input[[1,1]]
  
  # Create the R1 and R2 reads and store them in corresponding lists
  lists = create_R1_R2(sequence, overlap, read_length)

  R1_list = lists$R1_list
  R2_list = lists$R2_list
  
  ## FASTQ FORMAT
  # sequence identifier: @ <allele> <R1/R2> <nr> \n
  # sequence: <R1/R2_list[i]> \n
  # separator: + \n
  # quality scores: (replaced by read length) <length(R1/R2_list[i])> \n
  
  ## Write R1 to fastq file
  # Create filename and dirname, to combine into output location
  filename = paste("R1_o", overlap, "_l", read_length, ".fastq", sep = "")
  dirname = paste("./Output/")
  
  R1_output_file = paste(dirname, filename, sep = "") # "./Output/R1_o30_l151.fastq"
  
  # If file exists and is not empty, delete contents
  if (file.exists(R1_output_file)){
    if (file.size(R1_output_file) > 0 ){
      close( file( R1_output_file, open="w" ) )
    }
  }
  # for each R1 read
  for (i in 1:length(R1_list)){
    # Create the entries' elements
    id=paste("@",gsub(">","", allele), " R1 ", i,sep = "")
    # sequence is R1_list[i]
    # separator is "+"
    quality_score = paste("Read length (instead of quality score): ", nchar(R1_list[i]), sep = "") 
      # could be hardcoded "A","E" etc for all, but I chose to show the read length here
    
    # Create the full entry
    entry = paste (id, R1_list[i], "+", quality_score, sep = "\n")
    
    # Append entry to fastq file. 
    cat(entry, file = R1_output_file, append=TRUE, sep = "\n")
      # I've looked into existing packages for writing fastq files, 
      # but due to time constraints, I decided to create the files manually
  
  }
  
  ## Write R2 to fastq file. See comments above.
  filename = paste("R2_o", overlap, "_l", read_length, ".fastq", sep = "")
  dirname = paste("./Output/")
  
  R2_output_file = paste(dirname, filename, sep = "") # "./Output/R2_o30_l151.fastq"
  
  if (file.exists(R2_output_file)){
    if (file.size(R2_output_file) > 0 ){
      close( file( R2_output_file, open="w" ) )
    }
  }
  
  for (i in 1:length(R2_list)){
    id=paste("@",gsub(">","", allele), " R2 ", i,sep = "")
    # sequence is R2_list[i]
    # separator is "+"
    quality_score = paste("Read length (instead of quality score): ", nchar(R2_list[i]), sep = "") 

    entry = paste (id, R2_list[i], "+", quality_score, sep = "\n")
    
    cat(entry, file = R2_output_file, append=TRUE, sep = "\n")    
  }
}


#'Create corresponding lists of R1 and R2 reads.
#'R1 is created from the genomic file. R2 is the reverse strands of R1.
#'The length as well as the overlap of the reads are determined by user input.
#'
#'@usage create_R1_R2(sequence, overlap, read_length)
#'@param sequence The genomic sequence of the HLA-A allele.
#'@param overlap The numbers of bp the reads should overlap.
#'@param read_length The length in bp the reads should be.
#'@return Two lists containing the created R1 and R2 reads, respectively
create_R1_R2 = function(sequence, overlap, read_length){

# Create empty lists
R1_reads = list()
R2_reads = list()

# Calculate the length of the 'frame' based on the users' input
frame = read_length - overlap

count = 0
# for the # of times the frame occurs in the sequence
for (i in 1:ceiling(nchar(sequence)/frame)) {
  ## R1
  # Adjust start position
  if (count == 0)
    {stop = count + read_length} 
  else 
    {stop = count + (read_length-1)}
  
  # Get R1 reads in fragments of read_length bp
  R1 = substring(sequence, count, stop)
  
  # Add R1 read to R1 list
  R1_reads[[length(R1_reads) + 1]] = R1
  
  
  ## R2
  # Create R2 by reversing R1
    # Split R1
    split_R1 = strsplit(R1, NULL)[[1]]
    # Reverse and paste
    R2 = paste(rev(split_R1), collapse = "")
  
  # Add R2 read to R2 list  
  R2_reads[[length(R2_reads) + 1]] = R2

  ## Move onto the next frame, keeping the # of bp overlap
  if (count == 0) 
    {count = count + (frame+1)}
  else 
    {count = count + frame}
}

return_list = list(R1_list=R1_reads, R2_list=R2_reads)
return(return_list)

}

# ## Test read-write fastq packages
# R1_fastq = readFastq("./Output/", pattern = "R1")
# # Didn't work, because the quality score is not as it should be.
