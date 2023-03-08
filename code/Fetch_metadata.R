## Extract extra info table
# Guillermo Reales
# Last modified: 2023-03-07
# Download and modify our main table for use in our exploration

library(data.table)
library(googlesheets4)
# Specify the user:
options(gargle_oauth_email = "grealesm@gmail.com")
gs4_deauth()

message("Fetching metadata...")
metadata <- read_sheet("https://docs.google.com/spreadsheets/d/16B4ANehcS4psdAFReTBQUJLYXuf5_RrpmjSPaASa2Nw/edit?usp=sharing", sheet = 1)
metadata <- metadata[,c("File_ID", "First_Author", "Year", "Reference","Trait_ID_2.0", "Trait_long", "Trait_class", "N0", "N1", "N", "URL", "Chip", "Population", "Collection", "Retrieval_date", "Public", "Notes")]
metadata <- metadata[!is.na(metadata$File_ID),]
metadata$Trait <- sapply(strsplit(metadata$File_ID,"-"),`[`, 1)
metadata <- data.table(metadata)

date <- format(Sys.time(), format="%Y%m%d")
version  <- 1
metadataname  <- paste("Metadata_", date, "-v",version, ".tsv", sep="")


while(metadataname %in% dir()){
  version  <- version + 1
  metadataname  <- paste("Metadata_", date, "-v",version, ".tsv", sep="")
}

fwrite(metadata, paste0("../data/", metadataname), sep = "\t", na = NA)
message("Done!")
