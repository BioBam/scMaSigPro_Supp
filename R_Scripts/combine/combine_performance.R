##################################
## Author: Priyansh Srivastava ###
## Email: spriyansh29@gmail.com ##
## Script: ScMaSigPro ############
##################################

# Prefix
base_string <- "../scMaSigPro_supp_data/"
base_string_2 <- ""
dirPath <- paste0(base_string, "analysis_public_data/")
tabPath <- paste0(base_string, "tables/")
helpScriptsDir <- paste0(base_string_2, "R_Scripts/helper_function/")
figPath <- paste0(base_string, "figures/")
figPath_hd <- paste0(figPath, "hd/")
figPath_lr <- paste0(figPath, "lr/")

# Load all the tables
all_tabs <- list.files(tabPath,
  pattern = "Performance.Table.tsv",
  full.names = TRUE
)

# Set the excel file
excelFile <- paste0(tabPath, "Additional_Table_1_All_Performance_Measures_Results.xlsx")
if (file.exists(excelFile)) {
  file.remove(excelFile)
}
write.xlsx(as.data.frame(matrix(data = NA)), excelFile, sheetName = "dummy")

for (filePath in all_tabs) {
  # load data
  data <- read.table(filePath, header = T, sep = "\t")

  # Set filename
  sheetName <- str_split_i(str_remove_all(str_split_i(filePath, pattern = "_Performance", i = 1), paste0(tabPath, "|/0")), pattern = "\\_", i = 2)

  # All Clusters
  workbook <- loadWorkbook(excelFile)
  newSheet <- createSheet(workbook, sheetName = sheetName)
  addDataFrame(data, newSheet, row.names = F)
  saveWorkbook(workbook, excelFile)
}
