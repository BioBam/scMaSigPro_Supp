##################################
## Author: Priyansh Srivastava ###
## Email: spriyansh29@gmail.com ##
## Script: ScMaSigPro ############
##################################

# Load all the tables
all_tabs <- list.files("/supp_data/Tables",
  pattern = "Performance",
  full.names = TRUE
)

# Set xlsx
excelFile <- "/supp_data/Tables/Additional_Table_1_All_Performance_Measures_Results.xlsx"
file.remove(excelFile)
write.xlsx(as.data.frame(matrix(data = NA)), excelFile, sheetName = "dummy")

for (filePath in all_tabs) {
  # load data
  data <- read.table(filePath, header = T, sep = "\t")

  # Set filename
  sheetName <- str_remove(str_split_i(filePath, pattern = "\\.", i = 1), "Tables/0")

  # All Clusters
  workbook <- loadWorkbook(excelFile)
  newSheet <- createSheet(workbook, sheetName = sheetName)
  addDataFrame(data, newSheet, row.names = F)
  saveWorkbook(workbook, excelFile)
}
