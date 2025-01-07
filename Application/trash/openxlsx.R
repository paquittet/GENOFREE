
# LIEN UTILE : https://ycphs.github.io/openxlsx/articles/Introduction.html

# data <- readColumns(sheet,
#                     startColumn = 1,
#                     endColumn = 17,
#                     startRow = 9,
#                     endRow = NULL,
#                     header = T,
#                     as.data.frame = T)


library(openxlsx)

# Loading workbook
wb <- openxlsx::loadWorkbook(file = "C:/Users/quittet/Documents/marmot-quantitative-genetics/PEDIGREE_GENETIQUE/GENOSCREEN/prout.xlsx")

# Read workbook sheet as data.frame
data <- openxlsx::readWorkbook(xlsxFile = wb,
                       sheet = 2,  
                       startRow = 9,
                       colNames = T)

# Filter by .fsa
data_filtered <- data[data$Référence.du.fichier %in% "D02_1953_010.fsa",]

# Get rows for reconstruction
row.index <- rownames(data_filtered)


# random mdofiication
data_filtered[1, 2] <- "PROUT"

# Add modification to data
data[row.index, ] <- data_filtered


# Define a custom font style
font_style <- openxlsx::createStyle(
  fontColour = "#0000ff",  # Set font color (blue)
  textDecoration = "bold",  # Bold text
  halign = "center",
  valign = "center",
  border = "TopBottomLeftRight"
)


for (i in row.index) {
  # Write the data
  openxlsx::writeData(
    wb = wb,
    sheet = 2,
    x = t(unlist(data[i, ])),  # Write row as a horizontal vector
    startRow = i,
    colNames = FALSE
  )
  
  # Apply the font style to the written data
  openxlsx::addStyle(
    wb = wb,
    sheet = 2,
    style = font_style,
    rows = i,
    cols = 1:ncol(data),  # Apply the style to all columns in the row
    gridExpand = TRUE     # Expand the style across all specified cells
  )
}



# Save work book
openxlsx::saveWorkbook(wb,
             file = "C:/Users/quittet/Documents/marmot-quantitative-genetics/PEDIGREE_GENETIQUE/GENOSCREEN/prout2.xlsx", 
             overwrite = T)



