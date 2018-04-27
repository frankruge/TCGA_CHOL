library(shiny)
library(ggplot2)
library(ggsignif)
library(grid)
#setwd("TCGA_CHOL/")
all = read.table("TCGA_CHOL_DATA.tsv", header = T); head(all)
tissueType = as.vector(read.table("tissueType.tsv", header = F, sep = '\t')$V1); tissueType

get_gene = function(df, gene, group_factor, whisker_width = 0.15, outlier_color = "#FFFFFF00", dp_col="#00000001", jitter=T, jitpos=0.07, dp_size=4){
  data = as.data.frame(unlist(df[df$hgnc_symbol == gene,], use.names = F)[-c(1,2)])  # the grep version below is inferior and causes problems when i.e. MYB and MYB-AS1. the \\b 
  #data = as.data.frame(unlist(df[grep(paste('\\',gene,"\\b", sep=""), df$hgnc_symbol),], use.names = F)[-c(1,2)])   #doesn't handle the hyphen
  names(data) = "FPKM" 
  data$FPKM=as.numeric(as.character(data$FPKM))
  data$group=as.factor(group_factor)
  p=ggplot(data, aes(x=group, y=FPKM, group=group)) +
    ggtitle(gene) +
    theme(plot.title = element_text(face="bold", hjust=0.5) ) +
    stat_boxplot(geom ='errorbar', width=whisker_width) + 
    geom_boxplot(outlier.color = outlier_color ) +                                                   
    geom_signif(comparisons = list(c("Primary Tumor","Solid Tissue Normal")),
                map_signif_level = TRUE,textsize = 7, test='t.test') 
  if(jitter==T){p = p + geom_jitter( position=position_jitter(jitpos), color=dp_col, size=dp_size)}
  p
  return(p)
}

############################################################################################
#debug
options(shiny.sanitize.errors = TRUE)
#test the function
#get_gene(df = all, gene = "STAT5A", group_factor = tissueType, jitter = T,dp_col = "#22234221")
# Define UI for application that draws a histogram
ui <- pageWithSidebar(
  headerPanel(textOutput("gene_name")),
  sidebarPanel(
    selectizeInput(inputId = "gn", label="Gene Name", choices = sort(unique(all$hgnc_symbol)), multiple = FALSE, options = NULL),
    #selectInput(inputId = "dset", label = "Select set of genes", choices = c("Downregulated" = "down", "Upregulated" = "up")),
    selectInput(inputId = "dp_color", label = "data points color", choices = c("#22234221",colors()), selected = "#22234221"), #c("lightgrey", "grey", "dodgerblue", "black")),
    sliderInput(inputId = "dp_size", label = "data point size", min = 0, max = 5,  value = 3, step = 0.1),
    sliderInput(inputId = "whisker", label = "whisker", min = 0, max = 1,  value = 0.15, step = 0.01),
    sliderInput(inputId = "jitpos", label = "data point jitter", min = 0, max = 0.4 ,  value = 0.07, step = 0.01),
    checkboxInput("jitter", label = "jitter", value = TRUE)
  ),
  mainPanel(
    # Use imageOutput to place the image on the page
    imageOutput("myImage")
    
  )
)

server <- function(input, output, session){
  #output$gene_name = renderText({paste("GO terms with FDR < ", round(10^-input$logFDR, 6))})
  output$myImage <- renderImage({
    ct=rgb(red = col2rgb(input$dp_color)[1], green = col2rgb(input$dp_color)[2], blue = col2rgb(input$dp_color)[3], maxColorValue = 255)
    ###print(ct)
    # A temp file to save the output.
    # This file will be removed later by renderImage
    outfile <- tempfile(fileext = '.png')
    
    # Generate the PNG
    #png(outfile, width = 500, height = ifelse(test = input$logFDR >=1, (1/input$logFDR)*800, 800))
    png(outfile, width = 600, height = 600)
    grid.draw(get_gene(df = all, gene = input$gn , group_factor = tissueType, dp_col = input$dp_color, dp_size=input$dp_size,
                       whisker_width = input$whisker, jitter = input$jitter, jitpos = input$jitpos))
    dev.off()
    # Return a list containing the filename
    list(src = outfile,
         contentType = 'image/png',
         width = 600,
         height = 600,
         alt = "This is alternate text")
  }, deleteFile = TRUE)
}

shinyApp(ui = ui, server = server)