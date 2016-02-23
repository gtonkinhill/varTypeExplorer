library(shiny)
library(data.table)
library(VennDiagram)
library(ggplot2)

colours5 <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")

ui <- fluidPage(
    navbarPage("Type Explorer", id="nav",
                           
      tabPanel("Inputs",                         
        fileInput(inputId = 'otuFile', label = 'Choose OTU Matrix File',
              accept=c('text/csv', 
                       'text/comma-separated-values,text/plain', 
                       '.csv')),
    
        fileInput(inputId = 'isolateFile', label = 'Choose Isolate CSV File',
              accept=c('text/csv', 
                       'text/comma-separated-values,text/plain', 
                       '.csv'))
      ),
  
      tabPanel("Otu Matrix",
        dataTableOutput('otuTable')
      ),
  
      tabPanel("Isolate Data",
        dataTableOutput('isolateTable')
      ),
      
      tabPanel("Venn Diagrams",
        fluidPage(
          sidebarLayout(
            sidebarPanel(
              uiOutput('vennSetCheckBox1'),
              uiOutput('vennSetCheckBox2')
            ),
            mainPanel(
                imageOutput('vennPlot'),
                plotOutput('frequencyPlot')
            )
          )
        )
      )
   )
)

server <- function(input, output, session){
  options(shiny.maxRequestSize=100*1024^2) 
  
  otuData <- reactive({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    
    inFile <- input$otuFile
    
    if (is.null(inFile))
      return(NULL)
    
    fread(inFile$datapath)
  })
  
  isolateData <- reactive({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    
    inFile <- input$isolateFile
    
    if (is.null(inFile))
      return(NULL)
    
    fread(inFile$datapath)
  })
  
  output$otuTable <- renderDataTable({
    otuData()
  })
  
  output$isolateTable <- renderDataTable({
    isolateData()
  })
  
  output$vennSetCheckBox1 <- renderUI({
    options <- colnames(isolateData())
    radioButtons(inputId = 'vennCheckBox1', label="Choose isolate variable of interest",
                       choices = options, selected=options[2])
  })
  
  output$vennSetCheckBox2 <- renderUI({
    options <- unique(isolateData()[[input$vennCheckBox1]])
    checkboxGroupInput(inputId = 'vennCheckBox2', label="Choose sets to compare",
                       choices = options, selected=options[1])
  })
  
  
  output$vennPlot <- renderImage({
    #Generate venn image
    sets = list()
    i <- 1
    for (set in input$vennCheckBox2){
      isolates <- isolateData()[isolateData()[[input$vennCheckBox1]]==set,][[1]]
      setCoverage <- rowSums(otuData()[,isolates, with = FALSE])
      sets[[i]] <- otuData()[[1]][setCoverage>0]
      i <- i+1
    }
    names(sets) <- input$vennCheckBox2

    # This file will be removed later by renderImage
    outfile <- tempfile(fileext='.png')
    
    venn.diagram(sets, filename=outfile, imagetype="png",
                 fill=colours5[1:length(sets)],
                 width=session$clientData$output_vennPlot_width,
                 height=session$clientData$output_vennPlot_height,
                 resolution=72*session$clientData$pixelratio)
    
    list(src = outfile,
         contentType = 'image/png',
         width = session$clientData$output_vennPlot_width,
         height = session$clientData$output_vennPlot_height,
         alt = "Problem loading the image!")
    
  }, deleteFile = TRUE)
  
  
  output$frequencyPlot <- renderPlot({
    isolates <- isolateData()[isolateData()[[input$vennCheckBox1]] %in% input$vennCheckBox2,][[1]]
    otuFreqs <- data.frame(otuOccurenceCounts=rowSums(otuData()[,isolates, with = FALSE]))
    
    gg <- ggplot(otuFreqs, aes(x=otuOccurenceCounts))
    gg <- gg + geom_histogram(binwidth=1)
    gg <- gg + theme_bw()
    gg <- gg + scale_y_sqrt()
    gg
    })
  
}

shinyApp(ui = ui, server = server)

