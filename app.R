library(shiny)
library(data.table)
library(VennDiagram)

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
        uiOutput('vennSetCheckBox1'),
        uiOutput('vennSetCheckBox2'),
        imageOutput('vennPlot')
      )
   )
)

server <- function(input, output){
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
  
  output$vennSetCheckBox2 <- renderUI({
    options <- unique(isolateData()[[2]])
    checkboxGroupInput(inputId = 'vennCheckBox', label="Choose sets to compare",
                       choices = options)
  })
  
  output$vennSetCheckBox2 <- renderUI({
    options <- unique(isolateData()[[2]])
    checkboxGroupInput(inputId = 'vennCheckBox', label="Choose sets to compare",
                       choices = options)
  })
  
  
  output$vennPlot <- renderImage({
    #Generate venn image
    sets = list()
    i <- 1
    for (set in input$vennCheckBox){
      isolates <- isolateData()[isolateData()[[2]]==set,][[1]]
      setCoverage <- rowSums(otuData()[,isolates, with = FALSE])
      sets[[i]] <- otuData()[[1]][setCoverage>0]
      i <- i+1
    }
    names(sets) <- input$vennCheckBox

    # This file will be removed later by renderImage
    outfile <- tempfile(fileext='.png')
    
    venn.diagram(sets, filename=outfile, imagetype="png")
    list(src = outfile,
         contentType = 'image/png',
         width = 400,
         height = 300,
         alt = "This is alternate text")
  }, deleteFile = TRUE)
}

shinyApp(ui = ui, server = server)

