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
                tableOutput('summaryTable'),
                imageOutput('vennPlot'),
                plotOutput('frequencyPlot'),
                plotOutput('cummulativePLot')
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
  
  isolates <- reactive({isolateData()[isolateData()[[input$vennCheckBox1]] %in% input$vennCheckBox2,][[1]]})
  
  sets <- reactive({
    sets = list()
    i <- 1
    for (set in input$vennCheckBox2){
      isolates <- isolateData()[isolateData()[[input$vennCheckBox1]]==set,][[1]]
      setCoverage <- rowSums(otuData()[,isolates, with = FALSE])
      sets[[i]] <- otuData()[[1]][setCoverage>0]
      i <- i+1
    }
    names(sets) <- input$vennCheckBox2
    sets
  })
  
  otuFreqs <- reactive ({data.frame(otuOccurenceCounts=rowSums(otuData()[,isolates(), with = FALSE]))})
  
  
  output$summaryTable <- renderTable({
    
    print(sets())
    print(length(Reduce(intersect, sets())))
    inter <- length(Reduce(intersect, sets()))
    
    relate <- inter/mean(as.numeric(lapply(sets(), length)))
    
    data.frame(Total_Isolates=length(isolates()), Total_DBLa_Types=length(otuFreqs()),
               Largest_OTU=max(otuFreqs()), relatedness=relate)
  })
  
  output$vennPlot <- renderImage({
    #Generate venn image
    sets <- sets()

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
    isolates <- isolates()
    otuFreqs <- otuFreqs()
    
    gg <- ggplot(otuFreqs, aes(x=otuOccurenceCounts))
    gg <- gg + geom_histogram(binwidth=1)
    gg <- gg + theme_bw()
    gg <- gg + scale_y_sqrt()
    gg
    })
  
  output$cummulativePLot <- renderPlot({
    isolates <- isolates()
    otus <- otuData()[,isolates, with = FALSE]
    cummulative <- data.frame(number=1:ncol(otus), cummulative=colSums(t(apply(otus, 1, cumsum))>0))
  
    gg <- ggplot(data=cummulative, aes(x=number, y=cummulative))
    gg <- gg + geom_line()
    gg <- gg + geom_point()
    gg <- gg + theme_bw() 
    gg <- gg + xlab("Number of isolates") + ylab("Cumulative total number of DBLalpha types")
    gg
    })
  
}

shinyApp(ui = ui, server = server)

