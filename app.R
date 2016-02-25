library(shiny)
library(shinythemes)
library(data.table)
library(VennDiagram)
library(ggplot2)
library(dplyr)

library(leaflet)
library(RColorBrewer)
library(scales)
library(lattice)

colours5 <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")

ui <- fluidPage(theme = shinytheme("flatly"),
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
      
      tabPanel("Diagrams",
        fluidPage(
          sidebarLayout(
            sidebarPanel(
              uiOutput('vennSetCheckBox1'),
              uiOutput('vennSetCheckBox2'),
              radioButtons(inputId="graphType", label="Choose chart type", 
                           choices=c(
                             "Venn Diagram"="venn",
                             "Frequency Plot"="freq",
                             "Cummulative Diversity Plot"="cumDiv"),
                           selected=character(0))
            ),
            mainPanel(width = 8,
                tableOutput('summaryTable'),
                conditionalPanel(
                  condition = "input.graphType == 'venn'", imageOutput('vennPlot', width="100%")),
                conditionalPanel(
                  condition = "input.graphType == 'freq'", plotOutput('frequencyPlot', width="100%")),
                conditionalPanel(
                  condition = "input.graphType == 'cumDiv'", plotOutput('cummulativePLot', width="100%"))
                
            )
          )
        )
      ),
      
      tabPanel("Interactive Map",
        div(class="outer",
                   
          tags$head(
            # Include our custom CSS
            includeCSS("styles.css"),
            includeScript("gomap.js")
          ),
                   
          leafletOutput("map", width="100%", height="100%"),

          # Shiny versions prior to 0.11 should use class="modal" instead.
          conditionalPanel(
            condition = "clickData.clickedMarker != NULL",
              absolutePanel(id = "controls", class = "panel panel-default", fixed = TRUE,
                        draggable = TRUE, top = 60, left = "auto", right = 20, bottom = "auto",
                        width = 330, height = "auto",
                        
                        h2("Summary plots"),
    
                        plotOutput("sideBarFreqPlot", height = 200),
                        plotOutput("sideBarCumPlot", height = 250)
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
    
    inter <- length(Reduce(intersect, sets()))
    
    relate <- inter/mean(as.numeric(lapply(sets(), length)))
    
    data.frame(Total_Isolates=length(isolates()), Total_DBLa_Types=sum(otuFreqs()$otuOccurenceCounts>0),
               Largest_OTU=max(otuFreqs()), relatedness=relate)
  })
  
  output$vennPlot <- renderImage({
    #Generate venn image
    sets <- sets()

    # This file will be removed later by renderImage
    outfile <- tempfile(fileext='.png')
    
    venn.diagram(sets, filename=outfile, imagetype="png",
                 fill=colours5[1:length(sets)],
                 cex=1.5,
                 cat.cex=1.5,
                 margin=0.1,
                 fontface = "bold",
                 width=800, #session$clientData$output_vennPlot_width,
                 height=600, #session$clientData$output_vennPlot_height,
                 resolution=72*session$clientData$pixelratio)
    
    list(src = outfile,
         contentType = 'image/png',
         width = 800, #session$clientData$output_vennPlot_width,
         height = 600, #session$clientData$output_vennPlot_height,
         alt = "Problem loading the image!")
    
  }, deleteFile = TRUE)
  
  
  output$frequencyPlot <- renderPlot({
    isolates <- isolates()
    otuFreqs <- otuFreqs()
    
    gg <- ggplot(otuFreqs, aes(x=otuOccurenceCounts))
    gg <- gg + geom_histogram(binwidth=1)
    gg <- gg + theme_bw()
    gg <- gg + scale_y_sqrt()
    gg <- gg + xlab("Number of isolates the DBLa type is observed in")
    gg <- gg + ylab("Count of DBLa types")
    gg
    }, height = 600, width = 800)
  
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
    }, height = 600, width = 800)
  
  ## Interactive Map ###########################################
  
  locationData <- reactive({
    isolates <- isolateData()
    otus <- otuData()
    
    otuMatrix <- as.matrix(otus[, 2:ncol(otus), with=FALSE])
    rownames(otuMatrix) <-  otus[[1]]
    otuMatrix[otuMatrix>0] <- 1
    isoTotals <- apply(otuMatrix, 2, function(m){rownames(otuMatrix)[m>0]})
    
    isolates$Otus <- isoTotals[isolates[[1]]]
    
    isolates <- data.frame(isolates)
    
    locationData <- isolates %>%
      group_by(Location) %>%
      summarize(
        TotalIsolates = n(),
        TotalTypes = length(unique(unlist(Otus))),
        Dates = paste(unique(Date), collapse = ''),
        Longitude = unique(Longitude),
        Latitude = unique(Latitude)
      )
    
    locationData    
  })
  
  # Create the map
  output$map <- renderLeaflet({
    leaflet() %>%
      addTiles(
        urlTemplate = "//{s}.tiles.mapbox.com/v3/jcheng.map-5ebohr46/{z}/{x}/{y}.png",
        attribution = 'Maps by <a href="http://www.mapbox.com/">Mapbox</a>'
      ) %>%
      setView(lng = -93.85, lat = 37.45, zoom = 4)
  })
  
  observeEvent(input$map_bounds,{
        
        if (!is.null(input$otuFile) & !is.null(input$isolateFile)){

        #set circle radius to reflect the number of isolates
        radius <- locationData()$TotalIsolates / max(locationData()$TotalIsolates) * 150000
        
        #set circle colour to reflect the diversity at particular location
        colorData <- locationData()$TotalTypes/locationData()$TotalIsolates
        pal <- colorBin("Spectral", colorData, 7, pretty = FALSE)
#         
        leafletProxy("map", data = locationData()) %>%
          clearShapes() %>%
          addCircles(~Longitude, ~Latitude, radius=radius, layerId=~Location,
                     stroke=FALSE, fillOpacity=0.4, fillColor=pal(colorData)) %>%
          addLegend("bottomleft", pal=pal, values=colorData, title="Diversity measure #OTUs/#Isolates",
                    layerId="colorLegend")

        }

  })
  
  # A reactive expression that returns the set of zips that are
  # in bounds right now
#   locationsInBounds <- reactive({
#     if (is.null(input$map_bounds))
#       return(locationData()[FALSE,])
#     bounds <- input$map_bounds
#     latRng <- range(bounds$north, bounds$south)
#     lngRng <- range(bounds$east, bounds$west)
#     
#     subset(locationData(),
#            Latitude >= latRng[1] & Latitude <= latRng[2] &
#              Longitude >= lngRng[1] & Longitude <= lngRng[2])
#   })

  output$sideBarFreqPlot <- renderPlot({
    event <- input$map_shape_click
    if (is.null(event) || is.null(clickData$clickedMarker))
      return()
    
    isolates <- isolateData()[isolateData()[['Location']] %in% event$id,][[1]]
    otuFreqs <- data.frame(otuOccurenceCounts=rowSums(otuData()[,isolates, with = FALSE]))
        
    gg <- ggplot(otuFreqs, aes(x=otuOccurenceCounts))
    gg <- gg + geom_histogram(binwidth=1)
    gg <- gg + theme_bw()
    gg <- gg + scale_y_sqrt()
    gg <- gg + xlab("Number of isolates the DBLa type is observed in")
    gg <- gg + ylab("Count of DBLa types")
    gg
  })
  
  output$sideBarCumPlot <- renderPlot({
    event <- input$map_shape_click
    if (is.null(event) || is.null(clickData$clickedMarker))
      return()
    
    isolates <- isolateData()[isolateData()[['Location']] %in% event$id,][[1]]
    
    otus <- otuData()[,isolates, with = FALSE]
    cummulative <- data.frame(number=1:ncol(otus), cummulative=colSums(t(apply(otus, 1, cumsum))>0))
    
    gg <- ggplot(data=cummulative, aes(x=number, y=cummulative))
    gg <- gg + geom_line()
    gg <- gg + geom_point()
    gg <- gg + theme_bw() 
    gg <- gg + xlab("Number of isolates") + ylab("Cumulative total number of DBLalpha types")
    gg
  })
  
  # Show a popup at the given location
  showLocationPopup <- function(location, lat, lng) {
    print(location)
    
    selectedLocation <- locationData()[locationData()$Location == location,]
    content <- as.character(tagList(
      tags$h4("Location:"),
      tags$strong(HTML(sprintf("%s", selectedLocation$Location))), tags$br(),
      sprintf("Number of isolates: %s", as.integer(selectedLocation$TotalIsolates)), tags$br(),
      sprintf("Number of types: %s", as.integer(selectedLocation$TotalTypes)), tags$br(),
      sprintf("Experiment Dates: %s", selectedLocation$Dates)
    ))
    leafletProxy("map") %>% addPopups(lng, lat, content, layerId = location)
  }
  
  # When map is clicked, show a popup with city info
  observe({
    leafletProxy("map") %>% clearPopups()
    event <- input$map_shape_click
    if (is.null(event))
      return()
    
    isolate({
      print(event)
      print(event$id)
      print(event$lat)
      print(event$lng)
      showLocationPopup(event$id, event$lat, event$lng)
    })
  })

  clickData <- reactiveValues(clickedMarker=NULL)

  observeEvent(input$map_shape_click, {clickData$clickedMarker <- input$map_shape_click})

  observeEvent(input$map_click, {clickData$clickedMarker <- NULL})
  
}

shinyApp(ui = ui, server = server)

