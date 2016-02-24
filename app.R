library(shiny)
library(data.table)
library(VennDiagram)
library(ggplot2)
library(dplyr)

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
      
      tabPanel("Diagrams",
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
          absolutePanel(id = "controls", class = "panel panel-default", fixed = TRUE,
            draggable = TRUE, top = 60, left = "auto", right = 20, bottom = "auto",
            width = 330, height = "auto",
                                 
            h2("Summary plots"),
                                                   
#             plotOutput("histCentile", height = 200),
#             plotOutput("scatterCollegeIncome", height = 250)
          ),
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
        Dates = paste(unique(Date), collapse = '')
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
    
    #set circle radius to reflect the number of isolates
    radius <- locationData()$TotalIsolates / max(locationData()$TotalIsolates) * 30000
      
    #set circle colour to reflect the diversity at particular location
    colorData <- locationData()$TotalTypes/locationData()$TotalIsolates
    pal <- colorBin("Spectral", colorData, 7, pretty = FALSE)
    
    leafletProxy("map", data = locationData()) %>%
      clearShapes() %>%
      addCircles(~Longitude, ~Latitude, radius=radius, layerId=~location,
                 stroke=FALSE, fillOpacity=0.4, fillColor=pal(colorData)) %>%
      addLegend("bottomleft", pal=pal, values=colorData, title=colorBy,
                layerId="colorLegend")
  })
  
  # A reactive expression that returns the set of zips that are
  # in bounds right now
  locationsInBounds <- reactive({
    if (is.null(input$map_bounds))
      return(locationData()[FALSE,])
    bounds <- input$map_bounds
    latRng <- range(bounds$north, bounds$south)
    lngRng <- range(bounds$east, bounds$west)
    
    subset(locationData(),
           Latitude >= latRng[1] & Latitude <= latRng[2] &
             Longitude >= lngRng[1] & Longitude <= lngRng[2])
  })

#   output$histCentile <- renderPlot({
#     # If no zipcodes are in view, don't plot
#     if (nrow(zipsInBounds()) == 0)
#       return(NULL)
#     
#     hist(zipsInBounds()$centile,
#          breaks = centileBreaks,
#          main = "SuperZIP score (visible zips)",
#          xlab = "Percentile",
#          xlim = range(allzips$centile),
#          col = '#00DD00',
#          border = 'white')
#   })
  
#   output$scatterCollegeIncome <- renderPlot({
#     # If no zipcodes are in view, don't plot
#     if (nrow(zipsInBounds()) == 0)
#       return(NULL)
#     
#     print(xyplot(income ~ college, data = zipsInBounds(), xlim = range(allzips$college), ylim = range(allzips$income)))
#   })
  
  # Show a popup at the given location
  showLocationcodePopup <- function(location, lat, lng) {
    selectedLocation <- locationData()[locationData()$Location == location,]
    content <- as.character(tagList(
      tags$h4("Number of isolates:", as.integer(selectedLocation$TotalIsolates)),
      tags$strong(HTML(sprintf("%s",
                               selectedLocation$Location
      ))), tags$br(),
      sprintf("Number of isolates: %s", as.integer(selectedLocation$TotalIsolates)), tags$br(),
      sprintf("Number of types: %s%%", as.integer(selectedLocation$TotalTypes)), tags$br(),
      sprintf("Dates: %s", selectedLocation$Dates)
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
      showLocationPopup(event$id, event$lat, event$lng)
    })
  })
  
}

shinyApp(ui = ui, server = server)

