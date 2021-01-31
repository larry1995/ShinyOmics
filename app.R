
library(shiny)
library(shinydashboard)
library(rgl)
library(grimon)
library(shinyWidgets)
library(shiny)
library(visNetwork)
library(shinyHeatmaply)
library(ggplot2)
library(visNetwork)
library(igraph)
library(RColorBrewer)
library(shiny)
library(heatmaply)
library(shinyHeatmaply)

source('./scripts/make_RNAseq_longtable.R')
source('./scripts/get_ci.R')
source('./scripts/make_node_table.R')
source('./scripts/files_validator.R')

exptsheet<-read.csv('./data/exptsheet.csv', header=T, stringsAsFactors = F)
ui<-fluidPage(
  
  # Application title
  titlePanel( h4(HTML("<b>KBCommons</b><br/> Omics Studio--3D Volcano"))),
  
  # Main panel has 4 tab panels
  mainPanel(
    
    tabPanel('Network', 
             fluidRow(),
             fluidRow(
               column(width=6,
                      uiOutput('network_selector'),
                      selectInput('network_experiment', 'Select Experiment', unique(exptsheet$Experiment)),
                      uiOutput('network_time_selector'),
                      h4("Significant Gene (fold change >2)",style="color:	#FF2D00"),
                      h4("Significant Gene (fold change <0.5)",style="color:	#0059FF"),
                      h4("Selected Gene", style = "color:magenta")
               ),# /network selectors
               column(width=6,
                      uiOutput('networkx_selector'),
                      uiOutput('networky_selector'),
                      uiOutput('networkcolor_selector'),
                      
                      checkboxInput('xaxis_log_net', 'log-scale x axis', FALSE),
                      checkboxInput('yaxis_log_net', 'log-scale y axis', FALSE),
                      checkboxInput("showviolins_net", "Show summary with Mean+-95% CI", FALSE)
                      
               ) # /scatterplot selectors
             ),
             fluidRow(
               column(width=6,
                      visNetworkOutput("networkplot")
               ),
               column(width=6,
                      plotOutput('networkstatsplot',brush = brushOpts(id="networkstats_brush")),
                      fluidRow(
                        downloadButton('downloadPlot_Panel4_png', 'Download plot (png)'),
                        downloadButton('downloadPlot_Panel4_svg', 'Download plot (svg)'),
                        downloadButton('downloadPlot_Panel4_pdf', 'Download plot (pdf)')
                      )
               )
             ),
             
             
             fluidRow(
               downloadButton('panel4download', 'Download Table')
             ), #/fluidRow for table DL 
             fluidRow(
               dataTableOutput("brushedTable_netstats")
               
             ) #/fluidRow

         
      ) # end tabPanel 4: Network
  
    ) # /tabsetPanel
  ) # /mainPanel




#######Server
server<-function(input, output) {
  
  validation <- files_validator()
  for (msg in validation){
    showNotification(msg, type="error", duration=0)
  }
  if (length(validation)>0){
    isvalidated <- FALSE
  }else{
    isvalidated <- TRUE
    showNotification("File Validation Step Passed!", type="message", duration=0)
  }
  
  values <- reactiveValues()
  
  
  netfiles <- Sys.glob('./data/network/*_Edges.csv')
  netnames <- as.character(sapply(netfiles, function(x) gsub(".*/network/(.+)_Edges.csv", "\\1",x)))
  net_table <- data.frame('Name'=netnames,'File'=netfiles)
  #selector for network
  output$network_selector <- renderUI({
    selectInput('network_selector', 'Select Network', netnames)
  })
  
  #selector for timepoint for network
  output$network_time_selector<-renderUI({
    thisexpt <- input$network_experiment
    t <- unique(exptsheet$Time[exptsheet$Experiment==thisexpt])
    selectInput('networkdatatime', 'Time(min)', t)
  })
  
  #load network data
  selectedNetwork <- reactive({
    thisnet <- input$network_selector
    edgetablefile <- as.character(net_table$File[net_table$Name==thisnet])
    print(c(thisnet, edgetablefile))
    edges <- read.csv(edgetablefile, header=T, stringsAsFactors = F)
    #print(head(edges))
    edges$source <- as.character(edges$source)
    edges$target <- as.character(edges$target)
    nodes <- make_node_table(edges)
    #add from and to node ID's to the edgetable
    edges$from <- nodes$id[match(edges$source, nodes$label)]
    edges$to <- nodes$id[match(edges$target, nodes$label)]
    
    return(list('nodes'=nodes,'edges'=edges))
  })
  
  #load RNAseq data
  panel4observer <- observe({
    thisexpt <- input$network_experiment
    thistime <- input$networkdatatime
    RNAdata <- make_RNAseq_longtable(exptsheet, thisexpt)
    RNAdata <- RNAdata[RNAdata$Time==thistime,]
    metafile <- exptsheet$MetadataFile[exptsheet$Experiment==thisexpt][1]
    metadata <- read.csv(metafile, header=T, stringsAsFactors = F)
    metacols <- names(metadata)[names(metadata)!="Gene"]
    values$metacols_network <- metacols
    # merge RNAseq data with metadata
    RNAdata <- merge(RNAdata, metadata, by="Gene", all.x=T, all.y=F, sort=F)
    # make names of resulting Data Frame unique
    names(RNAdata) <- make.names(names(RNAdata), unique=T)
    # if there are selected genes, filter the data frame
    if (!is.null(values$geneselection)){
      RNAdata <- RNAdata[RNAdata$Gene %in% values$geneselection,]
    }
    #print(names(RNAdata))
    values$RNAdata_network <- RNAdata
  }, suspended = T)
  
  # network plot
  
  #######Update RenderCachedPlot(cache expr =values$RNAdata_network )
  output$networkplot <- renderCachedPlot({
    validate(need(!is.null(values$RNAdata_network) & !is.null(input$network_selector), message="Waiting for datasets to be loaded..."))
    rnadata <- values$RNAdata_network
    #print(head(rnadata))
    networkdata <- selectedNetwork()
    print("network data loaded")
    values$networkdata <- networkdata
    nodes <- networkdata$nodes
    edges <- networkdata$edges
    
    df <- merge(nodes, rnadata, by.x="label", by.y="Gene", all.x=T, all.y=F)
    
    ## NETWORK LABELS & COLORS
    df$group[df$Sig & !is.na(df$Sig) & 
               df$Value>1 & !is.na(df$Value)] = "upSIG"
    df$group[df$Sig & !is.na(df$Sig) & 
               df$Value< -1 & !is.na(df$Value)] = "downSIG"
    
    df <- unique(df)
    values$networkdf <- df
    
    visNetwork(df, edges, width = "100%") %>%
      visPhysics(stabilization=F) %>%
      visEdges(smooth=F, color="grey", width=0.3)  %>%
      visNodes(color = list(background="white", highlight="magenta", border="black")) %>%
      visIgraphLayout() %>%
      visOptions(nodesIdSelection = list(enabled=T, useLabels=T,
                                         style = 'width: 200px; height: 26px;
                                         background: #f8f8f8;
                                         color: black;
                                         border:none;
                                         outline:none;'),
                 highlightNearest = list(enabled =TRUE, degree = 1, hover = T))%>%
      visGroups(groupname = "upSIG", color = "red") %>%
      visGroups(groupname = "downSIG", color = "blue") %>% 
      visExport(type = "png", name = "Network",
                float = "left", label = "Save network (png)", background = "white", style= "") 
    
    
  },
  ########Using the value of reactive expression in cacheKeyExpr https://shiny.rstudio.com/articles/plot-caching.html$#######
  ##Use values() instead values$xxx$###
  cacheKeyExpr = { list(values()) }
  )
  
  
  output$networkx_selector <- renderUI({
    selectInput('networkstats_x', 'X axis variable', 
                c(values$metacols_network))
  })
  output$networky_selector <- renderUI({
    selectInput('networkstats_y', 'Y axis variable', 
                c('Degree', 'Betweenness', 'Eigencentrality'))
  })
  output$networkcolor_selector <- renderUI({
    selectInput('networkstats_col', 'Color variable', 
                c(values$metacols_network, 'Degree', 'Betweenness', 'Eigencentrality'))
  })
  
  get_networkstats_axis_vars <- reactive({
    xaxisvar <- input$networkstats_x
    yaxisvar <- input$networkstats_y
    colorvar <- input$networkstats_col
    return(list('xaxisvar'=xaxisvar,
                'yaxisvar'=yaxisvar,
                'colorvar'=colorvar))
  })
  
  # network stats scatter plot
  output$networkstatsplot <- renderPlot({
    validate(need(!is.null(values$networkdf) , message="Waiting for datasets to be loaded..."))
    df <- values$networkdf
    axisvars <- get_networkstats_axis_vars()
    
    myxaxis <- axisvars$xaxisvar
    myyaxis <- axisvars$yaxisvar
    mycolor <- axisvars$colorvar
    validate(need(!all(is.na(df[mycolor])) & !all(is.na(df[myxaxis])) & !all(is.na(df[myyaxis])), 
                  message="Looks like the network and data don't come from the same organism. Please make sure all data is coming from the same species/strain!"))
    
    p=ggplot(df, aes_string(x=myxaxis, y=myyaxis, color=mycolor))+
      geom_point(alpha=0.5)+theme_bw()
    
    if(input$xaxis_log_net){p = p+scale_x_log10()}
    if(input$yaxis_log_net){p = p+scale_y_log10()}
    if(input$showviolins_net){p = p+#geom_violin(trim=FALSE)+ 
      stat_summary(fun.data=get_ci, conf.int=0.95, B=100, 
                   geom="pointrange", color="red")}
    
    values$Panel4_scatter <- p
    
    return(p)
  })
  

}



shinyAPP(ui,server)