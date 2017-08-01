

# LIBRARY -----------------------------------------------------------------

library(shiny)
library(plotly)
library(shinythemes)
library(ggplot2)
library(DT)
library(pheatmap)
library(threejs)
library(sm)
library(RColorBrewer)
library(mclust)
library(reshape)

shinyServer(function(input, output) {
  set.seed(1)
  
  # FUNCTIONS ------------------------------------------------------------------
  
  options(shiny.maxRequestSize = 2000 * 1024 ^ 2)
  
  dataTables <- reactiveValues(
    log2cpm = NULL,
    tsne.data = NULL,
    featuredata = NULL,
    selectedDge = NULL
  )
  n_fun <- function(x) {
    return(data.frame(y = -.5, label = paste0(length(x), "\ncells")))
  }
  
  diffLRT = function(x, y, xmin = 1) {
    lrtX = bimodLikData(x)
    lrtY = bimodLikData(y)
    lrtZ = bimodLikData(c(x, y))
    lrt_diff = 2 * (lrtX + lrtY - lrtZ)
    return(pchisq(lrt_diff, 3, lower.tail = F))
  }
  
  bimodLikData = function(x, xmin = 0) {
    x1 = x[x <= xmin]
    x2 = x[x > xmin]
    xal = minmax(length(x2) / length(x),
                 min = 1e-5,
                 max = (1 - 1e-5))
    likA = length(x1) * log(1 - xal)
    mysd = sd(x2)
    if (length(x2) < 2) {
      mysd = 1
    }
    likB = length(x2) * log(xal) + sum(dnorm(x2, mean(x2), mysd, log = TRUE))
    return(likA + likB)
  }
  
  ainb = function(a, b) {
    a2 = a[a %in% b]
    return(a2)
  }
  
  minmax = function(data, min, max) {
    data2 = data
    data2[data2 > max] = max
    data2[data2 < min] = min
    return(data2)
  }
  set.ifnull = function(x, y) {
    if (is.null(x))
      return(y)
    return(x)
  }
  
  expMean = function(x) {
    return(log(mean(exp(x) - 1) + 1))
  }
  
  
  DiffExpTest = function(expression,
                         cells.1,
                         cells.2,
                         genes.use = NULL,
                         print.bar = TRUE) {
    genes.use = set.ifnull(genes.use, rownames(expression))
    p_val = unlist(lapply(genes.use, function(x)
      diffLRT(
        as.numeric(expression[x, cells.1]), as.numeric(expression[x, cells.2])
      )))
    to.return = data.frame(p_val, row.names = genes.use)
    return(to.return)
  }
  
  # INPUTS ------------------------------------------------------------------
  
  data <- reactive({
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    
    load(inFile$datapath)
    
    cat(stderr(), 'Loaded')
    
    dataTables$log2cpm <- log2cpm
    dataTables$tsne.data <- tsne.data
    featuredata <-
      featuredata[which(featuredata$Chromosome.Name %in% c(unlist(lapply(
        seq(1, 22, 1), toString
      )), c("X", "Y", "MT"))), ]
    featuredata$Associated.Gene.Name <-
      toupper(featuredata$Associated.Gene.Name)
    featuredata<-featuredata[rownames(log2cpm),]
    dataTables$featuredata <- featuredata
    dataTables$positiveCells <- NULL
    dataTables$positiveCellsAll <- NULL
    
    #cat(stderr(),grep('^T_',colnames(log2cpm)))
  })
  
  
  # RENDER UI  ------------------------------------------------------------------
  output$clusters <- renderUI({
    noOfClusters <- max(dataTables$tsne.data$dbCluster)
    selectInput(
      "cluster",
      label = "Cluster",
      choices = c(0:noOfClusters),
      selected = 1
    )
  })
  
  output$clusters1 <- renderUI({
    noOfClusters <- max(dataTables$tsne.data$dbCluster)
    selectInput(
      "clusters1",
      label = "Cluster",
      choices = c(0:noOfClusters),
      selected = 1
    )
  })
  
  output$clusters2 <- renderUI({
    noOfClusters <- max(dataTables$tsne.data$dbCluster)
    selectInput(
      "clusters2",
      label = "Cluster",
      choices = c(0:noOfClusters),
      selected = 1
    )
  })
  
  output$clusters3 <- renderUI({
    noOfClusters <- max(dataTables$tsne.data$dbCluster)
    selectInput(
      "clusters3",
      label = "Cluster",
      choices = c(0:noOfClusters),
      selected = 1
    )
  })
  
  output$clusters4 <- renderUI({
    noOfClusters <- max(dataTables$tsne.data$dbCluster)
    selectInput(
      "clusters4",
      label = "Cluster",
      choices = c(c('All'),c(0:noOfClusters)),
      selected = 1
    )
  })  
  
  
  # MAIN 3D PLOT ------------------------------------------------------------------
  output$tsne_main <- renderPlotly({
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    data()
    tsne.data <- as.data.frame(dataTables$tsne.data)
    #cat(stderr(),colnames(tsne.data)[1:5])
    tsne.data$dbCluster <- as.factor(tsne.data$dbCluster)
    
    p <-
      plot_ly(
        tsne.data,
        x = ~ V1,
        y = ~ V2,
        z = ~ V3,
        type = "scatter3d",
        color =  ~ dbCluster,
        hoverinfo = "text",
        text = paste('Cluster:', tsne.data$dbCluster),
        mode = 'markers',
        marker =
          list(
            line = list(width = 0),
            size = rep(10, nrow(tsne.data)),
            sizeref = 3
          )
      )
    layout(p)
    
  })
  
  # SUMMARY STATS ----------------------------------------------------------------
  
  output$summaryStats<-renderUI({
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    line1<-paste('No. of cells:', dim(dataTables$log2cpm)[2],sep='\t')
    line2<-paste('Median UMIs:', median(t(dataTables$log2cpm['ENSGUMI',])),sep='\t')
    line3<-paste('Median Genes:', median(t(dataTables$log2cpm['ENSGGENES',])),sep='\t')
    line4<-paste('No. of clusters:', max(dataTables$tsne.data$dbCluster),sep='\t')
    HTML(
      paste0("Summary statistics of this dataset:", '<br/>','<br/>',
        line1, '<br/>', 
             line2, '<br/>',
             line3, '<br/>',
             line4)
    )
  })
  
  # EXPLORE TAB 3D PLOT ------------------------------------------------------------------
  genes <- eventReactive(input$goButton, {
    genesin <- input$genesList
    genesin <- strsplit(genesin, ',')
    return(genesin)
  })
  
  v <- reactiveValues(doPlot = FALSE)
  observeEvent(input$goButton, {
    v$doPlot <- input$goButton
  })
  
  output$tsne_plt <- renderPlotly({
    if (v$doPlot == FALSE)
      return()
    
    isolate({
      geneid <- rownames(dataTables$featuredata[which(dataTables$featuredata$Associated.Gene.Name ==
                                                        toupper(input$gene_id)), ])[1]
      
      expression <- dataTables$log2cpm[geneid, ]
      cat(file = stderr(), rownames(expression))
      
      validate(need(
        is.na(sum(expression)) != TRUE,
        'Gene symbol incorrect or gene not expressed'
      ))
      
      tsne.data <- cbind(dataTables$tsne.data, t(expression))
      names(tsne.data)[names(tsne.data) == geneid] <- 'values'
      
      p <-
        plot_ly(
          tsne.data,
          x = ~ V1,
          y = ~ V2,
          z = ~ V3,
          type = "scatter3d",
          hoverinfo = "text",
          text = paste('Cluster:', tsne.data$dbCluster),
          mode = 'markers',
          marker = list(
            size = 2,
            line = list(width = 0),
            color =  ~ values,
            colors = 'Greens'
          )
        )
      layout(p, title = toupper(input$gene_id))
    })
  })
  # EXPLORE TAB CLUSTER PLOT ------------------------------------------------------------------
  output$clusterPlot <- renderPlot({
    if (v$doPlot == FALSE)
      return()
    
    isolate({
      geneid <- rownames(dataTables$featuredata[which(dataTables$featuredata$Associated.Gene.Name ==
                                                        toupper(input$gene_id)), ])[1]
      
      expression <- dataTables$log2cpm[geneid, ]
      
      validate(need(is.na(sum(expression)) != TRUE, ''))
      
      tsne.data <- cbind(dataTables$tsne.data, t(expression))
      names(tsne.data)[names(tsne.data) == geneid] <- 'values'
      
      subsetData <- subset(tsne.data, dbCluster == input$cluster)
      p1 <-
        ggplot(subsetData,
               aes_string(x = input$dimension_x, y = input$dimension_y)) +
        geom_point(aes_string(size = 2, color = 'values')) +
        geom_point(shape = 1,
                   size = 4,
                   colour = "black") +
        theme_bw() +
        theme(
          axis.text.x = element_text(
            angle = 90,
            size = 12,
            vjust = 0.5
          ),
          axis.text.y = element_text(size = 12),
          strip.text.x = element_text(size = 16),
          strip.text.y = element_text(size = 14),
          axis.title.x = element_text(face = "bold", size = 16),
          axis.title.y = element_text(face = "bold", size = 16),
          legend.position = "none"
        ) +
        ggtitle(paste(toupper(input$gene_id), input$cluster, sep = '-Cluster')) +
        scale_colour_gradient2(low = 'grey50', high = "red")
      p1
    })
  })
  # EXPLORE TAB VIOLIN PLOT ------------------------------------------------------------------
  output$gene_vio_plot <- renderPlot({
    if (v$doPlot == FALSE)
      return()
    
    isolate({
      geneid <- rownames(dataTables$featuredata[which(dataTables$featuredata$Associated.Gene.Name ==
                                                        toupper(input$gene_id)), ])[1]
      
      expression <- dataTables$log2cpm[geneid, ]
      
      validate(need(is.na(sum(expression)) != TRUE, ''))
      
      tsne.data <- cbind(dataTables$tsne.data, t(expression))
      names(tsne.data)[names(tsne.data) == geneid] <- 'values'
      #tsne.data<-subset(tsne.data,dbCluster!=0)
      
      p1 <-
        ggplot(tsne.data, aes(factor(dbCluster), values, fill = factor(dbCluster))) +
        geom_violin(scale = "width") +
        stat_summary(
          fun.y = median,
          geom = "point",
          size = 5,
          color = 'black'
        ) +
        stat_summary(fun.data = n_fun, geom = "text") +
        theme_bw() +
        theme(
          axis.text.x = element_text(
            angle = 90,
            size = 12,
            vjust = 0.5
          ),
          axis.text.y = element_text(size = 12),
          strip.text.x = element_text(size = 16),
          strip.text.y = element_text(size = 14),
          axis.title.x = element_text(face = "bold", size = 16),
          axis.title.y = element_text(face = "bold", size = 16),
          legend.position = "none"
        ) +
        xlab('Cluster') +
        ylab('Expression') +
        ggtitle(toupper(input$gene_id))
      p1
    })
  })
  # EXPLORE TABL DOWNLOAD SELECTED WITH BRUSH ------------------------------------------------------------------
  output$downloadExpression <- downloadHandler(
    filename = function() {
      paste(input$cluster, "Selected_Expression_table.csv", sep = '_')
    },
    content = function(file) {
      geneid <- rownames(dataTables$featuredata[which(dataTables$featuredata$Associated.Gene.Name ==
                                                        toupper(input$gene_id)), ])[1]
      
      expression <- dataTables$log2cpm[geneid, ]
      #cat(stderr(),colnames(expression)[1:5])
      tsne.data <- cbind(dataTables$tsne.data, t(expression))
      #cat(file=stderr(),grep('^T_',rownames(tsne.data)))
      
      names(tsne.data)[names(tsne.data) == geneid] <- 'values'
      
      #cat(file=stderr(),grep('^T_',rownames(tsne.data)))
      
      subsetData <- subset(tsne.data, dbCluster == input$cluster)
      #cat(file=stderr(),rownames(subsetData)[1:5])
      cells.names <- brushedPoints(subsetData, input$b1, allRows = T)
      #cat(file=stderr(),colnames(cells.names))
      cells <-
        rownames(subsetData[which(cells.names$selected_ == TRUE), ])
      #cat(file=stderr(),cells[1:5])
      
      if (length(cells) == 1) {
        subsetExpression <- dataTables$log2cpm[, cells]
        subsetExpression <-
          as.data.frame(subsetExpression, row.names = rownames(dataTables$log2cpm))
        colnames(subsetExpression) <- cells
        subsetExpression$Associated.Gene.Name <-
          dataTables$featuredata[rownames(subsetExpression), 'Associated.Gene.Name']
        write.csv(subsetExpression, file)
      }
      else{
        subsetExpression <- dataTables$log2cpm[, cells]
        #cat(stderr(),colnames(subsetExpression)[1:5])
        
        subsetExpression$Associated.Gene.Name <-
          dataTables$featuredata[rownames(subsetExpression), 'Associated.Gene.Name']
        #cat(stderr(),colnames(subsetExpression))
        write.csv(subsetExpression, file)
      }
    }
  )
  # EXPLORE TAB PANEL PLOT------------------------------------------------------------------
  vvvvvvvv <- reactiveValues(doPlot = FALSE)
  observeEvent(input$goButton8, {
    vvvvvvvv$doPlot <- input$goButton8
  })
  
  
  output$panelPlot <- renderPlot({
    if (vvvvvvvv$doPlot == FALSE)
      return()
    
    isolate({
      genesin <- input$panelplotids
      genesin <- toupper(genesin)
      genesin <- strsplit(genesin, ',')
      genesin<-genesin[[1]]
      
      cat(file=stderr(),length(genesin))
      par(mfrow=c(ceiling(length(genesin)/4),4), mai = c(0, 0., 0., 0.))
      rbPal <- colorRampPalette(c('#f0f0f0','red'))
      cat(file=stderr(),input$clusters4)
      
      if (input$clusters4 == 'All') 
        {
        for (i in 1:length(genesin)){
          Col <- rbPal(10)[
            as.numeric(
              cut(
                as.numeric(
                  dataTables$log2cpm[
                    rownames(dataTables$featuredata[which(dataTables$featuredata$Associated.Gene.Name==genesin[i]),])
                    ,]
                ),breaks = 10))]
          plot(dataTables$tsne.data[,input$dimension_x4],dataTables$tsne.data[,input$dimension_y4],col=Col,pch=16,axes = FALSE,frame.plot = TRUE, ann=FALSE)
          title(genesin[i],line=-1.2,adj = 0.05,cex.main=2)
          cat(file=stderr(),genesin[i])
        }
      }
        else{
          for (i in 1:length(genesin)){
            
            subsetTSNE <- subset(dataTables$tsne.data, dbCluster == input$clusters4)

            Col <- rbPal(10)[
              as.numeric(
                cut(
                  as.numeric(
                    dataTables$log2cpm[
                      rownames(dataTables$featuredata[which(dataTables$featuredata$Associated.Gene.Name==genesin[i]),])
                      ,]
                  ),breaks = 10))]
            
            names(Col)<-rownames(dataTables$tsne.data)
            plotCol<-Col[rownames(subsetTSNE)]
            plot(subsetTSNE[,input$dimension_x4],subsetTSNE[,input$dimension_y4],col=plotCol,pch=16,axes = FALSE,frame.plot = TRUE, ann=FALSE)
            title(genesin[i],line=-1.2,adj = 0.05,cex.main=2)
            cat(file=stderr(),input$clusters4)
        }
      }
    })
  })
  # CO-EXPRESSION HEATMAP ALL CLUSTERS ------------------------------------------------------------------
  vv <- reactiveValues(doPlot = FALSE)
  observeEvent(input$goButton1, {
    vv$doPlot <- input$goButton1
  })
  
  # output$heatmap<-renderPlotly({
  #
  #   if (vv$doPlot == FALSE) return()
  #
  #   isolate({
  #
  #     genesin<-input$heatmap_geneids
  #     genesin<-toupper(genesin)
  #     genesin<-strsplit(genesin,',')
  #
  #     map<-rownames(dataTables$featuredata[
  #       which(dataTables$featuredata$Associated.Gene.Name %in% genesin[[1]]),])
  #     cat(file=stderr(),map[1])
  #
  #     #geneid<-rownames(
  #      # dataTables$featuredata[which(dataTables$featuredata$Associated.Gene.Name==toupper(input$gene_id)),])[1]
  #
  #     expression<-dataTables$log2cpm[map,]
  #
  #     validate(
  #       need(is.na(sum(expression))!=TRUE,'Gene symbol incorrect or genes not expressed')
  #     )
  #
  #     tsne.data<-dataTables$tsne.data
  #     #names(tsne.data)[names(tsne.data) == geneid] <- 'values'
  #
  #     #p <- plot_ly(tsne.data, x = ~V1, y = ~V2, z = ~V3, type = "scatter3d",
  #     #             hoverinfo = "text",
  #     #             text=paste('Cluster:',tsne.data$dbCluster),
  #     #             mode = 'markers',marker=list(size=2,line=list(width=0),
  #     #                                          color=~values,colors='Greens')
  #
  #     #expression<-rbind(expression,t(tsne.data$dbCluster))
  #
  #     p<-plot_ly()%>%
  #       add_trace(z=as.matrix(expression),
  #                y=dataTables$featuredata[rownames(expression),'Associated.Gene.Name'],
  #                hoverinfo='text',
  #                text=paste(rownames(expression),colnames(expression)),
  #                type='heatmap',
  #                zmax=5,
  #                zmin=-5,
  #                colorscale='RdBu',
  #                xtype='scaled'
  #                )%>%
  #     layout(p)
  #   })
  # })
  
  output$heatmap <- renderPlot({
    if (vv$doPlot == FALSE)
      return()
    
    isolate({
      genesin <- input$heatmap_geneids
      genesin <- toupper(genesin)
      genesin <- strsplit(genesin, ',')
      
      map <- rownames(dataTables$featuredata[which(dataTables$featuredata$Associated.Gene.Name %in% genesin[[1]]), ])
      cat(file = stderr(), length(map))
      
      expression <- dataTables$log2cpm[map, ]
      
      validate(need(
        is.na(sum(expression)) != TRUE,
        'Gene symbol incorrect or genes not expressed'
      ))
      
      tsne.data <- dataTables$tsne.data
      tsne.data <- tsne.data[order(tsne.data$dbCluster), ]
      
      expression <- expression[, rownames(tsne.data)]
      expression <- expression[complete.cases(expression), ]
      
      annotation <- data.frame(factor(tsne.data$dbCluster))
      rownames(annotation) <- colnames(expression)
      colnames(annotation) <- c('Cluster')
      
      h <-
        pheatmap(
          as.matrix(expression),
          cluster_rows = TRUE,
          cluster_cols = FALSE,
          scale = 'row',
          fontsize_row = 10,
          labels_col = colnames(expression),
          labels_row = dataTables$featuredata[rownames(expression), 'Associated.Gene.Name'],
          show_rownames = TRUE,
          annotation_col = annotation,
          show_colnames = FALSE,
          annotation_legend = TRUE,
          breaks = seq(-6, 6, by = .12),
          colorRampPalette(rev(brewer.pal(
            n = 6, name =
              "RdBu"
          )))(100)
          
        )
      h
      
      # h3<-heatmap.3(as.matrix(expression),
      #               Colv = F,
      #               Rowv=F,
      #               na.color = "gray95",
      #               #dendrogram = "row",
      #               col = colorRampPalette(rev(brewer.pal(n = 6, name =
      #                                                       "RdBu")))(100),
      #               #trace = "row",
      #               tracecol = NULL,
      #               linecol = "gray80",
      #               KeyValueName = "Z-score",
      #               ColSideColors = as.matrix(annotation$Cluster),
      #               lhei = c(2,5),
      #               lwid = c(2,5),
      #               side.height.fraction = 0.5
      #               #labRow = dataTables$featuredata[rownames(expression),'Associated.Gene.Name']
      #
      #               )
      # h3
    })
  })
  
  # CO EXPRESSION TAB CLUSTER PLOT ------------------------------------------------------------------
  vvvvv <- reactiveValues(doPlot = FALSE)
  observeEvent(input$goButton5, {
    vvvvv$doPlot <- input$goButton5
  })
  
  output$clusterPlot2 <- renderPlot({
    if (vvvvv$doPlot == FALSE)
      return()
    
    isolate({
      geneid <- rownames(dataTables$featuredata[which(dataTables$featuredata$Associated.Gene.Name ==
                                                        toupper(input$gene_id_sch)), ])[1]
      
      expression <- dataTables$log2cpm[geneid, ]
      
      validate(need(is.na(sum(expression)) != TRUE, ''))
      
      tsne.data <- cbind(dataTables$tsne.data, t(expression))
      names(tsne.data)[names(tsne.data) == geneid] <- 'values'
      
      subsetData <- subset(tsne.data, dbCluster == input$clusters2)
      p1 <-
        ggplot(subsetData,
               aes_string(x = input$dimension_x2, y = input$dimension_y2)) +
        geom_point(aes_string(size = 2, color = 'values')) +
        geom_point(shape = 1,
                   size = 4,
                   colour = "black") +
        theme_bw() +
        theme(
          axis.text.x = element_text(
            angle = 90,
            size = 12,
            vjust = 0.5
          ),
          axis.text.y = element_text(size = 10),
          strip.text.x = element_text(size = 16),
          strip.text.y = element_text(size = 14),
          axis.title.x = element_text(face = "bold", size = 16),
          axis.title.y = element_text(face = "bold", size = 16),
          legend.position = "none"
        ) +
        ggtitle(paste(toupper(input$gene_id_sch), input$clusters2, sep =
                        '-Cluster')) +
        scale_colour_gradient2(low = 'grey50', high = "red")
      p1
    })
  })
  # CO EXPRESSION TAB SELECTED HEATMAP ------------------------------------------------------------------
  vvvvvv <- reactiveValues(doPlot = FALSE)
  observeEvent(input$goButton6, {
    vvvvvv$doPlot <- input$goButton6
  })
  
  output$selectedHeatmap <- renderPlot({
    if (vvvvvv$doPlot == FALSE)
      return()
    
    isolate({
      genesin <- input$heatmap_geneids2
      genesin <- toupper(genesin)
      genesin <- strsplit(genesin, ',')
      
      subsetData <-
        subset(dataTables$tsne.data, dbCluster == input$clusters2)
      cells.1 <- rownames(brushedPoints(subsetData, input$scb1))
      
      
      map <- rownames(dataTables$featuredata[which(dataTables$featuredata$Associated.Gene.Name %in% genesin[[1]]), ])
      #cat(file=stderr(),map[1])
      
      expression <- dataTables$log2cpm[map, cells.1]
      cat(file = stderr(), rownames(expression))
      
      expression <- expression[complete.cases(expression), ]
      cat(file = stderr(), rownames(expression))
      mColor <- max(expression)
      
      validate(need(
        is.na(sum(expression)) != TRUE,
        'Gene symbol incorrect or genes not expressed'
      ))
      
      #annotation<-data.frame(factor(tsne.data$dbCluster))
      #rownames(annotation)<-colnames(expression)
      #colnames(annotation)<-c('Cluster')
      
      h <-
        pheatmap(
          as.matrix(expression),
          cluster_rows = TRUE,
          cluster_cols = TRUE,
          scale = 'row',
          fontsize_row = 10,
          labels_col = colnames(expression),
          labels_row = dataTables$featuredata[rownames(expression), 'Associated.Gene.Name'],
          show_rownames = TRUE,
          show_colnames = FALSE,
          breaks = seq(-6, 6, by = .12),
          colorRampPalette(rev(brewer.pal(
            n = 6, name =
              "RdBu"
          )))(100)
          
        )
      h
    })
  })
  
  # CO EXPRESSION TAB ON/OFF PLOT ------------------------------------------------------------------
  vvvvvvv <- reactiveValues(doPlot = FALSE)
  observeEvent(input$goButton7, {
    vvvvvvv$doPlot <- input$goButton7
  })
  
  output$plotCoExpression <- renderPlot({
    if (vvvvvvv$doPlot == FALSE)
      return()
    
    isolate({
      genesin <- input$mclustids
      genesin <- toupper(genesin)
      genesin <- strsplit(genesin, ',')
      
      subsetData <-
        subset(dataTables$tsne.data, dbCluster == input$clusters3)
      cells.1 <- rownames(subsetData)
      
      
      map <- rownames(dataTables$featuredata[which(dataTables$featuredata$Associated.Gene.Name %in% genesin[[1]]), ])
      #cat(file=stderr(),map[1])
      
      expression <- dataTables$log2cpm[map, ]
      #cat(file=stderr(),rownames(expression))
      
      #expression<-expression[complete.cases(expression),]
      #cat(file=stderr(),rownames(expression))
      
      validate(need(
        is.na(sum(expression)) != TRUE,
        'Gene symbol incorrect or genes not expressed'
      ))
      
      bin <- expression
      bin[] <- 0
      
      for (i in 1:nrow(expression))
      {
        x <- Mclust(expression[i, ], G = 2)
        bin[i, ] <- x$classification
      }
      bin <- bin - 1
      allexprs <- apply(bin, 2, sum)
      plotexprs <- allexprs
      plotexprs[] <- 0
      plotexprs[allexprs >= length(rownames(bin))] <- 1
      dataTables$positiveCells <- allexprs >= length(rownames(bin))
      dataTables$positiveCellsAll <- plotexprs
      #save(subsetData,bin,allexprs,file='~/Desktop/test.Rds')
      #cat(file=stderr(),names(allexprs))
      
      mergeExprs <- plotexprs[rownames(subsetData)]
      #cat(file=stderr(),length(mergeExprs))
      
      subsetData$CoExpression <- mergeExprs
      #cat(file=stderr(),colnames(subsetData))
      
      p1 <-
        ggplot(subsetData,
               aes_string(x = input$dimension_x3, y = input$dimension_y3)) +
        geom_point(aes_string(size = 2, color = 'CoExpression')) +
        geom_point(shape = 1,
                   size = 4,
                   colour = "black") +
        theme_bw() +
        theme(
          axis.text.x = element_text(
            angle = 90,
            size = 12,
            vjust = 0.5
          ),
          axis.text.y = element_text(size = 12),
          strip.text.x = element_text(size = 16),
          strip.text.y = element_text(size = 14),
          axis.title.x = element_text(face = "bold", size = 16),
          axis.title.y = element_text(face = "bold", size = 16),
          legend.position = "none"
        ) +
        #ggtitle(paste(toupper(input$gene_id),input$cluster,sep='-Cluster'))+
        scale_colour_gradient2(low = 'grey50', high = "red")
      p1
    })
  })
  
  
  # ONOFF TAB DOWNLOAD POSITIVECELLS ------------------------------------------------------------------
  output$downloadExpressionOnOff <- downloadHandler(
    filename = function() {
      paste(input$clusters3, "PositiveCells.csv", sep = '_')
    },
    content = function(file) {
      cells <- dataTables$positiveCells
      #cat(file=stderr(),cells[1:5])
      
      if (length(cells) == 1) {
        subsetExpression <- dataTables$log2cpm[, cells]
        subsetExpression <-
          as.data.frame(subsetExpression, row.names = rownames(dataTables$log2cpm))
        colnames(subsetExpression) <- cells
        subsetExpression$Associated.Gene.Name <-
          dataTables$featuredata[rownames(subsetExpression), 'Associated.Gene.Name']
        write.csv(subsetExpression, file)
      }
      else{
        subsetExpression <- dataTables$log2cpm[, cells]
        #cat(stderr(),colnames(subsetExpression)[1:5])
        subsetExpression$Associated.Gene.Name <-
          dataTables$featuredata[rownames(subsetExpression), 'Associated.Gene.Name']
        #cat(stderr(),colnames(subsetExpression))
        write.csv(subsetExpression, file)
      }
    }
  )
  
  # ONOFF TAB RENDER TABLE ALL CELLS ------------------------------------------------------------------
  output$onOffTable <- DT::renderDataTable({
    if (vvvvvvv$doPlot == FALSE)
      return()
    
    isolate({
      merge <- dataTables$tsne.data
      merge$CoExpression <- dataTables$positiveCellsAll
      df <-
        as.data.frame(table(merge[, c('dbCluster', 'CoExpression')]))
      dfOut <- cast(df, dbCluster ~ CoExpression)
      colnames(dfOut) <- c("Cluster", 'OFF', 'ON')
      #save(merge,file='~/Desktop/test.Rds')
      rownames(dfOut) <- dfOut$Cluster
      dfOut['Sum', ] <- c('', sum(dfOut$OFF), sum(dfOut$ON))
      DT::datatable(dfOut)
      
    })
  })
  
  
  # SUBCLUSTER DGE PLOT1 ------------------------------------------------------------------
  vvv <- reactiveValues(doPlot = FALSE)
  observeEvent(input$goButton2, {
    vvv$doPlot <- input$goButton2
  })
  
  output$dge_plot1 <- renderPlot({
    if (vvv$doPlot == FALSE)
      return()
    
    isolate({
      tsne.data <- dataTables$tsne.data
      
      subsetData <- subset(tsne.data, dbCluster == input$clusters1)
      p1 <-
        ggplot(subsetData,
               aes_string(x = input$dimension_x1, y = input$dimension_y1)) +
        geom_point() +
        geom_point(shape = 1,
                   size = 4,
                   color = "black") +
        theme_bw() +
        theme(
          axis.text.x = element_text(
            angle = 90,
            size = 12,
            vjust = 0.5
          ),
          axis.text.y = element_text(size = 12),
          strip.text.x = element_text(size = 16),
          strip.text.y = element_text(size = 14),
          axis.title.x = element_text(face = "bold", size = 16),
          axis.title.y = element_text(face = "bold", size = 16),
          legend.position = "none"
        ) +
        ggtitle(input$clusters1)
      p1
    })
  })
  # SUBCLUSTER DGE PLOT2 ------------------------------------------------------------------
  
  output$dge_plot2 <- renderPlot({
    if (vvv$doPlot == FALSE)
      return()
    
    isolate({
      tsne.data <- dataTables$tsne.data
      
      subsetData <- subset(tsne.data, dbCluster == input$clusters1)
      p1 <-
        ggplot(subsetData,
               aes_string(x = input$dimension_x1, y = input$dimension_y1)) +
        geom_point() +
        geom_point(shape = 1,
                   size = 4,
                   color = "black") +
        theme_bw() +
        theme(
          axis.text.x = element_text(
            angle = 90,
            size = 12,
            vjust = 0.5
          ),
          axis.text.y = element_text(size = 12),
          strip.text.x = element_text(size = 16),
          strip.text.y = element_text(size = 14),
          axis.title.x = element_text(face = "bold", size = 16),
          axis.title.y = element_text(face = "bold", size = 16),
          legend.position = "none"
        ) +
        ggtitle(input$clusters1)
      p1
    })
  })
  # SUBCLUSTER DGE ANALYSIS ------------------------------------------------------------------
  
  vvvv <- reactiveValues(doPlot = FALSE)
  observeEvent(input$goButton3, {
    vvvv$doPlot <- input$goButton3
  })
  
  
  dge <- reactive({
    if (vvvv$doPlot == FALSE)
      return()
    
    isolate({
      subsetData <- subset(dataTables$tsne.data, dbCluster == input$clusters1)
      cells.1 <- rownames(brushedPoints(subsetData, input$db1))
      
      cat(file = stderr(), cells.1[1:5])
      
      cells.2 <- rownames(brushedPoints(subsetData, input$db2))
      cat(file = stderr(), cells.2[1:5])
      
      subsetExpression <- dataTables$log2cpm[, union(cells.1, cells.2)]
      
      genes.use <- rownames(subsetExpression)
      data.1 = apply(subsetExpression[genes.use, cells.1], 1, expMean)
      data.2 = apply(subsetExpression[genes.use, cells.2], 1, expMean)
      total.diff = (data.1 - data.2)
      
      genes.diff = names(which(abs(total.diff) > .2))
      genes.use = ainb(genes.use, genes.diff)
      
      toReturn <-
        DiffExpTest(subsetExpression, cells.1, cells.2, genes.use = genes.use)
      toReturn[, "avg_diff"] = total.diff[rownames(toReturn)]
      toReturn$Associated.Gene.Name <-
        dataTables$featuredata[rownames(toReturn), 'Associated.Gene.Name']
      dataTables$selectedDge <- toReturn
      return(toReturn)
      cat(stderr(), rownames(toReturn)[1:5])
      
    })
  })
  # SUBCLUSTER DGE OUTPUT TABLE ------------------------------------------------------------------
  
  output$dge <- DT::renderDataTable({
    if (vvvv$doPlot == FALSE)
      return()
    
    isolate({
      top.genes <- dge()
      top.genes$Associated.Gene.Name <-
        dataTables$featuredata[rownames(top.genes), 'Associated.Gene.Name']
      if (dim(top.genes)[1] > 1) {
        DT::datatable(top.genes,
                      options = list(
                        orderClasses = TRUE,
                        lengthMenu = c(10, 30, 50),
                        pageLength = 10
                      ))
      }
    })
  })
  
  # SUBCLUSTER DGE DOWNLOADS ------------------------------------------------------------------
  output$download_dge_table <- downloadHandler(
    filename = function() {
      paste("SubCluster", "DGE_table.csv", sep = '_')
    },
    content = function(file) {
      write.csv(dataTables$selectedDge, file)
    }
  )
})
# VIOPLOT_SOURCE ----------------------------------------------------------

my.vioplot <-
  function (x,
            ...,
            range = 1.5,
            h = NULL,
            ylim = NULL,
            names = NULL,
            horizontal = FALSE,
            col = "magenta",
            border = "black",
            lty = 1,
            lwd = 1,
            rectCol = "black",
            colMed = "white",
            pchMed = 19,
            at,
            add = FALSE,
            wex = 1,
            drawRect = TRUE)
  {
    datas <- list(x, ...)
    n <- length(datas)
    if (missing(at))
      at <- 1:n
    upper <- vector(mode = "numeric", length = n)
    lower <- vector(mode = "numeric", length = n)
    q1 <- vector(mode = "numeric", length = n)
    q3 <- vector(mode = "numeric", length = n)
    med <- vector(mode = "numeric", length = n)
    base <- vector(mode = "list", length = n)
    height <- vector(mode = "list", length = n)
    baserange <- c(Inf,-Inf)
    args <- list(display = "none")
    if (!(is.null(h)))
      args <- c(args, h = h)
    for (i in 1:n) {
      data <- datas[[i]]
      data.min <- min(data)
      data.max <- max(data)
      q1[i] <- quantile(data, 0.25)
      q3[i] <- quantile(data, 0.75)
      med[i] <- median(data)
      iqd <- q3[i] - q1[i]
      upper[i] <- min(q3[i] + range * iqd, data.max)
      lower[i] <- max(q1[i] - range * iqd, data.min)
      est.xlim <- c(min(lower[i], data.min), max(upper[i],
                                                 data.max))
      smout <- do.call("sm.density", c(list(data, xlim = est.xlim),
                                       args))
      hscale <- 0.4 / max(smout$estimate) * wex
      base[[i]] <- smout$eval.points
      height[[i]] <- smout$estimate * hscale
      t <- range(base[[i]])
      baserange[1] <- min(baserange[1], t[1])
      baserange[2] <- max(baserange[2], t[2])
    }
    if (!add) {
      xlim <- if (n == 1)
        at + c(-0.5, 0.5)
      else
        range(at) + min(diff(at)) / 2 * c(-1, 1)
      if (is.null(ylim)) {
        ylim <- baserange
      }
    }
    if (is.null(names)) {
      label <- 1:n
    }
    else {
      label <- names
    }
    boxwidth <- 0.05 * wex
    if (!add)
      plot.new()
    if (!horizontal) {
      if (!add) {
        plot.window(xlim = xlim, ylim = ylim)
        axis(2)
        axis(1, at = at, label = label)
      }
      box()
      for (i in 1:n) {
        polygon(
          c(at[i] - height[[i]], rev(at[i] + height[[i]])),
          c(base[[i]], rev(base[[i]])),
          col = col[i],
          border = border,
          lty = lty,
          lwd = lwd
        )
        if (drawRect) {
          lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd,
                lty = lty)
          rect(at[i] - boxwidth / 2, q1[i], at[i] + boxwidth / 2,
               q3[i], col = rectCol)
          points(at[i], med[i], pch = pchMed, col = colMed)
        }
      }
    }
    else {
      if (!add) {
        plot.window(xlim = ylim, ylim = xlim)
        axis(1)
        axis(2, at = at, label = label)
      }
      box()
      for (i in 1:n) {
        polygon(
          c(base[[i]], rev(base[[i]])),
          c(at[i] - height[[i]],
            rev(at[i] + height[[i]])),
          col = col,
          border = border,
          lty = lty,
          lwd = lwd
        )
        if (drawRect) {
          lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd,
                lty = lty)
          rect(q1[i], at[i] - boxwidth / 2, q3[i], at[i] +
                 boxwidth / 2, col = rectCol)
          points(med[i], at[i], pch = pchMed, col = colMed)
        }
      }
    }
    invisible(list(
      upper = upper,
      lower = lower,
      median = med,
      q1 = q1,
      q3 = q3
    ))
  }

# HEATMAP.3_SOURCE ------------------------------------------------------------------
heatmap.3 <- function(x,
                      Rowv = TRUE,
                      Colv = if (symm)
                        "Rowv"
                      else
                        TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both", "row", "column", "none"),
                      symm = FALSE,
                      scale = c("none", "row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv, "Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) ||
                        scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column", "row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5, 5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction = 0.3,
                      cexRow = 0.2 + 1 / log10(nr),
                      cexCol = 0.2 + 1 / log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) ||
                        symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      ColSideColorsSize = 1,
                      RowSideColorsSize = 1,
                      KeyValueName = "Value",
                      ...) {
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else
      return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x,
                      low = min(x),
                      high = max(x)) {
    x <- (x - low) / (high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else
    match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning(
      "Using scale=\"row\" or scale=\"column\" when breaks are",
      "specified can produce unpredictable results.",
      "Please consider using only one or the other."
    )
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else
        dedrogram <- "none"
      warning(
        "Discrepancy: Rowv is FALSE, while dendrogram is `",
        dendrogram,
        "'. Omitting row dendogram."
      )
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else
        dendrogram <- "none"
      warning(
        "Discrepancy: Colv is FALSE, while dendrogram is `",
        dendrogram,
        "'. Omitting column dendogram."
      )
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else
      colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else
        t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else
        t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else
    rownames(x)
  else
    labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else
    colnames(x)
  else
    labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else
      breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1,] + 1, c(NA, 1), lmat[2,] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei = c(lhei[1], side.height.fraction * ColSideColorsSize / 2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <-
        cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <-
        c(lwid[1], side.height.fraction * RowSideColorsSize / 2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat,
         widths = lwid,
         heights = lhei,
         respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)) {
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[, rowInd, drop = F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(
          1,
          0:(dim(rsc)[2] - 1) / max(1, (dim(rsc)[2] - 1)),
          rownames(RowSideColors),
          las = 2,
          tick = FALSE
        )
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    if (!is.matrix(ColSideColors)) {
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop = F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(
          2,
          0:(dim(csc)[2] - 1) / max(1, (dim(csc)[2] - 1)),
          colnames(ColSideColors),
          las = 2,
          tick = FALSE
        )
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else
    iy <- 1:nr
  image(
    1:nc,
    1:nr,
    x,
    xlim = 0.5 + c(0, nc),
    ylim = 0.5 + c(0, nr),
    axes = FALSE,
    xlab = "",
    ylab = "",
    col = col,
    breaks = breaks,
    ...
  )
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) {
    # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(
      1:nc,
      1:nr,
      mmat,
      axes = FALSE,
      xlab = "",
      ylab = "",
      col = na.color,
      add = TRUE
    )
  }
  axis(
    1,
    1:nc,
    labels = labCol,
    las = 2,
    line = -0.5,
    tick = 0,
    cex.axis = cexCol
  )
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(
    4,
    iy,
    labels = labRow,
    las = 2,
    line = -0.5,
    tick = 0,
    cex.axis = cexRow
  )
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep)
      rect(
        xleft = csep + 0.5,
        ybottom = rep(0, length(csep)),
        xright = csep + 0.5 + sepwidth[1],
        ytop = rep(ncol(x) + 1, csep),
        lty = 1,
        lwd = 1,
        col = sepcolor,
        border = sepcolor
      )
  if (!missing(rowsep))
    for (rsep in rowsep)
      rect(
        xleft = 0,
        ybottom = (ncol(x) + 1 - rsep) - 0.5,
        xright = nrow(x) + 1,
        ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2],
        lty = 1,
        lwd = 1,
        col = sepcolor,
        border = sepcolor
      )
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals,
               col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(
        x = xv,
        y = yv,
        lwd = 1,
        col = tracecol,
        type = "s"
      )
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline,
               col = linecol,
               lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i,] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(
        x = xv,
        y = yv,
        lwd = 1,
        col = tracecol,
        type = "s"
      )
    }
  }
  if (!missing(cellnote))
    text(
      x = c(row(cellnote)),
      y = c(col(cellnote)),
      labels = c(cellnote),
      col = notecol,
      cex = notecex
    )
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(
      ddr,
      horiz = TRUE,
      axes = FALSE,
      yaxs = "i",
      leaflab = "none"
    )
  }
  else
    plot.new()
  par(mar = c(0, 0, if (!is.null(main))
    5
    else
      0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc,
         axes = FALSE,
         xaxs = "i",
         leaflab = "none")
  }
  else
    plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(
      z = matrix(z, ncol = 1),
      col = col,
      breaks = tmpbreaks,
      xaxt = "n",
      yaxt = "n"
    )
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else
      mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x,
            dens$y / max(dens$y) * 0.95,
            col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y) / max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(
        hx,
        hy / max(hy) * 0.95,
        lwd = 1,
        type = "s",
        col = denscol
      )
      axis(2, at = pretty(hy) / max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else
      title("Color Key")
  }
  else
    plot.new()
  retval$colorTable <-
    data.frame(
      low = retval$breaks[-length(retval$breaks)],
      high = retval$breaks[-1],
      color = retval$col
    )
  invisible(retval)
}
