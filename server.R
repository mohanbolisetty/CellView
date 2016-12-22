
# LIBRARY -----------------------------------------------------------------

library(shiny)
library(plotly)
library(shinythemes)


shinyServer(function(input, output) {
  
  set.seed(1)
  
  # INPUTS ------------------------------------------------------------------
  
  #log2cpm<-read.csv('Data/Expression.csv',row.names=1,stringsAsFactors = F, as.is=T, check.names=F)
  
  #featuredata<-read.csv('Databases/HG19_v74_FeatureData.csv',row.names=1,stringsAsFactors = F, as.is=T)
  #tsne.data<-read.csv('Data/Filtered_log2cpm_TNSE_dbscan.csv',row.names=1,stringsAsFactors = F,as.is=T)
  
 #load('data.Rds')
  
  options(shiny.maxRequestSize = 100*1024^2)
  
  
  saved_plots_and_tables <- reactiveValues(log2cpm = NULL,
                                           tsne.data=NULL,
                                           featuredata=NULL)
  n_fun <- function(x){
    return(data.frame(y = -.5, label = paste0(length(x),"\ncells")))
  }
  
  data<-reactive({

      inFile <- input$file1
      
      if (is.null(inFile))
        return(NULL)
      
      
      load(inFile$datapath)
      
      #cat(stderr(),ls())
      
      saved_plots_and_tables$log2cpm<-log2cpm
      saved_plots_and_tables$tsne.data<-tsne.data
      featuredata$Associated.Gene.Name<-toupper(featuredata$Associated.Gene.Name)
      saved_plots_and_tables$featuredata<-featuredata
      
      #cat(stderr(),grep('^T_',colnames(log2cpm)))
  })
  
  
  
  output$clusters <- renderUI({
    noOfClusters <- max(saved_plots_and_tables$tsne.data$dbCluster)
    selectInput("cluster",label="Cluster", choices=c(0:noOfClusters),selected=1)
  })
  
  

  output$tsne_main<-renderPlotly({
      
      inFile <- input$file1
      
      if (is.null(inFile))
        return(NULL)
      
      data()
      tsne.data<-as.data.frame(saved_plots_and_tables$tsne.data)
      #cat(stderr(),colnames(tsne.data)[1:5])
      tsne.data$dbCluster<-as.factor(tsne.data$dbCluster)

      p <- plot_ly(tsne.data, x = ~V1, y = ~V2, z = ~V3, type = "scatter3d",color=~dbCluster,
                   hoverinfo = "text",
                   text=paste('Cluster:',tsne.data$dbCluster),
                   mode = 'markers',marker=list(size=2,line=list(width=0)
                                                )
      )      
      layout(p)
    
  })
  genes <- eventReactive(input$goButton, {
    genesin<-input$genesList
    genesin<-strsplit(genesin,',')
    return(genesin)
  })
  
  v <- reactiveValues(doPlot = FALSE)
  observeEvent(input$goButton, {
    v$doPlot <- input$goButton
  })
  
  output$tsne_plt<-renderPlotly({
    
    if (v$doPlot == FALSE) return()
    
    isolate({
      
      geneid<-rownames(
        saved_plots_and_tables$featuredata[which(saved_plots_and_tables$featuredata$Associated.Gene.Name==toupper(input$gene_id)),])[1]
      
      expression<-saved_plots_and_tables$log2cpm[geneid,]
      cat(file=stderr(),rownames(expression))
      
      validate(
        need(is.na(sum(expression))!=TRUE,'Gene symbol incorrect or gene not expressed')
        )
      
      tsne.data<-cbind(saved_plots_and_tables$tsne.data,t(expression))
      names(tsne.data)[names(tsne.data) == geneid] <- 'values'
      
      p <- plot_ly(tsne.data, x = ~V1, y = ~V2, z = ~V3, type = "scatter3d",
                   hoverinfo = "text",
                   text=paste('Cluster:',tsne.data$dbCluster),
                   mode = 'markers',marker=list(size=2,line=list(width=0),
                                                color=~values,colors='Greens')
      )      
      layout(p)
    })
  })
  
  output$clusterPlot<-renderPlot({
    
    if (v$doPlot == FALSE) return()
    
    isolate({
      
      geneid<-rownames(
        saved_plots_and_tables$featuredata[which(saved_plots_and_tables$featuredata$Associated.Gene.Name==toupper(input$gene_id)),])[1]
      
      expression<-saved_plots_and_tables$log2cpm[geneid,]
      
      validate(
        need(is.na(sum(expression))!=TRUE,'')
      )
      
      tsne.data<-cbind(saved_plots_and_tables$tsne.data,t(expression))
      names(tsne.data)[names(tsne.data) == geneid] <- 'values'
      
      subsetData<-subset(tsne.data,dbCluster==input$cluster)
      p1<-ggplot(subsetData,aes_string(x=input$dimension_x,y=input$dimension_y))+
        geom_point(aes_string(size=2,color='values'))+
        geom_point(shape = 1,size = 4,colour = "black")+
        theme_bw()+
        theme(axis.text.x=element_text(angle=90, size=12, vjust=0.5), 
              axis.text.y=element_text(size=12), strip.text.x = element_text(size=16), 
              strip.text.y = element_text(size=14), axis.title.x = element_text(face="bold", size=16),
              axis.title.y = element_text(face="bold", size=16),
              legend.position="none")+
        scale_colour_gradient2(low='grey50', high="red")
      p1
    })
  })
  
  output$gene_vio_plot<-renderPlot({
    
    if (v$doPlot == FALSE) return()
    
    isolate({
      
      geneid<-rownames(
        saved_plots_and_tables$featuredata[which(saved_plots_and_tables$featuredata$Associated.Gene.Name==toupper(input$gene_id)),])[1]
      
      expression<-saved_plots_and_tables$log2cpm[geneid,]
      
      validate(
        need(is.na(sum(expression))!=TRUE,'')
      )
      
      tsne.data<-cbind(saved_plots_and_tables$tsne.data,t(expression))
      names(tsne.data)[names(tsne.data) == geneid] <- 'values'
      #tsne.data<-subset(tsne.data,dbCluster!=0)
      
      p1<-p1<-ggplot(tsne.data,aes(factor(dbCluster),values,fill=factor(dbCluster)))+
        geom_violin(scale = "width")+
        stat_summary(fun.y=median, geom="point", size=5,color='black')+
        stat_summary(fun.data = n_fun, geom = "text")+
        theme_bw()+
        theme(axis.text.x=element_text(angle=90, size=12, vjust=0.5), 
              axis.text.y=element_text(size=12), strip.text.x = element_text(size=16), 
              strip.text.y = element_text(size=14), axis.title.x = element_text(face="bold", size=16),
              axis.title.y = element_text(face="bold", size=16),
              legend.position="none")+
        xlab('Cluster')+
        ylab('Expression')
      p1
    })
  })
  
  output$downloadExpression <- downloadHandler(
    
    filename = function() { paste(input$cluster,"Selected_Expression_table.csv",sep='_')},
    content = function(file) {
      
      geneid<-rownames(
        saved_plots_and_tables$featuredata[which(saved_plots_and_tables$featuredata$Associated.Gene.Name==toupper(input$gene_id)),])[1]
      
      expression<-saved_plots_and_tables$log2cpm[geneid,]
      #cat(stderr(),colnames(expression)[1:5])
      tsne.data<-cbind(saved_plots_and_tables$tsne.data,t(expression))
      #cat(file=stderr(),grep('^T_',rownames(tsne.data)))
      
      names(tsne.data)[names(tsne.data) == geneid] <- 'values'
      
      cat(file=stderr(),grep('^T_',rownames(tsne.data)))
      
      subsetData<-subset(tsne.data,dbCluster==input$cluster)
      #cat(file=stderr(),rownames(subsetData)[1:5])
      cells<-rownames(brushedPoints(subsetData,input$b1))
      cat(file=stderr(),cells[1:5])
      subsetExpression<-saved_plots_and_tables$log2cpm[,cells]
      #cat(stderr(),colnames(subsetExpression)[1:5])
      
      subsetExpression$Associated.Gene.Name<-saved_plots_and_tables$featuredata[rownames(subsetExpression),'Associated.Gene.Name']
      #cat(stderr(),colnames(subsetExpression))
      write.csv(subsetExpression, file)
    })
  
})
# VIOPLOT_SOURCE ----------------------------------------------------------

my.vioplot<-function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
                      horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
                      lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
                      at, add = FALSE, wex = 1, drawRect = TRUE) 
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
  baserange <- c(Inf, -Inf)
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
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
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
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              c(base[[i]], rev(base[[i]])), col = col[i], border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
              lty = lty)
        rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
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
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                              rev(at[i] + height[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
               boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}
