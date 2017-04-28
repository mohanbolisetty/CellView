library(shiny)
library(plotly)
library(shinythemes)
library(ggplot2)
library(DT)
library(edgeR)
library(pheatmap)
library(threejs)

jscode <- '
$(function() {
var $els = $("[data-proxy-click]");
$.each(
$els,
function(idx, el) {
var $el = $(el);
var $proxy = $("#" + $el.data("proxyClick"));
$el.keydown(function (e) {
if (e.keyCode == 13) {
$proxy.click();
}
});
}
);
});
'

shinyUI(
  navbarPage(
    theme = shinytheme("flatly"),
    a('SingleCellBiology',
      href = 'https://www.jax.org/research-and-faculty/tools/scientific-research-services/genome-tech-single-cell-biology/single-cell-biology',
      target = '_blank'),
    windowTitle = 'Single Cell Biology',
    
    tabPanel(
      'overview',
      fluidRow(div(h3('CellView'), align = 'center')),
      br(),
      fluidRow(div(
        h5(
          'This app is designed for exploratory data analysis of
          processed RNA-Seq data of single cell experiments.'
        ),
        align = 'center'
      )),
      br(),
      br(),
      fluidRow(column(
        5,
        offset = 4,
        fileInput(
          'file1',
          'Choose .Rds file to upload',
          accept = c(
            '.Rds',
            'text/comma-separated-values',
            'text/tab-separated-values',
            'text/plain',
            '.csv',
            '.tsv'
          )
        )
      )),
      fluidRow(column(5, offset = "3",
                      plotlyOutput('tsne_main')))
    ),
    
    navbarMenu(
      'Explore',
      
      tabPanel(
        'Expression',
        
        fluidRow(div(
          p(strong('\tInformation:')),
          tags$ul(
            tags$li(
              strong('Clustering'),
              ':Clustering was performed with t-SNE followed by identification using DBSCAN'
            ),
            tags$li(
              strong('Cluster 0'),
              ':Cells that cannot be assigned to any cluster'
            ),
            tags$li(
              strong('3D Plot'),
              ':Enter gene name to visualize expression in a single cell'
            ),
            tags$li(
              strong('2D Plot'),
              ':Pick a cluster, highlight cells of interest to download gene expression matrix'
            )
          )
        )),
        br(),
        br(),
        fluidRow(
          column(
            2,
            tags$head(tags$script(HTML(jscode))),
            tagAppendAttributes(
              textInput('gene_id', 'Enter gene', value = 'Il6'),
              `data-proxy-click` = "goButton"
            ),
            actionButton('goButton', 'Run')
          ),
          column(2,
                 uiOutput("clusters")),
          column(
            2,
            selectInput(
              'dimension_x',
              label = 'X',
              choice = c('V1', 'V2', 'V3'),
              selected = 'V1'
            )
          ),
          column(
            2,
            selectInput(
              'dimension_y',
              label = 'Y',
              choice = c('V1', 'V2', 'V3'),
              selected = 'V2'
            )
          ),
          column(
            2,
            div(
              align = "center",
              style = "margin-center:50px; margin-top:25px",
              downloadButton("downloadExpression", "Download Expression")
            )
          )
        ),
        br(),
        br(),
        br(),
        fluidRow(
          column(5, offset = 1,
                 plotlyOutput('tsne_plt')),
          column(5, offset = 0,
                 plotOutput('clusterPlot', brush = brushOpts(id =
                                                               'b1')))
        ),
        fluidRow(column(10, offset = 1,
                        plotOutput('gene_vio_plot')))
      )
    ),
    
    navbarMenu('Co-expression',
               
               tabPanel('AllClusters',
                        
                        fluidRow(
                          column(
                            2,
                            tags$head(tags$script(HTML(jscode))),
                            tagAppendAttributes(
                              textInput('heatmap_geneids', 'Comma seperated gene names', value = 'Il6,Cd3d,Genes,Umi'),
                              `data-proxy-click` = "goButton1"
                            ),
                            actionButton('goButton1', 'Run')
                          )
                        ),
                        
                        fluidRow(
                          column(10,offset=1,
                                 plotOutput('heatmap')
                          )
                        )
               ),
               
               tabPanel('Selected cells',
                        tags$ul(
                          tags$li(
                            strong('Subclustering'),
                            ':Select a group of cells in plot1 based based on a single gene expression. Enter multiple gene ids to assess the co-expression of genes in these cells'
                          )
                          
                        ),
                        fluidRow(
                          column(
                            2,
                            tags$head(tags$script(HTML(jscode))),
                            tagAppendAttributes(
                              textInput('gene_id_sch', 'Enter gene', value = 'Il6'),
                              `data-proxy-click` = "goButton5"
                            ),
                            actionButton('goButton5', 'Run')
                          ),
                          column(2,
                                 uiOutput("clusters2")),
                          column(
                            2,
                            selectInput(
                              'dimension_x2',
                              label = 'X',
                              choice = c('V1', 'V2', 'V3'),
                              selected = 'V1'
                            )
                          ),
                          column(
                            2,
                            selectInput(
                              'dimension_y2',
                              label = 'Y',
                              choice = c('V1', 'V2', 'V3'),
                              selected = 'V2'
                            )
                          )
                        ),
                        fluidRow(
                          column(5, offset = 1,
                                 plotOutput('clusterPlot2',brush = brushOpts(id =
                                                                               'scb1')
                                 )
                          )
                        ),
                        fluidRow(
                          column(
                            2,
                            tags$head(tags$script(HTML(jscode))),
                            tagAppendAttributes(
                              textInput('heatmap_geneids2', 'Comma seperated gene names', value = 'Il6,Cd3d'),
                              `data-proxy-click` = "goButton6"
                            ),
                            actionButton('goButton6', 'Run')
                          )
                        ),
                        fluidRow(
                          column(10, offset = 1,
                                 plotOutput('selectedHeatmap' )
                          )
                        )
               )
  ),
  navbarMenu(
    'Subcluster-analysis',
    tabPanel(
      'DGE Analysis',
      tags$ul(
        tags$li(
          strong('Subclustering'),
          ':Select a group of cells in plot1 and a different group of cells in plot2 for identifying differential features between these subclusters'
        )
        
      ),
      
      fluidRow(
        column(2,
               uiOutput("clusters1")),
        column(
          2,
          selectInput(
            'dimension_x1',
            label = 'X',
            choice = c('V1', 'V2', 'V3'),
            selected = 'V1'
          )
        ),
        column(
          2,
          selectInput(
            'dimension_y1',
            label = 'Y',
            choice = c('V1', 'V2', 'V3'),
            selected = 'V2'
          )
        ),
        column(2,
               tags$head(tags$script(HTML(jscode))),
               actionButton('goButton2', 'Plot'))
      ),
      
      fluidRow(column(
        4,
        plotOutput('dge_plot1', brush = brushOpts(id =
                                                    "db1"))
      ),
      column(
        4,
        plotOutput('dge_plot2', brush = brushOpts(id =
                                                    'db2'))
      ),
      column(2,
             tags$head(tags$script(HTML(jscode))),
             actionButton('goButton3', 'Differential'))
      ),
      
      fluidRow(
        h4('Top Differentially Expressed Genes', offset = 1),
        DT::dataTableOutput('dge')
      ),
      fluidRow(
        div(align = "right", style = "margin-right:15px; margin-bottom:10px",
            downloadButton("download_dge_table", "Download DGE Table"))
        
      )
    )
  ),
  
  navbarMenu(
    'Gene set analysis',
    tabPanel(
      'KEGG',
      tags$ul(
        tags$li(
          strong('KEGG Pathway'),
          ':Select a cluster to investigate significant upregulated KEGG pathways using GAGE'
        )
      )
    ),
    tabPanel(
      'GO'
    )
  )
)
)
