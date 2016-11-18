library(shiny)
library(DT)
library(plotly)
library(shinythemes)

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

shinyUI(navbarPage(theme = shinytheme("flatly"),
  a(
    'SingleCellBiology',
    href = 'https://www.jax.org/research-and-faculty/tools/scientific-research-services/genome-tech-single-cell-biology/single-cell-biology',
    target = '_blank'
  ),
  windowTitle = 'Single Cell Biology',
  
  tabPanel(
    'overview',
    fluidRow(div(h3('Cell View'), align = 'center')),
    br(),
    fluidRow(div(
      h5(
        'This app is designed for exploratory data analysis of
        processed RNA-Seq data of 10X single cell experiments.'
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
    fluidRow(column(
      5,offset="3",
      plotlyOutput('tsne_main')
    ))
  ),
  
  navbarMenu('Explore',
             
             tabPanel(
               'Expression',
               
               fluidRow(div(
                 p(strong('\tInformation:')),
                 tags$ul(
                   tags$li(strong('Clustering'), ':Clustering was performed with t-SNE followed by identification using DBSCAN'),
                   tags$li(strong('Cluster 0'), ':Cells that cannot be assigned to any cluster'),
                   tags$li(strong('3D Plot'), ':Enter gene name to visualize expression in a single cell'),
                   tags$li(strong('2D Plot'), ':Pick a cluster, highlight cells of interest to download gene expression matrix')
                 )
                 )
                 ),
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
                 column(
                   2,
                   uiOutput("clusters")
                 ),
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
                   div(align = "center", style = "margin-center:50px; margin-top:25px",
                       downloadButton("downloadExpression", "Download Expression"))
                 )
               ),
               
               fluidRow(
                 column(5,offset=1,
                        plotlyOutput('tsne_plt')
                 ),
                 column(5,offset=0,
                        plotOutput('clusterPlot',brush=brushOpts(id='b1'))
                 )     
               ),
               fluidRow(
                 column(10,offset=1,
                 plotOutput('gene_vio_plot')                 
               )
               )
             ))
))