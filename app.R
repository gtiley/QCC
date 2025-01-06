#Load dependencies
library(shiny)
library(bslib)
library(mathjaxr)
library(dplyr)
library(kableExtra)
#----------------#
#Pre-app computations
#
#----------------#

# Application UI
# Sliders for T1, T2, and gamma
# T1 [0,5]
# T2 [0,5]
# gamma [0,1]
# theta [0,0.1]

ui <- fluidPage(
    # Application title and appearance
    titlePanel(tags$b("Quartet Coalescent Calculator")),
    theme = bs_theme(preset = "journal"),
    withMathJax(),
    # Sidebar with a slider inputs
    sidebarLayout(
        sidebarPanel(
          sliderInput("T1",
                      "Coalescent unit 1:",
                      min = 0,
                      max = 5,
                      value = 1,
                      step = 0.1
          ),
          sliderInput("T2",
                      "Coalescent unit 2:",
                      min = 0,
                      max = 5,
                      value = 1,
                      step = 0.1
          ),
          sliderInput("gamma",
                      "Introgression probability:",
                      min = 0,
                      max = 1,
                      value = 0.1,
                      step = 0.01
          ),
          numericInput("theta",
                      "Nucleotide diversity:",
                      min = 0,
                      max = 1,
                      value = 0.01
          ),
          numericInput("mu",
                       "Mutation rate:",
                       min = 0,
                       max = 1,
                       value = 1e-8
          ),
          sliderInput("ploidy",
                      "Ploidy:",
                      min = 1,
                      max = 8,
                      value = 2,
                      step = 1
          )
        ),
    # Show a plot of the generated distribution
        mainPanel(
          card(
            plotOutput(outputId = "expectationPlots")
          ),
          card(
            tableOutput(outputId = "paramtable")
          ),
          card(
            uiOutput(outputId = "info")
          )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    #create plots of expectations for gene tree frequencies under MSC and NMSC
    source("./scripts/gene_tree_probs.R", local = TRUE)
    output$expectationPlots <- renderPlot({
      par(mfrow=c(1,2))
        
      par(bg=NA)
#      c1 <- p_s_msc(input$T1)
#      c2 <- p_d_msc(input$T1)
      curve(p_s_msc(x),from=0,to=5.0,xlab="",ylab="",xlim=c(0,5.0),ylim=c(0,1.0),bty='n')
      par(new=TRUE)
      curve(p_d_msc(x),from=0,to=5.0,xlab="",ylab="",xlim=c(0,5.0),ylim=c(0,1.0),bty='n',lty=2)
      segments(input$T1,0,input$T1,1,col="blue",lty=2)
      mtext(expression(paste(T[1])),side=1,cex=1.2,line=3)
      mtext("Probability",side=2,cex=1.2,line=3)
      
      par(bg=NA)
#      c3 <- p_s_nmsc(input$T2, input$T1, input$gamma)
#      c4 <- p_major_nmsc(input$T2, input$T1, input$gamma)
#      c5 <- p_minor_nmsc(input$T2, input$T1, input$gamma)
      curve(p_s_nmsc(x, input$T1, input$gamma),from=0,to=5.0,xlab="",ylab="",xlim=c(0,5.0),ylim=c(0,1.0),bty='n')
      par(new=TRUE)
      curve(p_major_nmsc(x, input$T1, input$gamma),from=0,to=5.0,xlab="",ylab="",xlim=c(0,5.0),ylim=c(0,1.0),bty='n', col="red")
      par(new=TRUE)
      curve(p_minor_nmsc(x, input$T1, input$gamma),from=0,to=5.0,xlab="",ylab="",xlim=c(0,5.0),ylim=c(0,1.0),bty='n',col="blue")
      segments(0,p_d1_msc(input$T1),5,p_d1_msc(input$T1),lty=2)
      segments(input$T2,0,input$T2,1,col="blue",lty=2)
      mtext(expression(paste(T[2])),side=1,cex=1.2,line=3)
    })
    
    #create a box of parameter values
    output$paramtable <- function (){
      ne = input$theta/((2 * input$ploidy) * input$mu)
      t1 = (input$ploidy * (input$theta/((2 * input$ploidy) * input$mu)) * input$T1)
      t2 = (input$ploidy * (input$theta/((2 * input$ploidy) * input$mu)) * input$T2)
      df <- data.frame(
        parameter <- c("Nucleotide Diversity","Mutation Rate", "Effective Population Size","Branch Length 1 Generations","Branch Length 2 Generations"),
        value <- c(input$theta, input$mu, ne, t1, t2)
      )
      knitr::kable(df, col.names = c("Parameter","Value"), digits = 16, format.args = list(scientific = TRUE)) %>% kable_styling("striped", full_width = F)
    }
    
    #create a box with some helpful interpretations of params
    output$info <- renderUI({
      theta = input$theta
      mu = input$mu
      ne = input$theta/(4 * input$mu)
      t1 = (2 * (input$theta/(4 * input$mu)) * input$T1)
      t2 = (2 * (input$theta/(4 * input$mu)) * input$T2)
      withMathJax(sprintf("The average pairwise nucleotide divergence of \\(\\theta = %s\\) is in substitutions per site and mutation rate of \\(\\mu = %s\\), is the number of mutations per site per generation. These parameters are used to calculate \\(T_{1}\\) and \\(T_{2}\\) in generations. The dashed vertical lines in plots show the user-selected values of \\(T_{1}\\) and \\(T_{2}\\).",
                          theta,
                          mu))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
