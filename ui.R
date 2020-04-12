#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

secuencias <- list.files()
secuencias <- secuencias[-c(2,5,6,7,8,9)]

library(shiny)
library(DT)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("ARN - COVID19 similitudes Perú"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      
      #selectInput("seq1", "Secuencia1", choices = unique(secuencias), selected = NULL, multiple = FALSE,
       #           selectize = TRUE, width = NULL, size = NULL),
      
      selectInput("seq2", "Escoja la Secuencia a comparar:", choices = unique(secuencias), selected = NULL, multiple = FALSE,
                  selectize = TRUE, width = NULL, size = NULL),
      
      
       sliderInput("bins",
                   "Number of bins:",
                   min = 1,
                   max = 200,
                   value = 30),
      
      h5("Herramienta para comparar el ARN del COVID19 peruano con el de otros países."),
      h5("Data obtenida de GISAID.org"),
      tags$a("Usted puede encontrar el código abierto AQUI", href="https://github.com/robertoruiz911/PERU_ARN_COVID"),
      h5("Contacto: Roberto Ruiz Icochea"),
      h5("mail: robertoruizicochea@gmail.com")
    ),
    
    
    
    # Show a plot of the generated distribution
    mainPanel(

      
        
                 
                 plotOutput("distPlot1"),
                 
                 radioButtons("desconocidos", h4("Quitar enes"), choices = c("si", "no"), selected = "no"),
                 
                 verbatimTextOutput("secuencia2_nombre"),
                 
                 dataTableOutput("acidos_diferentes")
                 
      
       
       
    )
  )
))
