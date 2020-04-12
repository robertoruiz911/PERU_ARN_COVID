
library(shiny)


secuencias <- list.files("~/Documents/corona/adn/corona/otros")

shinyUI(fluidPage(
  
  # Application title
  titlePanel("ARN - COVID19 similitudes Perú"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      
      #selectInput("seq1", "Secuencia1", choices = unique(secuencias), selected = NULL, multiple = FALSE,
       #           selectize = TRUE, width = NULL, size = NULL),
      
      selectInput("seq2", "Secuencia2", choices = unique(secuencias), selected = NULL, multiple = FALSE,
                  selectize = TRUE, width = NULL, size = NULL),
      
      
       sliderInput("bins",
                   "Number of bins:",
                   min = 1,
                   max = 200,
                   value = 30)
    
    ),
    
    
    
    # Show a plot of the generated distribution
    mainPanel(

      
        
                 
                 plotOutput("distPlot1"),
                 
                 radioButtons("desconocidos", h4("Quitar enes"), choices = c("si", "no"), selected = "no"),
                 
                 dataTableOutput("acidos_diferentes")
                
          
    )
  )
))
