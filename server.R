#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#




##############################################



library(shiny)
library(data.table)
library(DT)
# Define server logic required to draw a histogram

shinyServer(function(input, output) {
  
  
  output$distPlot1 <- renderPlot({
    
  #secuencia1 <- reactive({
    
    llave1 <- read.fasta("EPI_ISL_415787.fasta")
    
    codigo <- names(llave1)
    
    llave <- as.matrix(llave1[[codigo]], nrow = length(llave1[[codigo]]), ncol = 1)
    
    llave <- as.data.frame(llave)
    llave$V1 <- as.character(llave$V1)
    colnames(llave) <- "Secuencia 1"
    
    llave$clave <- as.numeric(as.character(row.names(llave)))
    
    
    ##########################################################
    
    cerradura1 <- read.fasta(paste("~/Documents/corona/adn/corona/otros/", input$seq2, sep =""))
    codigo <- names(cerradura1)
    
    cerradura <- as.matrix(cerradura1[[codigo]], nrow = length(cerradura1[[codigo]]), ncol = 1)
    cerradura <- as.data.frame(cerradura)
    cerradura$V1 <- as.character(cerradura$V1)
    colnames(cerradura) <- "Secuencia 2"
    
    
    mejor <- as.data.frame(t(c(length(iguales$similitud), length(diferentes$similitud))))
    mejor <- cbind(0, mejor)
    
    colnames(mejor) <- c("Calibracion", "Iguales", "Diferentes")
    
    for(i in c(-20:-1, 1:20)){
      
      cerradura$clave <- as.numeric(as.character(row.names(cerradura)))
      cerradura$clave <- cerradura$clave + i
      
      candado <- full_join(llave, cerradura, by= "clave")
      candado$`Secuencia 1` <- as.character(candado$`Secuencia 1`)
      candado$`Secuencia 2` <- as.character(candado$`Secuencia 2`)
      
      candado$similitud <- ifelse(candado$`Secuencia 1` == candado$`Secuencia 2`, 1, 0)
      
      
      ############################################################################################
      
      diferentes <- subset(candado, candado$similitud == 0)
      iguales <- subset(candado, candado$similitud == 1)
      
      mejor2 <- as.data.frame(t(c(length(iguales$similitud), length(diferentes$similitud))))
      mejor2 <- cbind(i, mejor2)
      
      colnames(mejor2) <- c("Calibracion", "Iguales", "Diferentes")
      
      mejor <- rbind(mejor, mejor2)
      
      
    }
    
    mejor2 <- subset(mejor, mejor$Iguales == max(mejor$Iguales))
    
    cerradura$clave <- as.numeric(as.character(row.names(cerradura)))
    
    cerradura$clave <- cerradura$clave + max(mejor2$Calibracion)
    
    
    ###########################################################
    
    
    candado <- full_join(llave, cerradura, by= "clave")
    
    candado
    
    
    candado$similitud <- ifelse(candado$`Secuencia 1` == candado$`Secuencia 2`, 1, 0)
    candado$similitud <- ifelse(is.na(candado$similitud), 0, candado$similitud)
    
    
 # })
  
  
   

    
   # candado <- secuencia1()
    
    diferentes <- subset(candado, candado$similitud == 0)
    iguales <- subset(candado, candado$similitud == 1)
    
    plotOverlappingHist(iguales$clave, diferentes$clave, breaks=input$bins)
    
  })
  
  
  
  
  output$acidos_diferentes <- renderDataTable({
    
    
    llave1 <- read.fasta("EPI_ISL_415787.fasta")
    
    codigo <- names(llave1)
    
    llave <- as.matrix(llave1[[codigo]], nrow = length(llave1[[codigo]]), ncol = 1)
    
    llave <- as.data.frame(llave)
    llave$V1 <- as.character(llave$V1)
    colnames(llave) <- "Secuencia 1"
    
    llave$clave <- as.numeric(as.character(row.names(llave)))
    
    
    ##########################################################
    
    cerradura1 <- read.fasta(paste("~/Documents/corona/adn/corona/otros/",input$seq2, sep = ""))
    codigo <- names(cerradura1)
    
    cerradura <- as.matrix(cerradura1[[codigo]], nrow = length(cerradura1[[codigo]]), ncol = 1)
    cerradura <- as.data.frame(cerradura)
    cerradura$V1 <- as.character(cerradura$V1)
    
    colnames(cerradura) <- "Secuencia 2"
    
    
    mejor <- as.data.frame(t(c(length(iguales$similitud), length(diferentes$similitud))))
    mejor <- cbind(0, mejor)
    
    colnames(mejor) <- c("Calibracion", "Iguales", "Diferentes")
    
    for(i in c(-20:-1, 1:20)){
      
      cerradura$clave <- as.numeric(as.character(row.names(cerradura)))
      cerradura$clave <- cerradura$clave + i
      
      candado <- full_join(llave, cerradura, by= "clave")
      candado$`Secuencia 1` <- as.character(candado$`Secuencia 1`)
      candado$`Secuencia 2` <- as.character(candado$`Secuencia 2`)
      
      candado$similitud <- ifelse(candado$`Secuencia 1` == candado$`Secuencia 2`, 1, 0)
      
      
      ############################################################################################
      
      diferentes <- subset(candado, candado$similitud == 0)
      iguales <- subset(candado, candado$similitud == 1)
      
      mejor2 <- as.data.frame(t(c(length(iguales$similitud), length(diferentes$similitud))))
      mejor2 <- cbind(i, mejor2)
      
      colnames(mejor2) <- c("Calibracion", "Iguales", "Diferentes")
      
      mejor <- rbind(mejor, mejor2)
      }
    
    mejor2 <- subset(mejor, mejor$Iguales == max(mejor$Iguales))
    
    cerradura$clave <- as.numeric(as.character(row.names(cerradura)))
    
    cerradura$clave <- cerradura$clave + max(mejor2$Calibracion)
    
    ###########################################################
    
    candado <- full_join(llave, cerradura, by = "clave")

    candado$similitud <- ifelse(candado$`Secuencia 1` == candado$`Secuencia 2`, 1, 0)
    candado$similitud <- ifelse(is.na(candado$similitud), 0, candado$similitud)
    
    diferentes <- subset(candado, candado$similitud == 0)
  
    diferentes$Posicion_Peru <- row.names(diferentes)
    diferentes$similitud <- NULL
    colnames(diferentes) <- c("ARN-Peru", "Posicion_calibrada", "ARN2", "Posicion_Peru")
    diferentes$Posicion_Pais2 <- diferentes$Posicion_calibrada - max(mejor2$Calibracion)
    diferentes <- diferentes[,c("Posicion_Peru", "Posicion_Pais2", "Posicion_calibrada", "ARN-Peru",  "ARN2" )]
    
    #diferentes <- ifelse(input$desconocidos == "si", subset(diferentes, !(diferentes$ARN2 %in% "n")), diferentes)
    
    diferentes <- as.data.table(diferentes)
    
    diferentes
    
  })
  
})
