#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

# Funcion -----------------------------------------------------------------


plotOverlappingHist <- function(a, b, colors=c("green","red","blue"),
                                breaks=NULL, xlim=NULL, ylim=NULL){
  
  ahist=NULL
  bhist=NULL
  
  if(!(is.null(breaks))){
    ahist=hist(a,breaks=breaks,plot=F)
    bhist=hist(b,breaks=breaks,plot=F)
  } else {
    ahist=hist(a,plot=F)
    bhist=hist(b,plot=F)
    
    dist = ahist$breaks[2]-ahist$breaks[1]
    breaks = seq(min(ahist$breaks,bhist$breaks),max(ahist$breaks,bhist$breaks),dist)
    
    ahist=hist(a,breaks=breaks,plot=F)
    bhist=hist(b,breaks=breaks,plot=F)
  }
  
  if(is.null(xlim)){
    xlim = c(min(ahist$breaks,bhist$breaks),max(ahist$breaks,bhist$breaks))
  }
  
  if(is.null(ylim)){
    ylim = c(0,max(ahist$counts,bhist$counts))
  }
  
  overlap = ahist
  for(i in 1:length(overlap$counts)){
    if(ahist$counts[i] > 0 & bhist$counts[i] > 0){
      overlap$counts[i] = min(ahist$counts[i],bhist$counts[i])
    } else {
      overlap$counts[i] = 0
    }
  }
  
  plot(ahist, xlim=xlim, ylim=ylim, col=colors[1], main = "Histograma de Similitudes", xlab = "Posicion", ylab="Numero de Acidos")
  plot(bhist, xlim=xlim, ylim=ylim, col=colors[2], add=T)
  plot(overlap, xlim=xlim, ylim=ylim, col=colors[3], add=T)
}





##############################################



library(shiny)
library(data.table)
library(DT)
library(seqinr)
library(dplyr)

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
    
    cerradura1 <- read.fasta(input$seq2)
    codigo <- names(cerradura1)
    
    cerradura <- as.matrix(cerradura1[[codigo]], nrow = length(cerradura1[[codigo]]), ncol = 1)
    cerradura <- as.data.frame(cerradura)
    cerradura$V1 <- as.character(cerradura$V1)
    colnames(cerradura) <- "Secuencia 2"
    cerradura$clave <- as.numeric(as.character(row.names(cerradura)))
    
    candado <- full_join(llave, cerradura, by= "clave")
    candado$`Secuencia 1` <- as.character(candado$`Secuencia 1`)
    candado$`Secuencia 2` <- as.character(candado$`Secuencia 2`)
    
    candado$similitud <- ifelse(candado$`Secuencia 1` == candado$`Secuencia 2`, 1, 0)
    
    iguales <- subset(candado, candado$similitud == 1)
    diferentes <- subset(candado, candado$similitud == 0)
    
    mejor <- as.data.frame(t(c(length(iguales$similitud), length(diferentes$similitud))))
    mejor <- cbind(0, mejor)
    
    colnames(mejor) <- c("Calibracion", "Iguales", "Diferentes")
    
    for(i in c(-20:20)){
      
      cerradura$clave <- as.numeric(as.character(row.names(cerradura)))
      cerradura$clave <- cerradura$clave + i
      
      candado <- full_join(llave, cerradura, by= "clave")
      candado$`Secuencia 1` <- as.character(candado$`Secuencia 1`)
      candado$`Secuencia 2` <- as.character(candado$`Secuencia 2`)
      
      candado$similitud <- ifelse(candado$`Secuencia 1` == candado$`Secuencia 2`, 1, 0)
      iguales <- subset(candado, candado$similitud == 1)
      diferentes <- subset(candado, candado$similitud == 0)
      
      
      ############################################################################################
    
      
      mejor2 <- as.data.frame(t(c(length(iguales$similitud), length(diferentes$similitud))))
      mejor2 <- cbind(i, mejor2)
      
      colnames(mejor2) <- c("Calibracion", "Iguales", "Diferentes")
      
      mejor <- rbind(mejor, mejor2)
      }
    
    dup <- duplicated(mejor)
    mejor <- subset(mejor, dup == FALSE)
    
    mejor2 <- subset(mejor, mejor$Iguales == max(mejor$Iguales))
    
    cerradura$clave <- as.numeric(as.character(row.names(cerradura)))
    
    cerradura$clave <- cerradura$clave + max(mejor2$Calibracion)
    
    
    ###########################################################
    
    
    candado <- full_join(llave, cerradura, by= "clave")
    
    candado$similitud <- ifelse(candado$`Secuencia 1` == candado$`Secuencia 2`, 1, 0)
    candado$similitud <- ifelse(is.na(candado$similitud), 0, candado$similitud)
    
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
    
    cerradura1 <- read.fasta(input$seq2)
    codigo <- names(cerradura1)
    
    cerradura <- as.matrix(cerradura1[[codigo]], nrow = length(cerradura1[[codigo]]), ncol = 1)
    cerradura <- as.data.frame(cerradura)
    cerradura$V1 <- as.character(cerradura$V1)
    colnames(cerradura) <- "Secuencia 2"
    cerradura$clave <- as.numeric(as.character(row.names(cerradura)))
    
    candado <- full_join(llave, cerradura, by= "clave")
    candado$`Secuencia 1` <- as.character(candado$`Secuencia 1`)
    candado$`Secuencia 2` <- as.character(candado$`Secuencia 2`)
    
    candado$similitud <- ifelse(candado$`Secuencia 1` == candado$`Secuencia 2`, 1, 0)
    
    iguales <- subset(candado, candado$similitud == 1)
    diferentes <- subset(candado, candado$similitud == 0)
    
    
    mejor <- as.data.frame(t(c(length(iguales$similitud), length(diferentes$similitud))))
    mejor <- cbind(0, mejor)
    
    colnames(mejor) <- c("Calibracion", "Iguales", "Diferentes")
    
    for(i in c(-20:20)){
      
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
    
    
    dup <- duplicated(mejor)
    mejor <- subset(mejor, dup == FALSE)
    
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
    
    diferentes <- datatable(diferentes)
    
    print(diferentes)
    
  })
  
  
  output$secuencia2_nombre <- renderText({
    
    cerradura1 <- read.fasta(input$seq2)
    paste("La secuencia 2 es:", names(cerradura1))
    
    
  })
  
  
})



