
library(plsgenomics)

data("leukemia")
leukaemia <- leukemia




plot_heatmap <- function(huge_gen, p){
  dimnames(huge_gen$omega) = list(c(sapply(1:p, function(x) paste("Gene", x)))
                                  ,sapply(1:p, function(x) paste("Gene", x)))
  d3heatmap(huge_gen$omega, scale = "column", dendrogram = "none",
            colors = "Greys")
}
#___________________________________________
transform_data_to_plot <- function(huge_gen, p, g){
 # group <- rep(1,p)
  each_ <-(floor(p/g)+1)
  min <- each_ * g - p  
  group <- rep(1:g, each = each_)[-((1:min)*(each_))]
  melted_matrix <- melt(huge_gen$omega)
  melted_matrix <- melted_matrix %>%
    filter(value > 0.1 & as.numeric(X1)<as.numeric(X2))
  names <- as.factor(sapply(1:p, function(x) paste("Gene", x)))
  node <- data.frame (ID = 1:p,names, group =group, size = 5)
  link <- data.frame(source = as.integer(melted_matrix$X1) - 1,
                     target = as.integer(melted_matrix$X2) - 1,
                     value = melted_matrix$value)
  return(list(g_Links = link, g_Nodes = node))
}
#___________________________________________
plot_network <- function(huge_gen, p, g){
  d_t_plot <- transform_data_to_plot(huge_gen, p, g)
  forceNetwork(Links = d_t_plot$g_Links, Nodes = d_t_plot$g_Nodes, Source = "source",
               Target = "target", Value = "value", NodeID = "names",
               Group = "group", opacity = 0.9)
}
#_________________________________________

transform_data_to_plot_emp <- function(out_huge, p , lambda_n, g){
  each_ <-(floor(p/g)+1)
  min <- each_ * g - p  
  group <- rep(1:g, each = each_)[-((1:min)*(each_))]
  ifelse(is.null(out_huge$icov), 
         out_huge_icov <- as(out_huge$path[[lambda_n]], "matrix"),
         out_huge_icov <- as(out_huge$icov[[lambda_n]], "matrix"))
  melted_matrix <- melt(out_huge_icov)
  melted_matrix <- melted_matrix %>%
    filter(value > 0.05 & as.numeric(X1)<as.numeric(X2))
  names <- as.factor(sapply(1:p, function(x) paste("Gene", x)))
  node <- data.frame (ID = 1:p,names, group =group, size = 5)
  link <- data.frame(source = as.integer(melted_matrix$X1) - 1,
                     target = as.integer(melted_matrix$X2) - 1,
                     value = melted_matrix$value)
  return(list(g_Links = link, g_Nodes = node))
}
#____________________________________
plot_network_emp <- function(out_huge, p, lambda_n, g){
  d_t_plot <- transform_data_to_plot_emp(out_huge, p, lambda_n, g)
  forceNetwork(Links = d_t_plot$g_Links, Nodes = d_t_plot$g_Nodes, Source = "source",
               Target = "target", Value = "value", NodeID = "names",
               Group = "group", opacity = 0.9)
}





#______________________________________________________________________________________part 3
nazwa<-function(data,od=1,do=10,czynazwy=TRUE){
  if (czynazwy==TRUE){
  if (data=="leukaemia") {  nazwa<-leukaemia$gene.names[od:do,3]}
  if (data=="SRBCT") {  nazwa<-SRBCT$gene.names[od:do,2]}  
  if (data=="Colon") {  nazwa<-Colon$gene.names[od:do]}}
  if (czynazwy==FALSE){nazwa<- NULL}
  return(nazwa)
}

dobor_modelu<-function(data,od=1, do=10, method="glasso",nlambda = 10){
  if (data=="leukaemia") {  L<-leukaemia$X[,od:do]}
  if (data=="SRBCT") {  L<-SRBCT$X[,od:do]}  
  if (data=="Colon") {  L<-Colon$X[,od:do]} 
  
  if (method=="glasso"){
    model<-huge(L,method="glasso",nlambda =nlambda)  }
  if (method=="mb")
    model<-huge(L,method="mb",nlambda =nlambda)  
  if (method=="clime")
    model<-sugm(L,method="clime",nlambda =nlambda)  
  if (method=="tiger")
    model<-sugm(L ,method="tiger")  
  return(model)
}

plot_heatmap1 <- function(data, lambda_n = 4, names = NULL){
  if(!is.null(data$icov)) 
    omega <- as(data$icov[[lambda_n]], "matrix")
  else if(!is.null(data$path)){
    omega <- as(data$path[[lambda_n]], "matrix")
    p <- dim(omega)[1]
    print(p)
    omega <- omega + diag(p)
  }
  else if(!is.null(data$omega))
    omega <- data$omega
  else
    message("Wrong imput: Try imput data type huge generator or output
            model type huge or flare")
  p <- dim(omega)[1]
  if(is.null(names))
    names <- sapply(1:p, function(x) paste("Gene", x))
  print(names)
  dimnames(omega) = list(names
                         ,names)
  d3heatmap(omega, scale = "column", dendrogram = "none",
            colors = "Greys", zoom = T)
}



plot_network1 <- function(data, g = 0, lambda_n = 24, names = NULL, k = 1){
  # Get the right structure
  if(!is.null(data$icov)) 
    omega <- as(data$icov[[lambda_n]], "matrix")
  else if(!is.null(data$path)){
    omega <- as(data$path[[lambda_n]], "matrix")
    p <- dim(omega)[1]
    omega <- omega + diag(p)
  }
  else if(!is.null(data$omega))
    omega <- data$omega
  else
    message("Wrong imput: Try imput data type huge generator or output
            model type huge or flare")
  # Przypisywanie etykietek grup
  p <- dim(omega)[1]
  if (g>0){
    each_ <-(floor(p/g)+1)
    min <- each_ * g - p  
    group <- rep(1:g, each = each_)[-((1:min)*(each_))]
    print(group)
  }
  #Prepare matrix of links
  melted_matrix <- melt(omega)
  melted_matrix <- melted_matrix %>%
    filter((value > 0.01 | value < -0.01) & as.numeric(X1) < as.numeric(X2))
  if(g==0){
    g1 <- graph_color(melted_matrix$X1, melted_matrix$X2, names = 1:p, k = k)
    group <-  as.data.frame(g1[,3])
    group <- group$value2
    print(as.vector(group))
  }
  if(is.null(names))
    names <- as.factor(sapply(1:p, function(x) paste("Gene", x)))
  #Assign nodes and link
  node <- data.frame (ID = 1:p,names, group =group, size = 5)
  link <- data.frame(source = as.integer(melted_matrix$X1) - 1,
                     target = as.integer(melted_matrix$X2) - 1,
                     value = melted_matrix$value)
  #Plot the graph
  forceNetwork(Links = link, Nodes = node, Source = "source",
               Target = "target", Value = "value", NodeID = "names",
               Group = "group", opacity = 0.9, zoom = T)
}


X1 <-c(1,2,3,4,5 ,5, 5, 7)
X2 <-c(0,1,2,2,1,3,4, 4)
library(dplyr)
tmp <- table(c(X1,X2))
tmp1 <- data.frame(name = as.numeric(names(tmp)), value = as.numeric(tmp))
df <- data.frame(name = 1:7)
qq <- df %>% left_join(tmp1, by = "name")


graph_color <- function(X1, X2, names, k = 3){
  cum_table <- table(c(X1,X2))
  group <- data.frame(name = as.numeric(names(cum_table)), value = as.numeric(cum_table))
  df <- data_frame(name = names)
  df <- df %>% left_join(group, by = "name")
  df$value[is.na(df$value)] <- 0 
  df$value2 <- round(df$value/k,0)
  return(df)
}


graph_color(X1,X2,1:10,2)
