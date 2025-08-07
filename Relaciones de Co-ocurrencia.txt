# Cargar librerías necesarias
library(dplyr)
library(tidyr)
library(igraph)
library(readr)
library(Hmisc)  # Para calcular correlaciones con p-valores


# PREPARAR DATOS MICROBIOTA
# -------------------------
# Convertir Patologia a factor
datos_microbiota$Patologia <- as.factor(datos_microbiota$Patologia)

# Filtrar solo variables numéricas (géneros bacterianos en este caso) de las columnas 6 a 281
datos_num <- datos_microbiota[, 6:281]
datos_num <- datos_num[, sapply(datos_num, is.numeric)]  # Solo numéricas

# Añadir Patologia nuevamente para clasificación
datos_num$Patologia <- datos_microbiota$Patologia  

# Filtrar solo bacterias presentes en al menos 20% de las muestras, aquí lo fui cambiando (el valor de 0.3), para que me apareciesen más o menos bacterias en el gráfico
keep_taxa <- colSums(datos_num[, -ncol(datos_num)] > 0.3) >= (0.2 * nrow(datos_num))
datos_filtrados <- datos_num[, c(keep_taxa, TRUE)]  # Mantiene la variable categórica Patologia

# Separar en dos matrices: Controles y Fibromialgia, para que me represente sus correlaciones por separado
control_data <- datos_filtrados %>% filter(Patologia == "Control") %>% select(-Patologia)
fibromialgia_data <- datos_filtrados %>% filter(Patologia == "Fibromialgia") %>% select(-Patologia)


# FUNCIONES PARA CALCULAR CORRELACIÓN Y P-VALORES Y TODO ESTE FOLLÓN
# ------------------------------------------------------------------
calc_correlation <- function(data) {
  cor_results <- rcorr(as.matrix(data), type = "spearman")  # Spearman
  cor_matrix <- cor_results$r   # Matriz de correlaciones
  p_matrix <- cor_results$P     # Matriz de p-valores
  
  return(list(cor = cor_matrix, pval = p_matrix))
}

control_cor <- calc_correlation(control_data)
fibro_cor <- calc_correlation(fibromialgia_data)


# CONVERTIR LA MATRIZ (cor_matrix) A FORMATO LARGO Y FILTRAR SIGNIFICATIVOS IC 95%
# --------------------------------------------------------------------------------
format_cor_data <- function(cor_matrix, p_matrix, threshold = 0.05) {
  cor_long <- as.data.frame(as.table(cor_matrix)) %>%
    rename(node1 = Var1, node2 = Var2, correlation = Freq) %>%
    filter(node1 != node2)  # Eliminar auto-correlaciones
  
  p_long <- as.data.frame(as.table(p_matrix)) %>%
    rename(node1 = Var1, node2 = Var2, p_value = Freq) %>%
    filter(node1 != node2)
  
  # Unir correlaciones con sus p-valores
  cor_long <- left_join(cor_long, p_long, by = c("node1", "node2"))
  
  # Filtrar correlaciones significativas
  cor_long <- cor_long %>% filter(p_value < threshold)
  
  return(cor_long)
}


# Formatear datos de Control y Fibromialgia
control_network <- format_cor_data(control_cor$cor, control_cor$pval)
fibro_network <- format_cor_data(fibro_cor$cor, fibro_cor$pval)



# EXPORTAR DATOS PARA CoNet (Cytoscape)
# -------------------------------------
write_csv(control_network, "control_network.csv")
write_csv(fibro_network, "fibro_network.csv")

print("Datos exportados para CoNet en Cytoscape: 'control_network.csv' y 'fibro_network.csv'.")


# GRAFICAR LA RED BACTERIANA ASÍ CON SUS COLORITOS (igraph)
# ---------------------------------------------------------
plot_network <- function(network_data, title_text) {
  # Filtrar solo correlaciones fuertes (ajustar umbral si es necesario)
  threshold <- 0.5
  network_data <- network_data %>% filter(abs(correlation) > threshold)
  
  # Crear la red de correlación
  g <- graph_from_data_frame(network_data, directed = TRUE)
  
  # Asignar colores según si la correlación es positiva o negativa
  E(g)$color <- ifelse(E(g)$correlation > 0, "green", "red")  # verde para positiva, roja para negativa
  
  # Asignar el grosor de las aristas según la fuerza de la correlación
  E(g)$width <- abs(E(g)$correlation) * 1  # Aumenta el grosor en función de la correlación, aquí lo dejo en 1 para que no se hagan gorditas las flechas, ya aumenta el nodo, no hace falta todo grande también que queda feo
  
  # Asignar tamaño de nodos según el número de conexiones (grado del nodo)
  V(g)$size <- degree(g) * 1  # Tamaño proporcional al grado (número de conexiones), lo dejo en 1 porque sino hace pelotas de voley playa el colega
  
  # Graficar
  plot(g, vertex.label.cex = 0.7, edge.arrow.size = 0.05, 
       main = title_text, 
       vertex.label.dist = 1.5, vertex.size2 = 15, 
       edge.width = E(g)$width, 
       vertex.color = "lightyellow", edge.curved = 0.1)
}

# Graficar redes
plot_network(control_network, "Red de Co-ocurrencia - Control")  #Con esta graficamos la red de mb de los controles
plot_network(fibro_network, "Red de Co-ocurrencia - Fibromialgia")  #Con esta igual pero de FM
