### Canonical Correlation Analysis (CCA)

# Instalar y cargar la librería necesaria
install.packages("vegan")  # Si no está instalada
library(vegan)
library(car)


# 1️⃣ Preparar los datos ------------------------------------------------------
datos_continuos_todos$Patologia <- as.factor(datos_continuos_todos$Patologia)
datos_microbiota$Patologia <- as.factor(datos_microbiota$Patologia)

datos_num_continuos <- datos_continuos_todos[, 6:119]
datos_num_microb <- datos_microbiota[, 6:281]

# 🛑 Asegurar que las muestras coincidan en ambos conjuntos
comunes <- intersect(rownames(datos_microbiota), rownames(datos_continuos_todos))
datos_num_microb <- datos_microbiota[comunes, 6:281]
datos_num_continuos <- datos_continuos_todos[comunes, 6:119]

# 🛑 Eliminar variables constantes (varianza 0)
datos_num_microb <- datos_num_microb[, apply(datos_num_microb, 2, var, na.rm = TRUE) > 0]
datos_num_continuos <- datos_num_continuos[, apply(datos_num_continuos, 2, var, na.rm = TRUE) > 0]

# 🛑 Eliminar NA
datos_num_microb <- na.omit(datos_num_microb)
datos_num_continuos <- na.omit(datos_num_continuos)

# 2️⃣ Función para eliminar colinealidad (correlación superior a 0.85) -----------
eliminar_colinealidad_cor <- function(data, threshold = 0.85) {
  cor_matrix <- cor(data, use = "pairwise.complete.obs")
  cor_matrix[lower.tri(cor_matrix, diag = TRUE)] <- NA  # Evitar duplicados
  variables_a_eliminar <- c()
  
  # Buscamos correlaciones altas y eliminamos variables
  while (TRUE) {
    pares_cor <- which(abs(cor_matrix) > threshold, arr.ind = TRUE)
    if (nrow(pares_cor) == 0) break  # Si no hay pares altamente correlacionados, salir
    
    # Elegir la variable con más correlaciones altas
    vars_contadas <- table(c(rownames(cor_matrix)[pares_cor[, 1]], colnames(cor_matrix)[pares_cor[, 2]]))
    variable_a_eliminar <- names(vars_contadas)[which.max(vars_contadas)]
    
    if (!variable_a_eliminar %in% colnames(data)) break  # Evitar errores si ya se eliminó
    
    variables_a_eliminar <- c(variables_a_eliminar, variable_a_eliminar)
    cor_matrix <- cor_matrix[!rownames(cor_matrix) %in% variable_a_eliminar, 
                             !colnames(cor_matrix) %in% variable_a_eliminar]
  }
  
  print(paste(" Eliminando variables colineales:", paste(variables_a_eliminar, collapse = ", ")))
  return(data[, !colnames(data) %in% variables_a_eliminar])
}

# Eliminar colinealidad en los datos numéricos continuos
datos_num_continuos <- eliminar_colinealidad_cor(datos_num_continuos)


# 3️⃣ Realizar el análisis CCA -------------------------------------------------
modelo_cca <- cca(datos_num_microb ~ ., data = datos_num_continuos)

# 4️⃣ Resumen del modelo ------------------------------------------------------
summary(modelo_cca)

# 5️⃣ Visualizar resultados ---------------------------------------------------
plot(modelo_cca, scaling = 2, main = "CCA - Relación entre microbiota y variables clínicas")
text(modelo_cca, display = "bp", col = "red", cex = 1.2)
text(modelo_cca, display = "species", col = "blue", cex = 0.8)

# 6️⃣ Evaluar la significancia del modelo -----------------------------------
anova(modelo_cca, by = "axis", permutations = 999)
