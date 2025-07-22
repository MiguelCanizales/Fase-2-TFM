#-------------------------------------------------------------------------------
# PASO 0: INSTALAR Y CARGAR LIBRERÍAS
#-------------------------------------------------------------------------------.

if (!require("terra")) install.packages("terra")
if (!require("writexl")) install.packages("writexl") # Paquete para escribir archivos Excel
if (!require("utils")) install.packages("utils") # Para la barra de progreso
if (!require("officer")) install.packages("officer") # NUEVO: Para crear documentos de Word

library(terra)   # Librería moderna y de alto rendimiento para análisis espacial (rásters, vectores).
library(writexl) # Librería ligera y rápida para exportar datos a formato .xlsx sin dependencias.
library(utils)   # Librería base de R que contiene utilidades como la barra de progreso.
library(officer) # NUEVO: Permite la creación y manipulación de documentos .docx.

#-------------------------------------------------------------------------------
# PASO 1: CONFIGURACIÓN GENERAL
#-------------------------------------------------------------------------------
# En esta sección se definen todos los parámetros que controlan la ejecución del script.
# Modificar estos valores permite cambiar el comportamiento del proceso sin alterar la lógica principal.

# --- Directorio de Salida ---
# Define la carpeta raíz donde se guardará toda la estructura de carpetas y archivos generados.
# ¡IMPORTANTE! Asegúrate de que esta ruta sea válida en tu sistema.
base_output_dir <- "C:/Users/skull/Downloads/tfm1_Lunes" 
if (!dir.exists(base_output_dir)) {
  # Si el directorio no existe, el script lo creará recursivamente.
  dir.create(base_output_dir, recursive = TRUE)
}

# --- Configuración del Archivo de Registro (Log) ---
# Se crea un fichero de texto para guardar un registro detallado de toda la ejecución,
# incluyendo mensajes de estado, diagnósticos, advertencias y posibles errores.
log_file_path <- file.path(base_output_dir, "registro_del_proceso.log")
# Se abre una conexión al fichero en modo "append" (añadir al final).
log_file_connection <- file(log_file_path, open = "a")
# Se redirige toda la salida de la consola (tanto mensajes estándar como errores) a este fichero.
sink(log_file_connection, type = "output")
sink(log_file_connection, type = "message")

# Se escribe una cabecera en el log para marcar el inicio de una nueva ejecución.
cat(paste("\n\n==========================================================\n"))
cat(paste("--- INICIO DEL PROCESO:", Sys.time(), "---\n"))
cat(paste("==========================================================\n"))

# --- Parámetros de los MDE y la Simulación ---
# Dimensiones de los MDE sintéticos en píxeles (ancho y alto).
pixels <- 200

# Lista de posiciones solares a simular.
sun_positions <- data.frame(
  azimuth   = c(135, 180, 315),
  elevation = c(45, 75, 25)
)

# --- PARÁMETRO DE MUESTREO ACTUALIZADO ---
n_puntos_muestreo_intra_celda <- 9 # (3x3) Puntos para promediar el ángulo.

#-------------------------------------------------------------------------------
# PASO 2: FUNCIONES Y DEFINICIONES
#-------------------------------------------------------------------------------
# Esta sección contiene todas las funciones personalizadas y las estructuras de datos
# que se usarán en el script.

# --- Funciones para generar los MDE sintéticos ---
create_linear_inclined_plane <- function(n) { outer(1:n, 1:n, FUN = function(x, y) 0.5 * x + 0.3 * y) }
create_linear_wave <- function(n) { x <- 1:n; matrix(50 * sin(x / 20), nrow = n, ncol = n, byrow = TRUE) }
create_linear_ramp <- function(n) { x <- 1:n; matrix(0.8 * x, nrow=n, ncol=n, byrow=TRUE) }

create_quadratic_paraboloid <- function(n) { x <- seq(-1, 1, length.out = n); y <- seq(-1, 1, length.out = n); outer(x, y, FUN = function(x, y) 50 * x^2 + 30 * y^2) }
create_quadratic_peak <- function(n) { x <- seq(-3, 3, length.out = n); y <- seq(-3, 3, length.out = n); 150 * outer(x, y, FUN = function(x,y) exp(-(x^2+y^2))) }
create_quadratic_saddle <- function(n) { x <- seq(-1, 1, length.out = n); y <- seq(-1, 1, length.out = n); 100 * outer(x, y, FUN = function(x,y) x^2 - y^2) }

create_cubic_complex <- function(n) { x <- seq(-2, 2, length.out = n); y <- seq(-2, 2, length.out = n); mde <- outer(x, y, FUN = function(x, y) x^3 - 3*x + y^2); 50 * (mde - min(mde)) / (max(mde) - min(mde)) }
create_cubic_folds <- function(n) { x <- seq(-5, 5, length.out = n); y <- seq(-5, 5, length.out = n); 20 * outer(x,y, FUN = function(x,y) sin(x) + 0.05*y^3) }
create_cubic_asymmetric_valley <- function(n) { x <- seq(-3, 3, length.out=n); y <- seq(-3, 3, length.out=n); 10 * outer(x,y, FUN=function(x,y) 0.1*x^3 + 0.5*y^2) }

mde_generators <- list(
  lineal = list(
    lineal_plano_inclinado = create_linear_inclined_plane,
    lineal_ondulado = create_linear_wave,
    lineal_rampa = create_linear_ramp
  ),
  cuadratico = list(
    cuadratico_paraboloide = create_quadratic_paraboloid,
    cuadratico_pico = create_quadratic_peak,
    cuadratico_silla_de_montar = create_quadratic_saddle
  ),
  cubico = list(
    cubico_complejo = create_cubic_complex,
    cubico_pliegues = create_cubic_folds,
    cubico_valle_asimetrico = create_cubic_asymmetric_valley
  )
)

equations_text <- list(
  lineal = list(
    lineal_plano_inclinado = "Z = 0.5*x + 0.3*y",
    lineal_ondulado = "Z = 50 * sin(x / 20)",
    lineal_rampa = "Z = 0.8*x"
  ),
  cuadratico = list(
    cuadratico_paraboloide = "Z = 50*x^2 + 30*y^2",
    cuadratico_pico = "Z = 150 * exp(-(x^2+y^2))",
    cuadratico_silla_de_montar = "Z = 100 * (x^2 - y^2)"
  ),
  cubico = list(
    cubico_complejo = "Z = x^3 - 3*x + y^2 (escalada al rango [0, 50])",
    cubico_pliegues = "Z = 20 * (sin(x) + 0.05*y^3)",
    cubico_valle_asimetrico = "Z = 10 * (0.1*x^3 + 0.5*y^2)"
  )
)

# --- FUNCIONES DE CÁLCULO ADAPTADAS DE LA FASE 2 ---
ajustar_superficie_modelo <- function(puntos_xyz, tipo_modelo) {
  puntos_entrenamiento_requeridos <- c(lineal=3, cuadratico=6, cubico=10)
  if (nrow(puntos_xyz) < puntos_entrenamiento_requeridos[tipo_modelo]) stop("Puntos insuficientes.")
  
  puntos_xyz_centrados <- puntos_xyz
  mean_x <- mean(puntos_xyz_centrados$x)
  mean_y <- mean(puntos_xyz_centrados$y)
  
  puntos_xyz_centrados$x_c <- puntos_xyz_centrados$x - mean_x
  puntos_xyz_centrados$y_c <- puntos_xyz_centrados$y - mean_y
  
  # Se usan polinomios ortogonales para mayor estabilidad numérica
  if (tipo_modelo == "lineal") {
    formula_str <- "z ~ x_c + y_c"
  } else if (tipo_modelo == "cuadratico") {
    formula_str <- "z ~ poly(x_c, y_c, degree = 2, raw = FALSE)"
  } else if (tipo_modelo == "cubico") {
    formula_str <- "z ~ poly(x_c, y_c, degree = 3, raw = FALSE)"
  } else {
    stop("Modelo no reconocido.")
  }
  
  modelo_lm <- lm(as.formula(formula_str), data = puntos_xyz_centrados)
  attr(modelo_lm, "mean_x") <- mean_x; attr(modelo_lm, "mean_y") <- mean_y
  return(modelo_lm)
}

predecir_z_superficie <- function(nuevos_puntos_xy, modelo_ajustado_lm) {
  mean_x <- attr(modelo_ajustado_lm, "mean_x"); mean_y <- attr(modelo_ajustado_lm, "mean_y")
  data_prediccion <- data.frame(x_c = nuevos_puntos_xy$x - mean_x, y_c = nuevos_puntos_xy$y - mean_y)
  predict(modelo_ajustado_lm, newdata = data_prediccion)
}

calcular_cos_angulo_promedio_celda_sintetico <- function(puntos_xyz_celda, modelo_seleccionado_str, vector_solar, centro_xy, resolucion_xy, n_muestras = 9) {
  
  puntos_requeridos <- c(lineal = 3, cuadratico = 6, cubico = 10)
  if (is.null(puntos_xyz_celda) || is.na(modelo_seleccionado_str) || nrow(puntos_xyz_celda) < puntos_requeridos[modelo_seleccionado_str]) {
    return(NA_real_)
  }
  
  tryCatch({
    modelo_final_lm <- ajustar_superficie_modelo(puntos_xyz_celda, modelo_seleccionado_str)
    
    if (modelo_seleccionado_str == "lineal") {
      coefs <- coef(modelo_final_lm)
      dz_dx_c <- if ("x_c" %in% names(coefs)) coefs["x_c"] else 0
      dz_dy_c <- if ("y_c" %in% names(coefs)) coefs["y_c"] else 0
      normal_vec <- c(-dz_dx_c, -dz_dy_c, 1)
    } else {
      lado_grid <- sqrt(n_muestras)
      
      paso_x <- resolucion_xy[1] / lado_grid
      paso_y <- resolucion_xy[2] / lado_grid
      offset_x <- seq(-resolucion_xy[1]/2 + paso_x/2, resolucion_xy[1]/2 - paso_x/2, length.out = lado_grid)
      offset_y <- seq(-resolucion_xy[2]/2 + paso_y/2, resolucion_xy[2]/2 - paso_y/2, length.out = lado_grid)
      grid_muestreo <- expand.grid(x = centro_xy[1] + offset_x, y = centro_xy[2] + offset_y)
      
      delta <- 0.01
      
      grid_dx <- grid_muestreo; grid_dx$x <- grid_dx$x + delta
      grid_dy <- grid_muestreo; grid_dy$y <- grid_dy$y + delta
      
      z_orig <- predecir_z_superficie(grid_muestreo, modelo_final_lm)
      z_dx   <- predecir_z_superficie(grid_dx, modelo_final_lm)
      z_dy   <- predecir_z_superficie(grid_dy, modelo_final_lm)
      
      dz_dx <- (z_dx - z_orig) / delta
      dz_dy <- (z_dy - z_orig) / delta
      
      normales_mat <- cbind(-dz_dx, -dz_dy, 1)
      
      magnitudes <- sqrt(rowSums(normales_mat^2))
      magnitudes[magnitudes == 0] <- 1 
      normales_norm_mat <- sweep(normales_mat, 1, magnitudes, "/")
      
      cosenos <- normales_norm_mat %*% vector_solar
      
      return(mean(pmax(pmin(cosenos, 1), -1), na.rm = TRUE))
    }
    
    norm_magnitude <- sqrt(sum(normal_vec^2, na.rm = TRUE))
    if (norm_magnitude == 0 || is.na(norm_magnitude)) normal_vec <- c(0, 0, 1)
    else normal_vec <- normal_vec / norm_magnitude
    
    cos_theta <- sum(normal_vec * vector_solar)
    return(max(min(cos_theta, 1), -1))
    
  }, error = function(e) {
    warning(paste("Error calculando ángulo para celda:", e$message))
    return(NA_real_)
  })
}

calcular_vector_solar_rad <- function(azimuth_grados, elevation_grados) {
  az_rad <- azimuth_grados * pi / 180
  el_rad <- elevation_grados * pi / 180
  x_sol <- sin(az_rad) * cos(el_rad)
  y_sol <- cos(az_rad) * cos(el_rad)
  z_sol <- sin(el_rad)
  v_sol <- c(x_sol, y_sol, z_sol)
  return(v_sol / sqrt(sum(v_sol^2)))
}

#-------------------------------------------------------------------------------
# PASO 3: BUCLE PRINCIPAL DE PROCESAMIENTO Y EXPORTACIÓN
#-------------------------------------------------------------------------------

# --- Inicializar la barra de progreso ---
total_mdes <- sum(sapply(mde_generators, length))
cat("\nIniciando la generación de MDEs y modelos de iluminación...\n")
pb <- txtProgressBar(min = 0, max = total_mdes, style = 3, width = 50, char = "█")
progress_counter <- 0

# Bucle externo: recorre cada TIPO de MDE (lineal, cuadratico, cubico)
for (mde_type_name in names(mde_generators)) {
  
  # --- Crear la carpeta principal para la categoría de MDE ---
  mde_type_dir <- file.path(base_output_dir, mde_type_name)
  if (!dir.exists(mde_type_dir)) dir.create(mde_type_dir)
  
  # --- NUEVO: Crear un documento Word con las ecuaciones para esta categoría ---
  cat(paste("\n--- Generando documento Word con ecuaciones para la categoría:", mde_type_name, "---\n"))
  tryCatch({
    doc <- read_docx()
    doc <- body_add_par(doc, paste("Ecuaciones para Modelos de Tipo:", mde_type_name), style = "heading 1")
    variations_in_type <- mde_generators[[mde_type_name]]
    equations_in_type <- equations_text[[mde_type_name]]
    for (variation_name in names(variations_in_type)) {
      doc <- body_add_par(doc, gsub("_", " ", variation_name), style = "heading 2")
      doc <- body_add_par(doc, paste("Ecuación:", equations_in_type[[variation_name]]), style = "Normal")
      doc <- body_add_par(doc, "")
    }
    docx_filename <- file.path(mde_type_dir, paste0("ECUACIONES_", mde_type_name, ".docx"))
    print(doc, target = docx_filename)
  }, error = function(e) {
    cat(paste("      ADVERTENCIA: No se pudo crear el documento Word. Error:", e$message, "\n"))
  })
  
  mde_type_list <- mde_generators[[mde_type_name]]
  
  # Bucle intermedio: recorre cada VARIACIÓN dentro del tipo
  for (mde_variation_name in names(mde_type_list)) {
    
    # --- Crear la estructura de carpetas ---
    mde_variation_dir <- file.path(mde_type_dir, mde_variation_name)
    if (!dir.exists(mde_variation_dir)) dir.create(mde_variation_dir, recursive = TRUE)
    
    cat(paste("\n\n--- Procesando MDE:", mde_variation_name, "-> Guardando en:", mde_variation_dir, "---\n"))
    
    # --- Generar y guardar MDE base y Patrón ---
    mde_generator_function <- mde_type_list[[mde_variation_name]]
    mde_matrix <- mde_generator_function(pixels)
    mde_raster <- rast(mde_matrix, crs="local")
    writeRaster(mde_raster, file.path(mde_variation_dir, paste0("MDE_", mde_variation_name, ".tif")), overwrite = TRUE)
    
    cat("      Generando fichero patrón de referencia...\n")
    patron_raster <- aggregate(mde_raster, fact = 20, fun = "mean")
    patron_filename <- file.path(mde_variation_dir, paste0("PATRON_", mde_variation_name, ".tif"))
    writeRaster(patron_raster, patron_filename, overwrite = TRUE)
    
    # Extraer información geométrica del patrón una sola vez
    info_celdas <- as.data.frame(patron_raster, xy=TRUE)
    names(info_celdas) <- c("x_centro", "y_centro", "valor_patron")
    resolucion_patron <- res(patron_raster)
    
    # Bucle interno: recorre cada posición solar
    for (j in 1:nrow(sun_positions)) {
      
      azimuth <- sun_positions$azimuth[j]
      elevation <- sun_positions$elevation[j]
      
      illum_condition_folder_name <- paste0("Azimut_", azimuth, "_Elevacion_", elevation)
      illum_condition_dir <- file.path(mde_variation_dir, illum_condition_folder_name)
      if (!dir.exists(illum_condition_dir)) dir.create(illum_condition_dir)
      
      miramon_dir <- file.path(illum_condition_dir, "miramon")
      if (!dir.exists(miramon_dir)) dir.create(miramon_dir)
      
      cat(paste("    Calculando para Azimut:", azimuth, "y Altura:", elevation, "\n"))
      
      # --- ALGORITMO CON LA NUEVA METODOLOGÍA ---
      cat("      Ajustando superficie local para cada celda del patrón...\n")
      
      vector_solar <- calcular_vector_solar_rad(azimuth, elevation)
      illum_final_raster <- rast(patron_raster)
      values(illum_final_raster) <- NA
      
      for (k in 1:ncell(patron_raster)) {
        cell_extent <- ext(patron_raster, k)
        high_res_subset <- crop(mde_raster, cell_extent)
        points_df <- as.data.frame(high_res_subset, xy=TRUE)
        colnames(points_df) <- c("x", "y", "z")
        
        centro_actual_xy <- as.numeric(info_celdas[k, c("x_centro", "y_centro")])
        
        # Llamada a la nueva función que implementa el muestreo
        cos_i <- calcular_cos_angulo_promedio_celda_sintetico(
          points_df, 
          mde_type_name, 
          vector_solar, 
          centro_actual_xy,
          resolucion_patron,
          n_puntos_muestreo_intra_celda
        )
        
        # Se aplica la fórmula de escalado a todos los valores de cos(i)
        valor_final <- round(cos_i * 100 + 100)
        
        illum_final_raster[k] <- valor_final
      }
      
      minmax_vals <- minmax(illum_final_raster)
      cat(paste("      Diagnóstico: Valores de Iluminación (min/max) =", round(minmax_vals[1], 2), "/", round(minmax_vals[2], 2), "\n"))
      
      mde_suffix <- substr(toupper(mde_type_name), 1, 3)
      illum_tif_filename <- paste0("COSi_", mde_suffix, "_E", elevation, "_A", azimuth, ".tif")
      writeRaster(illum_final_raster, file.path(illum_condition_dir, illum_tif_filename), overwrite = TRUE, datatype = "INT1U")
      
      cat("      Generando archivo Excel para comparación...\n")
      
      df_patron <- as.data.frame(patron_raster, xy=TRUE)
      names(df_patron) <- c("Columna", "Fila", "Valor_Patron")
      df_illum <- as.data.frame(illum_final_raster, xy=TRUE)
      names(df_illum) <- c("Columna", "Fila", "Coseno_R")
      
      df_merged <- merge(df_patron, df_illum, by = c("Columna", "Fila"))
      df_merged$Pixel_ID <- paste0("F", 1:nrow(df_merged))
      df_merged$Coseno_Miramon_Esperado <- NA
      df_merged$Delta_cos_i <- NA
      
      df_comparacion <- df_merged[, c("Pixel_ID", "Fila", "Columna", "Valor_Patron", "Coseno_R", "Coseno_Miramon_Esperado", "Delta_cos_i")]
      df_comparacion <- df_comparacion[order(df_comparacion$Fila, df_comparacion$Columna), ]
      
      num_pixels <- nrow(df_comparacion)
      formula_mae <- paste0("=AVERAGE(ABS(G2:G", num_pixels + 1, "))")
      formula_rmse <- paste0("=SQRT(SUMSQ(G2:G", num_pixels + 1, ")/COUNT(G2:G", num_pixels + 1, "))")
      
      metrics_rows <- data.frame(
        Pixel_ID = c("", "Métricas de Validación", "MAE (Mean Absolute Error)", "RMSE (Root Mean Square Error)"),
        Fila = c("", "-------------------", "", ""), Columna = c("", "", "", ""), Valor_Patron = c("", "", "", ""),
        Coseno_R = c("", "", "", ""), Coseno_Miramon_Esperado = c("", "", "Valor Calculado:", "Valor Calculado:"),
        Delta_cos_i = c("", "", formula_mae, formula_rmse)
      )
      
      df_comparacion_char <- as.data.frame(lapply(df_comparacion, as.character))
      names(metrics_rows) <- names(df_comparacion_char)
      df_final_sheet <- rbind(df_comparacion_char, metrics_rows)
      
      df_matriz_mde <- as.data.frame(mde_matrix)
      
      sheets <- list( "Tabla_Comparacion_y_Metricas" = df_final_sheet, "Matriz_MDE_200x200" = df_matriz_mde )
      
      excel_filename <- file.path(illum_condition_dir, paste0("COMPARACION_", mde_variation_name, ".xlsx"))
      write_xlsx(sheets, excel_filename)
    }
    
    progress_counter <- progress_counter + 1
    setTxtProgressBar(pb, progress_counter)
  }
}

# --- Cerrar la barra de progreso ---
close(pb)

cat(paste("\n\n==========================================================\n"))
cat(paste("--- PROCESO COMPLETADO:", Sys.time(), "---\n"))
cat(paste("--- Revisa la carpeta:", base_output_dir, "---\n"))
cat(paste("==========================================================\n"))

#-------------------------------------------------------------------------------
# PASO 4: LIMPIEZA FINAL
#-------------------------------------------------------------------------------
# --- Restaurar la salida a la consola ---
sink(type = "message")
sink()
close(log_file_connection)

# Imprime un mensaje final en la consola (esto ya no se guardará en el log)
message("Proceso finalizado. El registro se ha guardado en: ", log_file_path)
