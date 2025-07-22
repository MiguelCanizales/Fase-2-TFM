# --- Script R Integral para Fase 2 TFM: Método Autoadaptativo ---

# --- 0. Configuración Inicial y Carga de Librerías ---
cat("INFO: Cargando librerías...\n")
library(terra)
library(dplyr)
library(purrr)
library(stats)
library(logger)
library(knitr)
library(rmarkdown)
library(ggplot2)
library(DT)
library(viridisLite)
library(utils)
library(tidyr)
library(future)
library(furrr)

# --- CONFIGURACIÓN DE TERRA PARA MEJOR RENDIMIENTO ---
terraOptions(memfrac = 0.8) 

# --- Parámetros Globales ---
cat("INFO: Definiendo parámetros globales...\n")
ruta_mde_tif <- "C:/Users/Miguel C/Downloads/Fase 2 TFM/MDE_8/MDE_ICGC_2m.tif"
area_prueba_pix_definida <- NULL 

# Usaremos el patrón de 10m. Factor = 10m / 2m = 5
factor_agregacion_patron <- 5

num_repeticiones_mc <- 10
solar_azimuth_grados <- 260
solar_elevation_grados <- 60

n_muestras_intracelda <- 16

ruta_salida_base <- "C:/Users/Miguel C/Downloads/Fase 2 TFM/MDE_8/19_Prueba_Final/"
dir.create(ruta_salida_base, showWarnings = FALSE, recursive = TRUE)

# Ficheros de salida
fichero_log_txt <- file.path(ruta_salida_base, "registro_proceso_R.txt")
fichero_resultados_csv <- file.path(ruta_salida_base, "resultados_celdas_R.csv")
fichero_mapa_modelos_tif <- file.path(ruta_salida_base, "mapa_modelos_seleccionados.tif")
fichero_mapa_angulos_autoadapt_tif <- file.path(ruta_salida_base, "mapa_angulos_autoadaptativo_scaled.tif")
fichero_mapa_angulos_global_lineal_tif <- file.path(ruta_salida_base, "mapa_angulos_local_lineal_scaled.tif")
fichero_mapa_angulos_global_cuadratico_tif <- file.path(ruta_salida_base, "mapa_angulos_local_cuadratico_scaled.tif")
fichero_mapa_angulos_global_cubico_tif <- file.path(ruta_salida_base, "mapa_angulos_local_cubico_scaled.tif")
fichero_patron_tif <- file.path(ruta_salida_base, "mapa_patron_celdas_alineado_10m_extent.tif") 
fichero_mapa_mascara_lineal_tif <- file.path(ruta_salida_base, "mapa_mascara_lineal.tif")
fichero_mapa_mascara_cuadratico_tif <- file.path(ruta_salida_base, "mapa_mascara_cuadratico.tif")
fichero_mapa_mascara_cubico_tif <- file.path(ruta_salida_base, "mapa_mascara_cubico.tif")
fichero_informe_rmd <- file.path(ruta_salida_base, "informe_resultados_tfm_fase2_R.Rmd")
fichero_informe_html <- file.path(ruta_salida_base, "informe_resultados_tfm_fase2_R.html")

# --- 1. Funciones de Registro (Logging) ---
if (exists("log_remove_appender") && "package:logger" %in% search() ) {
  try(log_remove_appender(1), silent = TRUE)
}
log_appender(appender_file(fichero_log_txt))
log_threshold(INFO)
log_info("--- INICIO DEL PROCESO TFM FASE 2 (VERSIÓN OPTIMIZADA Y MEJORADA) ---")

# --- 2. Funciones de Carga y Creación de Patrón ---
cargar_mde <- function(ruta_mde_tif_func, area_prueba_pix_func = NULL) {
  log_info(paste("Intentando cargar MDE desde:", ruta_mde_tif_func))
  tryCatch({
    if (!file.exists(ruta_mde_tif_func)) {
      log_fatal(paste("El fichero MDE no existe en la ruta:", ruta_mde_tif_func))
      stop(paste("El fichero MDE no existe en la ruta:", ruta_mde_tif_func))
    }
    mde_r_completo <- rast(ruta_mde_tif_func)
    log_info(paste("MDE completo cargado. Dimensiones (filas,cols):", nrow(mde_r_completo), "x", ncol(mde_r_completo)))
    
    if (!is.null(area_prueba_pix_func)) {
      log_info(paste("Recortando MDE a un área de prueba de", area_prueba_pix_func[1], "x", area_prueba_pix_func[2], "píxeles."))
      ext_recorte <- ext(mde_r_completo, 1, area_prueba_pix_func[2], 1, area_prueba_pix_func[1])
      mde_r_completo <- crop(mde_r_completo, ext_recorte)
      log_info(paste("MDE de prueba. Dimensiones (filas,cols):", nrow(mde_r_completo), "x", ncol(mde_r_completo)))
    }
    return(mde_r_completo)
  }, error = function(e) {
    log_fatal(paste("Error fatal al cargar MDE:", e$message))
    stop("Fallo al cargar MDE.")
  })
}

crear_patron_raster <- function(mde_raster_func, factor_agregacion_func, nombre_fichero_salida_func) {
  log_info(paste("Creando RASTER patrón con método optimizado (Factor:", factor_agregacion_func, ")"))
  tryCatch({
    patron_recortado <- aggregate(mde_raster_func, fact = factor_agregacion_func, fun = "mean", na.rm = TRUE)
    plantilla_final <- rast(
      extent = ext(mde_raster_func),
      resolution = res(mde_raster_func) * factor_agregacion_func,
      crs = crs(mde_raster_func)
    )
    patron_final_alineado <- resample(patron_recortado, plantilla_final, method = "near")
    values(patron_final_alineado) <- 1:ncell(patron_final_alineado)
    names(patron_final_alineado) <- "id_celda_patron"
    writeRaster(patron_final_alineado, nombre_fichero_salida_func, overwrite=TRUE, datatype="INT4S")
    log_info(paste("Mapa de patrón con EXTENSIÓN EXACTA guardado en:", nombre_fichero_salida_func))
    return(patron_final_alineado)
  }, error = function(e) {
    log_error(paste("Error al crear el ráster patrón:", e$message))
    stop("Fallo al crear el ráster patrón.")
  })
}

# --- 3 y 4. Funciones de Ajuste, Error y Procesamiento ---
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

calcular_rmse <- function(z_reales, z_predichas) {
  if(length(z_reales)!=length(z_predichas) || length(z_reales)==0) return(NA_real_)
  sqrt(mean((z_reales - z_predichas)^2, na.rm = TRUE))
}

calcular_mae <- function(z_reales, z_predichas) {
  if(length(z_reales)!=length(z_predichas) || length(z_reales)==0) return(NA_real_)
  mean(abs(z_reales - z_predichas), na.rm = TRUE)
}

procesar_grupo_de_puntos <- function(puntos_celda_actual, num_repeticiones) {
  id_celda <- puntos_celda_actual$id_celda_patron[1]
  n_puntos_totales_celda <- nrow(puntos_celda_actual)
  puntos_entrenamiento_fijos <- c(lineal=3, cuadratico=6, cubico=10)
  min_puntos_test <- 1
  
  if (n_puntos_totales_celda < (puntos_entrenamiento_fijos["lineal"] + min_puntos_test)) {
    return(tibble::tibble(id_celda_patron = id_celda, estado_celda = "descartada", razon_descarte = "puntos_insuficientes_inicial"))
  }
  
  resultados_rmse_celda <- list()
  resultados_mae_celda <- list()
  
  for (tipo_modelo_actual in c("lineal", "cuadratico", "cubico")) {
    n_entrenamiento <- puntos_entrenamiento_fijos[tipo_modelo_actual]
    if (n_puntos_totales_celda < (n_entrenamiento + min_puntos_test)) next
    
    errores_rmse_iter <- numeric(num_repeticiones)
    errores_mae_iter <- numeric(num_repeticiones)
    
    for (i in 1:num_repeticiones) {
      tryCatch({
        indices_entrenamiento <- sample(1:n_puntos_totales_celda, n_entrenamiento)
        puntos_entrenamiento <- puntos_celda_actual[indices_entrenamiento, ]
        puntos_prueba <- puntos_celda_actual[-indices_entrenamiento, ]
        
        modelo_ajustado <- ajustar_superficie_modelo(puntos_entrenamiento, tipo_modelo_actual)
        z_predichas_prueba <- predecir_z_superficie(puntos_prueba, modelo_ajustado)
        
        errores_rmse_iter[i] <- calcular_rmse(puntos_prueba$z, z_predichas_prueba)
        errores_mae_iter[i] <- calcular_mae(puntos_prueba$z, z_predichas_prueba)
        
      }, error = function(e) { 
        errores_rmse_iter[i] <- NA_real_ 
        errores_mae_iter[i] <- NA_real_
      })
    }
    
    if (any(!is.na(errores_rmse_iter))) {
      resultados_rmse_celda[[tipo_modelo_actual]] <- mean(errores_rmse_iter, na.rm = TRUE)
      resultados_mae_celda[[tipo_modelo_actual]] <- mean(errores_mae_iter, na.rm = TRUE)
    }
  }
  
  if (length(resultados_rmse_celda) == 0) {
    return(tibble::tibble(id_celda_patron = id_celda, estado_celda = "descartada", razon_descarte = "ningun_modelo_valido"))
  }
  
  modelo_seleccionado <- names(which.min(unlist(resultados_rmse_celda)))
  rmse_mejor_modelo <- min(unlist(resultados_rmse_celda), na.rm = TRUE)
  mae_mejor_modelo <- resultados_mae_celda[[modelo_seleccionado]]
  
  return(tibble::tibble(
    id_celda_patron = id_celda,
    estado_celda = "ajustada", 
    razon_descarte = NA_character_,
    modelo_seleccionado = modelo_seleccionado, 
    rmse_mejor_modelo = rmse_mejor_modelo,
    mae_mejor_modelo = mae_mejor_modelo
  ))
}

# --- 5. Funciones de Cálculo de Ángulos (MODIFICADAS) ---
calcular_vector_solar_rad <- function(azimuth_grados, elevation_grados) {
  az_rad <- azimuth_grados * pi / 180
  el_rad <- elevation_grados * pi / 180
  x_sol <- sin(az_rad) * cos(el_rad)
  y_sol <- cos(az_rad) * cos(el_rad)
  z_sol <- sin(el_rad)
  v_sol <- c(x_sol, y_sol, z_sol)
  return(v_sol / sqrt(sum(v_sol^2)))
}

# --- NUEVA FUNCIÓN UNIFICADA Y MEJORADA ---
# Se utiliza diferenciación numérica con predict() para un cálculo más robusto de la normal.
calcular_cos_angulo_promedio_celda <- function(puntos_xyz_celda, modelo_seleccionado_str, vector_solar, centro_xy, resolucion_xy, n_muestras = 9) {
  
  puntos_requeridos <- c(lineal = 3, cuadratico = 6, cubico = 10)
  if (is.null(puntos_xyz_celda) || is.na(modelo_seleccionado_str) || nrow(puntos_xyz_celda) < puntos_requeridos[modelo_seleccionado_str]) {
    return(NA_real_)
  }
  
  tryCatch({
    # Se usa el centro de la media de los datos para el ajuste del modelo
    modelo_final_lm <- ajustar_superficie_modelo(puntos_xyz_celda, modelo_seleccionado_str)
    
    if (modelo_seleccionado_str == "lineal") {
      # Para un plano, la normal es constante y se puede derivar de los coeficientes.
      coefs <- coef(modelo_final_lm)
      dz_dx_c <- if ("x_c" %in% names(coefs)) coefs["x_c"] else 0
      dz_dy_c <- if ("y_c" %in% names(coefs)) coefs["y_c"] else 0
      normal_vec <- c(-dz_dx_c, -dz_dy_c, 1)
    } else {
      # Para superficies curvas, usamos diferenciación numérica en una malla geométrica regular.
      lado_grid <- sqrt(n_muestras)
      
      paso_x <- resolucion_xy[1] / lado_grid
      paso_y <- resolucion_xy[2] / lado_grid
      offset_x <- seq(-resolucion_xy[1]/2 + paso_x/2, resolucion_xy[1]/2 - paso_x/2, length.out = lado_grid)
      offset_y <- seq(-resolucion_xy[2]/2 + paso_y/2, resolucion_xy[2]/2 - paso_y/2, length.out = lado_grid)
      grid_muestreo <- expand.grid(x = centro_xy[1] + offset_x, y = centro_xy[2] + offset_y)
      
      delta <- 0.01 # Pequeño diferencial para el cálculo numérico de la derivada
      
      grid_dx <- grid_muestreo; grid_dx$x <- grid_dx$x + delta
      grid_dy <- grid_muestreo; grid_dy$y <- grid_dy$y + delta
      
      # Predecir Z para los tres sets de puntos
      z_orig <- predecir_z_superficie(grid_muestreo, modelo_final_lm)
      z_dx   <- predecir_z_superficie(grid_dx, modelo_final_lm)
      z_dy   <- predecir_z_superficie(grid_dy, modelo_final_lm)
      
      # Calcular derivadas parciales
      dz_dx <- (z_dx - z_orig) / delta
      dz_dy <- (z_dy - z_orig) / delta
      
      normales_mat <- cbind(-dz_dx, -dz_dy, 1)
      
      magnitudes <- sqrt(rowSums(normales_mat^2))
      magnitudes[magnitudes == 0] <- 1 
      normales_norm_mat <- sweep(normales_mat, 1, magnitudes, "/")
      
      cosenos <- normales_norm_mat %*% vector_solar
      
      return(mean(pmax(pmin(cosenos, 1), -1), na.rm = TRUE))
    }
    
    # Este bloque solo se ejecuta para el caso lineal
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


# ===================================================================
# --- 6. BUCLE PRINCIPAL DE PROCESAMIENTO (VERSIÓN PARALELIZADA) ---
# ===================================================================
log_info("--- Sección 6: Bucle Principal Optimizado y en Paralelo ---")

mde_activo <- cargar_mde(ruta_mde_tif, area_prueba_pix_definida)
patron_raster <- crear_patron_raster(mde_activo, factor_agregacion_patron, fichero_patron_tif)
log_info("MDE y RASTER patrón preparados.")

log_info("Preparando datos para el análisis (operación única)...")
patron_upsampled <- terra::resample(patron_raster, mde_activo, method = "near")
names(mde_activo) <- "z"
pila_combinada <- c(mde_activo, patron_upsampled)
datos_extraidos_df <- terra::as.data.frame(pila_combinada, xy = TRUE, na.rm = TRUE)
log_info(paste("Se extrajeron", nrow(datos_extraidos_df), "puntos del MDE para el análisis."))
rm(pila_combinada, patron_upsampled)

# --- CONFIGURACIÓN PARA EJECUTAR EN PARALELO ---
workers_a_usar <- max(1, future::availableCores() - 1)
plan(multisession, workers = workers_a_usar)
log_info(paste("Iniciando procesamiento en paralelo con", workers_a_usar, "núcleos."))

# Dividimos el dataframe en una lista, donde cada elemento es el dataframe de una celda
lista_de_celdas <- datos_extraidos_df %>%
  group_by(id_celda_patron) %>%
  group_split()

# --- 6.1 Selección de modelos en paralelo ---
log_info("Fase 6.1: Seleccionando el mejor modelo para cada celda en paralelo...")
df_resultados_final <- future_map_dfr(lista_de_celdas, ~ {
  procesar_grupo_de_puntos(.x, num_repeticiones_mc)
}, .progress = TRUE, .options = furrr_options(seed = TRUE, globals = c("procesar_grupo_de_puntos", "ajustar_superficie_modelo", "predecir_z_superficie", "calcular_rmse", "calcular_mae", "num_repeticiones_mc")))
log_info("Selección de modelos completada.")


# =============================================================================
# --- 7. POST-PROCESAMIENTO: GENERACIÓN DE MAPAS Y RESULTADOS FINALES ---
# =============================================================================
log_info("--- Sección 7: Post-Procesamiento y Generación de Salidas (en paralelo) ---")

# --- 7.1 Cálculo de Ángulos de Incidencia (Autoadaptativo) ---
log_info("Fase 7.1: Calculando cosenos para el método autoadaptativo en paralelo...")
vector_solar <- calcular_vector_solar_rad(solar_azimuth_grados, solar_elevation_grados)

info_celdas <- as.data.frame(patron_raster, xy = TRUE)
names(info_celdas) <- c("x_centro", "y_centro", "id_celda_patron")
resolucion_patron <- res(patron_raster)

celdas_ajustadas <- df_resultados_final %>% 
  filter(estado_celda == "ajustada") %>%
  select(id_celda_patron, modelo_seleccionado)

# Se procesan solo las celdas que fueron ajustadas exitosamente
cosenos_calculados <- future_map_dbl(celdas_ajustadas$id_celda_patron, function(id_actual) {
  modelo_actual <- celdas_ajustadas$modelo_seleccionado[celdas_ajustadas$id_celda_patron == id_actual]
  puntos_actuales <- datos_extraidos_df[datos_extraidos_df$id_celda_patron == id_actual, ]
  centro_actual <- info_celdas[info_celdas$id_celda_patron == id_actual, c("x_centro", "y_centro")]
  
  calcular_cos_angulo_promedio_celda(puntos_actuales, modelo_actual, vector_solar, as.numeric(centro_actual), resolucion_patron, n_muestras_intracelda)
}, .progress = TRUE, .options = furrr_options(seed = TRUE, globals = c("datos_extraidos_df", "celdas_ajustadas", "info_celdas", "resolucion_patron", "calcular_cos_angulo_promedio_celda", "ajustar_superficie_modelo", "predecir_z_superficie", "vector_solar", "n_muestras_intracelda")))

celdas_ajustadas$cos_autoadapt <- cosenos_calculados

df_resultados_completo <- df_resultados_final %>%
  left_join(celdas_ajustadas, by = c("id_celda_patron", "modelo_seleccionado"))


# --- 7.1.bis: GENERACIÓN DE MAPAS DE ÁNGULOS CON MODELOS LOCALES FIJOS ---
log_info("--- Sección 7.1.bis: Generando mapas con modelos locales fijos para validación (en paralelo) ---")

plantilla_mapa <- rast(patron_raster)
values(plantilla_mapa) <- NA

# Obtenemos la lista de IDs de celdas que tienen suficientes puntos para ser procesadas
ids_celdas_validas <- unique(df_resultados_final$id_celda_patron)

for (modelo_fijo in c("lineal", "cuadratico", "cubico")) {
  log_info(paste("Procesando con modelo local fijo:", modelo_fijo))
  tryCatch({
    
    cosenos_modelo_fijo <- future_map_dbl(ids_celdas_validas, function(id_actual) {
      puntos_actuales <- datos_extraidos_df[datos_extraidos_df$id_celda_patron == id_actual, ]
      centro_actual <- info_celdas[info_celdas$id_celda_patron == id_actual, c("x_centro", "y_centro")]
      
      calcular_cos_angulo_promedio_celda(puntos_actuales, modelo_fijo, vector_solar, as.numeric(centro_actual), resolucion_patron, n_muestras_intracelda)
    }, .progress = TRUE, .options = furrr_options(seed = TRUE, globals = c("datos_extraidos_df", "info_celdas", "resolucion_patron", "calcular_cos_angulo_promedio_celda", "ajustar_superficie_modelo", "predecir_z_superficie", "vector_solar", "n_muestras_intracelda", "modelo_fijo")))
    
    df_cosenos_modelo_fijo <- tibble::tibble(
      id_celda_patron = ids_celdas_validas,
      cos_angulo = cosenos_modelo_fijo
    )
    
    fichero_salida_actual <- switch(modelo_fijo,
                                    "lineal"     = fichero_mapa_angulos_global_lineal_tif,
                                    "cuadratico" = fichero_mapa_angulos_global_cuadratico_tif,
                                    "cubico"     = fichero_mapa_angulos_global_cubico_tif
    )
    
    mapa_angulos_fijo <- plantilla_mapa
    celdas_validas <- df_cosenos_modelo_fijo %>% filter(!is.na(cos_angulo))
    if(nrow(celdas_validas) > 0) {
      valores_escalados <- round(celdas_validas$cos_angulo * 100 + 100)
      mapa_angulos_fijo[celdas_validas$id_celda_patron] <- valores_escalados
    }
    
    writeRaster(mapa_angulos_fijo, fichero_salida_actual, overwrite=TRUE, datatype="INT1U")
    log_info(paste("Mapa de ángulos (Local Fijo -", modelo_fijo, ") guardado en:", fichero_salida_actual))
  }, error = function(e) {
    log_error(paste("Fallo al procesar el modelo local fijo", modelo_fijo, ":", e$message))
  })
}

# --- Volvemos al modo secuencial al finalizar todo el procesamiento pesado ---
plan(sequential) 
log_info("Procesamiento en paralelo de todas las celdas completado.")

# --- 7.2 Creación de Mapas de Salida (Autoadaptativo) ---
log_info("Creando mapas de salida del método autoadaptativo...")

mapa_angulos_autoadapt <- plantilla_mapa
celdas_validas_autoadapt <- df_resultados_completo %>% filter(!is.na(cos_autoadapt))
if(nrow(celdas_validas_autoadapt) > 0) {
  mapa_angulos_autoadapt[celdas_validas_autoadapt$id_celda_patron] <- round(celdas_validas_autoadapt$cos_autoadapt * 100 + 100)
}
writeRaster(mapa_angulos_autoadapt, fichero_mapa_angulos_autoadapt_tif, overwrite=TRUE, datatype="INT1U")
log_info("Mapa de ángulos (Autoadaptativo) guardado.")

mapa_modelos_num <- plantilla_mapa
modelo_codigos <- c("lineal"=1L, "cuadratico"=2L, "cubico"=3L)
df_codigos <- df_resultados_completo %>% 
  mutate(codigo = modelo_codigos[modelo_seleccionado]) %>%
  filter(!is.na(codigo))
if(nrow(df_codigos) > 0) {
  mapa_modelos_num[df_codigos$id_celda_patron] <- df_codigos$codigo
}
writeRaster(mapa_modelos_num, fichero_mapa_modelos_tif, overwrite=TRUE, datatype="INT2S")
log_info("Mapa de modelos seleccionados guardado.")

for (modelo_actual in c("lineal", "cuadratico", "cubico")) {
  fichero_salida <- file.path(ruta_salida_base, paste0("mapa_mascara_", modelo_actual, ".tif"))
  mapa_mascara_actual <- plantilla_mapa
  ids_celdas <- df_resultados_completo %>%
    filter(modelo_seleccionado == modelo_actual, !is.na(modelo_seleccionado)) %>%
    pull(id_celda_patron)
  if (length(ids_celdas) > 0) {
    mapa_mascara_actual[ids_celdas] <- 1
  }
  writeRaster(mapa_mascara_actual, fichero_salida, overwrite = TRUE, datatype = "INT1U")
  log_info(paste("Mapa de máscara para modelo '", modelo_actual, "' guardado."))
}

# --- 7.3 Guardado de la tabla de resultados ---
log_info("Guardando tabla de resultados final en CSV...")
write.csv2(df_resultados_completo, file = fichero_resultados_csv, row.names = FALSE, na = "")
log_info(paste("Tabla de resultados completa guardada en:", fichero_resultados_csv))

# --- 7.4 Generación del Informe en R Markdown ---
log_info("Generando informe final en R Markdown...")
informe_env <- new.env(parent = globalenv())
informe_env$df_resultados_final <- df_resultados_completo
informe_env$mapa_modelos_num <- mapa_modelos_num
informe_env$mapa_angulos_cos_autoadapt <- mapa_angulos_autoadapt

rmd_content_integral <- '
---
title: "Informe de Resultados TFM Fase 2 (Patrón 10m): Método Autoadaptativo Mejorado"
author: "Miguel C. (Generado con R)"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: cerulean
    code_folding: hide
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, fig.align = "center", dpi=150, fig.width=9, fig.height=6)
library(dplyr); library(knitr); library(ggplot2); library(DT); library(terra); library(viridisLite)
```

## Resumen General del Procesamiento
```{r resumen_vars, include=FALSE}
num_ajustadas <- sum(df_resultados_final$estado_celda == "ajustada", na.rm = TRUE)
num_descartadas <- sum(df_resultados_final$estado_celda == "descartada", na.rm = TRUE)
```
* **Número total de celdas patrón analizadas:** `r nrow(df_resultados_final)`
* **Número de celdas ajustadas exitosamente:** `r num_ajustadas`
* **Número de celdas descartadas:** `r num_descartadas`

## Detalles por Celda Ajustada
```{r tabla_ajustadas, echo=FALSE}
# Ahora la tabla puede mostrar ambas métricas de error
df_resultados_final %>%
  filter(estado_celda == "ajustada") %>%
  select(id_celda_patron, modelo_seleccionado, rmse_mejor_modelo, mae_mejor_modelo) %>%
  mutate(
    rmse_mejor_modelo = round(rmse_mejor_modelo, 3),
    mae_mejor_modelo = round(mae_mejor_modelo, 3)
  ) %>%
  DT::datatable(caption = "Tabla con los detalles de las celdas que fueron ajustadas exitosamente.", options = list(pageLength = 10, scrollX=TRUE, autoWidth = TRUE))
```

## Mapas de Salida Principales
### Mapa de Modelos Seleccionados (Autoadaptativo)
```{r mapa_modelos_plot, echo=FALSE, fig.cap="Mapa de Modelos Seleccionados (Autoadaptativo)", fig.width=8, fig.height=7}
try({
  mapa_modelos <- as.factor(mapa_modelos_num)
  levels(mapa_modelos) <- data.frame(id = c(1,2,3), label = c("Lineal","Cuadrático","Cúbico"))
  colores_modelos <- c("Lineal" = "#440154FF", "Cuadrático" = "#21908CFF", "Cúbico" = "#FDE725FF")

  plot(
    mapa_modelos, type = "classes", col = colores_modelos,
    main = "Mapa de Modelos Seleccionados", mar = c(3.1, 3.1, 2.1, 11.1),
    plg = list(cex = 0.9)
  )
}, silent=FALSE)
```

### Mapa de Cosenos de Ángulos de Incidencia (Método Autoadaptativo)
```{r mapa_angulos_plot, echo=FALSE, fig.cap="Valores Escalados de Ángulo de Incidencia (Autoadaptativo)", fig.width=8, fig.height=7}
try({
  plot(
      mapa_angulos_cos_autoadapt, main="Valores Escalados de Ángulo de Incidencia", 
      col=viridis(100), mar=c(3.1, 3.1, 2.1, 4.1),
      plg=list(title="Valor Escalado\n(0-200)", cex=0.8)
  )
}, silent = FALSE)
```
'

tryCatch({
  writeLines(rmd_content_integral, fichero_informe_rmd)
  log_info(paste("Fichero R Markdown generado en:", fichero_informe_rmd))
  
  rmarkdown::render(
    fichero_informe_rmd,
    output_file = basename(fichero_informe_html),
    output_dir = dirname(fichero_informe_html),
    envir = informe_env,
    quiet = FALSE
  )
  log_info(paste("Informe HTML generado en:", fichero_informe_html))
  
}, error = function(e_report) {
  log_error(paste("Error al generar el informe R Markdown:", e_report$message))
  log_error(paste("Detalles del error de R Markdown:", format(e_report)))
  cat("ADVERTENCIA: Asegúrate de tener 'rmarkdown' y 'pandoc' instalados y configurados.\n")
})

log_info("--- FIN DEL PROCESO FINAL Y COMPLETO ---")
