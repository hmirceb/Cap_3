#############################################################################################
#### Script para descargar datos de forma automática del Atlas Climático de Aragon (ACA) ####
#############################################################################################

# ACA_year: Disponibilidad de datos desde 1950 hasta 2020
# ACA_coord: Disponibilidad de datos desde 1950 hasta 2020 
# ACA_promedios: Disponibilidad de datos desde 1990 hasta 2020.

# Rásters de 1km2 en UTM con EPSG:25830

# Función para descargar datos de precipitación, tmax o tmin por año 
# en formato raster 
ACA_year <- function(inicio = NULL, 
         final = NULL, 
         x = c('pcp', 'tmin', 'tmax'),
         dest = NULL){
  x <- match.arg(x, choices = c('pcp', 'tmin', 'tmax'))
  message('El destino por defecto de las descargas es su carpeta de trabajo actual')
  message('Indique una de las siguientes variables para descargar:
         Precipitación = pcp, 
         Temperatura máxima = tmax,
         Temperatura mínima = tmin')
  
  if(is.null(inicio) & is.null(final)){
    stop('Es necesario indicar una fechas')
  }
  
  # Si falta uno de los años, usamos el otro como fecha de inicio y final
  if(is.null(inicio) | is.null(final)){
  inicio <- final <-  c(inicio, final)[which(!is.null(c(inicio, final)))]
  }
  
  for (year in inicio:final) {
    # Creamos la url
    url <- paste0('https://idearagon.aragon.es/datosdescarga/descarga.php?file=CartoTema/Medioambiente/Clima/2020/',
    x, '/', x, '_', as.character(year), '.zip')
    # Descargamos el archivo
    if(is.null(dest)){
      download.file(url = url, 
                    destfile = basename(url), 
                    method="curl", extra="-k")
    }else{
      download.file(url = url, 
                    destfile = paste0(dest,
                                     basename(url)), 
                    method="curl", extra="-k")
    }
  }
}

# Funcion para descargar datos de precipitación, tmin y tmax usando
# coordenadas UTM en (EPSG:25830) dentro de Aragón
ACA_coord <- function(lon = NULL, 
                            lat = NULL,
                            dest = NULL,
                      quiet = T){
  if(isFALSE(quiet)){message('El destino por defecto de las descargas es su carpeta de trabajo actual')}
  
  if(is.null(lon) | is.null(lat)){
    stop('Es necesario indicar latitud y longitud (EPSG:25830)')
  }
  if(length(lon) != length(lat)){
    stop('Cada latitud tiene que tener una longitud y viceversa')
  }
  if(!is.numeric(lon) | !is.numeric(lat)){
    stop('Latitud y longitud tienen que ser números')
  }
  # IDE Aragón utiliza cuadrículas de 1km por lo que hay que redondear
  # las coordenadas cambiando los tres últimos digitos por 500
  if(endsWith(as.character(lon), '500')){lon <- as.integer(lon)}else{ lon <- as.integer(as.numeric(
    paste0(substr(as.character(round(lon, 0)), 
                 start = 1, 
                 stop = nchar(as.character(round(lon, 0)))-3), 
          '500')), 0)}
  if(endsWith(as.character(lat), '500')){lat <- as.integer(lat)}else{lat <- as.integer(as.numeric(
    paste0(substr(as.character(round(lat, 0)), 
                 start = 1, 
                 stop = nchar(as.character(round(lat, 0)))-3), 
          '500')), 0)}
  
  lonlats <- matrix(c(lon, lat), ncol = 2, byrow = F)
  for (pixel in seq_len(dim(lonlats)[1])) {
    url <- paste0('https://idearagon.aragon.es/datosdescarga/descarga.php?file=CartoTema/Medioambiente/Clima/pixel_hco/',
                 lonlats[pixel, 1], 
                 '_', 
                 lonlats[pixel, 2],
                 '.zip')
    
    if(is.null(dest)){
      download.file(url = url, 
                    destfile = basename(url), 
                    method="curl", extra="-k",
                    quiet = quiet)
    }else{
      download.file(url = url, 
                    destfile = paste0(dest,
                                      basename(url)), 
                    method="curl", extra="-k",
                    quiet = quiet)
    }
  }
}

# Función para descargar los valores promedio de precipitación, tmin, tmax,
# y tmed para los periodos mensual, estacional y anual.
# Mensual y estacional devuelven un zip con una capa raster correspondiente
# al mes o estacion. Anual devuelve una sola capa raster. Todas ellas
# con UTM EPSG:25830
ACA_promedios <- function(periodo = c('mensual', 'estacional', 'anual'),
                          x = c('pcp', 'tmin', 'tmed', 'tmax'),
                            dest = NULL){
  x <- tolower(x)
  x <- match.arg(x, choices = c('pcp', 'tmin', 'tmed', 'tmax'))
  periodo <- tolower(periodo)
  periodo <- match.arg(periodo, choices = c('mensual', 'estacional', 'anual'))
  
  if(periodo == 'mensual'){periodo <- 'mensuales'}
  if(periodo == 'estacional'){periodo <- 'estacionales'}
  
  message('El destino por defecto de las descargas es su carpeta de trabajo actual')
  message('Indique una combinación de las siguientes variables y períodos para descargar:
         Variables:
         Precipitación = pcp, 
         Temperatura máxima = tmax,
         Temperatura mínima = tmin,
          Temperatura media = tmed
          Períodos:
          Mensual,
          Estacional,
          Anual')

    # Creamos la url
  if(periodo == 'anual'){
    url <- paste0('https://idearagon.aragon.es/datosdescarga/descarga.php?file=CartoTema/Medioambiente/Clima/2020/promedios/',
                  x, '/', x, '_media_', periodo, '.tif')
  }else{
    url <- paste0('https://idearagon.aragon.es/datosdescarga/descarga.php?file=CartoTema/Medioambiente/Clima/2020/promedios/',
                  x, '/', x, '_medias_', periodo, '.zip')
  }
    # Descargamos el archivo
    if(is.null(dest)){
      download.file(url = url, 
                    destfile = basename(url), 
                    method="curl", extra="-k")
    }else{
      download.file(url = url, 
                    destfile = paste0(dest,
                                      basename(url)), 
                    method="curl", extra="-k")
  }
}
