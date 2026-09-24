## Definir la función para calcular la regresión logística
calcular_prevalencia <- function(dataset, coeficientes) {
  # Verificar que los coeficientes incluyen al menos el intercepto
  if (!("Intercept" %in% names(coeficientes))) {
    stop("El vector de coeficientes debe incluir 'Intercepto' como el término independiente.")
  }
  
  # Extraer el término independiente (intercepto)
  intercepto <- coeficientes["Intercept"]
  
  # Calcular el logit inicial con el intercepto
  dataset$Logit <- intercepto 
  
  # Iterar sobre los coeficientes restantes y multiplicarlos por las variables del dataset
  for (variable in setdiff(names(coeficientes), "Intercept")) {
    if (variable %in% names(dataset)) {
      dataset$Logit <- dataset$Logit + dataset[[variable]] * coeficientes[variable]
    } else {
      warning(paste("La variable", variable, "no está en el dataset y será ignorada."))
    }
  }
  
  # Aplicar la función sigmoide para obtener la probabilidad
  dataset$Probabilidad <- 1 / (1 + exp(-dataset$Logit))
  
  return(dataset)
}
