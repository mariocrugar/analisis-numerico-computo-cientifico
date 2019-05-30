##########################################################
# -------------------------------------------------------#
# Trabajo_MNO                                            #
# Equipo 6                                               #
#                                                        #
# Portafolio de Media-Varianza                           #
# Utilizaci?n del m?todo de programaci?n cuadr?tica      #
#--------------------------------------------------------#
##########################################################

# Portafolio sin posiciones cortas
# Funci?n Portafolio_Optimo

#Los inputs son el rendimiento de los activos y el rendimiento establecido por el inversionista
Portfolio_Optimo=function(Rendimientos_Activos, Rend_Obj)
{ # Obtenemos el n?mero de activos
  n= ncol(Rendimientos_Activos)
  # Determinamos la matriz de covarianzas de los rendimientos
  Dmat = cov(Rendimientos_Activos)
  # Hacemos que el vector d sea de puros ceros para que -d^Tb=0
  dvec = rep(0, times=n)
  # Definimos la matriz de restricciones
  Amat = t(rbind(Rendimiento=colMeans(Rendimientos_Activos),Presupuesto=rep(1, n),LongOnly=diag(n)))
  # Definimos lo que debe cumplirse: Que el rendimiento sea igual al rendimiento objetivo y que los pesos sumen 1
  bvec = c(Rendimiento=Rend_Obj, Presupuesto=1, LongOnly=rep(0, times=n))
  # Le decimos a la funci?n que s?lo las dos primeras restricci?nes son de igualdad
  meq = 2
  # 2 Optimize Weights:
  # La funci?n <solve.QP> imiplementa el m?todo dual de Golfarb & Idnani (1982) para resolver problemas
  # de programaci?n cuadr?tica (quadratic programming) de la forma 
  # min (-d' * b + 1/2 * b' D * b) con la restricci?n A' * b >= b_0
  # En la que Dmat es la funci?n objetivo
  # Amat son las restricciones
  # bvec es el equivalente a b_0
  # Pesos optimizados:
  Portfolio = solve.QP(Dmat, dvec, Amat, bvec, meq)
  Pesos = Portfolio$solution
  names(Pesos) = colnames(Rendimientos_Activos)
  # Rendimiento:
  list(Pesos = 100*Pesos, Riesgo = Portfolio$value, Rendimiento = Rend_Obj)
}

# Cargar las librer?as necesarias
library(quadprog)
args(solve.QP)
library(fBasics)
# Obtener los rendimeintos
Rendimientos_Activos=100 * LPP2005REC[, 1:6]
Rendimientos_Activos
# Verificar los rendimientos
names(Rendimientos_Activos)
# Proponemos como rendimiento objetivo la media de todos los activos del portafolio:
Rend_Obj=mean(colMeans(Rendimientos_Activos))
Rend_Obj
# Optimizamos:
Portafolio=Portfolio_Optimo(Rendimientos_Activos, Rend_Obj)
Portafolio
# Obtenemos los pesos del portafolio
Pesos = Portafolio$Pesos
Pesos
# Revisamos que la suma de los pesos da 100%
sum(Pesos)
# Revisamos que el rendimiento ponderado coincide con el rendimiento objetivo
Rendimiento_Ponderado = Pesos %*% colMeans(Rendimientos_Activos)
Rendimiento_Ponderado

#---------------------------------------------------------------------
# Problema 2
# Se define la funci?n del Portafolio de Tangencia.
# Dicha funci?n incluye una funci?n que calcula el ?ndice de Sharp
#---------------------------------------------------------------------

Portafolio_Tangencia=function(Rendimientos_Activos, Tasa_Libre_Riesgo=0)
{
  # Funci?n del ?ndice de Sharpe:
  Indice_Sharpe=function(x, Rendimientos_Activos, Tasa_Libre_Riesgo)
  {
    # Primero calculamos el ?ndice de Sharpe mediante el rendimiento objetivo, la tasa libre de riesgo y el rendimiento de los activos:
    # Se define el rendimiento objetivo
    Rend_Obj = x
    # Se obtienen los pesos del portafolio de m?nima varianza
    Portafolio=Portfolio_Optimo(Rendimientos_Activos, Rend_Obj)
    Pesos = Portafolio$Pesos
    # Se define el riesgo objetivo
    Riesgo_Obj = sqrt( Pesos %*% cov(Rendimientos_Activos) %*% Pesos )[[1]]
    # Se obtiene el valor del ?ndice (Cociente) de Sharpe y sus atributos
    Cociente = (Rend_Obj - Tasa_Libre_Riesgo)/Riesgo_Obj
    attr(Cociente, "Pesos") <- Pesos
    attr(Cociente, "Riesgo_Obj") <- Riesgo_Obj
    Cociente
  }
  # Segundo: Optimizamos el portafolio de tangencia
  # N?mero de activos
  n = ncol(Rendimientos_Activos)
  # rendimiento de cada uno de los activos
  miu = colMeans(Rendimientos_Activos)
  # Matriz de Varianzas y Covarianzas de los Activos Financieros
  Cov = cov(Rendimientos_Activos)
  # Optimizaci?n del portafolio de tangencia
  Port_Tangencia=optimize(f=Indice_Sharpe, interval=range(miu), maximum=TRUE,
    Rendimientos_Activos=Rendimientos_Activos, Tasa_Libre_Riesgo=Tasa_Libre_Riesgo)
  # Tercero: Caracteristicas del portafolio de tangencia:
  # Rendimiento del portafolio de tangencia
  Rendimiento_PT = Port_Tangencia$maximum
  # Riesgo del portafolio de tangencia
  Riesgo_PT = attr(Port_Tangencia$objective, "Riesgo_Obj")
  # Pesos del portafolio de tangencia
  Pesos_PT = attr(Port_Tangencia$objective, "Pesos")
  # ?ndice de Sharpe del portafolio de tangencia
  Indice_Sharpe = Indice_Sharpe(Rendimiento_PT, Rendimientos_Activos, Tasa_Libre_Riesgo)[[1]]
  # Valor de los endimientos
  list(Indice_Sharpe=Indice_Sharpe,Riesgo_PT=Riesgo_PT, Rendimiento_PT=Rendimiento_PT, Pesos_PT=Pesos_PT)
}

#-----------------------------------------------------------------
# Ejemplo:
#-----------------------------------------------------------------
# Cargar la libreria

library(fBasics)

# Rendimientos de los activos
Rendimientos_Activos=100 * LPP2005REC[, 1:6]

# Portafolio de Tangencia
Portafolio_Tangencia(Rendimientos_Activos, Tasa_Libre_Riesgo = 0)

###################################################################
#------------------------------------------------------------------
# Portafolio de m?nima varianza y la frontera eficiente

Pesos_FE = function(Rendimientos_Activos, Rend_Obj)
{
  n = ncol(Rendimientos_Activos)
  portfolio = solve.QP( Dmat = cov(Rendimientos_Activos),dvec = rep(0, times=n),Amat = t(rbind(Rendimiento=colMeans(Rendimientos_Activos),
                   Presupuesto=rep(1, n), LongOnly=diag(n))),
    bvec = c(Rendimiento=Rend_Obj, Presupuesto=1, LongOnly=rep(0, times=n)), meq=2)
  Pesos_FE = portfolio$solution
  Pesos_FE
}


Frontera_Eficiente = function(Rendimientos_Activos, N_Portafs=20)
{
  # N?mero de activos
  n = ncol(Rendimientos_Activos)
  # rendimiento de cada uno de los activos
  miu = colMeans(Rendimientos_Activos)
  # Rendimientos objetivos
  Rendimientos_Objetivo=seq(min(miu), max(miu), length=N_Portafs)
  # Pesos optimizados
  Pesos = rep(0, n)
  Pesos[which.min(miu)] = 1
  for (i in 2:(N_Portafs-1)) {
    Pesos_Nuevos=Pesos_FE(Rendimientos_Activos, Rendimientos_Objetivo[i])
    Pesos = rbind(Pesos, Pesos_Nuevos)
  }
  
  Pesos_Nuevos = rep(0, n)
  Pesos_Nuevos[which.max(miu)] = 1
  Pesos = rbind(Pesos, Pesos_Nuevos)
  Pesos = round(Pesos, 4)
  colnames(Pesos) = colnames(Rendimientos_Activos)
  rownames(Pesos) = 1:N_Portafs
  # Obtenemos los pesos:
  Pesos
}

# Aplicaci?n:
library(fBasics)
Rendimientos_Activos=100 * LPP2005REC[, 1:6]

# Pesos de la Frontera Eficiente
Pesos = Frontera_Eficiente(Rendimientos_Activos, N_Portafs = 20)

# Imprimir los pesos de la FE
print(Pesos)

# Definir las medias de los rendimientos
miu = colMeans(Rendimientos_Activos)
Rendimientos_Objetivo = seq(min(miu), max(miu), length = nrow(Pesos))
Riesgo_Objet = NULL
for (i in 1:nrow(Pesos)) {
  Nuevo_Riesgo_Obj = sqrt(Pesos[i, ] %*% cov(Rendimientos_Activos) %*% Pesos[i, ])
  Riesgo_Objet = c(Riesgo_Objet, Nuevo_Riesgo_Obj)
}

# Graficar la Frontera Eficiente
plot(Riesgo_Objet, Rendimientos_Objetivo, pch = 19, col="Blue")
title(main = "Frontera Eficiente")
