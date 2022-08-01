# Capitulo 9: Diseño de la compensacion

## Descripcion

Se diseña el compensador adaptivo basado en modelos desarrollado por Fermandois et al. (2020)
Corresponde a una modificacion de los archivos encontrados en https://github.com/FermandoisLab/RobustAMBCvRTHS/


## Archivos

### AMB_1_Modelo_inicial.m (Matlab script)

Script contiene la informacion del modelo inicial del compensador, los parametros de simulacion, terremotos, estructuras de referencia.

### AMB_2_Calibracion.m (Matlab script)

Contiene la informacion de las plantas, terremotos y subestructuras numericas de calibracion.
 
### AMB_4_Analisis_sensibilidad.m (Matlab script)

Analiza el indicador R2 con diferentes valores de las ganancias adaptivas.

### AMB_5_optimizacion.m (Matlab script)

Corre la optimizacion con enjambre de particulas para encontrar las ganancias optimas. 

### AMB_6_R2function.m (Matlab function)

Corresponde a la funcion para calcular los valores de R2 en N simulaciones.

### AMB_9_calibration_simulation.slx (Simulink model)

Modelo de Simulink para las simulaciones de calibracion.


## Carpetas

### Resultados

Carpeta donde se guardan los resultados de la calibracion
 

 
