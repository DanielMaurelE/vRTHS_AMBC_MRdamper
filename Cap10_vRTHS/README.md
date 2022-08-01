# Capitulo 10: Implementacion de la RTHS virtual

## Descripcion

Se implementa compensacion dinamica adaptiva basada en modelos en ensayos de RTHS virtuales, con amortiguadores magneto-reologicos como especimen de ensayo.
Se consideran 3 estructuras de referencia (Caso I, II y III), con los terremotos de El Centro, Kobe y del Maule.


## Archivos

### vRTHS_Simulaciones.m (Matlab script)

Script que permite obtener los resultados de las diferentes vRTHS con las diferentes estructuras y terremotos de entrada.

### vRTHS.m (Matlab function)

Contiene la informacion para correr las simulaciones hibridas virtuales.

### indicadores.m (Matlab function)

Contiene una funcion con la que se calculan los diferentes indicadores de la RTHS. Ademas, tiene la opcion de graficar los diferentes indicadores en tiempo-historia.

### vRTHS_Referencia.slx (Simulink model)

Contiene el modelo para las simulaciones del sistema de referencia con amortiguador MR y el modelo sin control.

### vRTHS_Compensacion.slx (Simulink model)

Contiene los modelos para las simulaciones del sistema sin compensacion y con compensacion.


## Carpetas

### Modelos

Carpeta que contiene los modelos tanto del amortiguador MR como del actuador.

### Resultados

Carpeta donde se guardan los resultados.

### Terremotos

Carpeta que contiene los registros sismicos.

 
