##################################
#
# TRABAJO R
# Realizado por: Andrea López Guirado
#
#################################


################
#PREGUNTA 1
################
#Lo primero que hay que hacer es subir los datos que están en el campus virtual
datos<-read.table("datos-trabajoR.txt",header=TRUE, sep="\t")
datos

#Ahora vemos las dimensiones de los datos:

head(datos)
#cuando lo ejecutamos podemos ver las 6 primeras filas de nuestros datos

summary(datos)
#este comando nos proporciona un resumen de los tres quartiles, la media, la mediana, maximos y mínimos de los datos del tratamiento, variable1 y variable2

dim(datos)
#esta función es empleada para ver las dimensiones de los datos, cuantas columnas y filas tiene: en este caso nos dice que tiene 50 líneas y 3 columnas

str(datos)
#esta función nos ayuda a ver los datos de una forma más compacta

#viendo los resultados podemos comprobar que hay 2 variables y 5 tratamientos


################
#PREGUNTA 2
################
#Procedemos a hacer un boxplot, mostrando la variable 1 en rojo y la variable 2 en azul
#usamos este gráfico de bigotes y cajas para representar nuestros datos tanto de la variable 1 como de la variable 2 enfrentados con el tratamiento (es el factor de las variables)
#el comando "col" sirve para poner color a nuestras variables y "main" para ponerle un título 
boxplot(Variable1~Tratamiento, data=datos, col=c("red"), main="Variable 1/Tratamiento - Boxplot")
boxplot(Variable2~Tratamiento, data=datos, col=c("blue"), main="Variable 2/Tratamiento - Boxplot")


################
#PREGUNTA 3
################
#Ahora hacemos un gráfico de dispersión con las dos variables. Cada tratamiento va de un color distinto
#usamos plot() para crear nuestro gráfico, "col = datos$Tratamiento" para que cada un de los 5 tratamiento estuviera en un color diferente, "xlab" para poner la variable1 en el eje X y el ylab para poner la variable2 en el eje Y
#por otro lado, usamos "main" para ponerle título y "pch=1" para representar los puntos 
plot(datos$Variable1, datos$Variable2, col = datos$Tratamiento, pch = 1, xlab = "Variable1", ylab = "Variable2", main = "Gráfico de Dispersión")


################
#PREGUNTA 4
################
#Le ponemos la leyenda al gráfico de dispersión en el margen inferior derecho utilizando el comando legend()
#como queremos que nuestra leyenda este abajo a la derecha, escribimos: "bottomright"
legend("bottomright", legend = levels(factor(datos$Tratamiento)), col = unique(datos$Tratamiento), pch = 1, title = "Tratamiento")


################
#PREGUNTA 5
################
#Hacemos un histograma manteniendo el color de cada una de las variables (Variable1 roja y Variable 2 azul)
hist(datos$Variable1, col = "red", main = "Histograma Variable1", xlab = "Variable1", ylab = "Frecuencia")
hist(datos$Variable2, col = "blue", main = "Histograma Variable2", xlab = "Variable2", ylab = "Frecuencia")


################
#PREGUNTA 6
################
#Creamos un factor en la columna tratamiento y lo guardamos en una variable
datos$tratamiento_factor <- factor(datos$Tratamiento)

#comprobamos que lo hemos hecho bien
head(datos$tratamiento_factor)
datos$tratamiento_factor

################
#PREGUNTA 7
################
#Ahora calculamos la media y la desviación estándar para cada tratamiento gracias a la función aggregate()
#calculamos la media y desviación estándar para Variable1 
media_variable1 <- aggregate(Variable1 ~ datos$tratamiento_factor, datos, mean)
desviacion_variable1 <- aggregate(Variable1 ~ datos$tratamiento_factor, datos, sd)

#vemos los resultados para la Variable1 con print()
print("Media de Variable1 por Tratamiento:")
print(media_variable1)
print("Desviación Estándar de Variable1 por Tratamiento:")
print(desviacion_variable1)

#hacemos lo mismo para la media y desviación estándar de la Variable2 
media_variable2 <- aggregate(Variable2 ~ datos$tratamiento_factor, datos, mean)
desviacion_variable2 <- aggregate(Variable2 ~ datos$tratamiento_factor, datos, sd)

#para poder ver los resultados para Variable2 hacemos lo mismo
print("Media de Variable2 por Tratamiento:")
print(media_variable2)
print("Desviación Estándar de Variable2 por Tratamiento:")
print(desviacion_variable2)


################
#PREGUNTA 8
################
#Para averiguar cuántos elementos tiene cada tratamiento, ejecutamos este comando con el factor table() 
elementos_por_tratamiento <- table(datos$tratamiento_factor)
elementos_por_tratamiento
#podemos ver que cada tratamiento (1-5) tiene 10 elementos 


################
#PREGUNTA 9
################
#Ahora extraemos los datos para el tratamiento 1 y el tratamiento 4 y los guardamos cada uno en una variable diferente

#para los datos para el tratamiento 1
datos_tratamiento_1 <- datos[datos$tratamiento_factor== "1", ]
datos_tratamiento_1

#para los datos para el tratamiento 4 y los guardo en otra variable
datos_tratamiento_4 <- datos[datos$tratamiento_factor== "4", ]
datos_tratamiento_4


################
#PREGUNTA 10
################
#Nuestra hipótesis nula es que las medias de tratamiento 1 y tratamiento 4 para la Variable 1 son iguales


#¿Puedes comprobarlo? Para ello, necesitarás comprobar primero si los datos se distribuyen de forma normal
#usamos el shapiro.text() para ver si mis datos siguen una distribución normal

# Prueba de normalidad para tratamiento 1
shapiro.test(datos$Variable1[datos$tratamiento_factor == 1])
#tiene un p value de 0.06434

# Prueba de normalidad para tratamiento 4
shapiro.test(datos$Variable1[datos$tratamiento_factor == 4])
#tiene un p value de 0.1564

#como obervamos nuestros p values no son menores de 0.05, indicando que no podemos rechazar nuestra hipótesis nula (H0=las medias de tratamiento 1 y tratamiento 4 para la Variable 1 son iguales)
#gracias a esto sabemos que nuestras muestras no siguen una distribución normal


#En función del resultado de la prueba de normalidad, ¿qué test usarías? 
#usaríamos la prueba de t-student si nuestras variables son iguales. En este caso no lo son, por ello usamos el test de Wilcoxon-Mann-Whitney
#como no siguen una distribución normal, usamos el test de Wilcoxon-Mann-Whitney para tratamiento 1 y tratamiento 4
wilcox.test(datos_tratamiento_1$Variable1, datos_tratamiento_4$Variable1,exact=FALSE)
#este test lo usamos para ver si el tratamiento 1 y 4 en base a la variable 1 son iguales o no
#una realizado llegamos a la conclusión de que no son iguales debido a que el pvalue obtenido es menor de 0,05

#teniendo todo esto en cuenta, podemos rechazar la hipótesis nula y decimos que los tratamientos no son iguales estadísticamente


#En general, asumimos que las muestras son independientes, pero ¿son sus varianzas iguales? 
#uso el comando var.test() para comparar las varianzas según el pvalue
var.test(datos_tratamiento_1$Variable1, datos_tratamiento_4$Variable1)
#nos da como resultado p-value = 4.595e-07. Si el valor de pvalue es mayor de 0.05 las varianzas son iguales, sin embargo si es menor de 0.05 no lo son

#teniendo todo esto en cuenta, podemos decir que las varianzas son estadísticamente diferentes.


