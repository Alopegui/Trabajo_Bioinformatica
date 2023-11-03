#############################################################################
#
# PRACTICA 1
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa
#
## ENTREGA EL 8 OCTUBRE 23:59
## Adjuntar en la entrega el PDF final y el archivo con los genes
#
##############################################################################

# Instalar RCurl

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")

# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/GSE5583_data", followlocation = TRUE)
data = as.matrix(read.table (text = url, row.names = 1, header = T))

# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
dim(data) #esto mide las dimensiones
head(data)
tail(data)

# Hacemos un primer histograma para explorar los datos
hist(data)

# Transformamos los datos con un logaritmo 
# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve?
data_log=log2(data) 
#usamos los logaritmos para trasformar la gráfica y hacerla más bonita y ver mejor los datos, tiene una distribución normal
hist(data_log)

# Hacemos un boxplot con los datos transformados. ¿Qué significan los parámetros que hemos empleado?
# ¿Qué es un boxplot?
boxplot(data_log) 
#usamos este tipo de gráfico de barras de cajas y bigotes para visualizar nuestros datos de otra forma
boxplot(data_log, col=c("blue","blue","blue","orange","orange","orange"),main="GSE5583-boxplots", las=2) 
#ponemos colores (col), lo nombramos con un título (main) y ponemos los ejes en vertical (las)

# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación
# de los valores de expresión. ¿Es correcta la separación?
#agrupamos los datos por clasificaciones, en este caso con wildtype y knockout, la separación es correcta ya que vemos 3 de cada grupo
hc = hclust(as.dist(1-cor(data_log)))
plot(hc, main="GSE5583 - Hierarchical Clustering")

#######################################
# Análisis de Expresión Diferencial 
#######################################

# Primero separamos las dos condiciones. ¿Qué tipo de datos has generado?
#hemos separado los datos de los wildtype de los knockout por medio de una matriz
wt <- data[,1:3]
ko <- data[,4:6]
class(wt) #tipo de datos que tenemos
head(wt) #para ver los encabezados

# Calcula las medias de las muestras para cada condición. Usa apply
#para calcular todas las medias (en este caso) de toda la tabla de golpe
wt.mean = apply(wt, 1, mean)
wt.mean
head(wt.mean) 

ko.mean = apply(ko, 1, mean)
ko.mean
head(ko.mean)

# ¿Cuál es la media más alta?
max(wt.mean)
max(ko.mean)
#la media es más alta en knockout

# Ahora hacemos un scatter plot (gráfico de dispersión)
#un gráfico de dispersión es un tipo de gráfico de puntos; compara las medias en el eje x e y; enfrenta una lista a otra lista
plot(ko.mean ~ wt.mean)
plot(ko.mean ~ wt.mean, xlab = "WT", ylab = "KO", main = "GSE5583 - Gráfico de dispersión")

# Añadir una línea diagonal con abline
#necesiitamos una línea de regresión con pendiente (y=ax+b; y=x)
abline(0, 1, col= "red")
#ahora añadimos una línea horizontal
abline(h=2,col="blue")
abline(v=5,col="green")

# ¿Eres capaz de añadirle un grid?
#para añadir cuadrícula a la gráfica
grid()

# Calculamos la diferencia entre las medias de las condiciones
#wt y ko son nuestras variables; restamos dos matrices y obtenemos una tabla
diff.mean = wt.mean - ko.mean

# Hacemos un histograma de las diferencias de medias
hist(diff.mean,col="pink")

# Calculamos la significancia estadística con un t-test.
#hay que ver si la diferencia es significativa o no; por cada fila hacemos un t-test
#Primero crea una lista vacía para guardar los p-values
#los p-values son lo que nos interesa del t-test; hacemos una lista: 1 para cada gen, hay 12488 genes y valores
# Segundo crea una lista vacía para guardar las estadísticas del test.
# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué?
#porque manipulamos los datos y pueden no ser fiables, hay que usar los datos tal cual
# ¿Cuántas valores tiene cada muestra?
#cada condición tiene 3 valores

pvalue = NULL
tstat = NULL
for(i in 1 : nrow(data)) { #Para cada gen 
	x=wt[i,] # gene wt número 1 (va a tener tres valores)
	y=ko[i,] # gene ko número i (va a tener tres valores)
	
	# Hacemos el test (guardamos el resultado en la variable t, es un vector)
	t = t.test(x,y)

	# Añadimos el p-value a la lista
	pvalue[i]=t$p.value
	# Añadimos las estadísticas a la lista
	tstat[i] = t$statistic
}

head(pvalue)

length(pvalue)

# Ahora comprobamos que hemos hecho TODOS los cálculos

# Hacemos un histograma de los p-values.
hist(pvalue)

# ¿Qué pasa si le ponemos con una transformación de -log10?
#para modificar la gráfica y ver los valores que nos interesan, cambiando la esccala y reflejar más pa proporción de la distribución
hist(-log10(pvalue), col="purple")

# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística
#diferencia de medias en contra del log10; los pvalue significativos están en la parte de arriba
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano")

# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico?
diff.mean_cutoff = 2
pvalue_cutoff = 0.01
abline(v = diff.mean_cutoff, col = "blue", lwd = 3)
abline(h = -log10(pvalue_cutoff), col= "green", lwd =3)
#los valores significativos están por encima de la línea verde y los que se separan po la línea azul son los reprimidos y los sobreexpresados

# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold)
filter_by_diff.mean = abs(diff.mean) >= diff.mean_cutoff
dim(data[filter_by_diff.mean, ])
#11859  6

# Ahora el filtro de p-value
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(data[filter_by_pvalue, ])
#426   6

# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios?
filter_combined = filter_by_diff.mean & filter_by_pvalue
filtered = data[filter_combined,]
dim(filtered)
head(filtered)
#ya tengo el filtro combinado y los datos metidos de los genes significativos
#hay 426 genes que cumplen ambos criterios

# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano #2")
points (diff.mean[filter_combined], -log10(pvalue[filter_combined]), col = "red")

# Ahora vamos a marcar los que estarían sobreexpresados (negativo, rojo) y reprimidos (positivo, azul). ¿Por qué parece que están al revés?
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano #3")
points (diff.mean[filter_combined & diff.mean < 0],
	-log10(pvalue [filter_combined & diff.mean < 0]), col = "red")
points (diff.mean[filter_combined & diff.mean > 0],
	-log10(pvalue [filter_combined & diff.mean > 0]), col = "blue")
#porque hemos calculado la diferencia de medidas al revés, por lo que los reprimidos estarán en el lado positivo y los sobreexpresados en el negativo

# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 
# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap?
# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors
#lo hacemos con los datos significativos
heatmap(filtered)
rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv = as.dendrogram(hclust(as.dist(1-cor(filtered))))
heatmap(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,labRow=FALSE)
#vemos una diferencia notable de expresión, los de notable expresión estan en rojo

# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer
#if (!requireNamespace("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install(c("gplots","RColorBrewer"))
install.packages("gplots")		
install.packages("RColorBrewer")	

library(gplots)

# Hacemos nuestro heatmap


# Lo guardamos en un archivo PDF


# Guardamos los genes diferencialmente expresados y filtrados en un fichero
write.table (filtered, "GSE5583_DE.txt", sep = "\t",
	quote = FALSE)
