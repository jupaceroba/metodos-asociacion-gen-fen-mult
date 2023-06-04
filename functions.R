## Reducción de dimensionalidad
dim_reduc_chisq <- function(datos, n1 = NULL, n2 = NULL){
    n1 <- ifelse(is.null(n1), 1, n1)
    n2 <- ifelse(is.null(n2), ncol(datos), n2)
    datos = datos[,n1:n2]
    n = length(n1:n2)
    matriz_ind=matrix(ncol=n, nrow = n)
    combinations = combn(1:ncol(datos), 2)    
    for(k in 1:ncol(combinations)){
        i = combinations[1, k]
        j = combinations[2, k]
        chi_value = 1-chisq.test(datos[,i], datos[,j])$p.value
        matriz_ind[i,j]=chi_value
        matriz_ind[j,i]=chi_value
    }
    
    diag(matriz_ind) =1

    colnames(matriz_ind)=colnames(datos)
    rownames(matriz_ind)=colnames(datos)
    datosNA=ifelse(matriz_ind<=1,matriz_ind,NA)

    return(matriz_ind)

}


# Metodología a aplicar
metodo<-function(X,G,Y){
  P=diag(length(X))-X%*%solve(t(X)%*%X)%*%t(X)
  raiz=NULL

  for (i in 1:ncol(G)) {
    raiz[i]=sqrt(t(G[,i])%*%P%*%G[,i])
  }

  G_gorro=matrix(ncol = ncol(G),nrow = nrow(G))
  for (i in 1:ncol(G)) {
    G_gorro[,i]=(P%*%G[,i])/raiz[i]
  }

  Beta=matrix(ncol=ncol(G), nrow = ncol(Y))

  for (i in 1:ncol(Y)) {
    Beta[i,]=t(Y[,i])%*%G_gorro
  }

  sigma2=matrix(ncol=ncol(G), nrow = ncol(Y))

  for (i in 1:ncol(G)) {
    for (j  in 1:ncol(Y)) {
      sigma2[j,i]=norm((P%*%Y[,j]-G_gorro[,i]*(Beta[j,i])))^2/(nrow(G)-3)
    }
  }

  Z=Beta/sigma2

  ##Matríz de covarianzas 
  matcov=matrix(ncol = nrow(Z), nrow = nrow(Z))
  for (i in 1:nrow(Z)) {
    for (j in 1:nrow(Z)) {
      matcov[i,j]=(t(G[,i])%*%G[,j])/sqrt((t(G[,i])%*%G[,i])*(t(G[,j])%*%G[,j]))
    }
  }
  re=list(m.cov= matcov,Z.score = Z)
  return(re)
}

draw.correlated.binary<-function(nrep,d,mean.vec,corr.mat){
  if((nrep<1)|(floor(nrep)!=nrep)){
    # Si el numero de repeticiones es menor que uno o
    #  floor devuelve el numero entero m?s grande no mayor que el numero que da
    # Por ejemplo floor(-3.5) = -4 o floor(2.5)=2
    stop("El numero de repeticiones debe ser un n?mero entreo cuyo valor
          sea al menos 1!\n")}
  if((d<2)|(floor(d)!=d)){
    stop("La dimension debe ser un n?mero entero cuyo valor sea al menos de 2!\n")}
  if((max(mean.vec)>=1)|(min(mean.vec)<=0)){
    # Vector de medias debe estar compuesto por proporciones 
    stop("Las medias deben estar entre 0 y 1!\n")}
  if(length(mean.vec)!=d){
    stop("El vector de medias esta mal especifiacdo, la dimesion es incorrecta!\n")}
  if((ncol(corr.mat)!=d)|(nrow(corr.mat)!=d)){
    stop("La matriz de correlaci?n esta mal especificada, la dimensi?n es incorrecta!\n")}
  if(sum(corr.mat!=t(corr.mat))>0){
    stop("La matriz de correlaci?n no es simetrica!\n")}
  if(sum(diag(corr.mat)!=rep(1,d))>0){
    stop("No todos los elementos de la diagonal de la matriz de correlaci?n son 1!\n")}
  if((max(corr.mat)>1)|(min(corr.mat)<0)){
    stop("Las correlaciones deben estar entre 0 y 1!\n")}
  
  alpha<-matrix(0,d,d) ; cor.limit<-matrix(0,d,d)
  for (i in 1:d){  # Calcula el valor maximo que puede tomar las correlaciones entre las variables
    for (j in 1:d){ # Formula (2.4)
      cor.limit[i,j]<-min(sqrt(mean.vec[j]*(1-mean.vec[i]))/(mean.vec[i]*(1-mean.vec[j])),
                          sqrt((mean.vec[i]*(1-mean.vec[j]))/(mean.vec[j]*(1-mean.vec[i]))))}} 
  
  if(sum(cor.limit>=corr.mat)<d^2){
    stop("Las correlaciones est?n m?s all? de sus l?mites superiores impuestos por las expectativas")}
  
  for (i in 1:d){for (j in 1:d){ #Formula (2.3)
    alpha[i,j]<-log(1+corr.mat[i,j]*sqrt((1-mean.vec[i])* (1-mean.vec[j])/(mean.vec[i]*mean.vec[j])))}}
  #Valores de alphaij para las variables X poisson 
  
  beta<-matrix(0,d,d*d)
  
  summ<-1 ; counter<-0 
  while (summ>0){
    counter<-counter+1 # Contador
    minloc<-min.loc.finder(alpha,d) # Localizacion del minimo global, indices para Sl
    w<-matrix(1,d,d) 
    my.min<-apply(matrix(alpha[,-minloc],d,d-length(unique(minloc))),2,min) 
    
    # Unique quita los valores duplicados 
    # Calcula el valor minimo para cada columna sin tener en cuenta las columnas correspondiente a las
    # coordenadas de localizacion del minimo calculado anteriormente 
    
    
    if (length(my.min)==1){
      w[,-minloc][my.min==0]<-0
      w[-minloc,][my.min==0]<-0}
    
    if (length(my.min)>1){
      w[,-minloc][,my.min==0]<-0
      w[-minloc,][my.min==0,]<-0
      w[alpha==0]<-0}
    
    for (i in 1:d){
      beta[i,counter]<- alpha[minloc[1],minloc[2]]*1*((minloc[1]==i)|(minloc[2]==i)|(sum(w[,i])==d))}
    
    alpha<-alpha-alpha[minloc[1],minloc[2]]*w
    summ<-sum(alpha)} 

  #hacer el while hasta que la diferencia entre los alpha sea cero
  
  tbeta<-t(beta)
  w<-(tbeta!=0) # ve si los valores de los betas son diferentes de cero
  x<-matrix(0,nrep,d); y<-matrix(0,nrep,d) # Matriz que almacenara los valores Poisson correlacionados
  pois<-numeric(nrow(tbeta)) # Vector numerico del tama?o de las filas de beta traspuesta
  sump<-numeric(d) #Vector del tama?o de las variables
  
  for (k in 1:nrep){
    for (j in 1:nrow(tbeta)){ 
      pois[j]<-rpois(1,max(tbeta[j,]))} #Genero un valor Poisson con el parametro beta seleccionado
    for (i in 1:d){
      sump[i]<-sum(pois*w[,i])} # Almacena las suma solo para los Poisson que son verdaderos
    x[k,]<-sump} # Calculo de las variables Y convolucionadas
  y[x==0]<-1 # Donde sea verdadero ponga un 1
  y[x!=0]<-0 # Donde sea vardadero ponga un 0
  y} # Para lo que en las copias es Z

min.loc.finder<-function(my.mat,d){ # Funci?n que localiza el minimo de los alpha
  w<-is.matrix(my.mat) # Verfica si la matriz corresponda a la clase que se necesita
  
  if (w==F){stop("Esto no es una matriz!\n")} # Si la matriz no es de la clase "matrix" arroja error
  
  if (nrow(my.mat)!=ncol(my.mat)){  # Verifica si la matriz es cuadrada
    stop("Esta no es una matriz cuadrada!\n")}
  
  n<-nrow(my.mat) # Numero de filas o observaciones
  my.vec<-as.vector(t(my.mat)) # Vuelve vector la matrix agrada
  my.vec[my.vec==0]<-999
  my.index<-min((1:length(my.vec))[my.vec==min(my.vec)]) #Celda donde esta el minimo valor
  row.index<-floor((my.index-1)/n)+1 # Numero de la fila donde se encuentra el valor minimo
  col.index<-my.index-d*floor((my.index-1)/n) #Numero de la columna donde sen encuentra el valor minimo
  c(row.index,col.index)} # Coordenada dentro de la tabla, donde esta el valor minimo

create_var<-function(n){
  a=diag(x=1, n, n)    
  for (i in 1:ncol(a)) {
    for (j in 1:ncol(a)) {
    if((i==1)&(j==1)){
        a[i,]=a[,i]=runif(ncol(a), 0.7, 0.9)
    }
      
    }}
  for (i in 2:ncol(a)) {
      if(i!=1){
        a[i,2:ncol(a)]=a[2:ncol(a),i]=runif(ncol(a)-1, 0, 0.1)
      }
    }
  for (i in 1:ncol(a)) {
    for (j in 1:ncol(a)) {
      if(i==j){
        a[i,j]=1
      }else{
        a
      }
    }
  }
  return(a)
  }

random_norm<-function(x){
    y=NULL
    fenopapa=read_excel("data/Accessiones_Papa.xlsx")
   fenotipo_anto=fenopapa[,c("Taxa","Delphinidin",	"Cyanidin",	"Petunidin","Pelargonidin",	"Peonidin")]
   fenotipo_anto2=na.omit(fenotipo_anto)
    for (i in 1:length(x)) {
       if(x[i]==1){
        y[i]=rnorm(sum(x[i]==1),mean = mean(fenotipo_anto2$Delphinidin), sd = sd(fenotipo_anto2$Delphinidin))
       }else{
        y[i]=rnorm(sum(x[i]==0),mean = mean(fenotipo_anto2$Peonidin), sd = sd(fenotipo_anto2$Peonidin))
       }
    }
    return(y)
}

random_norm_cov<-function(x,a, b,c){
    y=NULL
    for (i in 1:length(x)) {
       if(x[i]==1){
        y[i]=rnorm(sum(x[i]==1),mean = a, sd = c)
       }else{
        y[i]=rnorm(sum(x[i]==0),mean = b, sd = c)
       }
    }
    return(y)
}

ET_function<-function(matcov, Z){
  valprop=eigen(matcov)$values[1]
  vecprop=eigen(matcov)$vectors[,1]
  i=1
  ET=NULL
  for (i in 1:ncol(Z)){
    ET[i]=((t(Z[,i])%*%vecprop)^2)/valprop
  }
  ET_Pval=pchisq(ET, df = 1)
  re=list(ET.test=ET, ET.pvalue=ET_Pval)
  return(re)
}

OT_function<-function(matcov, Z){
  OT=NULL
  for (i in 1:ncol(Z)) {
    OT[i]=t(Z[,i])%*%solve(matcov)%*%(Z[,i])
  }

  OT_Pval=pchisq(OT, df=6)
  re=list(OT.test=OT, OT.pvalue=OT_Pval)
  return(re)
}


AT_function<-function(ET, OT){
   p=NULL
   AT=matrix(ncol = length(ET), nrow=1000)
   p=seq(0,1,length.out=1000)
   for (i in 1:nrow(AT)) {
      AT[i,]=(1-p[i])*(OT-ET)+p[i]*ET
   }
   AT_Pval=pchisq(AT, df=5)
   AT_Pval_min=NULL
   for (i in 1:ncol(AT)) {
      AT_Pval_min[i]=min(AT_Pval[,i])
}
re=list(AT.test=AT, AT.pvalue=AT_Pval_min)
return(re)
}

AT_function2<-function(ET, OT, frecuencias){
   p=NULL
   AT=matrix(ncol = length(ET), nrow=1000)
   p=seq(0,1,length.out=1000)
   for (i in 1:nrow(AT)) {
      AT[i,]=(1-p[i])*(OT*frecuencias-ET*frecuencias)+p[i]*ET*frecuencias
   }
   AT_Pval=pchisq(AT, df=5)
   AT_Pval_min=NULL
   for (i in 1:ncol(AT)) {
      AT_Pval_min[i]=min(AT_Pval[,i])
}
re=list(AT.test=AT, AT.pvalue=AT_Pval_min)
return(re)
}

AT_function3<-function(ET, OT, freq2){
   p=NULL
   AT=matrix(ncol = length(ET), nrow=1000)
   p=seq(0,1,length.out=1000)
   for (i in 1:nrow(AT)) {
      AT[i,]=(1-p[i])*(OT/freq2-ET/freq2)+p[i]*ET/freq2
   }
   AT_Pval=pchisq(AT, df=5)
   AT_Pval_min=NULL
   for (i in 1:ncol(AT)) {
      AT_Pval_min[i]=min(AT_Pval[,i])
}
re=list(AT.test=AT, AT.pvalue=AT_Pval_min)
return(re)
}

simulacion_completa<-function(nsnps.cor, ndata, nsnps.no.cor, datos=NULL){
   a=create_var(nsnps.cor)
   datos2=draw.correlated.binary(ndata, ncol(a), rep(0.5,ncol(a)), abs(a))
   n1=nsnps.no.cor
   #snps2=draw.correlated.binary(1000, n1, rep(0.5,n1), diag(x=1, ncol = n1, nrow = n1))
   fenopapa=read_excel("Accessiones_Papa.xlsx")
   fenotipo_anto=fenopapa[,c("Taxa","Delphinidin",	"Cyanidin",	"Petunidin","Pelargonidin",	"Peonidin")]
   fenotipo_anto2=na.omit(fenotipo_anto)

   probs=c(seq(0.05, 0.35, length.out=n1))
   random_binom=matrix(ncol=n1, nrow = nrow(datos2))
   for (i in 1:n1) {
      random_binom[,i]=rbinom(nrow(datos2),1, probs[i])
   }
   snps2=random_binom

   y_continua=random_norm(datos2[,1])
   y_continua2=random_norm_cov(datos2[,1],a=2, b=3,c=1)
   y_continua3=random_norm_cov(datos2[,1],a=2, b=4.5,c=2)
   y_continua4=random_norm_cov(datos2[,1],a=1.5, b=3.5,c=2)
   y_continua5=random_norm_cov(datos2[,1],a=2.5, b=5.5,c=2)
   x_continua=random_norm_cov(datos2[,1],a=3, b=5,c=2)
   Y=cbind(y_continua, y_continua2,y_continua3,y_continua4,y_continua5,datos2[,1])
   X=rep(x_continua)
   G=cbind(datos2[,2:ncol(a)],snps2)

   medias_columna=colMeans(G)
   freq_alelo_raro=ifelse(medias_columna>=0.5, 1-medias_columna, medias_columna)
   frecuencia_raro=freq_alelo_raro*nrow(G)

   colnames(G)=seq(1:ncol(G))
   Z=metodo(X,G,Y)$Z.score
   matcov=metodo(X,G,Y)$m.cov
   ET=ET_function(matcov, Z)$ET.test
   OT=OT_function(matcov, Z)$OT.test

   ET_Pval=ET_function(matcov, Z)$ET.pvalue
   OT_Pval=OT_function(matcov, Z)$OT.pvalue
   AT_Pval_min=AT_function(ET, OT)$AT.pvalue
   AT_Pval_min2=AT_function2(ET, OT, freq_alelo_raro)$AT.pvalue
   AT_Pval_min3=AT_function3(ET, OT, frecuencia_raro)$AT.pvalue
   result=tibble(snp=colnames(G),ET_Pval,OT_Pval,AT_Pval_min,AT_Pval_min2,AT_Pval_min3)
   return(result)
}

potencia=function(data,alpha,M, variable){
    h1=sum(data$H_0)
    A=nrow(data)
    lista=c()
    for (i in 1:nrow(data)){
        if((data$H_0[i]==0)&(data[[variable]][i]>alpha)){
            lista[i]=1
        } 
        else {
           lista[i]=0
        } 
}
U_n=sum(lista)
power=1-(1/A)*(1/h1)*U_n
return(power)
}

error_tipo_I=function(data,alpha,M, variable){
    A=nrow(data)
    lista=c()
    for (i in 1:nrow(data)){
        if((data$H_0[i]==1)&(data[[variable]][i]<=alpha)){
            lista[i]=1
        } 
        else {
           lista[i]=0
        }  
}
V_n=sum(lista)
fder=1-(1/A)*V_n
return(fder)
}


error_tipo_I_combinado=function(data,M, variable){
    alpha=c(0.0000000000000000001,0.000000000000000001, 0.00000000000000001, 0.0000000000000001, 0.000000000000001,0.00000000000001,0.0000000000001,0.000000000001, 0.00000000001, 0.0000000001,0.000000001,0.00000001,0.0000001, 0.000001,0.00001, 0.0001, 0.001, seq(0.01, 0.5, by=0.01) )
    potencias=c()
    for (i in 1:length(alpha)) {
        potencias[i]=error_tipo_I(data=data,alpha= alpha[i], M=M,variable)
    }
return(potencias)
}

potencia_combinada=function(data,M, variable){
    alpha=c(0.0000000000000000001,0.000000000000000001, 0.00000000000000001, 0.0000000000000001, 0.000000000000001,0.00000000000001,0.0000000000001,0.000000000001, 0.00000000001, 0.0000000001,0.000000001,0.00000001,0.0000001, 0.000001,0.00001, 0.0001, 0.001, seq(0.01, 0.5, by=0.01) )
    potencias=c()
    for (i in 1:length(alpha)) {
        potencias[i]=potencia(data=data,alpha= alpha[i], M=M,variable)
    }
return(potencias)
}


Data.frame_pvalores<-function(X, G, Y, frecuencias, freq2){
   Z=metodo(X,G,Y)$Z.score
   matcov=metodo(X,G,Y)$m.cov
   ET=ET_function(matcov, Z)$ET.test
   OT=OT_function(matcov, Z)$OT.test

   ET_Pval=ET_function(matcov, Z)$ET.pvalue
   OT_Pval=OT_function(matcov, Z)$OT.pvalue
   AT_Pval_min=AT_function(ET, OT)$AT.pvalue
   AT_Pval_min2=AT_function2(ET, OT, frecuencias)$AT.pvalue
   AT_Pval_min3=AT_function3(ET, OT, freq2)$AT.pvalue
   result=tibble(snp=colnames(G),ET_Pval,OT_Pval,AT_Pval_min,AT_Pval_min2,AT_Pval_min3)
   return(result)
}