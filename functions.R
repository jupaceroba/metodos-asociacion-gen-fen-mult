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