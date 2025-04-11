# Main function iDIP.phylo
iDIP.phylo =
  function(abun,struc,tree){
    phyloData <- newick2phylog(tree)
    Temp <- as.matrix(abun[names(phyloData$leaves), ])
    nodenames = c(names(phyloData$leaves),names(phyloData$nodes));
    M = matrix(0,nrow=length(phyloData$leaves),ncol=length(nodenames),dimnames=list(names(phyloData$leaves),nodenames))
    for(i in 1:length(phyloData$leaves)){M[i,][unlist(phyloData$paths[i])]=rep(1,length(unlist(phyloData$paths[i])))}
    pA=matrix(0,ncol=ncol(abun),nrow=length(nodenames),dimnames=list(nodenames,colnames(abun)))
    for(i in 1:ncol(abun)){pA[,i]=Temp[,i]%*%M;}
    pB=c(phyloData$leaves,phyloData$nodes)
    n=sum(abun);N=ncol(abun);
    ga=rowSums(pA);
    gp=ga/n;TT=sum(gp*pB);
    G=sum(-pB[gp>0]*gp[gp>0]/TT*log(gp[gp>0]/TT))
    PD=sum(pB[gp>0]);
    H=nrow(struc);
    A=numeric(H-1);W=numeric(H-1);B=numeric(H-1);
    Diff=numeric(H-1);Prop=numeric(H-1);
    wi=colSums(abun)/n;
    W[H-1]=-sum(wi[wi>0]*log(wi[wi>0]));
    pi=sapply(1:N,function(k) pA[,k]/sum(abun[,k]))
    Ai=sapply(1:N,function(k) -sum(pB[pi[,k]>0]*pi[,k][pi[,k]>0]/TT*log(pi[,k][pi[,k]>0]/TT)))
    A[H-1]=sum(wi*Ai);
    if(H>2){
      for(i in 2:(H-1)){
        I=unique(struc[i,]);NN=length(I);
        pi=matrix(0,ncol=NN,nrow=nrow(pA));ni=numeric(NN);
        for(j in 1:NN){
          II=which(struc[i,]==I[j]);
          if(length(II)==1) {pi[,j]=pA[,II]/sum(abun[,II]);ni[j]=sum(abun[,II]);
          }else{pi[,j]=rowSums(pA[,II])/sum(abun[,II]);ni[j]=sum(abun[,II])}
        }
        #pi=sapply(1:NN,function(k) ai[,k]/sum(ai[,k]));
        wi=ni/sum(ni);
        W[i-1]=-sum(wi*log(wi))
        Ai=sapply(1:NN,function(k) -sum(pB[pi[,k]>0]*pi[,k][pi[,k]>0]/TT*log(pi[,k][pi[,k]>0]/TT)))
        A[i-1]=sum(wi*Ai);
      }
    }
    total=G-A[H-1];
    Diff[1]=(G-A[1])/W[1];
    Prop[1]=(G-A[1])/total;
    B[1]=exp(G)/exp(A[1]);
    if(H>2){14
      for(i in 2:(H-1)){
        Diff[i]=(A[i-1]-A[i])/(W[i]-W[i-1]);
        Prop[i]=(A[i-1]-A[i])/total;
        B[i]=exp(A[i-1])/exp(A[i]);
      }}
    #Gamma=exp(G)/TT;Alpha=exp(A)/TT;Diff=Diff;Prop=Prop;
    Gamma=exp(G); Alpha=exp(A); Diff=Diff; Prop=Prop;
    out=matrix(c(PD,TT,Gamma,Alpha,B,Prop,Diff),ncol=1)
    #out1=iDIP(abun,struc);
    #out=cbind(out1,out2);
    rownames(out) <- 
      c(
        paste("Faith's PD"),
        paste("mean_T"),
        paste0("PD_gamma"),
        paste0("PD_alpha.", (H-1):1),
        paste0("PD_beta.", (H-1):1),
        paste0("PD_prop.", (H-1):1),
        paste0("PD_diff.",(H-1):1)
      )
    return(out)
  }
