"composition" <- function(M){
    if(is.vector(M)){
        M <- rbind(M)
        rownames(M) <- NULL
    }	
    
    jj <- is_valid_composition(M)
    if(isTRUE(jj)){	      
        class(M) <- "composition"  # This is the only occurence of class(.) <- "composition"
        return(M)
    } else {
        stop(jj)
    }
}	      

"as.composition" <- function(x){
    x <- rbind(x)
    rownames(x) <- NULL
    composition(sweep(x,1,rowSums(x),"/")) # much MUCH faster than t(apply(ta, 1, function(x) x/sum(x)))
}

"is_valid_composition" <- function(x,tol=1e-10){
    if (!is.matrix(x)){	
        return("x is not a matrix")
    } else if (any(is.nan(x))){
        return("x has NaN elements")
    } else if (any(is.na(x))){
        return("x has NA elements")
    } else if (any(x<0)){
        return("x has negative elements")
    } else if (any(abs(rowSums(x)-1) > tol)){
        return("x has rowSums not equal to 1 to within tolerance")
    } else {    # everything is OK;
        return(TRUE)
    }
}

"[<-.composition" <- function(x,i,j,...){
    stop("only complete rows of a composition may be changed (the columns are not independent)")
} 

"CVA" <- function(x){   # CVA == compositional variation array
    x <- as.composition(x)
    
    D <- ncol(x)
    d <- D-1
    out <- matrix(0,D,D)
    dimnames(out) <- list(means=colnames(x),variances=colnames(x))
    
    for(i in seq_len(d)){
        for(j in seq(from=i+1,to=D,by=1)){
            out[i,j] <-  var(log(x[,i]/x[,j]))   # variances; tau[i,j]
            out[j,i] <- mean(log(x[,i]/x[,j]))   # means; xi[i,j]; note swapped indices on LHS
        }
    }
    return(out)
}

"mean.composition" <- function(x){
    colMeans(as.composition(x))
}

"summary.composition" <- function(x){
    x <- as.composition(x)
    return(list(means=mean.composition(x),CVA=CVA(x)))
}

"print.composition" <- function(x){
    print(unclass(drop(x)))
}

"print.summary.composition" <- function(x){
print(x)
}




B <- function(a,w,log=TRUE){   # normalizing factor for weighted Dirichlet
    out <- (
        -lgamma(sum(a))
        +sum(lgamma(a))
        +sum(a*log(w))
        )
    if(log){
        return(out)
    } else {
        return(exp(out))
    }
}

dwd_single <- function(q, a, w, log=FALSE){  # PDF for a single observation of q's
    q <- q/sum(q)
    out <- (
        -B(a,w,log=TRUE)
        +sum((a-1)*log(q))
        -sum(a)*log(sum(q/w))
        )
    if(log){
        return(out)
    } else {
        return(exp(out))
    }
}
    
dwd <- function(q, a, w, log=FALSE){    # PDF for a matrix 'q' with rows being independent
    q <- q/rowSums(q)  #silently normalize so rowSums(q)==1
    qbar <- rowSums(sweep(q,2,w,"/"))
    out <- (
        -nrow(q)*B(a,w,log=TRUE)
        +sum(tcrossprod(a-1,log(q)))
        -sum(log(qbar)*sum(a))
        )
    if(log){
        return(out)
    } else {
        return(exp(out))
    }
}

p_to_q <- function(p,w){   # convert from p to q using weights w
    if(is.vector(p)){
        return(p*w/sum(p*w))
    } else if(is.matrix(p)){
        jj <- sweep(p, 2, w, "*")
        return(jj/rowSums(jj))
    } else {
        stop('p must be a vector or a matrix')
    }
}
        
q_to_p <- function(q,w){p_to_q(q,1/w)}  # inverse transform from q to p

pbar <- function(p,w){sum(p*w)}
qbar <- function(p,w){pbar(q,1/w)}

rwd <- function(n,a,w){
    la <- length(a)
    jj <- matrix(rgamma(n*la, a, 1), ncol=la,  byrow=TRUE)
    jj <- jj/rowSums(jj)

    jj <- sweep(jj, 2, w, "*")
    jj/rowSums(jj)
}

OO  <- c(1,5,6,         1.6,1.65,1.3)

f <- function(y){
    a <- y[1:3]
    w <- y[4:6]
    return(-dwd(x,a,w,log=TRUE))
}
    
Hessian <- function(q,a,w){
    N <- nrow(q)
    n <- length(a)
    iW <- diag(1/w)
    P <- q_to_p(q,w)
    bit1 <- N*(diag(psigamma(a))-matrix(1,n,n)*psigamma(sum(a)))
    bit2 <- N*iW - iW %*% t(P) %*% matrix(1,N,length(a))
    bit3 <- iW %*% crossprod(P) %*% iW
    return(    list(bit1,bit2,bit3))
}
