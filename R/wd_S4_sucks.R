setClass("composition", representation = representation(x="matrix") )
extends("composition","matrix")
.as_matrix <- function(from){ return(from@x) }   # this is the only occurence of "@"

"composition" <- function(from){ new("composition",x=from) }    # This is the only occurence of new()

setAs("composition", "matrix", .as_matrix)   #coerces from a matrix to a compostion
setMethod("as.matrix",signature(x="composition"), function(x){as(x,"matrix")})

setAs( "matrix", "composition", composition)   #coerces from a matrix to a compostion
setMethod("as.matrix",signature(x="composition"), function(x){as(x,"composition")})

as.composition <- function(x){ #user-friendly version; coerces to matrix
    if(is.vector(x)){
        x <- rbind(x)
        rownames(x) <- NULL
    }
    return(composition(x))
}
 
setReplaceMethod("[","composition",function(x){stop("not yet implemented")})

setMethod("[","composition",function(x){stop("niiiiiiiot yet implemented")})
    
".comp_valid" <- function(object){
    x <- .as_matrix(object)


        if (!is.matrix(x)){
        return("x is not a matrix")
    } else if (any(is.nan(x))){
        return("x has NaN elements")
    } else if (any(is.na(x))){
        return("x has NA elements")
    } else if (any(x<0)){
        return("x has negative elements")
    } else if (any(abs(rowSums(x)-1) > 1e-7)){
        return("x has rowSums out of tolerance")
    } else {    # everything is OK;
        return(TRUE)
    }
}
setValidity("composition", .comp_valid)

C <- function(M){ sweep(M,1,rowSums(M),"/") }

dot <- function(M1,M2){
    if(nrow(M1)==1){
        return(C(sweep(M2,1,M1,"*")))
    } else if (nrow(M2)==1){
        return(C(sweep(M1,1,M2,"*")))
    } else {
        return(C(M1*M2))
    }
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
x <- rwd(23,a=OO[1:3],w=OO[4:6])
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
