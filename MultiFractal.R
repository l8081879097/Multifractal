# Coded by Jierui LI on July 29,2022

# ts: time series data (no less than 0, or should be normalized)
# q: maximum moment
# step: step between moment array
# nbox: number of box
# t: box size array (based on the length of ts and nbox)

multi_fractal <- function(ts,q,step,nbox){
  len_ts <- length(ts)
  t <- array(NA,dim = nbox)
  t[nbox] <- len_ts # should include the box containing all mass
  for (x in 1:(nbox-1)){
    t[nbox-x] <- len_ts/2^x
  } # tr_array
  len_t <- length(t) # t
  qmin <- -q
  qmax <- q
  
  Mqt <- array(NA,dim = c(len_ts,len_t,length(seq(qmin,qmax,step)))) # matrix to store mass values
  for (q in 1:length(seq(qmin,qmax,step))){ # loop for each moment
    for (j in 1:len_t){ # loop for each scale (box size)
      # time resolution
      T <- t[j]
      # interval amount
      n <- len_ts%/%T
      # group ts
      w_mat <- matrix(NA,nrow = n,ncol = 4) # weight matrix
      tag <- seq(1,len_ts,T)[-1]
      m_tag <- floor(tag)
      w_mat[1,1:4] <- c(1,tag[1]-m_tag[1],1,m_tag[1])
      w_mat[n,1:4] <- c(m_tag[n-1]+1-tag[n-1],1,m_tag[n-1],len_ts)
      if (n > 2){
        for (z in 2:(n-1)) {
          w_mat[z,c(1,3,4)] <- c(1-w_mat[z-1,2],m_tag[z-1],m_tag[z])
          w_mat[z,2] <- T-w_mat[z,1]-(w_mat[z,4]-w_mat[z,3]-1)
        }
      }
      for (i in 1:n){ # loop for each interval
        if (is.integer(T)){
          Mqt[i,j,q] <- sum(ts[((i-1)*T+1):(i*T)],na.rm = T)
        }else{
          Mqt[i,j,q] <- as.numeric(sum(ts[(w_mat[i,3]+1):(w_mat[i,4]-1)])+ts[w_mat[i,3]]*w_mat[i,1]+ts[w_mat[i,4]]*w_mat[i,2])
        }  
        if (Mqt[i,j,q]==0){
          Mqt[i,j,q] <- 0
        } else {
          Mqt[i,j,q] <- Mqt[i,j,q]^(seq(qmin,qmax,step))[q]
        }
      }
    }
  }
  # linear regression
  M <- matrix(NA,nrow = length(seq(qmin,qmax,step)),ncol = len_t)
  rm_q <- array(1,dim = length(seq(qmin,qmax,step)))
  for (q in 1:length(seq(qmin,qmax,step))){
    for (j in 1:len_t){
      M[q,j] <- sum(Mqt[,j,q],na.rm = T)
    }
    if (any(is.na(M[q,]))){
      rm_q[q] <- 0
    }
  }
  tau <- array(NA,dim = length(seq(qmin,qmax,step)))
  for (qod in 1:length(seq(qmin,qmax,step))){
    if (rm_q[qod]==0){
      tau[qod] <- NA
    } else {
      tau[qod] <- lm(log(M[qod,])~log(t))$coefficients[2] # reason for deleting resolution 1
    }
  }
  fitmod <- loess(tau~seq(1,length(seq(qmin,qmax,step)),1),span=1, family="gaussian")
  tau_raw <- tau
  tau <- predict(fitmod, seq(1,length(seq(qmin,qmax,step)),1))
  
  # multifractal parameters
  alpha <- array(NA,dim = length(tau)-1)
  for (alod in 1:length(alpha)){
    alpha[alod] <- lm(tau[alod:(alod+1)]~seq(qmin,qmax,step)[alod:(alod+1)])$coefficients[2]
  }
  tau_new <- array(NA,dim = length(tau)-1)
  for (tau_len in 1:(length(tau)-1)){
    tau_new[tau_len] <- (tau[tau_len]+tau[tau_len+1])/2
  }
  q_new <- array(NA,dim = length(seq(qmin,qmax,step))-1)
  for (seqq in 1:(length(seq(qmin,qmax,step))-1)){
    q_new[seqq] <- mean(seq(qmin,qmax,step)[seqq:(seqq+1)])
  }
  f <- q_new*alpha-tau_new
  out <- cbind(alpha,f)
}
