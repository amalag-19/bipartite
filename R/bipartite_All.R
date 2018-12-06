########################################################################################################
########################################################################################################
########################### ALL TESTED BIPARTITE CODE FOR DIFFERENT CASES ##############################
########################################################################################################
########################################################################################################
## All models below consider static undirected density networks for Bipartite with proportional odds model for discrete weighted rating networks
library(Rcpp)
file_path<-"/Users/Amal/Box Sync/PSU/Fall 2018/Main_Research/Network Models/Project 5 (Bipartite)/code/bipartite_All/"

########################################################################################################
########################################################################################################
## Model 7 features (works):
## 1) Block structure over intercept
## 2) Proportional odds model for ratings based on my own derivations

wrapper_bipartite_model7<-function(net_adjacency,net_rating,nclust1,nclust2,thres=10^(-6),n_iter_min=50,n_iter_max=500,theta_init,delta_0_init,sim_indicator,theta_true=NA, delta_0_true=NA,K_true_U=NA,K_true_P=NA,cluster_ids_true_U=NA,cluster_ids_true_P=NA,R=NA){
  ## This is the combined final updated and working code for EM undirected case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(quadprog)
  library(combinat)
  library(MASS)
  library(Rcpp)
  library(Matrix)
  library(gtools)
  #sourceCpp(file = paste0(file_path,"model_7.cpp"))
  library(biSBM)
  #################################################
  ## Defining a function to update variational parameters gamma using quadratic program solver. Input: 
  ## gamma.curr is a N*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
  ## estimates, theta.curr is a T_grid*K array for current theta estimates. Given Network is assumed to be an N*N*T_data array
  solve_QP_wrapper<-function(quad_lin_coeff,constraint_matrix,constraint_vector,node_ID,gamma.curr){
    try_QP<-try(gamma_next_vec<-solve.QP((-2*diag(as.vector(quad_lin_coeff[node_ID,,1]))),as.vector(quad_lin_coeff[node_ID,,2]),t(constraint_matrix),constraint_vector)$solution,silent = TRUE)
    if('numeric' %in% class(try_QP)){
      gamma_next_vec<-try_QP
    }
    else{
      gamma_next_vec<-gamma.curr[node_ID,]
      print("red flag")
      print(paste0("Node ID for red flag is ", node_ID))
    }
    return(gamma_next_vec)
  }
  
  gamma_U.update.wrapper<-function(gamma_U.curr,gamma_P.curr,pi_U.curr,theta.curr, logitInv_delta_0.curr, net_adjacency, net_rating, N1,N2,K1,K2,R){
    gamma_U.next<-matrix(NA_real_,N1,K1)
    constraint_matrix<-matrix(NA_real_,2+K1,K1)
    constraint_matrix[1,]<-rep(1,K1)
    constraint_matrix[2,]<-rep(-1,K1)
    constraint_matrix[3:(K1+2),]<-diag(K1)
    constraint_vector<-c(1,-1,rep(0,K1))
    
    quad_lin_coeff<-gamma_U_update_stat_undir(gamma_U=gamma_U.curr, gamma_P=gamma_P.curr, pi_U=pi_U.curr, theta=theta.curr,logitInv_delta_0=logitInv_delta_0.curr,net_adjacency=net_adjacency, net_rating = net_rating, N1=N1, N2=N2, K1=K1, K2=K2, R=R)
    
    #print(quad_lin_coeff)
    for (i in 1:N1){
      if(sum(is.nan(quad_lin_coeff[i,,1]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,1])),1]<--Inf
      }
      if(sum(is.nan(quad_lin_coeff[i,,2]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,2])),2]<-Inf
      }
      if(sum(quad_lin_coeff[i,,1]<=-10^8)==0){
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma_U.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_U.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K1)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma_U.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_U.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K1){
          gamma_U.next[i,]<-gamma.curr[i,]
        }
      }else if((sum(quad_lin_coeff[i,,1]<=-10^8)>0)&&(sum(quad_lin_coeff[i,,1]<=-10^8)<K1)){
        quad_lin_coeff[i,which(quad_lin_coeff[i,,1]<=-10^8),1]<--(10^10)
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma_U.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_U.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K1)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma_U.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_U.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K1){
          gamma_U.next[i,]<-gamma_U.curr[i,]
        }
      }else if(sum(quad_lin_coeff[i,,1]<=-10^8)==K1){
        gamma_U.next[i,]<-gamma_U.curr[i,]
      }
      #print(i)
    }
    #print('gammaU')
    #print(quad_lin_coeff)
    # for (i in 1:N1){
    #   gamma_U.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    #   print(i)
    # }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma_U.next<-t(apply(X = gamma_U.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma_U.next)
  }
  
  gamma_P.update.wrapper<-function(gamma_U.curr,gamma_P.curr,pi_P.curr,theta.curr, logitInv_delta_0.curr, net_adjacency, net_rating, N1,N2,K1,K2,R){
    gamma_P.next<-matrix(NA_real_,N2,K2)
    constraint_matrix<-matrix(NA_real_,2+K2,K2)
    constraint_matrix[1,]<-rep(1,K2)
    constraint_matrix[2,]<-rep(-1,K2)
    constraint_matrix[3:(K2+2),]<-diag(K2)
    constraint_vector<-c(1,-1,rep(0,K2))
    
    quad_lin_coeff<-gamma_P_update_stat_undir(gamma_U=gamma_U.curr, gamma_P=gamma_P.curr, pi_P=pi_P.curr, theta=theta.curr, logitInv_delta_0=logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating = net_rating, N1=N1, N2=N2, K1=K1, K2=K2, R=R)
    #print('gammaP')
    #print(quad_lin_coeff)
    
    #print(quad_lin_coeff)
    for (i in 1:N2){
      if(sum(is.nan(quad_lin_coeff[i,,1]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,1])),1]<--Inf
      }
      if(sum(is.nan(quad_lin_coeff[i,,2]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,2])),2]<-Inf
      }
      if(sum(quad_lin_coeff[i,,1]<=-10^8)==0){
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma_P.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_P.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K2)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma_P.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_P.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K2){
          gamma_U.next[i,]<-gamma.curr[i,]
        }
      }else if((sum(quad_lin_coeff[i,,1]<=-10^8)>0)&&(sum(quad_lin_coeff[i,,1]<=-10^8)<K2)){
        quad_lin_coeff[i,which(quad_lin_coeff[i,,1]<=-10^8),1]<--(10^10)
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma_P.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_P.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K2)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma_P.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_P.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K2){
          gamma_P.next[i,]<-gamma_P.curr[i,]
        }
      }else if(sum(quad_lin_coeff[i,,1]<=-10^8)==K2){
        gamma_P.next[i,]<-gamma_P.curr[i,]
      }
      #print(i)
    }
    
    # for (i in 1:N2){
    #   gamma_P.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    #   print(i)
    # }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma_P.next<-t(apply(X = gamma_P.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma_P.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of pi
  pi.update<-function(gamma.curr){
    pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## normalization of pi
    pi.next<-pi.next/sum(pi.next)
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions K1,K2
  theta.update<-function(theta.curr,gamma_U,gamma_P,net_adjacency,N1,N2,K1,K2){
    theta.next<-matrix(NA_real_,K1,K2)
    gradient<-grad_bipartite_stat_undir_theta(theta=theta.curr, gamma_U=gamma_U, gamma_P=gamma_P, net_adjacency=net_adjacency, N1=N1, N2=N2, K1=K1, K2=K2)
    hess<-hess_bipartite_stat_undir_theta(theta=theta.curr, gamma_U=gamma_U, gamma_P=gamma_P, N1=N1, N2=N2, K1=K1, K2=K2)
    for(k in 1:K1){
      for(l in 1:K2){
        theta.next[k,l]<-theta.curr[k,l]-gradient[k,l]/hess[k,l]
      }
    }
    return(theta.next)
  }
  
  #################################################
  ## Defining the update function of full proportional odds nuisance parameters delta_0
  delta_0.update<-function(gamma_U,gamma_P,delta_0,logitInv_delta_0,logitInv_phi,net_adjacency,net_rating_statistic,N1,N2,K1,K2,R){
    delta_0.next<-array(NA_real_,dim=c(K1,K2,(R-1)))
    
    ## Getting the estimated cluster memberships of all nodes in two sets from variational E step
    #clust_est_U<-as.vector(apply(X = gamma_U,MARGIN = 1,FUN = which.max))
    #clust_est_P<-as.vector(apply(X = gamma_P,MARGIN = 1,FUN = which.max))
    
    #ratings_list<-tie_clust_partition(clust_est_U=clust_est_U, clust_est_P=clust_est_P, net_adjacency=net_adjacency, net_rating=net_rating, N1=N1, N2=N2, K1=K1, K2=K2)
    
    gradient_delta_0<-grad_bipartite_stat_undir_delta_0(gamma_U=gamma_U, gamma_P=gamma_P, delta_0=delta_0, logitInv_delta_0=logitInv_delta_0, net_adjacency=net_adjacency, net_rating=net_rating, N1=N1, N2=N2, K1=K1, K2=K2, R=R)
    
    for (k in 1:K1){
      for (l in 1:K2){
        hess_mat<-hess_bipartite_stat_undir_delta_0(k=k-1, l=l-1, gamma_U=gamma_U, gamma_P=gamma_P, delta_0=delta_0, logitInv_delta_0=logitInv_delta_0, net_adjacency=net_adjacency, net_rating=net_rating, N1=N1, N2=N2, R=R)
        delta_0.next[k,l,]<-sort(delta_0[k,l,]-as.vector((solve(hess_mat+(10^(-6)*diag(R-1)))%*%gradient_delta_0[k,l,])))
        #delta_0.next[k,l,]<-sort(delta_0[k,l,]-as.vector((solve(hess_mat)%*%gradient_delta_0[k,l,])))
      }
    }
    #print(delta_0.next)
    return(delta_0.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,net_adjacency,net_rating,K1,K2,n_iter,thres,n_iter_min,n_iter_max,R){
    ## Assuming a network adjacency matrix of dimensions N1*N2 with rows corresponding to users and columns corresponding to products and each entry is an indicator whether the user rated the product. This framework can be generalized.
    N1<-dim(net_adjacency)[1]
    N2<-dim(net_adjacency)[2]
    
    ## initializing the arrays for parameters
    ## Initializing the cluster memberships gamma
    gamma_U<-array(NA_real_,dim=c(N1,K1,n_iter)) ## assuming there are K1 clusters in user set
    gamma_P<-array(NA_real_,dim=c(N2,K2,n_iter)) ## assuming there are K2 clusters in product set
    
    ## Initializing the mixture proportions pi
    pi_U<-matrix(NA_real_,K1,n_iter)
    pi_P<-matrix(NA_real_,K2,n_iter)
    
    ## Initializing the network parameters theta
    theta<-array(NA_real_,dim=c(K1,K2,n_iter))
    delta_0<-array(NA_real_,dim=c(K1,K2,R-1,n_iter))
    
    gamma_U[,,1]<-start[[1]]
    gamma_P[,,1]<-start[[2]]
    pi_U[,1]<-start[[3]]
    pi_P[,1]<-start[[4]]
    theta[,,1]<-start[[5]]
    delta_0[,,,1]<-start[[6]]
    
    #print(theta[,,1])
    #print(omega[,1])
    #print(mu[,,,1])
    
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_curr<-10^10
    time_per_iter<-c()
    while(((error>thres)&(iter_index<n_iter_max))|(iter_index<n_iter_min)){ 
      #while((error>thres|iter_index<200)&(iter_index<250)&(iter_index==2|error<10)){ 
      #while(iter_index<21){ 
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      
      ## Initializing the G arrays
      #Gvals_cube_0<-array(NA_real_,dim=c(K1,K2,R))
      #Gvals_cube_1<-array(NA_real_,dim=c(K1,K2,R))
      #Gvals_cube_2<-array(NA_real_,dim=c(K1,K2,R))
      
      ## Updating the G arrays
      #print(mu[,,,iter_index-1])
      #if (R==1){mu[,,,iter_index-1]=array(mu[,,,iter_index-1],dim=c(2,2,1))}
      #Gvals_cube_0<-get_Gval_0(omega=omega[,,iter_index-1], mu=mu[,,,iter_index-1],K1=K1,K2=K2,R=R)
      #Gvals_cube_1<-get_Gval_1(omega=omega[,,iter_index-1], mu=mu[,,,iter_index-1],K1=K1,K2=K2,R=R)
      #Gvals_cube_2<-get_Gval_2(omega=omega[,,iter_index-1], mu=mu[,,,iter_index-1],K1=K1,K2=K2,R=R)
      #print('Gvals_cube_0')
      #print(Gvals_cube_0)
      #Gvals_mat_0=apply(X = Gvals_cube_0,MARGIN = c(1,2),sum)
      #Gvals_mat_1=apply(X = Gvals_cube_1,MARGIN = c(1,2),sum)
      #Gvals_mat_2=apply(X = Gvals_cube_2,MARGIN = c(1,2),sum)
      
      #log_Gval_0 <- matrix(NA_real_,K1,K2)
      #log_Gval_0 <- get_log_Gval_0(omega=omega[,,iter_index-1], mu=mu[,,,iter_index-1],K1=K1,K2=K2,R=R)
      
      ## Calculating the logit inverse of delta_0
      logitInv_delta_0.curr<-logitInv_delta_0_cal(delta_0=delta_0[,,,iter_index-1], K1=K1, K2=K2, R=R)
      
      ## Transforming logit inverse of delta_0 to phi
      #phi.curr<-phi_cal(logitInv_delta_0=logitInv_delta_0.curr, K1=K1, K2=K2, R=R)
      
      ## Calculating the logit inverse of phi
      #logitInv_phi.curr<-logitInv_phi_cal(phi=phi.curr, K1=K1, K2=K2, R=R)
      
      ## Updating gamma_U, the mixed memberships for users
      gamma_U[,,iter_index]<-gamma_U.update.wrapper(gamma_U.curr=gamma_U[,,iter_index-1], gamma_P.curr=gamma_P[,,iter_index-1], pi_U.curr=pi_U[,iter_index-1], theta.curr=theta[,,iter_index-1], logitInv_delta_0.curr=logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating=net_rating, N1=N1,N2=N2,K1=K1,K2=K2,R=R)
      #print('gamma_U[,,iter_index]')
      # print(gamma_U[,,iter_index])
      ## Updating the mixture proportions pi_U
      pi_U[,iter_index]<-pi.update(gamma.curr=gamma_U[,,iter_index])
      #print('pi_U[,,iter_index]')
      #print(pi_U[,,iter_index])
      ## Updating gamma_P, the mixed memberships for products
      gamma_P[,,iter_index]<-gamma_P.update.wrapper(gamma_U.curr=gamma_U[,,iter_index],gamma_P.curr=gamma_P[,,iter_index-1],pi_P.curr=pi_P[,iter_index-1], theta.curr=theta[,,iter_index-1], logitInv_delta_0.curr=logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating=net_rating, N1=N1,N2=N2,K1=K1,K2=K2,R=R)
      # print('gamma_P[,,iter_index]')
      #print(gamma_P[,,iter_index])
      ## Updating the mixture proportions pi_P
      pi_P[,iter_index]<-pi.update(gamma.curr=gamma_P[,,iter_index])
      #print('pi_P[,,iter_index]')
      #print(pi_P[,,iter_index])
      ## Updating the network parameters theta
      theta[,,iter_index]<-theta.update(theta.curr=theta[,,iter_index-1],gamma_U=gamma_U[,,iter_index], gamma_P=gamma_P[,,iter_index], net_adjacency=net_adjacency, N1=N1, N2=N2, K1=K1, K2=K2)
      #print('theta')
      #print(theta[,,iter_index])
      
      ## Updating the omega matrix
      delta_0[,,,iter_index]<-delta_0.update(gamma_U=gamma_U[,,iter_index], gamma_P=gamma_P[,,iter_index], delta_0=delta_0[,,,iter_index-1],logitInv_delta_0 = logitInv_delta_0.curr,logitInv_phi = logitInv_phi.curr, net_adjacency=net_adjacency, net_rating_statistic = net_rating_statistic, N1=N1, N2=N2, K1=K1, K2=K2, R=R)
      #print('delta_0')
      #print(delta_0[,,,iter_index])
      #Sys.sleep(1.5)
      
      ELBO_prev<-ELBO_curr
      ELBO_curr<-ELBO_conv_bipartite_stat_undir(gamma_U=gamma_U[,,iter_index], gamma_P=gamma_P[,,iter_index], pi_U=pi_U[,iter_index], pi_P=pi_P[,iter_index], theta = theta[,,iter_index], logitInv_delta_0=logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating=net_rating, N1=N1, N2=N2, K1=K1, K2=K2)
      
      error<-abs(((ELBO_prev-ELBO_curr)/ELBO_prev))
      #print('error')
      #print(error)
      
      #print("probs")
      #print(Gvals_cube_0)
      # print(apply(X = Gvals_cube_0,MARGIN = c(1,2),FUN = function(x){
      #   
      #   return(x/sum(x))
      # }))
      
      if((iter_index%%1)==0){
        #print(theta[,,iter_index])
        print('iter_index')
        print(iter_index)
        print('error')
        print(error)
        print(proc.time()-ptm)
      }
      ptm_diff<-proc.time()-ptm
      time_per_iter<-c(time_per_iter,ptm_diff[3])
      
      #print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    names(time_per_iter)<-NULL
    total_iterations<-iter_index-1
    return(list(gamma_U,gamma_P,pi_U,pi_P,theta,delta_0,time_per_iter,total_iterations))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining Model Selection functions based on converged paramters
  
  #################################################
  ## Defining a function to calculate the estimate of complete log likelihood
  comp_loglik_cal<-function(cluster_ids_est_U = NA, cluster_ids_est_P = NA, pi_U = NA, pi_P = NA, theta, Gvals_cube_0, Gvals_mat_0, net_adjacency, net_rating, N1, N2, K1, K2, R=R){
    if((K1!=1)&(K2!=1)){
      ## 1st term
      t1<-0
      for(i in 1:N1){
        for(j in 1:N2){
          cluster_id_i<-cluster_ids_est_U[i]
          cluster_id_j<-cluster_ids_est_P[j]
          exp_val<-exp(theta[cluster_id_i,cluster_id_j])
          if(net_adjacency[i,j]==0){
            t1<-t1-(log(1+exp_val))
          }else{
            t1<-t1+(theta[cluster_id_i,cluster_id_j]-log(1+exp_val))+log(Gvals_cube_0[cluster_id_i,cluster_id_j,net_rating[i,j]])-log(Gvals_mat_0[cluster_id_i,cluster_id_j])
          }
          
        }
      }
      ## 2nd term
      t2<-0
      for (i in 1:N1){
        t2<-t2+log(pi_U[cluster_ids_est_U[i]])
      }
      t3<-0
      for (j in 1:N2){
        t3<-t3+log(pi_P[cluster_ids_est_P[j]])
      }
      
      comp_val<-t1+t2+t3
    }
    return(comp_val)
  }
  
  #################################################
  ## Defining a function to calculate the integrated classification likelihood
  ICL_cal<-function(cluster_ids_est_U, cluster_ids_est_P, pi_U, pi_P, theta, Gvals_cube_0, Gvals_mat_0, net_adjacency, net_rating, N1, N2, K1, K2, R){
    t1<-comp_loglik_cal(cluster_ids_est_U = cluster_ids_est_U, cluster_ids_est_P = cluster_ids_est_P, pi_U = pi_U, pi_P = pi_P, theta = theta, Gvals_cube_0 = Gvals_cube_0, Gvals_mat_0 = Gvals_mat_0, net_adjacency = net_adjacency, net_rating = net_rating, N1 = N1, N2 = N2, K1 = K1, K2 = K2, R=R)
    ## Penalization over pi
    t2<-(K1-1)*log(N1)+(K2-1)*log(N2)
    ## Penalization over theta
    t3<-(K1*K2)*log(N1*N2)
    ## Calculating the number of edges
    N0_val<-N0_cal(net_adjacency=net_adjacency, N1=N1, N2=N2)
    ## Penalization over omega, mu
    t4<-(K1*K2*(R-1))*log(N0_val)
    ICL_val<-t1-t2-t3-t4
    return(ICL_val)
  }
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI_cal<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  #################################################
  ## RASE functions
  RASE_cal<-function(param_est, param_true){
    if(is.vector(param_est)!=1){
      RASE_val<-sqrt(sum((param_est-param_true)^2)/prod(dim(param_est)))
    }else{RASE_val<-sqrt(sum((param_est-param_true)^2)/length(param_est))}
    return(RASE_val)
  }
  
  ########################################################################################################
  ## Defining the functions to evaluate clustering performance using normalized mutual information.
  probab_clust_cal<-function(cluster_ids,N,K){
    probab_vec<-rep(NA_real_,K)
    for(k in 1:K){
      probab_vec[k]<-sum(cluster_ids==k)/N
    }
    return(probab_vec)
  }
  NMI_cal<-function(cluster_ids_true_U,cluster_ids_true_P,K_true_U,K_true_P,cluster_ids_est_U,cluster_ids_est_P,K1,K2,N1,N2){
    probab_vec_true_U<-probab_clust_cal(cluster_ids = cluster_ids_true_U,N = N1+N2,K = K_true_U)
    probab_vec_true_P<-probab_clust_cal(cluster_ids = cluster_ids_true_P,N = N1+N2,K = K_true_P)
    probab_vec_true<-c(probab_vec_true_U,probab_vec_true_P)
    total_entropy_true<-0
    for(k in 1:length(probab_vec_true)){
      if(probab_vec_true[k]!=0){
        total_entropy_true<-total_entropy_true-(probab_vec_true[k]*log(probab_vec_true[k],base = 2))
      }
    }
    
    probab_vec_est_U<-probab_clust_cal(cluster_ids = cluster_ids_est_U,N = N1+N2,K = K1)
    probab_vec_est_P<-probab_clust_cal(cluster_ids = cluster_ids_est_P,N = N1+N2,K = K2)
    probab_vec_est<-c(probab_vec_est_U,probab_vec_est_P)
    total_entropy_est<-0
    for(k in 1:length(probab_vec_est)){
      if(probab_vec_est[k]!=0){
        total_entropy_est<-total_entropy_est-(probab_vec_est[k]*log(probab_vec_est[k],base = 2))
      }
    }
    
    clust_ids_true_total<-c(cluster_ids_true_U,(cluster_ids_true_P+K_true_U))
    clust_ids_est_total<-c(cluster_ids_est_U,(cluster_ids_est_P+K1))
    
    conditional_entropy<-rep(NA_real_,(K1+K2))
    for(k in 1:(K1+K2)){
      probab_conditional_vec<-rep(0,(K_true_U+K_true_P))
      for(m in 1:(K_true_U+K_true_P)){
        if(sum(clust_ids_est_total==k)!=0){
          probab_conditional_vec[m]<-sum((clust_ids_true_total==m)&(clust_ids_est_total==k))/sum(clust_ids_est_total==k)
        }
      }
      plogp_conditional_vec<-rep(0,length(probab_conditional_vec))
      for(m in 1:length(probab_conditional_vec)){
        if(probab_conditional_vec[m]!=0){
          plogp_conditional_vec[m]<-probab_conditional_vec[m]*log(x = probab_conditional_vec[m],base = 2)
        }
      }
      
      conditional_entropy[k]<--probab_vec_est[k]*sum(plogp_conditional_vec)
    }
    
    mutual_info<-total_entropy_true-sum(conditional_entropy)
    NMI_val<-(2*mutual_info)/(total_entropy_true+total_entropy_est)
    return(NMI_val)
  }
  
  ########################################################################################################
  ## Calculating the recovery rate
  recovery_rate_cal <- function(cluster_ids_est, cluster_ids_true){
    all_perm<-permutations(n=4,r=2,repeats.allowed = T)
    clust_perm_mat<-matrix(NA_integer_,4,4)
    for(i in 1:nrow(all_perm)){
      clust_perm_mat[all_perm[i,1],all_perm[i,2]] <- sum((cluster_ids_true==all_perm[i,1])&(cluster_ids_est==all_perm[i,2]))
    }
    recovery_rate_val<-sum(apply(X = clust_perm_mat, MARGIN = 1, FUN = max))/length(cluster_ids_est)
    return(recovery_rate_val)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N1<-dim(net_adjacency)[1] ## Number of nodes from the network
  N2<-dim(net_adjacency)[2]
  ## Total time points in the network for which we have data in the form of adjacency matrix
  K1<-nclust1 ## Defining the number of clusters
  K2<-nclust2
  K <-K1*K2
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  start[[1]]<-matrix(rep(1/K1,N1*K1),N1,K1)
  start[[2]]<-matrix(rep(1/K2,N2*K2),N2,K2)
  start[[3]]<-rep(1/K1,K1)
  start[[4]]<-rep(1/K2,K2)
  start[[5]]<-theta_init
  start[[6]]<-delta_0_init
  if(K1!=1&K2!=1){ param<-iterator(start=start, net_adjacency = net_adjacency, net_rating=net_rating, K1=K1, K2=K2, n_iter=1000, thres=thres, n_iter_min=n_iter_min, n_iter_max=n_iter_max, R=R)}
  
  #################################################
  ## extracting the coverged parameter values and calculating ICL
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    temp<-is.na(param[[1]][1,1,n_iter])
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  #print(param)
  
  if(K1!=1&K2!=1){
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,,n_last]
    param_converge[[3]]<-param[[3]][,n_last]
    param_converge[[4]]<-param[[4]][,n_last]
    param_converge[[5]]<-param[[5]][,,n_last]
    param_converge[[6]]<-param[[6]][,,,n_last]
    output_list<-param_converge
    cluster_ids_est_U<-as.vector(apply(X = matrix(1:N1),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[1]][x,])
      return(cluster_id)
    }))
    cluster_ids_est_P<-as.vector(apply(X = matrix(1:N2),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[2]][x,])
      return(cluster_id)
    }))
    ICL_val<-NA
    #ICL_val<-ICL_cal(cluster_ids_est_U = cluster_ids_est_U, cluster_ids_est_P = cluster_ids_est_P, pi_U=param_converge[[3]], pi_P=param_converge[[4]] ,theta = param_converge[[5]], Gvals_cube_0 = Gvals_cube_0, Gvals_mat_0 = Gvals_mat_0, net_adjacency=net_adjacency, net_rating=net_rating, N1 = N1, N2 = N2, K1 = K1, K2 = K2, R = R)
    if(sim_indicator==1){
      RI_val_U<-RI_cal(cluster_ids_est = cluster_ids_est_U,cluster_ids_true = cluster_ids_true_U)
      RI_val_P<-RI_cal(cluster_ids_est = cluster_ids_est_P,cluster_ids_true = cluster_ids_true_P)
      NMI_val<-NMI_cal(cluster_ids_true_U = cluster_ids_true_U,cluster_ids_true_P = cluster_ids_true_P,K_true_U = K_true_U,K_true_P = K_true_P,cluster_ids_est_U = cluster_ids_est_U,cluster_ids_est_P = cluster_ids_est_P, K1 = K1, K2 = K2, N1 = N1, N2 = N2)
      Recovery_rate_val_U<-recovery_rate_cal(cluster_ids_est = cluster_ids_est_U,cluster_ids_true = cluster_ids_true_U)
      Recovery_rate_val_P<-recovery_rate_cal(cluster_ids_est = cluster_ids_est_P,cluster_ids_true = cluster_ids_true_P)
      if((K1==K_true_U)&(K2==K_true_P)){
        K1_permute_mat<-do.call(rbind,permn(1:K1))
        K2_permute_mat<-do.call(rbind,permn(1:K2))
        
        RASE_theta_mat<-matrix(NA_real_,nrow(K1_permute_mat),nrow(K2_permute_mat))
        for (k in 1:nrow(K1_permute_mat)){
          for (l in 1:nrow(K2_permute_mat)){
            theta_true_rotated<-as(as.integer(K1_permute_mat[k,]), "pMatrix")%*%theta_true%*%t(as(as.integer(K2_permute_mat[l,]), "pMatrix"))
            RASE_theta_mat[k,l]<-RASE_cal(param_est = param_converge[[5]], param_true = theta_true_rotated)
          }
        }
        
        RASE_theta_val<-min(RASE_theta_mat)
        rotate_indices<-which(RASE_theta_mat == min(RASE_theta_mat), arr.ind = TRUE)
        
        delta_0_true_rotated<-array(NA_real_,dim=c(K1,K2,(R-1)))
        for (r in 1:(R-1)){
          delta_0_true_rotated[,,r]<-as(as.integer(K1_permute_mat[rotate_indices[1],]), "pMatrix")%*%delta_0_true[,,r]%*%t(as(as.integer(K2_permute_mat[rotate_indices[2],]), "pMatrix"))
          
        }
        RASE_delta_0_val<-RASE_cal(param_est = param_converge[[6]], param_true = delta_0_true_rotated)
        
        output_list<-list("parameters_converged"=param_converge,"Set1_estimated_cluster_IDs"=cluster_ids_est_U,"Set2_estimated_cluster_IDs"=cluster_ids_est_P,"ICL"=ICL_val,"Set1_Rand_Index"=RI_val_U,"Set2_Rand_Index"=RI_val_P,"Normalized_Mutual_Info"=NMI_val,"Set1_Recovery_Rate"=Recovery_rate_val_U, "Set2_Recovery_Rate"=Recovery_rate_val_P,"RASE_theta"=RASE_theta_val,"RASE_delta_0"=RASE_delta_0_val,"time_per_iter_vec"=param[[7]], "total_iterations"=param[[8]])
      }else{
        output_list<-list("parameters_converged"=param_converge,"Set1_estimated_cluster_IDs"=cluster_ids_est_U,"Set2_estimated_cluster_IDs"=cluster_ids_est_P,"ICL"=ICL_val,"Set1_Rand_Index"=RI_val_U,"Set2_Rand_Index"=RI_val_P,"Normalized_Mutual_Info"=NMI_val,"Set1_Recovery_Rate"=Recovery_rate_val_U, "Set2_Recovery_Rate"=Recovery_rate_val_P,"time_per_iter_vec"=param[[7]], "total_iterations"=param[[8]])
      }
    }else{output_list<-list("parameters_converged"=param_converge,"Set1_estimated_cluster_IDs"=cluster_ids_est_U,"Set2_estimated_cluster_IDs"=cluster_ids_est_P,"ICL"=ICL_val,"time_per_iter_vec"=param[[7]], "total_iterations"=param[[8]])}
    #}
  }
  return(output_list)
}

simulate_network_bipartite_model7<-function(N1,N2,pi_U,pi_P,theta,delta_0,K1,K2,R){
  sourceCpp(file = paste0(file_path,"model_7.cpp"))
  #library(Matrix)
  net_adjacency<-matrix(0,N1,N2)
  net_rating<-matrix(0,N1,N2)
  #net_adjacency<-Matrix(0,N1,N2,sparse = T)
  #net_rating<-Matrix(0,N1,N2,sparse = T)
  if((K1==1)&(K2==1)){
    clust_id_U_sampling<-rmultinom(n = N1, size = 1, prob = pi_U)
    clust_id_P_sampling<-rmultinom(n = N2, size = 1, prob = pi_P)
    clust_ids_U<-apply(X = clust_id_U_sampling,MARGIN = 2,FUN = function(x){
      y<-which(x==1)
      return(y)
    })
    clust_ids_P<-apply(X = clust_id_P_sampling,MARGIN = 2,FUN = function(x){
      y<-which(x==1)
      return(y)
    })
    for (i in 1:N1){
      for (j in 1:N2){
        exp_val<-exp(theta)
        probab<-exp_val/(1+exp_val)
        edge_sample<-rbinom(n = 1, size = 1, prob = probab)
        if(edge_sample==0){
          network[i,j]<-0
        }else if(edge_sample==1){
          G_vec<-rep(NA_real_,R)
          for (r in 1:R){
            G_vec[r]<-exp(mu[r]*omega-(mu[r])^2/2)
          }
          probab_discrete<-G_vec/sum(G_vec)
          network[i,j] <- which(as.vector(rmultinom(n = 1, size = 1, prob = probab_discrete))==1)
        }
      }
    }
  }
  else{
    clust_id_U_sampling<-rmultinom(n = N1, size = 1, prob = pi_U)
    clust_id_P_sampling<-rmultinom(n = N2, size = 1, prob = pi_P)
    clust_ids_U<-apply(X = clust_id_U_sampling,MARGIN = 2,FUN = function(x){
      y<-which(x==1)
      return(y)
    })
    clust_ids_P<-apply(X = clust_id_P_sampling,MARGIN = 2,FUN = function(x){
      y<-which(x==1)
      return(y)
    })
    
    logitInv_delta_0<-logitInv_delta_0_cal(delta_0=delta_0, K1=K1, K2=K2, R=R)
    probab_block_ratings<-array(NA_real_,dim=c(K1,K2,R))
    for(k in 1:K1){
      for(l in 1:K2){
        for(r in 1:R){
          if(r==1){
            probab_block_ratings[k,l,r]<-logitInv_delta_0[k,l,r]
          }else if(r==R){
            probab_block_ratings[k,l,r]<-1-logitInv_delta_0[k,l,(R-1)]
          }else{
            probab_block_ratings[k,l,r]<-logitInv_delta_0[k,l,r]-logitInv_delta_0[k,l,(r-1)]
          }
        }
      }
    }
    
    for (i in 1:N1){
      for (j in 1:N2){
        exp_val<-exp(theta[clust_ids_U[i],clust_ids_P[j]])
        probab<-exp_val/(1+exp_val)
        edge_sample<-rbinom(n = 1, size = 1, prob = probab)
        if(edge_sample==1){
          net_adjacency[i,j]<-1
          net_rating[i,j] <- which(as.vector(rmultinom(n = 1, size = 1, prob = probab_block_ratings[clust_ids_U[i],clust_ids_P[j],]))==1)
        }
      }
    }
  }
  return(list(net_adjacency,net_rating,clust_ids_U,clust_ids_P))
}

########################################################################################################
########################################################################################################
## Model 8 features (works):
## 1) Block structure over intercept
## 2) Proportional odds model for ratings based on my own derivations
## 3) Segmentation level (K*1) Covariate parameters in outer edge part

wrapper_bipartite_model8<-function(net_adjacency,net_rating,nclust1,nclust2,thres=10^(-6),n_iter_min=50,n_iter_max=500,theta_init,beta_u_init,beta_p_init,delta_0_init,cov_u_outer,cov_p_outer,sim_indicator,theta_true=NA,beta_u_true=NA,beta_p_true=NA, delta_0_true=NA,K_true_U=NA,K_true_P=NA,cluster_ids_true_U=NA,cluster_ids_true_P=NA,R=NA){
  ## This is the combined final updated and working code for EM undirected case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(quadprog)
  library(combinat)
  library(MASS)
  library(Rcpp)
  library(Matrix)
  library(gtools)
  #sourceCpp(file = paste0(file_path,"model_8.cpp"))
  library(biSBMcovNoSEM)
  #################################################
  ## Defining a function to update variational parameters gamma using quadratic program solver. Input:
  ## gamma.curr is a N*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
  ## estimates, theta.curr is a T_grid*K array for current theta estimates. Given Network is assumed to be an N*N*T_data array
  solve_QP_wrapper<-function(quad_lin_coeff,constraint_matrix,constraint_vector,node_ID,gamma.curr){
    try_QP<-try(gamma_next_vec<-solve.QP((-2*diag(as.vector(quad_lin_coeff[node_ID,,1]))),as.vector(quad_lin_coeff[node_ID,,2]),t(constraint_matrix),constraint_vector)$solution,silent = TRUE)
    if('numeric' %in% class(try_QP)){
      gamma_next_vec<-try_QP
    }
    else{
      gamma_next_vec<-gamma.curr[node_ID,]
      print("red flag")
      print(paste0("Node ID for red flag is ", node_ID))
    }
    return(gamma_next_vec)
  }
  
  gamma_U.update.wrapper<-function(gamma_U.curr,gamma_P.curr,pi_U.curr,theta.curr, logitInv_delta_0.curr, net_adjacency, net_rating, cov_beta_u, cov_beta_p, N1,N2,K1,K2,R){
    gamma_U.next<-matrix(NA_real_,N1,K1)
    constraint_matrix<-matrix(NA_real_,2+K1,K1)
    constraint_matrix[1,]<-rep(1,K1)
    constraint_matrix[2,]<-rep(-1,K1)
    constraint_matrix[3:(K1+2),]<-diag(K1)
    constraint_vector<-c(1,-1,rep(0,K1))
    
    quad_lin_coeff<-gamma_U_update_stat_undir(gamma_U=gamma_U.curr, gamma_P=gamma_P.curr, pi_U=pi_U.curr, theta=theta.curr, logitInv_delta_0=logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating=net_rating, cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p,  N1=N1, N2=N2, K1=K1, K2=K2, R=R)
    
    #print(quad_lin_coeff)
    for (i in 1:N1){
      if(sum(is.nan(quad_lin_coeff[i,,1]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,1])),1]<--Inf
      }
      if(sum(is.nan(quad_lin_coeff[i,,2]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,2])),2]<-Inf
      }
      if(sum(quad_lin_coeff[i,,1]<=-10^8)==0){
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma_U.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_U.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K1)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma_U.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_U.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K1){
          gamma_U.next[i,]<-gamma.curr[i,]
        }
      }else if((sum(quad_lin_coeff[i,,1]<=-10^8)>0)&&(sum(quad_lin_coeff[i,,1]<=-10^8)<K1)){
        quad_lin_coeff[i,which(quad_lin_coeff[i,,1]<=-10^8),1]<--(10^10)
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma_U.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_U.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K1)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma_U.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_U.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K1){
          gamma_U.next[i,]<-gamma_U.curr[i,]
        }
      }else if(sum(quad_lin_coeff[i,,1]<=-10^8)==K1){
        gamma_U.next[i,]<-gamma_U.curr[i,]
      }
      #print(i)
    }
    #print('gammaU')
    #print(quad_lin_coeff)
    # for (i in 1:N1){
    #   gamma_U.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    #   print(i)
    # }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma_U.next<-t(apply(X = gamma_U.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma_U.next)
  }
  
  gamma_P.update.wrapper<-function(gamma_U.curr,gamma_P.curr,pi_P.curr,theta.curr, logitInv_delta_0.curr, net_adjacency, net_rating, cov_beta_u, cov_beta_p, N1,N2,K1,K2,R){
    gamma_P.next<-matrix(NA_real_,N2,K2)
    constraint_matrix<-matrix(NA_real_,2+K2,K2)
    constraint_matrix[1,]<-rep(1,K2)
    constraint_matrix[2,]<-rep(-1,K2)
    constraint_matrix[3:(K2+2),]<-diag(K2)
    constraint_vector<-c(1,-1,rep(0,K2))
    
    quad_lin_coeff<-gamma_P_update_stat_undir(gamma_U=gamma_U.curr, gamma_P=gamma_P.curr, pi_P=pi_P.curr, theta=theta.curr, logitInv_delta_0=logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating=net_rating, cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, N1=N1, N2=N2, K1=K1, K2=K2, R=R)
    #print('gammaP')
    #print(quad_lin_coeff)
    
    #print(quad_lin_coeff)
    for (i in 1:N2){
      if(sum(is.nan(quad_lin_coeff[i,,1]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,1])),1]<--Inf
      }
      if(sum(is.nan(quad_lin_coeff[i,,2]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,2])),2]<-Inf
      }
      if(sum(quad_lin_coeff[i,,1]<=-10^8)==0){
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma_P.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_P.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K2)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma_P.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_P.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K2){
          gamma_U.next[i,]<-gamma.curr[i,]
        }
      }else if((sum(quad_lin_coeff[i,,1]<=-10^8)>0)&&(sum(quad_lin_coeff[i,,1]<=-10^8)<K2)){
        quad_lin_coeff[i,which(quad_lin_coeff[i,,1]<=-10^8),1]<--(10^10)
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma_P.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_P.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K2)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma_P.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_P.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K2){
          gamma_P.next[i,]<-gamma_P.curr[i,]
        }
      }else if(sum(quad_lin_coeff[i,,1]<=-10^8)==K2){
        gamma_P.next[i,]<-gamma_P.curr[i,]
      }
      #print(i)
    }
    
    # for (i in 1:N2){
    #   gamma_P.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    #   print(i)
    # }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma_P.next<-t(apply(X = gamma_P.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma_P.next)
  }
  
  #################################################
  ## Defining a function to update variational parameters gamma using quadratic program solver. Input: 
  ## gamma.curr is a N*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
  ## estimates, theta.curr is a T_grid*K array for current theta estimates. Given Network is assumed to be an N*N*T_data array
  # gamma_U.update.wrapper<-function(gamma_U.curr,gamma_P.curr,pi_U.curr,theta.curr,logitInv_delta_0.curr,net_adjacency,net_rating,N1,N2,K1,K2,R){
  #   gamma_U.next<-matrix(NA_real_,N1,K1)
  #   constraint_matrix<-matrix(NA_real_,2+K1,K1)
  #   constraint_matrix[1,]<-rep(1,K1)
  #   constraint_matrix[2,]<-rep(-1,K1)
  #   constraint_matrix[3:(K1+2),]<-diag(K1)
  #   constraint_vector<-c(1,-1,rep(0,K1))
  #   
  #   quad_lin_coeff<-gamma_U_update_stat_undir(gamma_U=gamma_U.curr, gamma_P=gamma_P.curr, pi_U=pi_U.curr, theta=theta.curr, beta_u=beta_u.curr, beta_p=beta_p.curr, logitInv_delta_0=logitInv_delta_0.curr,net_adjacency=net_adjacency, net_rating = net_rating, N1=N1, N2=N2, K1=K1, K2=K2, R=R)
  #   #print('gammaU')
  #   #print(quad_lin_coeff)
  #   for (i in 1:N1){
  #     gamma_U.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
  #   }
  #   
  #   ## normalizing gamma_i. deal with case outside (0,1) later
  #   gamma_U.next<-t(apply(X = gamma_U.next,MARGIN = 1,FUN = function(x){
  #     x_norm<-x/sum(x)
  #     return(x_norm)
  #   }))
  #   return(gamma_U.next)
  # }
  # 
  # gamma_P.update.wrapper<-function(gamma_U.curr,gamma_P.curr,pi_P.curr,theta.curr,beta_u.curr, beta_p.curr,logitInv_delta_0.curr,net_adjacency,net_rating,N1,N2,K1,K2,R){
  #   gamma_P.next<-matrix(NA_real_,N2,K2)
  #   constraint_matrix<-matrix(NA_real_,2+K2,K2)
  #   constraint_matrix[1,]<-rep(1,K2)
  #   constraint_matrix[2,]<-rep(-1,K2)
  #   constraint_matrix[3:(K2+2),]<-diag(K2)
  #   constraint_vector<-c(1,-1,rep(0,K2))
  #   
  #   quad_lin_coeff<-gamma_P_update_stat_undir(gamma_U=gamma_U.curr, gamma_P=gamma_P.curr, pi_P=pi_P.curr, theta=theta.curr, beta_u=beta_u.curr, beta_p=beta_p.curr, logitInv_delta_0=logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating = net_rating, N1=N1, N2=N2, K1=K1, K2=K2, R=R)
  #   #print('gammaP')
  #   #print(quad_lin_coeff)
  #   for (i in 1:N2){
  #     gamma_P.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
  #   }
  #   
  #   ## normalizing gamma_i. deal with case outside (0,1) later
  #   gamma_P.next<-t(apply(X = gamma_P.next,MARGIN = 1,FUN = function(x){
  #     x_norm<-x/sum(x)
  #     return(x_norm)
  #   }))
  #   return(gamma_P.next)
  # }
  
  #################################################
  ## Defining a function to update K*1 vector of pi
  pi.update<-function(gamma.curr){
    pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## normalization of pi
    pi.next<-pi.next/sum(pi.next)
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions K1,K2
  theta.update<-function(theta.curr,gamma_U,gamma_P,net_adjacency, cov_beta_u, cov_beta_p, N1,N2,K1,K2){
    theta.next<-matrix(NA_real_,K1,K2)
    gradient<-grad_bipartite_stat_undir_theta(theta=theta.curr, gamma_U=gamma_U, gamma_P=gamma_P, net_adjacency=net_adjacency,  cov_beta_u = cov_beta_u, cov_beta_p = cov_beta_p, N1=N1, N2=N2, K1=K1, K2=K2)
    hess<-hess_bipartite_stat_undir_theta(theta=theta.curr, gamma_U=gamma_U, gamma_P=gamma_P, cov_beta_u = cov_beta_u, cov_beta_p = cov_beta_p, N1=N1, N2=N2, K1=K1, K2=K2)
    for(k in 1:K1){
      for(l in 1:K2){
        theta.next[k,l]<-theta.curr[k,l]-gradient[k,l]/hess[k,l]
      }
    }
    return(theta.next)
  }
  
  #################################################
  ## Defining the update function of full vector beta_u with dimensions K1*C1_outer
  beta_u.update<-function(beta_u.curr,gamma_U,gamma_P,theta,net_adjacency, cov_beta_u, cov_beta_p, cov_u_outer, N1,N2,K1,K2,C1_outer){
    beta_u.next<-array(NA_real_,dim=c(K1,C1_outer))
    gradient<-grad_bipartite_stat_undir_beta_u(gamma_U=gamma_U, gamma_P=gamma_P, theta=theta, net_adjacency=net_adjacency, cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, cov_u_outer = cov_u_outer,  N1=N1, N2=N2, K1=K1, K2=K2, C1_outer=C1_outer)
    hess<-hess_bipartite_stat_undir_beta_u(gamma_U=gamma_U, gamma_P=gamma_P, theta=theta, cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, cov_u_outer = cov_u_outer, N1=N1, N2=N2, K1=K1, K2=K2, C1_outer=C1_outer)
    for(c_u in 1:C1_outer){
      for(k in 1:K1){
        beta_u.next[k,c_u]<-beta_u.curr[k,c_u]-(gradient[k,c_u]/hess[k,c_u])
      }
    }
    return(beta_u.next)
  }
  
  #################################################
  ## Defining the update function of full vector beta_p with dimensions K2*C2_outer
  beta_p.update<-function(beta_p.curr,gamma_U,gamma_P,theta,net_adjacency, cov_beta_u, cov_beta_p, cov_p_outer, N1,N2,K1,K2,C2_outer){
    beta_p.next<-array(NA_real_,dim=c(K2,C2_outer))
    gradient<-grad_bipartite_stat_undir_beta_p(gamma_U=gamma_U, gamma_P=gamma_P, theta=theta, net_adjacency=net_adjacency, cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, cov_p_outer = cov_p_outer, N1=N1, N2=N2, K1=K1, K2=K2, C2_outer=C2_outer)
    hess<-hess_bipartite_stat_undir_beta_p(gamma_U=gamma_U, gamma_P=gamma_P, theta=theta, cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, cov_p_outer = cov_p_outer, N1=N1, N2=N2, K1=K1, K2=K2, C2_outer=C2_outer)
    for(c_p in 1:C2_outer){
      for(l in 1:K2){
        beta_p.next[l,c_p]<-beta_p.curr[l,c_p]-(gradient[l,c_p]/hess[l,c_p])
      }
    }
    return(beta_p.next)
  }
  
  #################################################
  ## Defining the update function of full proportional odds nuisance parameters delta_0
  delta_0.update<-function(gamma_U,gamma_P,delta_0,logitInv_delta_0,logitInv_phi,net_adjacency,net_rating_statistic,N1,N2,K1,K2,R){
    delta_0.next<-array(NA_real_,dim=c(K1,K2,(R-1)))
    
    ## Getting the estimated cluster memberships of all nodes in two sets from variational E step
    #clust_est_U<-as.vector(apply(X = gamma_U,MARGIN = 1,FUN = which.max))
    #clust_est_P<-as.vector(apply(X = gamma_P,MARGIN = 1,FUN = which.max))
    
    #ratings_list<-tie_clust_partition(clust_est_U=clust_est_U, clust_est_P=clust_est_P, net_adjacency=net_adjacency, net_rating=net_rating, N1=N1, N2=N2, K1=K1, K2=K2)
    
    gradient_delta_0<-grad_bipartite_stat_undir_delta_0(gamma_U=gamma_U, gamma_P=gamma_P, delta_0=delta_0, logitInv_delta_0=logitInv_delta_0, net_adjacency=net_adjacency, net_rating=net_rating, N1=N1, N2=N2, K1=K1, K2=K2, R=R)
    
    for (k in 1:K1){
      for (l in 1:K2){
        hess_mat<-hess_bipartite_stat_undir_delta_0(k=k-1, l=l-1, gamma_U=gamma_U, gamma_P=gamma_P, delta_0=delta_0, logitInv_delta_0=logitInv_delta_0, net_adjacency=net_adjacency, net_rating=net_rating, N1=N1, N2=N2, R=R)
        delta_0.next[k,l,]<-sort(delta_0[k,l,]-as.vector((solve(hess_mat+(10^(-6)*diag(R-1)))%*%gradient_delta_0[k,l,])))
        #delta_0.next[k,l,]<-sort(delta_0[k,l,]-as.vector((solve(hess_mat)%*%gradient_delta_0[k,l,])))
      }
    }
    #print(delta_0.next)
    return(delta_0.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,net_adjacency,net_rating,K1,K2,n_iter,thres,n_iter_min,n_iter_max,R,cov_u_outer,cov_p_outer){
    ## Assuming a network adjacency matrix of dimensions N1*N2 with rows corresponding to users and columns corresponding to products and each entry is an indicator whether the user rated the product. This framework can be generalized.
    N1<-dim(net_adjacency)[1]
    N2<-dim(net_adjacency)[2]
    C1_outer<-dim(cov_u_outer)[2]
    C2_outer<-dim(cov_p_outer)[2]
    
    ## initializing the arrays for parameters
    ## Initializing the cluster memberships gamma
    gamma_U<-array(NA_real_,dim=c(N1,K1,n_iter)) ## assuming there are K1 clusters in user set
    gamma_P<-array(NA_real_,dim=c(N2,K2,n_iter)) ## assuming there are K2 clusters in product set
    
    ## Initializing the mixture proportions pi
    pi_U<-matrix(NA_real_,K1,n_iter)
    pi_P<-matrix(NA_real_,K2,n_iter)
    
    ## Initializing the network parameters theta
    theta<-array(NA_real_,dim=c(K1,K2,n_iter))
    
    ## Initializing the covariate parameters beta_u and beta_p
    beta_u<-array(NA_real_,dim=c(K1,C1_outer,n_iter))
    beta_p<-array(NA_real_,dim=c(K2,C2_outer,n_iter))
    
    ## Initializing the proportional odds model intercept parameters
    delta_0<-array(NA_real_,dim=c(K1,K2,R-1,n_iter))
    
    gamma_U[,,1]<-start[[1]]
    gamma_P[,,1]<-start[[2]]
    pi_U[,1]<-start[[3]]
    pi_P[,1]<-start[[4]]
    theta[,,1]<-start[[5]]
    beta_u[,,1]<-start[[6]]
    beta_p[,,1]<-start[[7]]
    delta_0[,,,1]<-start[[8]]
    
    #print(theta[,,1])
    #print(omega[,1])
    #print(mu[,,,1])
    
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_curr<-10^10
    time_per_iter<-c()
    while(((error>thres)&(iter_index<n_iter_max))|(iter_index<n_iter_min)){ 
      #while((error>thres|iter_index<200)&(iter_index<250)&(iter_index==2|error<10)){ 
      #while(iter_index<21){ 
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      
      ## Calculating cov_beta_u and cov_beta_p
      cov_beta_u<-cov_beta_mul(cov=cov_u_outer,beta=t(beta_u[,,iter_index-1]))
      cov_beta_p<-cov_beta_mul(cov=cov_p_outer,beta=t(beta_p[,,iter_index-1]))
      
      ## Calculating the logit inverse of delta_0
      logitInv_delta_0.curr<-logitInv_delta_0_cal(delta_0=delta_0[,,,iter_index-1], K1=K1, K2=K2, R=R)
      
      ## Transforming logit inverse of delta_0 to phi
      #phi.curr<-phi_cal(logitInv_delta_0=logitInv_delta_0.curr, K1=K1, K2=K2, R=R)
      
      ## Calculating the logit inverse of phi
      #logitInv_phi.curr<-logitInv_phi_cal(phi=phi.curr, K1=K1, K2=K2, R=R)
      
      ## Updating gamma_U, the mixed memberships for users
      gamma_U[,,iter_index]<-gamma_U.update.wrapper(gamma_U.curr=gamma_U[,,iter_index-1], gamma_P.curr=gamma_P[,,iter_index-1], pi_U.curr=pi_U[,iter_index-1], theta.curr=theta[,,iter_index-1], cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, logitInv_delta_0.curr=logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating=net_rating, N1=N1,N2=N2,K1=K1,K2=K2,R=R)
      ## Updating the mixture proportions pi_U
      pi_U[,iter_index]<-pi.update(gamma.curr=gamma_U[,,iter_index])
      
      ## Updating gamma_P, the mixed memberships for products
      gamma_P[,,iter_index]<-gamma_P.update.wrapper(gamma_U.curr=gamma_U[,,iter_index],gamma_P.curr=gamma_P[,,iter_index-1],pi_P.curr=pi_P[,iter_index-1], theta.curr=theta[,,iter_index-1], cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, logitInv_delta_0.curr=logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating=net_rating, N1=N1,N2=N2,K1=K1,K2=K2,R=R)
      ## Updating the mixture proportions pi_P
      pi_P[,iter_index]<-pi.update(gamma.curr=gamma_P[,,iter_index])
     
      ## Updating the network parameters theta
      theta[,,iter_index]<-theta.update(theta.curr=theta[,,iter_index-1],gamma_U=gamma_U[,,iter_index], gamma_P=gamma_P[,,iter_index], net_adjacency=net_adjacency, cov_beta_u = cov_beta_u, cov_beta_p = cov_beta_p, N1=N1, N2=N2, K1=K1, K2=K2)
      #print('theta')
      #print(theta[,,iter_index])
      
      ## Updating the covariate parameters beta_u
      beta_u[,,iter_index]<-beta_u.update(beta_u.curr = matrix(beta_u[,,iter_index-1],K1,C1_outer),gamma_U=gamma_U[,,iter_index], gamma_P=gamma_P[,,iter_index], theta=theta[,,iter_index], net_adjacency=net_adjacency, cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, cov_u_outer = cov_u_outer, N1=N1, N2=N2, K1=K1, K2=K2, C1_outer= C1_outer)
      #print('beta_u and beta_p')
      #print(beta_u[,,iter_index])
      
      ## Updating the covariate parameters beta_p
      beta_p[,,iter_index]<-beta_p.update(beta_p.curr = matrix(beta_p[,,iter_index-1],K2,C2_outer),gamma_U=gamma_U[,,iter_index], gamma_P=gamma_P[,,iter_index], theta=theta[,,iter_index], net_adjacency=net_adjacency, cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, cov_p_outer = cov_p_outer, N1=N1, N2=N2, K1=K1, K2=K2, C2_outer=C2_outer)
      #print(beta_p[,,iter_index])
      
      ## Updating the omega matrix
      delta_0[,,,iter_index]<-delta_0.update(gamma_U=gamma_U[,,iter_index], gamma_P=gamma_P[,,iter_index], delta_0=delta_0[,,,iter_index-1],logitInv_delta_0 = logitInv_delta_0.curr,logitInv_phi = logitInv_phi.curr, net_adjacency=net_adjacency, net_rating_statistic = net_rating_statistic, N1=N1, N2=N2, K1=K1, K2=K2, R=R)
      #print('delta_0')
      #print(delta_0[,,,iter_index])
      #Sys.sleep(1.5)
      
      ELBO_prev<-ELBO_curr
      ELBO_curr<-ELBO_conv_bipartite_stat_undir(gamma_U=gamma_U[,,iter_index], gamma_P=gamma_P[,,iter_index], pi_U=pi_U[,iter_index], pi_P=pi_P[,iter_index], theta = theta[,,iter_index], logitInv_delta_0=logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating=net_rating, cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, N1=N1, N2=N2, K1=K1, K2=K2)
      
      error<-abs(((ELBO_prev-ELBO_curr)/ELBO_prev))
      
      if((iter_index%%1)==0){
        #print(theta[,,iter_index])
        print('iter_index')
        print(iter_index)
        print('error')
        print(error)
        print(proc.time()-ptm)
      }
      ptm_diff<-proc.time()-ptm
      time_per_iter<-c(time_per_iter,ptm_diff[3])
      
      #print("probs")
      #print(Gvals_cube_0)
      # print(apply(X = Gvals_cube_0,MARGIN = c(1,2),FUN = function(x){
      #   
      #   return(x/sum(x))
      # }))
      print('iter_index')
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    names(time_per_iter)<-NULL
    total_iterations<-iter_index-1
    return(list(gamma_U,gamma_P,pi_U,pi_P,theta,beta_u,beta_p,delta_0,time_per_iter,total_iterations))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining Model Selection functions based on converged paramters
  
  #################################################
  ## Defining a function to calculate the estimate of complete log likelihood
  comp_loglik_cal<-function(cluster_ids_est_U = NA, cluster_ids_est_P = NA, pi_U = NA, pi_P = NA, theta, beta_u, beta_p, delta_0, net_adjacency, net_rating, cov_u_outer, cov_p_outer, N1, N2, K1, K2, R){
    if((K1!=1)&(K2!=1)){
      logitInv_delta_0<-logitInv_delta_0_cal(delta_0=delta_0, K1=K1, K2=K2, R=R)
      probab_block_ratings<-array(NA_real_,dim=c(K1,K2,R))
      for(k in 1:K1){
        for(l in 1:K2){
          for(r in 1:R){
            if(r==1){
              probab_block_ratings[k,l,r]<-logitInv_delta_0[k,l,r]
            }else if(r==R){
              probab_block_ratings[k,l,r]<-1-logitInv_delta_0[k,l,(R-1)]
            }else{
              probab_block_ratings[k,l,r]<-logitInv_delta_0[k,l,r]-logitInv_delta_0[k,l,(r-1)]
            }
          }
        }
      }
      
      ## 1st term
      t1<-0
      for(i in 1:N1){
        for(j in 1:N2){
          cluster_id_i<-cluster_ids_est_U[i]
          cluster_id_j<-cluster_ids_est_P[j]
          predictor<-theta[cluster_id_i,cluster_id_j]+sum(cov_u_outer[i,]*beta_u[cluster_id_i,])+sum(cov_p_outer[j,]*beta_p[cluster_id_j,])
          exp_val<-exp(predictor)
          if(net_adjacency[i,j]==0){
            t1<-t1-(log(1+exp_val))
          }else{
            for(r in 1:R){
              if(net_rating[i,j]==r){
                t1<-t1+(predictor-log(1+exp_val)+log(probab_block_ratings[cluster_id_i,cluster_id_j,r]))
              }
            }
          }
        }
      }
      ## 2nd term
      t2<-0
      for (i in 1:N1){
        t2<-t2+log(pi_U[cluster_ids_est_U[i]])
      }
      t3<-0
      for (j in 1:N2){
        t3<-t3+log(pi_P[cluster_ids_est_P[j]])
      }
      
      comp_val<-t1+t2+t3
    }
    return(comp_val)
  }
  
  #################################################
  ## Defining a function to calculate the integrated classification likelihood
  ICL_cal<-function(cluster_ids_est_U, cluster_ids_est_P, pi_U, pi_P, theta, beta_u, beta_p, delta_0, net_adjacency, net_rating, cov_u_outer, cov_p_outer, N1, N2, K1, K2, R){
    comp_loglik_val<-comp_loglik_cal(cluster_ids_est_U = cluster_ids_est_U, cluster_ids_est_P = cluster_ids_est_P, pi_U = pi_U, pi_P = pi_P, theta = theta, beta_u = beta_u, beta_p = beta_p, delta_0 = delta_0, cov_u_outer = cov_u_outer, cov_p_outer = cov_p_outer, net_adjacency = net_adjacency, net_rating = net_rating, N1 = N1, N2 = N2, K1 = K1, K2 = K2, R=R)
    ## Penalization over pi
    t2<-(K1-1)*log(N1)+(K2-1)*log(N2)
    ## Penalization over theta
    t3<-(K1*K2)*log(N1*N2)
    ## Penalization over beta_u
    C1_outer<-dim(cov_u_outer)[2]
    t4<-(K1*C1_outer)*log(N1*N2)
    ## Penalization over beta_p
    C2_outer<-dim(cov_p_outer)[2]
    t5<-(K2*C2_outer)*log(N1*N2)
    ## Calculating the number of edges
    N0_val<-N0_cal(net_adjacency=net_adjacency, N1=N1, N2=N2)
    ## Penalization over omega, mu
    t6<-(K1*K2*(R-1))*log(N0_val)
    
    ## all penalty terms together:
    penalty<-t2+t3+t4+t5+t6
    ICL_val<-comp_loglik_val-penalty
    ICL_vec<-c(comp_loglik_val,penalty,ICL_val)
    return(ICL_vec)
  }
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI_cal<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  #################################################
  ## RASE functions
  RASE_cal<-function(param_est, param_true){
    if(is.vector(param_est)!=1){
      RASE_val<-sqrt(sum((param_est-param_true)^2)/prod(dim(param_est)))
    }else{RASE_val<-sqrt(sum((param_est-param_true)^2)/length(param_est))}
    return(RASE_val)
  }
  
  ########################################################################################################
  ## Defining the functions to evaluate clustering performance using normalized mutual information.
  probab_clust_cal<-function(cluster_ids,N,K){
    probab_vec<-rep(NA_real_,K)
    for(k in 1:K){
      probab_vec[k]<-sum(cluster_ids==k)/N
    }
    return(probab_vec)
  }
  
  ########################################################################################################
  ## Normalized Mutual Info
  NMI_cal<-function(cluster_ids_true_U,cluster_ids_true_P,K_true_U,K_true_P,cluster_ids_est_U,cluster_ids_est_P,K1,K2,N1,N2){
    probab_vec_true_U<-probab_clust_cal(cluster_ids = cluster_ids_true_U,N = N1+N2,K = K_true_U)
    probab_vec_true_P<-probab_clust_cal(cluster_ids = cluster_ids_true_P,N = N1+N2,K = K_true_P)
    probab_vec_true<-c(probab_vec_true_U,probab_vec_true_P)
    total_entropy_true<-0
    for(k in 1:length(probab_vec_true)){
      if(probab_vec_true[k]!=0){
        total_entropy_true<-total_entropy_true-(probab_vec_true[k]*log(probab_vec_true[k],base = 2))
      }
    }
    
    probab_vec_est_U<-probab_clust_cal(cluster_ids = cluster_ids_est_U,N = N1+N2,K = K1)
    probab_vec_est_P<-probab_clust_cal(cluster_ids = cluster_ids_est_P,N = N1+N2,K = K2)
    probab_vec_est<-c(probab_vec_est_U,probab_vec_est_P)
    total_entropy_est<-0
    for(k in 1:length(probab_vec_est)){
      if(probab_vec_est[k]!=0){
        total_entropy_est<-total_entropy_est-(probab_vec_est[k]*log(probab_vec_est[k],base = 2))
      }
    }
    
    clust_ids_true_total<-c(cluster_ids_true_U,(cluster_ids_true_P+K_true_U))
    clust_ids_est_total<-c(cluster_ids_est_U,(cluster_ids_est_P+K1))
    
    conditional_entropy<-rep(NA_real_,(K1+K2))
    for(k in 1:(K1+K2)){
      probab_conditional_vec<-rep(0,(K_true_U+K_true_P))
      for(m in 1:(K_true_U+K_true_P)){
        if(sum(clust_ids_est_total==k)!=0){
          probab_conditional_vec[m]<-sum((clust_ids_true_total==m)&(clust_ids_est_total==k))/sum(clust_ids_est_total==k)
        }
      }
      plogp_conditional_vec<-rep(0,length(probab_conditional_vec))
      for(m in 1:length(probab_conditional_vec)){
        if(probab_conditional_vec[m]!=0){
          plogp_conditional_vec[m]<-probab_conditional_vec[m]*log(x = probab_conditional_vec[m],base = 2)
        }
      }
      
      conditional_entropy[k]<--probab_vec_est[k]*sum(plogp_conditional_vec)
    }
    
    mutual_info<-total_entropy_true-sum(conditional_entropy)
    NMI_val<-(2*mutual_info)/(total_entropy_true+total_entropy_est)
    return(NMI_val)
  }
  
  ########################################################################################################
  ## Calculating the recovery rate
  recovery_rate_cal <- function(cluster_ids_est, cluster_ids_true){
    all_perm<-permutations(n=4,r=2,repeats.allowed = T)
    clust_perm_mat<-matrix(NA_integer_,4,4)
    for(i in 1:nrow(all_perm)){
      clust_perm_mat[all_perm[i,1],all_perm[i,2]] <- sum((cluster_ids_true==all_perm[i,1])&(cluster_ids_est==all_perm[i,2]))
    }
    recovery_rate_val<-sum(apply(X = clust_perm_mat, MARGIN = 1, FUN = max))/length(cluster_ids_est)
    return(recovery_rate_val)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N1<-dim(net_adjacency)[1] ## Number of nodes from the network
  N2<-dim(net_adjacency)[2]
  ## Total time points in the network for which we have data in the form of adjacency matrix
  K1<-nclust1 ## Defining the number of clusters
  K2<-nclust2
  K <-K1*K2
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  start[[1]]<-matrix(rep(1/K1,N1*K1),N1,K1)
  start[[2]]<-matrix(rep(1/K2,N2*K2),N2,K2)
  start[[3]]<-rep(1/K1,K1)
  start[[4]]<-rep(1/K2,K2)
  start[[5]]<-theta_init
  start[[6]]<-beta_u_init
  start[[7]]<-beta_p_init
  start[[8]]<-delta_0_init
  if(K1!=1&K2!=1){ param<-iterator(start=start, net_adjacency = net_adjacency, net_rating=net_rating, K1=K1, K2=K2, n_iter=1000, thres=thres, n_iter_min=n_iter_min, n_iter_max=n_iter_max, R=R, cov_u_outer=cov_u_outer, cov_p_outer=cov_p_outer)}
  if(K1==1&K2!=1){ param<-iterator_K1(start=start, network=sim.net, K1=K1, K2=K2, n_iter=1000, thres=thres,R=R)}
  if(K1!=1&K2==1){ param<-iterator_K2(start=start, network=sim.net, K1=K1, K2=K2, n_iter=1000, thres=thres,R=R)}
  if(K1==1&K2==1){ param<-iterator_K1K2(start=start, network=sim.net, K1=K1, K2=K2, n_iter=1000, thres=thres,R=R)}
  
  #################################################
  ## extracting the coverged parameter values and calculating ICL
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    temp<-is.na(param[[1]][1,1,n_iter])
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  #print(param)
  
  if(K1!=1&K2!=1){
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,,n_last]
    param_converge[[3]]<-param[[3]][,n_last]
    param_converge[[4]]<-param[[4]][,n_last]
    param_converge[[5]]<-param[[5]][,,n_last]
    param_converge[[6]]<-param[[6]][,,n_last]
    param_converge[[7]]<-param[[7]][,,n_last]
    param_converge[[8]]<-param[[8]][,,,n_last]
    output_list<-param_converge
    cluster_ids_est_U<-as.vector(apply(X = matrix(1:N1),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[1]][x,])
      return(cluster_id)
    }))
    cluster_ids_est_P<-as.vector(apply(X = matrix(1:N2),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[2]][x,])
      return(cluster_id)
    }))
    ICL_vec<-NA
    #ICL_vec<-ICL_cal(cluster_ids_est_U = cluster_ids_est_U, cluster_ids_est_P = cluster_ids_est_P, pi_U=param_converge[[3]], pi_P=param_converge[[4]] ,theta = param_converge[[5]], beta_u = as.matrix(param_converge[[6]]), beta_p = as.matrix(param_converge[[7]]), delta_0 = param_converge[[8]], net_adjacency=net_adjacency, net_rating=net_rating, cov_u_outer = cov_u_outer, cov_p_outer = cov_p_outer, N1 = N1, N2 = N2, K1 = K1, K2 = K2, R = R)
    if(sim_indicator==1){
      RI_val_U<-RI_cal(cluster_ids_est = cluster_ids_est_U,cluster_ids_true = cluster_ids_true_U)
      RI_val_P<-RI_cal(cluster_ids_est = cluster_ids_est_P,cluster_ids_true = cluster_ids_true_P)
      NMI_val<-NMI_cal(cluster_ids_true_U = cluster_ids_true_U,cluster_ids_true_P = cluster_ids_true_P,K_true_U = K_true_U,K_true_P = K_true_P,cluster_ids_est_U = cluster_ids_est_U,cluster_ids_est_P = cluster_ids_est_P, K1 = K1, K2 = K2, N1 = N1, N2 = N2)
      Recovery_rate_val_U<-recovery_rate_cal(cluster_ids_est = cluster_ids_est_U,cluster_ids_true = cluster_ids_true_U)
      Recovery_rate_val_P<-recovery_rate_cal(cluster_ids_est = cluster_ids_est_P,cluster_ids_true = cluster_ids_true_P)
      if((K1==K_true_U)&(K2==K_true_P)){
        K1_permute_mat<-do.call(rbind,permn(1:K1))
        K2_permute_mat<-do.call(rbind,permn(1:K2))
        
        RASE_beta_u_vec<-rep(NA_real_,nrow(K1_permute_mat))
        for (k in 1:nrow(K1_permute_mat)){
          beta_u_est<-as.vector(param_converge[[6]])[K1_permute_mat[k,]]
          RASE_beta_u_vec[k]<-RASE_cal(param_est = beta_u_est, param_true = as.vector(beta_u_true))
        }
        permute_true_id_beta_u<-which.min(RASE_beta_u_vec)
        RASE_beta_u_val<-RASE_beta_u_vec[permute_true_id_beta_u]
        
        RASE_beta_p_vec<-rep(NA_real_,nrow(K2_permute_mat))
        for (k in 1:nrow(K2_permute_mat)){
          beta_p_est<-as.vector(param_converge[[7]])[K2_permute_mat[k,]]
          RASE_beta_p_vec[k]<-RASE_cal(param_est = beta_p_est, param_true = as.vector(beta_p_true))
        }
        permute_true_id_beta_p<-which.min(RASE_beta_p_vec)
        RASE_beta_p_val<-RASE_beta_p_vec[permute_true_id_beta_p]
        
        RASE_theta_mat<-matrix(NA_real_,nrow(K1_permute_mat),nrow(K2_permute_mat))
        for (k in 1:nrow(K1_permute_mat)){
          for (l in 1:nrow(K2_permute_mat)){
            theta_true_rotated<-as(as.integer(K1_permute_mat[k,]), "pMatrix")%*%theta_true%*%t(as(as.integer(K2_permute_mat[l,]), "pMatrix"))
            RASE_theta_mat[k,l]<-RASE_cal(param_est = param_converge[[5]], param_true = theta_true_rotated)
          }
        }
        
        RASE_theta_val<-min(RASE_theta_mat)
        rotate_indices<-which(RASE_theta_mat == min(RASE_theta_mat), arr.ind = TRUE)
        
        ## RASE for delta_0
        # RASE_delta_0_vec<-rep(NA_real_,R-1)
        # for (r in 1:(R-1)){
        #   delta_0_true_rotated<-as(as.integer(K1_permute_mat[rotate_indices[1],]), "pMatrix")%*%delta_0_true[,,r]%*%t(as(as.integer(K2_permute_mat[rotate_indices[2],]), "pMatrix"))
        #   RASE_delta_0_vec[r]<-RASE_cal(param_est = param_converge[[8]][,,r], param_true = delta_0_true_rotated)
        # }
        
        delta_0_true_rotated<-array(NA_real_,dim=c(K1,K2,(R-1)))
        for (r in 1:(R-1)){
          delta_0_true_rotated[,,r]<-as(as.integer(K1_permute_mat[rotate_indices[1],]), "pMatrix")%*%delta_0_true[,,r]%*%t(as(as.integer(K2_permute_mat[rotate_indices[2],]), "pMatrix"))
          
        }
        RASE_delta_0_val<-RASE_cal(param_est = param_converge[[8]], param_true = delta_0_true_rotated)
        
        output_list<-list("parameters_converged"=param_converge,"Set1_estimated_cluster_IDs"=cluster_ids_est_U,"Set2_estimated_cluster_IDs"=cluster_ids_est_P,"ICL"=ICL_vec,"Set1_Rand_Index"=RI_val_U,"Set2_Rand_Index"=RI_val_P,"Normalized_Mutual_Info"=NMI_val, "Set1_Recovery_Rate"=Recovery_rate_val_U, "Set2_Recovery_Rate"=Recovery_rate_val_P, "RASE_theta"=RASE_theta_val, "RASE_beta_u"=RASE_beta_u_val, "RASE_beta_p"=RASE_beta_p_val, "RASE_delta_0"=RASE_delta_0_val,"time_per_iter_vec"=param[[9]],"total_iterations"=param[[10]])
      }else{
        output_list<-list("parameters_converged"=param_converge,"Set1_estimated_cluster_IDs"=cluster_ids_est_U,"Set2_estimated_cluster_IDs"=cluster_ids_est_P,"ICL"=ICL_vec,"Set1_Rand_Index"=RI_val_U,"Set2_Rand_Index"=RI_val_P,"Normalized_Mutual_Info"=NMI_val, "Set1_Recovery_Rate"=Recovery_rate_val_U, "Set2_Recovery_Rate"=Recovery_rate_val_P, "time_per_iter_vec"=param[[9]],"total_iterations"=param[[10]])
      }
    }else{output_list<-list("parameters_converged"=param_converge,"Set1_estimated_cluster_IDs"=cluster_ids_est_U,"Set2_estimated_cluster_IDs"=cluster_ids_est_P,"ICL"=ICL_vec,"time_per_iter_vec"=param[[9]],"total_iterations"=param[[10]])}
    #}
  }
  return(output_list)
}

simulate_network_bipartite_model8<-function(N1,N2,pi_U,pi_P,theta,beta_u,beta_p,delta_0,cov_u_outer,cov_p_outer,K1,K2,R){
  library(Rcpp)
  sourceCpp(file = paste0(file_path,"model_8.cpp"))
  #library(Matrix)
  net_adjacency<-matrix(0,N1,N2)
  net_rating<-matrix(0,N1,N2)
  #net_adjacency<-Matrix(0,N1,N2,sparse = T)
  #net_rating<-Matrix(0,N1,N2,sparse = T)
  if((K1==1)&(K2==1)){
    clust_id_U_sampling<-rmultinom(n = N1, size = 1, prob = pi_U)
    clust_id_P_sampling<-rmultinom(n = N2, size = 1, prob = pi_P)
    clust_ids_U<-apply(X = clust_id_U_sampling,MARGIN = 2,FUN = function(x){
      y<-which(x==1)
      return(y)
    })
    clust_ids_P<-apply(X = clust_id_P_sampling,MARGIN = 2,FUN = function(x){
      y<-which(x==1)
      return(y)
    })
    for (i in 1:N1){
      for (j in 1:N2){
        exp_val<-exp(theta)
        probab<-exp_val/(1+exp_val)
        edge_sample<-rbinom(n = 1, size = 1, prob = probab)
        if(edge_sample==0){
          network[i,j]<-0
        }else if(edge_sample==1){
          G_vec<-rep(NA_real_,R)
          for (r in 1:R){
            G_vec[r]<-exp(mu[r]*omega-(mu[r])^2/2)
          }
          probab_discrete<-G_vec/sum(G_vec)
          network[i,j] <- which(as.vector(rmultinom(n = 1, size = 1, prob = probab_discrete))==1)
        }
      }
    }
  }
  else{
    clust_id_U_sampling<-rmultinom(n = N1, size = 1, prob = pi_U)
    clust_id_P_sampling<-rmultinom(n = N2, size = 1, prob = pi_P)
    clust_ids_U<-apply(X = clust_id_U_sampling,MARGIN = 2,FUN = function(x){
      y<-which(x==1)
      return(y)
    })
    clust_ids_P<-apply(X = clust_id_P_sampling,MARGIN = 2,FUN = function(x){
      y<-which(x==1)
      return(y)
    })
    
    logitInv_delta_0<-logitInv_delta_0_cal(delta_0=delta_0, K1=K1, K2=K2, R=R)
    probab_block_ratings<-array(NA_real_,dim=c(K1,K2,R))
    for(k in 1:K1){
      for(l in 1:K2){
        for(r in 1:R){
          if(r==1){
            probab_block_ratings[k,l,r]<-logitInv_delta_0[k,l,r]
          }else if(r==R){
            probab_block_ratings[k,l,r]<-1-logitInv_delta_0[k,l,(R-1)]
          }else{
            probab_block_ratings[k,l,r]<-logitInv_delta_0[k,l,r]-logitInv_delta_0[k,l,(r-1)]
          }
        }
      }
    }
    
    for (i in 1:N1){
      for (j in 1:N2){
        exp_val<-exp(theta[clust_ids_U[i],clust_ids_P[j]]+sum(cov_u_outer[i,]*beta_u[clust_ids_U[i],])+sum(cov_p_outer[j,]*beta_p[clust_ids_P[j],]))
        probab<-exp_val/(1+exp_val)
        edge_sample<-rbinom(n = 1, size = 1, prob = probab)
        if(edge_sample==1){
          net_adjacency[i,j]<-1
          net_rating[i,j] <- which(as.vector(rmultinom(n = 1, size = 1, prob = probab_block_ratings[clust_ids_U[i],clust_ids_P[j],]))==1)
        }
      }
    }
  }
  return(list(net_adjacency,net_rating,clust_ids_U,clust_ids_P))
}

########################################################################################################
########################################################################################################
## Model 10 features (This works!):
## 1) Block structure over intercept
## 2) Proportional odds model for ratings based on my own derivations
## 3) Segmentation level (K*1) Covariate parameters in outer edge part
## 4) Stochastic over K

wrapper_bipartite_model10<-function(net_adjacency,net_rating,nclust1,nclust2,thres=10^(-6),n_iter_min=50,n_iter_max=500,theta_init,beta_u_init,beta_p_init,delta_0_init,cov_u_outer,cov_p_outer,sim_indicator,theta_true=NA,beta_u_true=NA,beta_p_true=NA, delta_0_true=NA,K_true_U=NA,K_true_P=NA,cluster_ids_true_U=NA,cluster_ids_true_P=NA,R=NA){
  ## This is the combined final updated and working code for EM undirected case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(quadprog)
  library(combinat)
  library(MASS)
  library(Rcpp)
  library(Matrix)
  library(gtools)
  #sourceCpp(file = paste0(file_path,"model_10.cpp"))
  library(biSBMcovSEM)
  #################################################
  ## Defining a function to update variational parameters gamma using quadratic program solver. Input:
  ## gamma.curr is a N*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
  ## estimates, theta.curr is a T_grid*K array for current theta estimates. Given Network is assumed to be an N*N*T_data array
  solve_QP_wrapper<-function(quad_lin_coeff,constraint_matrix,constraint_vector,node_ID,gamma.curr){
    try_QP<-try(gamma_next_vec<-solve.QP((-2*diag(as.vector(quad_lin_coeff[node_ID,,1]))),as.vector(quad_lin_coeff[node_ID,,2]),t(constraint_matrix),constraint_vector)$solution,silent = TRUE)
    if('numeric' %in% class(try_QP)){
      gamma_next_vec<-try_QP
    }
    else{
      gamma_next_vec<-gamma.curr[node_ID,]
      print("red flag")
      print(paste0("Node ID for red flag is ", node_ID))
    }
    return(gamma_next_vec)
  }
  
  gamma_U.update.wrapper<-function(gamma_U.curr,gamma_P.curr,pi_U.curr,theta.curr, logitInv_delta_0.curr, net_adjacency, net_rating, cov_beta_u, cov_beta_p, N1,N2,K1,K2,R){
    gamma_U.next<-matrix(NA_real_,N1,K1)
    constraint_matrix<-matrix(NA_real_,2+K1,K1)
    constraint_matrix[1,]<-rep(1,K1)
    constraint_matrix[2,]<-rep(-1,K1)
    constraint_matrix[3:(K1+2),]<-diag(K1)
    constraint_vector<-c(1,-1,rep(0,K1))
    
    quad_lin_coeff<-gamma_U_update_stat_undir(gamma_U=gamma_U.curr, gamma_P=gamma_P.curr, pi_U=pi_U.curr, theta=theta.curr, logitInv_delta_0=logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating=net_rating, cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p,  N1=N1, N2=N2, K1=K1, K2=K2, R=R)
    
    #print(quad_lin_coeff)
    for (i in 1:N1){
      if(sum(is.nan(quad_lin_coeff[i,,1]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,1])),1]<--Inf
      }
      if(sum(is.nan(quad_lin_coeff[i,,2]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,2])),2]<-Inf
      }
      if(sum(quad_lin_coeff[i,,1]<=-10^8)==0){
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma_U.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_U.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K1)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma_U.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_U.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K1){
          gamma_U.next[i,]<-gamma.curr[i,]
        }
      }else if((sum(quad_lin_coeff[i,,1]<=-10^8)>0)&&(sum(quad_lin_coeff[i,,1]<=-10^8)<K1)){
        quad_lin_coeff[i,which(quad_lin_coeff[i,,1]<=-10^8),1]<--(10^10)
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma_U.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_U.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K1)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma_U.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_U.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K1){
          gamma_U.next[i,]<-gamma_U.curr[i,]
        }
      }else if(sum(quad_lin_coeff[i,,1]<=-10^8)==K1){
        gamma_U.next[i,]<-gamma_U.curr[i,]
      }
      #print(i)
    }
    #print('gammaU')
    #print(quad_lin_coeff)
    # for (i in 1:N1){
    #   gamma_U.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    #   print(i)
    # }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma_U.next<-t(apply(X = gamma_U.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma_U.next)
  }
  
  gamma_P.update.wrapper<-function(gamma_U.curr,gamma_P.curr,pi_P.curr,theta.curr, logitInv_delta_0.curr, net_adjacency, net_rating, cov_beta_u, cov_beta_p, N1,N2,K1,K2,R){
    gamma_P.next<-matrix(NA_real_,N2,K2)
    constraint_matrix<-matrix(NA_real_,2+K2,K2)
    constraint_matrix[1,]<-rep(1,K2)
    constraint_matrix[2,]<-rep(-1,K2)
    constraint_matrix[3:(K2+2),]<-diag(K2)
    constraint_vector<-c(1,-1,rep(0,K2))
    
    quad_lin_coeff<-gamma_P_update_stat_undir(gamma_U=gamma_U.curr, gamma_P=gamma_P.curr, pi_P=pi_P.curr, theta=theta.curr, logitInv_delta_0=logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating=net_rating, cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, N1=N1, N2=N2, K1=K1, K2=K2, R=R)
    #print('gammaP')
    #print(quad_lin_coeff)
    
    #print(quad_lin_coeff)
    for (i in 1:N2){
      if(sum(is.nan(quad_lin_coeff[i,,1]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,1])),1]<--Inf
      }
      if(sum(is.nan(quad_lin_coeff[i,,2]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,2])),2]<-Inf
      }
      if(sum(quad_lin_coeff[i,,1]<=-10^8)==0){
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma_P.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_P.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K2)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma_P.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_P.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K2){
          gamma_U.next[i,]<-gamma.curr[i,]
        }
      }else if((sum(quad_lin_coeff[i,,1]<=-10^8)>0)&&(sum(quad_lin_coeff[i,,1]<=-10^8)<K2)){
        quad_lin_coeff[i,which(quad_lin_coeff[i,,1]<=-10^8),1]<--(10^10)
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma_P.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_P.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K2)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma_P.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma_P.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K2){
          gamma_P.next[i,]<-gamma_P.curr[i,]
        }
      }else if(sum(quad_lin_coeff[i,,1]<=-10^8)==K2){
        gamma_P.next[i,]<-gamma_P.curr[i,]
      }
      #print(i)
    }
    
    # for (i in 1:N2){
    #   gamma_P.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    #   print(i)
    # }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma_P.next<-t(apply(X = gamma_P.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma_P.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of pi
  pi.update<-function(gamma.curr){
    pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## normalization of pi
    pi.next<-pi.next/sum(pi.next)
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions K1,K2
  theta.update<-function(theta.curr,cluster_ids_U, cluster_ids_P,net_adjacency, cov_beta_u, cov_beta_p, N1,N2,K1,K2){
    theta.next<-matrix(NA_real_,K1,K2)
    gradient<-grad_bipartite_stat_undir_theta(theta=theta.curr, cluster_ids_U=cluster_ids_U, cluster_ids_P=cluster_ids_P, net_adjacency=net_adjacency,  cov_beta_u = cov_beta_u, cov_beta_p = cov_beta_p, N1=N1, N2=N2, K1=K1, K2=K2)
    hess<-hess_bipartite_stat_undir_theta(theta=theta.curr, cluster_ids_U=cluster_ids_U, cluster_ids_P=cluster_ids_P, cov_beta_u = cov_beta_u, cov_beta_p = cov_beta_p, N1=N1, N2=N2, K1=K1, K2=K2)
    for(k in 1:K1){
      for(l in 1:K2){
        theta.next[k,l]<-theta.curr[k,l]-gradient[k,l]/hess[k,l]
      }
    }
    return(theta.next)
  }
  
  #################################################
  ## Defining the update function of full vector beta_u with dimensions K1*C1_outer
  beta_u.update<-function(beta_u.curr,cluster_ids_U, cluster_ids_P,theta,net_adjacency, cov_beta_u, cov_beta_p, cov_u_outer, N1,N2,K1,K2,C1_outer){
    beta_u.next<-array(NA_real_,dim=c(K1,C1_outer))
    gradient<-grad_bipartite_stat_undir_beta_u(cluster_ids_U=cluster_ids_U, cluster_ids_P=cluster_ids_P, theta=theta, net_adjacency=net_adjacency, cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, cov_u_outer = cov_u_outer,  N1=N1, N2=N2, K1=K1, K2=K2, C1_outer=C1_outer)
    hess<-hess_bipartite_stat_undir_beta_u(cluster_ids_U=cluster_ids_U, cluster_ids_P=cluster_ids_P, theta=theta, cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, cov_u_outer = cov_u_outer, N1=N1, N2=N2, K1=K1, K2=K2, C1_outer=C1_outer)
    for(c_u in 1:C1_outer){
      for(k in 1:K1){
        beta_u.next[k,c_u]<-beta_u.curr[k,c_u]-(gradient[k,c_u]/hess[k,c_u])
      }
    }
    return(beta_u.next)
  }
  
  #################################################
  ## Defining the update function of full vector beta_p with dimensions K2*C2_outer
  beta_p.update<-function(beta_p.curr,cluster_ids_U,cluster_ids_P,theta,net_adjacency, cov_beta_u, cov_beta_p, cov_p_outer, N1,N2,K1,K2,C2_outer){
    beta_p.next<-array(NA_real_,dim=c(K2,C2_outer))
    gradient<-grad_bipartite_stat_undir_beta_p(cluster_ids_U=cluster_ids_U, cluster_ids_P=cluster_ids_P, theta=theta, net_adjacency=net_adjacency, cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, cov_p_outer = cov_p_outer, N1=N1, N2=N2, K1=K1, K2=K2, C2_outer=C2_outer)
    hess<-hess_bipartite_stat_undir_beta_p(cluster_ids_U=cluster_ids_U, cluster_ids_P=cluster_ids_P, theta=theta, cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, cov_p_outer = cov_p_outer, N1=N1, N2=N2, K1=K1, K2=K2, C2_outer=C2_outer)
    for(c_p in 1:C2_outer){
      for(l in 1:K2){
        beta_p.next[l,c_p]<-beta_p.curr[l,c_p]-(gradient[l,c_p]/hess[l,c_p])
      }
    }
    return(beta_p.next)
  }
  
  #################################################
  ## Defining the update function of full proportional odds nuisance parameters delta_0
  delta_0.update<-function(cluster_ids_U,cluster_ids_P,delta_0,logitInv_delta_0,net_adjacency,net_rating_statistic,N1,N2,K1,K2,R){
    delta_0.next<-array(NA_real_,dim=c(K1,K2,(R-1)))
    gradient_delta_0<-grad_bipartite_stat_undir_delta_0(cluster_ids_U=cluster_ids_U, cluster_ids_P=cluster_ids_P, delta_0=delta_0, logitInv_delta_0=logitInv_delta_0, net_adjacency=net_adjacency, net_rating=net_rating, N1=N1, N2=N2, K1=K1, K2=K2, R=R)
    for (k in 1:K1){
      for (l in 1:K2){
        hess_mat<-hess_bipartite_stat_undir_delta_0(k=k-1, l=l-1, cluster_ids_U=cluster_ids_U, cluster_ids_P=cluster_ids_P, delta_0=delta_0, logitInv_delta_0=logitInv_delta_0, net_adjacency=net_adjacency, net_rating=net_rating, N1=N1, N2=N2, R=R)
        delta_0.next[k,l,]<-sort(delta_0[k,l,]-as.vector((solve(hess_mat+(10^(-6)*diag(R-1)))%*%gradient_delta_0[k,l,])))
        #delta_0.next[k,l,]<-sort(delta_0[k,l,]-as.vector((solve(hess_mat)%*%gradient_delta_0[k,l,])))
      }
    }
    #print(delta_0.next)
    return(delta_0.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,net_adjacency,net_rating,K1,K2,n_iter,thres,n_iter_min,n_iter_max,R,cov_u_outer,cov_p_outer){
    ## Assuming a network adjacency matrix of dimensions N1*N2 with rows corresponding to users and columns corresponding to products and each entry is an indicator whether the user rated the product. This framework can be generalized.
    N1<-dim(net_adjacency)[1]
    N2<-dim(net_adjacency)[2]
    C1_outer<-dim(cov_u_outer)[2]
    C2_outer<-dim(cov_p_outer)[2]
    
    ## initializing the arrays for parameters
    ## Initializing the cluster memberships gamma
    gamma_U<-array(NA_real_,dim=c(N1,K1,n_iter)) ## assuming there are K1 clusters in user set
    gamma_P<-array(NA_real_,dim=c(N2,K2,n_iter)) ## assuming there are K2 clusters in product set
    
    ## Initializing the mixture proportions pi
    pi_U<-matrix(NA_real_,K1,n_iter)
    pi_P<-matrix(NA_real_,K2,n_iter)
    
    ## Initializing the network parameters theta
    theta<-array(NA_real_,dim=c(K1,K2,n_iter))
    
    ## Initializing the covariate parameters beta_u and beta_p
    beta_u<-array(NA_real_,dim=c(K1,C1_outer,n_iter))
    beta_p<-array(NA_real_,dim=c(K2,C2_outer,n_iter))
    
    ## Initializing the proportional odds model intercept parameters
    delta_0<-array(NA_real_,dim=c(K1,K2,R-1,n_iter))
    
    gamma_U[,,1]<-start[[1]]
    gamma_P[,,1]<-start[[2]]
    pi_U[,1]<-start[[3]]
    pi_P[,1]<-start[[4]]
    theta[,,1]<-start[[5]]
    beta_u[,,1]<-start[[6]]
    beta_p[,,1]<-start[[7]]
    delta_0[,,,1]<-start[[8]]
    
    #print(theta[,,1])
    #print(omega[,1])
    #print(mu[,,,1])
    
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_curr<-10^10
    time_per_iter<-c()
    #iter_index<-iter_index+1
    while(((error>thres)&(iter_index<n_iter_max))|(iter_index<n_iter_min)){ 
      #while((error>thres|iter_index<200)&(iter_index<250)&(iter_index==2|error<10)){ 
      #while(iter_index<21){ 
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      
      ## Calculating cov_beta_u and cov_beta_p
      cov_beta_u<-cov_beta_mul(cov=cov_u_outer,beta=t(beta_u[,,iter_index-1]))
      cov_beta_p<-cov_beta_mul(cov=cov_p_outer,beta=t(beta_p[,,iter_index-1]))
      
      ## Calculating the logit inverse of delta_0
      logitInv_delta_0.curr<-logitInv_delta_0_cal(delta_0=delta_0[,,,iter_index-1], K1=K1, K2=K2, R=R)
      
      ## Transforming logit inverse of delta_0 to phi
      #phi.curr<-phi_cal(logitInv_delta_0=logitInv_delta_0.curr, K1=K1, K2=K2, R=R)
      
      ## Calculating the logit inverse of phi
      #logitInv_phi.curr<-logitInv_phi_cal(phi=phi.curr, K1=K1, K2=K2, R=R)
      
      ## Updating gamma_U, the mixed memberships for users
      gamma_U[,,iter_index]<-gamma_U.update.wrapper(gamma_U.curr=gamma_U[,,iter_index-1], gamma_P.curr=gamma_P[,,iter_index-1], pi_U.curr=pi_U[,iter_index-1], theta.curr=theta[,,iter_index-1], cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, logitInv_delta_0.curr=logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating=net_rating, N1=N1,N2=N2,K1=K1,K2=K2,R=R)
      ## Updating the mixture proportions pi_U
      pi_U[,iter_index]<-pi.update(gamma.curr=gamma_U[,,iter_index])
      
      ## Stochastic Step: Sampling the cluster memberships for each user node
      cluster_ids_U<-as.vector(apply(X = matrix(1:N1),MARGIN = 1,FUN = function(x){
        y<-as.vector(rmultinom(n = 1,size = 1,prob = gamma_U[x,,iter_index]))
        cluster_id<-which(y==1)
        return(cluster_id)
      }))
      
      ## Updating gamma_P, the mixed memberships for products
      gamma_P[,,iter_index]<-gamma_P.update.wrapper(gamma_U.curr=gamma_U[,,iter_index],gamma_P.curr=gamma_P[,,iter_index-1],pi_P.curr=pi_P[,iter_index-1], theta.curr=theta[,,iter_index-1], cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, logitInv_delta_0.curr=logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating=net_rating, N1=N1,N2=N2,K1=K1,K2=K2,R=R)
      ## Updating the mixture proportions pi_P
      pi_P[,iter_index]<-pi.update(gamma.curr=gamma_P[,,iter_index])
      
      ## Stochastic Step: Sampling the cluster memberships for each user node
      cluster_ids_P<-as.vector(apply(X = matrix(1:N2),MARGIN = 1,FUN = function(x){
        y<-as.vector(rmultinom(n = 1,size = 1,prob = gamma_P[x,,iter_index]))
        cluster_id<-which(y==1)
        return(cluster_id)
      }))
      
      ## Updating the network parameters theta
      theta[,,iter_index]<-theta.update(theta.curr=theta[,,iter_index-1],cluster_ids_U=cluster_ids_U, cluster_ids_P=cluster_ids_P, net_adjacency=net_adjacency, cov_beta_u = cov_beta_u, cov_beta_p = cov_beta_p, N1=N1, N2=N2, K1=K1, K2=K2)
      #print('theta')
      #print(theta[,,iter_index])
      
      ## Updating the covariate parameters beta_u
      beta_u[,,iter_index]<-beta_u.update(beta_u.curr = matrix(beta_u[,,iter_index-1],K1,C1_outer),cluster_ids_U=cluster_ids_U, cluster_ids_P=cluster_ids_P, theta=theta[,,iter_index], net_adjacency=net_adjacency, cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, cov_u_outer = cov_u_outer, N1=N1, N2=N2, K1=K1, K2=K2, C1_outer= C1_outer)
      #print('beta_u and beta_p')
      #print(beta_u[,,iter_index])
      
      ## Updating the covariate parameters beta_p
      beta_p[,,iter_index]<-beta_p.update(beta_p.curr = matrix(beta_p[,,iter_index-1],K2,C2_outer),cluster_ids_U=cluster_ids_U, cluster_ids_P=cluster_ids_P, theta=theta[,,iter_index], net_adjacency=net_adjacency, cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, cov_p_outer = cov_p_outer, N1=N1, N2=N2, K1=K1, K2=K2, C2_outer=C2_outer)
      #print(beta_p[,,iter_index])
      
      ## Updating the omega matrix
      delta_0[,,,iter_index]<-delta_0.update(cluster_ids_U=cluster_ids_U, cluster_ids_P=cluster_ids_P, delta_0=delta_0[,,,iter_index-1],logitInv_delta_0 = logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating_statistic = net_rating_statistic, N1=N1, N2=N2, K1=K1, K2=K2, R=R)
      #print('delta_0')
      #print(delta_0[,,,iter_index])
      #Sys.sleep(1.5)
      
      ELBO_prev<-ELBO_curr
      ELBO_curr<-ELBO_conv_bipartite_stat_undir(gamma_U=gamma_U[,,iter_index], gamma_P=gamma_P[,,iter_index], pi_U=pi_U[,iter_index], pi_P=pi_P[,iter_index], theta = theta[,,iter_index], logitInv_delta_0=logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating=net_rating, cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, N1=N1, N2=N2, K1=K1, K2=K2)
      
      error<-abs(((ELBO_prev-ELBO_curr)/ELBO_prev))
      
      # if(!is.nan(abs(((ELBO_prev-ELBO_curr)/ELBO_prev)))){
      #   error<-abs(((ELBO_prev-ELBO_curr)/ELBO_prev))
      # }
      
      #print('error')
      #print(error)
      
      #print("probs")
      #print(Gvals_cube_0)
      # print(apply(X = Gvals_cube_0,MARGIN = c(1,2),FUN = function(x){
      #   
      #   return(x/sum(x))
      # }))
      
      if((iter_index%%1)==0){
        #print(theta[,,iter_index])
        print('iter_index')
        print(iter_index)
        print('error')
        print(error)
        print(proc.time()-ptm)
      }
      ptm_diff<-proc.time()-ptm
      time_per_iter<-c(time_per_iter,ptm_diff[3])
      iter_index<-iter_index+1
    }
    names(time_per_iter)<-NULL
    total_iterations<-iter_index-1
    return(list(gamma_U,gamma_P,pi_U,pi_P,theta,beta_u,beta_p,delta_0,time_per_iter,total_iterations))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining Model Selection functions based on converged paramters
  
  #################################################
  ## Defining a function to calculate the estimate of complete log likelihood
  comp_loglik_cal<-function(cluster_ids_est_U = NA, cluster_ids_est_P = NA, pi_U = NA, pi_P = NA, theta, beta_u, beta_p, delta_0, net_adjacency, net_rating, cov_u_outer, cov_p_outer, N1, N2, K1, K2, R){
    if((K1!=1)&(K2!=1)){
      logitInv_delta_0<-logitInv_delta_0_cal(delta_0=delta_0, K1=K1, K2=K2, R=R)
      probab_block_ratings<-array(NA_real_,dim=c(K1,K2,R))
      for(k in 1:K1){
        for(l in 1:K2){
          for(r in 1:R){
            if(r==1){
              probab_block_ratings[k,l,r]<-logitInv_delta_0[k,l,r]
            }else if(r==R){
              probab_block_ratings[k,l,r]<-1-logitInv_delta_0[k,l,(R-1)]
            }else{
              probab_block_ratings[k,l,r]<-logitInv_delta_0[k,l,r]-logitInv_delta_0[k,l,(r-1)]
            }
          }
        }
      }
      
      ## 1st term
      t1<-0
      for(i in 1:N1){
        for(j in 1:N2){
          cluster_id_i<-cluster_ids_est_U[i]
          cluster_id_j<-cluster_ids_est_P[j]
          predictor<-theta[cluster_id_i,cluster_id_j]+sum(cov_u_outer[i,]*beta_u[cluster_id_i,])+sum(cov_p_outer[j,]*beta_p[cluster_id_j,])
          exp_val<-exp(predictor)
          if(net_adjacency[i,j]==0){
            t1<-t1-(log(1+exp_val))
          }else{
            for(r in 1:R){
              if(net_rating[i,j]==r){
                t1<-t1+(predictor-log(1+exp_val)+log(probab_block_ratings[cluster_id_i,cluster_id_j,r]))
              }
            }
          }
        }
      }
      ## 2nd term
      t2<-0
      for (i in 1:N1){
        t2<-t2+log(pi_U[cluster_ids_est_U[i]])
      }
      t3<-0
      for (j in 1:N2){
        t3<-t3+log(pi_P[cluster_ids_est_P[j]])
      }
      
      comp_val<-t1+t2+t3
    }
    return(comp_val)
  }
  
  #################################################
  ## Defining a function to calculate the integrated classification likelihood
  ICL_cal<-function(cluster_ids_est_U, cluster_ids_est_P, pi_U, pi_P, theta, beta_u, beta_p, delta_0, net_adjacency, net_rating, cov_u_outer, cov_p_outer, N1, N2, K1, K2, R){
    comp_loglik_val<-comp_loglik_cal(cluster_ids_est_U = cluster_ids_est_U, cluster_ids_est_P = cluster_ids_est_P, pi_U = pi_U, pi_P = pi_P, theta = theta, beta_u = beta_u, beta_p = beta_p, delta_0 = delta_0, cov_u_outer = cov_u_outer, cov_p_outer = cov_p_outer, net_adjacency = net_adjacency, net_rating = net_rating, N1 = N1, N2 = N2, K1 = K1, K2 = K2, R=R)
    ## Penalization over pi
    t2<-(K1-1)*log(N1)+(K2-1)*log(N2)
    ## Penalization over theta
    t3<-(K1*K2)*log(N1*N2)
    ## Penalization over beta_u
    C1_outer<-dim(cov_u_outer)[2]
    t4<-(K1*C1_outer)*log(N1*N2)
    ## Penalization over beta_p
    C2_outer<-dim(cov_p_outer)[2]
    t5<-(K2*C2_outer)*log(N1*N2)
    ## Calculating the number of edges
    N0_val<-N0_cal(net_adjacency=net_adjacency, N1=N1, N2=N2)
    ## Penalization over omega, mu
    t6<-(K1*K2*(R-1))*log(N0_val)
    
    ## all penalty terms together:
    penalty<-(t2+t3+t4+t5+t6)/2
    ICL_val<-comp_loglik_val-penalty
    ICL_vec<-c(comp_loglik_val,penalty,ICL_val)
    return(ICL_vec)
  }
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI_cal<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  #################################################
  ## RASE functions
  RASE_cal<-function(param_est, param_true){
    if(is.vector(param_est)!=1){
      RASE_val<-sqrt(sum((param_est-param_true)^2)/prod(dim(param_est)))
    }else{RASE_val<-sqrt(sum((param_est-param_true)^2)/length(param_est))}
    return(RASE_val)
  }
  
  ########################################################################################################
  ## Defining the functions to evaluate clustering performance using normalized mutual information.
  probab_clust_cal<-function(cluster_ids,N,K){
    probab_vec<-rep(NA_real_,K)
    for(k in 1:K){
      probab_vec[k]<-sum(cluster_ids==k)/N
    }
    return(probab_vec)
  }
  NMI_cal<-function(cluster_ids_true_U,cluster_ids_true_P,K_true_U,K_true_P,cluster_ids_est_U,cluster_ids_est_P,K1,K2,N1,N2){
    probab_vec_true_U<-probab_clust_cal(cluster_ids = cluster_ids_true_U,N = N1+N2,K = K_true_U)
    probab_vec_true_P<-probab_clust_cal(cluster_ids = cluster_ids_true_P,N = N1+N2,K = K_true_P)
    probab_vec_true<-c(probab_vec_true_U,probab_vec_true_P)
    total_entropy_true<-0
    for(k in 1:length(probab_vec_true)){
      if(probab_vec_true[k]!=0){
        total_entropy_true<-total_entropy_true-(probab_vec_true[k]*log(probab_vec_true[k],base = 2))
      }
    }
    
    probab_vec_est_U<-probab_clust_cal(cluster_ids = cluster_ids_est_U,N = N1+N2,K = K1)
    probab_vec_est_P<-probab_clust_cal(cluster_ids = cluster_ids_est_P,N = N1+N2,K = K2)
    probab_vec_est<-c(probab_vec_est_U,probab_vec_est_P)
    total_entropy_est<-0
    for(k in 1:length(probab_vec_est)){
      if(probab_vec_est[k]!=0){
        total_entropy_est<-total_entropy_est-(probab_vec_est[k]*log(probab_vec_est[k],base = 2))
      }
    }
    
    clust_ids_true_total<-c(cluster_ids_true_U,(cluster_ids_true_P+K_true_U))
    clust_ids_est_total<-c(cluster_ids_est_U,(cluster_ids_est_P+K1))
    
    conditional_entropy<-rep(NA_real_,(K1+K2))
    for(k in 1:(K1+K2)){
      probab_conditional_vec<-rep(0,(K_true_U+K_true_P))
      for(m in 1:(K_true_U+K_true_P)){
        if(sum(clust_ids_est_total==k)!=0){
          probab_conditional_vec[m]<-sum((clust_ids_true_total==m)&(clust_ids_est_total==k))/sum(clust_ids_est_total==k)
        }
      }
      plogp_conditional_vec<-rep(0,length(probab_conditional_vec))
      for(m in 1:length(probab_conditional_vec)){
        if(probab_conditional_vec[m]!=0){
          plogp_conditional_vec[m]<-probab_conditional_vec[m]*log(x = probab_conditional_vec[m],base = 2)
        }
      }
      
      conditional_entropy[k]<--probab_vec_est[k]*sum(plogp_conditional_vec)
    }
    
    mutual_info<-total_entropy_true-sum(conditional_entropy)
    NMI_val<-(2*mutual_info)/(total_entropy_true+total_entropy_est)
    return(NMI_val)
  }
  
  ########################################################################################################
  ## Calculating the recovery rate
  recovery_rate_cal <- function(cluster_ids_est, cluster_ids_true){
    all_perm<-permutations(n=4,r=2,repeats.allowed = T)
    clust_perm_mat<-matrix(NA_integer_,4,4)
    for(i in 1:nrow(all_perm)){
      clust_perm_mat[all_perm[i,1],all_perm[i,2]] <- sum((cluster_ids_true==all_perm[i,1])&(cluster_ids_est==all_perm[i,2]))
    }
    recovery_rate_val<-sum(apply(X = clust_perm_mat, MARGIN = 1, FUN = max))/length(cluster_ids_est)
    return(recovery_rate_val)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N1<-dim(net_adjacency)[1] ## Number of nodes from the network
  N2<-dim(net_adjacency)[2]
  ## Total time points in the network for which we have data in the form of adjacency matrix
  K1<-nclust1 ## Defining the number of clusters
  K2<-nclust2
  K <-K1*K2
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  start[[1]]<-matrix(rep(1/K1,N1*K1),N1,K1)
  start[[2]]<-matrix(rep(1/K2,N2*K2),N2,K2)
  start[[3]]<-rep(1/K1,K1)
  start[[4]]<-rep(1/K2,K2)
  start[[5]]<-theta_init
  start[[6]]<-beta_u_init
  start[[7]]<-beta_p_init
  start[[8]]<-delta_0_init
  if(K1!=1&K2!=1){ param<-iterator(start=start, net_adjacency = net_adjacency, net_rating=net_rating, K1=K1, K2=K2, n_iter=2000, thres=thres, n_iter_min = n_iter_min, n_iter_max = n_iter_max, R=R, cov_u_outer=cov_u_outer, cov_p_outer=cov_p_outer)}
  if(K1==1&K2!=1){ param<-iterator_K1(start=start, network=sim.net, K1=K1, K2=K2, n_iter=1000, thres=thres,R=R)}
  if(K1!=1&K2==1){ param<-iterator_K2(start=start, network=sim.net, K1=K1, K2=K2, n_iter=1000, thres=thres,R=R)}
  if(K1==1&K2==1){ param<-iterator_K1K2(start=start, network=sim.net, K1=K1, K2=K2, n_iter=1000, thres=thres,R=R)}
  
  #################################################
  ## extracting the coverged parameter values and calculating ICL
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    temp<-is.na(param[[1]][1,1,n_iter])
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  #print(param)
  
  if(K1!=1&K2!=1){
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,,n_last]
    param_converge[[3]]<-param[[3]][,n_last]
    param_converge[[4]]<-param[[4]][,n_last]
    param_converge[[5]]<-param[[5]][,,n_last]
    param_converge[[6]]<-param[[6]][,,n_last]
    param_converge[[7]]<-param[[7]][,,n_last]
    param_converge[[8]]<-param[[8]][,,,n_last]
    output_list<-param_converge
    cluster_ids_est_U<-as.vector(apply(X = matrix(1:N1),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[1]][x,])
      return(cluster_id)
    }))
    cluster_ids_est_P<-as.vector(apply(X = matrix(1:N2),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[2]][x,])
      return(cluster_id)
    }))
    #ICL_val<-NA
    ICL_vec<-ICL_cal(cluster_ids_est_U = cluster_ids_est_U, cluster_ids_est_P = cluster_ids_est_P, pi_U=param_converge[[3]], pi_P=param_converge[[4]] ,theta = param_converge[[5]], beta_u = as.matrix(param_converge[[6]]), beta_p = as.matrix(param_converge[[7]]), delta_0 = param_converge[[8]], net_adjacency=net_adjacency, net_rating=net_rating, cov_u_outer = cov_u_outer, cov_p_outer = cov_p_outer, N1 = N1, N2 = N2, K1 = K1, K2 = K2, R = R)
    ## Calculating cov_beta_u and cov_beta_p
    cov_beta_u<-cov_beta_mul(cov=cov_u_outer,beta=t(as.matrix(param_converge[[6]])))
    cov_beta_p<-cov_beta_mul(cov=cov_p_outer,beta=t(as.matrix(param_converge[[7]])))
    logitInv_delta_0.curr<-logitInv_delta_0_cal(delta_0=param_converge[[8]], K1=K1, K2=K2, R=R)
    ELBO_val<-ELBO_conv_bipartite_stat_undir(gamma_U=param_converge[[1]], gamma_P=param_converge[[2]], pi_U=param_converge[[3]], pi_P=param_converge[[4]], theta = param_converge[[5]], logitInv_delta_0=logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating=net_rating, cov_beta_u=cov_beta_u, cov_beta_p=cov_beta_p, N1=N1, N2=N2, K1=K1, K2=K2)
    ELBO_modsel_vec<-c(ELBO_val,ICL_vec[2],ELBO_val-ICL_vec[2])
    if(sim_indicator==1){
      RI_val_U<-RI_cal(cluster_ids_est = cluster_ids_est_U,cluster_ids_true = cluster_ids_true_U)
      RI_val_P<-RI_cal(cluster_ids_est = cluster_ids_est_P,cluster_ids_true = cluster_ids_true_P)
      NMI_val<-NMI_cal(cluster_ids_true_U = cluster_ids_true_U,cluster_ids_true_P = cluster_ids_true_P,K_true_U = K_true_U,K_true_P = K_true_P,cluster_ids_est_U = cluster_ids_est_U,cluster_ids_est_P = cluster_ids_est_P, K1 = K1, K2 = K2, N1 = N1, N2 = N2)
      Recovery_rate_val_U<-recovery_rate_cal(cluster_ids_est = cluster_ids_est_U,cluster_ids_true = cluster_ids_true_U)
      Recovery_rate_val_P<-recovery_rate_cal(cluster_ids_est = cluster_ids_est_P,cluster_ids_true = cluster_ids_true_P)
      if((K1==K_true_U)&(K2==K_true_P)){
        K1_permute_mat<-do.call(rbind,permn(1:K1))
        K2_permute_mat<-do.call(rbind,permn(1:K2))
        
        RASE_beta_u_vec<-rep(NA_real_,nrow(K1_permute_mat))
        for (k in 1:nrow(K1_permute_mat)){
          beta_u_est<-as.vector(param_converge[[6]])[K1_permute_mat[k,]]
          RASE_beta_u_vec[k]<-RASE_cal(param_est = beta_u_est, param_true = as.vector(beta_u_true))
        }
        permute_true_id_beta_u<-which.min(RASE_beta_u_vec)
        RASE_beta_u_val<-RASE_beta_u_vec[permute_true_id_beta_u]
        
        RASE_beta_p_vec<-rep(NA_real_,nrow(K2_permute_mat))
        for (k in 1:nrow(K2_permute_mat)){
          beta_p_est<-as.vector(param_converge[[7]])[K2_permute_mat[k,]]
          RASE_beta_p_vec[k]<-RASE_cal(param_est = beta_p_est, param_true = as.vector(beta_p_true))
        }
        permute_true_id_beta_p<-which.min(RASE_beta_p_vec)
        RASE_beta_p_val<-RASE_beta_p_vec[permute_true_id_beta_p]
        
        RASE_theta_mat<-matrix(NA_real_,nrow(K1_permute_mat),nrow(K2_permute_mat))
        for (k in 1:nrow(K1_permute_mat)){
          for (l in 1:nrow(K2_permute_mat)){
            theta_true_rotated<-as(as.integer(K1_permute_mat[k,]), "pMatrix")%*%theta_true%*%t(as(as.integer(K2_permute_mat[l,]), "pMatrix"))
            RASE_theta_mat[k,l]<-RASE_cal(param_est = param_converge[[5]], param_true = theta_true_rotated)
          }
        }
        
        RASE_theta_val<-min(RASE_theta_mat)
        rotate_indices<-which(RASE_theta_mat == min(RASE_theta_mat), arr.ind = TRUE)
        
        ## RASE for delta_0
        # RASE_delta_0_vec<-rep(NA_real_,R-1)
        # for (r in 1:(R-1)){
        #   delta_0_true_rotated<-as(as.integer(K1_permute_mat[rotate_indices[1],]), "pMatrix")%*%delta_0_true[,,r]%*%t(as(as.integer(K2_permute_mat[rotate_indices[2],]), "pMatrix"))
        #   RASE_delta_0_vec[r]<-RASE_cal(param_est = param_converge[[8]][,,r], param_true = delta_0_true_rotated)
        # }
        
        delta_0_true_rotated<-array(NA_real_,dim=c(K1,K2,(R-1)))
        for (r in 1:(R-1)){
          delta_0_true_rotated[,,r]<-as(as.integer(K1_permute_mat[rotate_indices[1],]), "pMatrix")%*%delta_0_true[,,r]%*%t(as(as.integer(K2_permute_mat[rotate_indices[2],]), "pMatrix"))
          
        }
        RASE_delta_0_val<-RASE_cal(param_est = param_converge[[8]], param_true = delta_0_true_rotated)
        
        output_list<-list("parameters_converged"=param_converge,"Set1_estimated_cluster_IDs"=cluster_ids_est_U,"Set2_estimated_cluster_IDs"=cluster_ids_est_P,"ICL"=ICL_vec,"ELBO_modsel"=ELBO_modsel_vec,"Set1_Rand_Index"=RI_val_U,"Set2_Rand_Index"=RI_val_P,"Normalized_Mutual_Info"=NMI_val, "Set1_Recovery_Rate"=Recovery_rate_val_U, "Set2_Recovery_Rate"=Recovery_rate_val_P, "RASE_theta"=RASE_theta_val, "RASE_beta_u"=RASE_beta_u_val, "RASE_beta_p"=RASE_beta_p_val, "RASE_delta_0"=RASE_delta_0_val,"time_per_iter_vec"=param[[9]], "total_iterations"=param[[10]])
      }else{
        output_list<-list("parameters_converged"=param_converge,"Set1_estimated_cluster_IDs"=cluster_ids_est_U,"Set2_estimated_cluster_IDs"=cluster_ids_est_P,"ICL"=ICL_vec,"ELBO_modsel"=ELBO_modsel_vec,"Set1_Rand_Index"=RI_val_U,"Set2_Rand_Index"=RI_val_P,"Normalized_Mutual_Info"=NMI_val, "Set1_Recovery_Rate"=Recovery_rate_val_U, "Set2_Recovery_Rate"=Recovery_rate_val_P, "time_per_iter_vec"=param[[9]], "total_iterations"=param[[10]])
      }
    }else{output_list<-list("parameters_converged"=param_converge,"Set1_estimated_cluster_IDs"=cluster_ids_est_U,"Set2_estimated_cluster_IDs"=cluster_ids_est_P,"ICL"=ICL_vec,"ELBO_modsel"=ELBO_modsel_vec,"time_per_iter_vec"=param[[9]], "total_iterations"=param[[10]])}
    #}
  }
  return(output_list)
}

simulate_network_bipartite_model10<-function(N1,N2,pi_U,pi_P,theta,beta_u,beta_p,delta_0,cov_u_outer,cov_p_outer,K1,K2,R){
  library(Rcpp)
  sourceCpp(file = paste0(file_path,"model_10.cpp"))
  #library(Matrix)
  net_adjacency<-matrix(0,N1,N2)
  net_rating<-matrix(0,N1,N2)
  #net_adjacency<-Matrix(0,N1,N2,sparse = T)
  #net_rating<-Matrix(0,N1,N2,sparse = T)
  if((K1==1)&(K2==1)){
    clust_id_U_sampling<-rmultinom(n = N1, size = 1, prob = pi_U)
    clust_id_P_sampling<-rmultinom(n = N2, size = 1, prob = pi_P)
    clust_ids_U<-apply(X = clust_id_U_sampling,MARGIN = 2,FUN = function(x){
      y<-which(x==1)
      return(y)
    })
    clust_ids_P<-apply(X = clust_id_P_sampling,MARGIN = 2,FUN = function(x){
      y<-which(x==1)
      return(y)
    })
    for (i in 1:N1){
      for (j in 1:N2){
        exp_val<-exp(theta)
        probab<-exp_val/(1+exp_val)
        edge_sample<-rbinom(n = 1, size = 1, prob = probab)
        if(edge_sample==0){
          network[i,j]<-0
        }else if(edge_sample==1){
          G_vec<-rep(NA_real_,R)
          for (r in 1:R){
            G_vec[r]<-exp(mu[r]*omega-(mu[r])^2/2)
          }
          probab_discrete<-G_vec/sum(G_vec)
          network[i,j] <- which(as.vector(rmultinom(n = 1, size = 1, prob = probab_discrete))==1)
        }
      }
    }
  }
  else{
    clust_id_U_sampling<-rmultinom(n = N1, size = 1, prob = pi_U)
    clust_id_P_sampling<-rmultinom(n = N2, size = 1, prob = pi_P)
    clust_ids_U<-apply(X = clust_id_U_sampling,MARGIN = 2,FUN = function(x){
      y<-which(x==1)
      return(y)
    })
    clust_ids_P<-apply(X = clust_id_P_sampling,MARGIN = 2,FUN = function(x){
      y<-which(x==1)
      return(y)
    })
    
    logitInv_delta_0<-logitInv_delta_0_cal(delta_0=delta_0, K1=K1, K2=K2, R=R)
    probab_block_ratings<-array(NA_real_,dim=c(K1,K2,R))
    for(k in 1:K1){
      for(l in 1:K2){
        for(r in 1:R){
          if(r==1){
            probab_block_ratings[k,l,r]<-logitInv_delta_0[k,l,r]
          }else if(r==R){
            probab_block_ratings[k,l,r]<-1-logitInv_delta_0[k,l,(R-1)]
          }else{
            probab_block_ratings[k,l,r]<-logitInv_delta_0[k,l,r]-logitInv_delta_0[k,l,(r-1)]
          }
        }
      }
    }
    
    for (i in 1:N1){
      for (j in 1:N2){
        exp_val<-exp(theta[clust_ids_U[i],clust_ids_P[j]]+sum(cov_u_outer[i,]*beta_u[clust_ids_U[i],])+sum(cov_p_outer[j,]*beta_p[clust_ids_P[j],]))
        probab<-exp_val/(1+exp_val)
        edge_sample<-rbinom(n = 1, size = 1, prob = probab)
        if(edge_sample==1){
          net_adjacency[i,j]<-1
          net_rating[i,j] <- which(as.vector(rmultinom(n = 1, size = 1, prob = probab_block_ratings[clust_ids_U[i],clust_ids_P[j],]))==1)
        }
      }
    }
  }
  return(list(net_adjacency,net_rating,clust_ids_U,clust_ids_P))
}

########################################################################################################
########################################################################################################
## Model 11 features (This works!):
## 1) Combined (reviewers and products) single mode network
## 2) Block structure over intercept
## 3) Proportional odds model for ratings based on my own derivations

wrapper_bipartite_model11<-function(net_adjacency,net_rating,nclust,thres=10^(-6),n_iter_min=50,n_iter_max=500,theta_init,delta_0_init,sim_indicator,theta_true=NA, delta_0_true=NA,K_true=NA,cluster_ids_true_U=NA, cluster_ids_true_P=NA, N1, N2, K1, K2, R=NA){
  ## This is the combined final updated and working code for EM undirected case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(quadprog)
  library(combinat)
  library(MASS)
  library(Rcpp)
  library(Matrix)
  library(gtools)
  #sourceCpp(file = paste0(file_path,"model_11.cpp"))
  library(SBMcov)
  #################################################
  ## Defining a function to update variational parameters gamma using quadratic program solver. Input:
  ## gamma.curr is a N*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
  ## estimates, theta.curr is a T_grid*K array for current theta estimates. Given Network is assumed to be an N*N*T_data array
  solve_QP_wrapper<-function(quad_lin_coeff,constraint_matrix,constraint_vector,node_ID,gamma.curr){
    try_QP<-try(gamma_next_vec<-solve.QP((-2*diag(as.vector(quad_lin_coeff[node_ID,,1]))),as.vector(quad_lin_coeff[node_ID,,2]),t(constraint_matrix),constraint_vector)$solution,silent = TRUE)
    if('numeric' %in% class(try_QP)){
      gamma_next_vec<-try_QP
    }
    else{
      gamma_next_vec<-gamma.curr[node_ID,]
      print("red flag")
      print(paste0("Node ID for red flag is ", node_ID))
    }
    return(gamma_next_vec)
  }
  
  gamma.update.wrapper<-function(gamma.curr,pi.curr,theta.curr, logitInv_delta_0.curr, net_adjacency, net_rating, N,K,R){
    gamma.next<-matrix(NA_real_,N,K)
    constraint_matrix<-matrix(NA_real_,2+K,K)
    constraint_matrix[1,]<-rep(1,K)
    constraint_matrix[2,]<-rep(-1,K)
    constraint_matrix[3:(K+2),]<-diag(K)
    constraint_vector<-c(1,-1,rep(0,K))
    
    quad_lin_coeff<-gamma_update_stat_undir(gamma=gamma.curr, pi=pi.curr, theta=theta.curr, logitInv_delta_0=logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating=net_rating,  N=N, K=K, R=R)
    
    #print(quad_lin_coeff)
    for (i in 1:N){
      if(sum(is.nan(quad_lin_coeff[i,,1]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,1])),1]<--Inf
      }
      if(sum(is.nan(quad_lin_coeff[i,,2]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,2])),2]<-Inf
      }
      if(sum(quad_lin_coeff[i,,1]<=-10^8)==0){
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K){
          gamma.next[i,]<-gamma.curr[i,]
        }
      }else if((sum(quad_lin_coeff[i,,1]<=-10^8)>0)&&(sum(quad_lin_coeff[i,,1]<=-10^8)<K)){
        quad_lin_coeff[i,which(quad_lin_coeff[i,,1]<=-10^8),1]<--(10^10)
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K){
          gamma.next[i,]<-gamma.curr[i,]
        }
      }else if(sum(quad_lin_coeff[i,,1]<=-10^8)==K){
        gamma.next[i,]<-gamma.curr[i,]
      }
      #print(i)
    }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma.next<-t(apply(X = gamma.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of pi
  pi.update<-function(gamma.curr){
    pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## normalization of pi
    pi.next<-pi.next/sum(pi.next)
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions K1,K2
  theta.update<-function(theta.curr,gamma,net_adjacency,N,K){
    theta.next<-matrix(NA_real_,K,K)
    gradient<-grad_bipartite_stat_undir_theta(theta=theta.curr, gamma=gamma, net_adjacency=net_adjacency, N=N, K=K)
    hess<-hess_bipartite_stat_undir_theta(theta=theta.curr, gamma=gamma, N=N, K=K)
    for(k in 1:K){
      for(l in 1:K){
        theta.next[k,l]<-theta.curr[k,l]-gradient[k,l]/hess[k,l]
      }
    }
    return(theta.next)
  }
  
  #################################################
  ## Defining the update function of full proportional odds nuisance parameters delta_0
  delta_0.update<-function(gamma,delta_0,logitInv_delta_0,net_adjacency,net_rating,N,K,R){
    delta_0.next<-array(NA_real_,dim=c(K,K,(R-1)))
    gradient_delta_0<-grad_bipartite_stat_undir_delta_0(gamma=gamma, delta_0=delta_0, logitInv_delta_0=logitInv_delta_0, net_adjacency=net_adjacency, net_rating=net_rating, N=N, K=K, R=R)
    for (k in 1:K){
      for (l in 1:K){
        hess_mat<-hess_bipartite_stat_undir_delta_0(k=k-1, l=l-1, gamma=gamma, delta_0=delta_0, logitInv_delta_0=logitInv_delta_0, net_adjacency=net_adjacency, net_rating=net_rating, N=N, R=R)
        
        try_delta_0_update<-try(delta_0.next[k,l,]<-sort(delta_0[k,l,]-as.vector((solve(hess_mat+(10^(-6)*diag(R-1)))%*%gradient_delta_0[k,l,]))))
        if('numeric' %in% class(try_delta_0_update)){
          delta_0.next[k,l,]<-try_delta_0_update
        }
        else{
          delta_0.next[k,l,]<-delta_0[k,l,]
        }
        
        #delta_0.next[k,l,]<-sort(delta_0[k,l,]-as.vector((solve(hess_mat+(10^(-6)*diag(R-1)))%*%gradient_delta_0[k,l,])))
        #delta_0.next[k,l,]<-sort(delta_0[k,l,]-as.vector((solve(hess_mat)%*%gradient_delta_0[k,l,])))
      }
    }
    #print(delta_0.next)
    return(delta_0.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,net_adjacency,net_rating,K,n_iter,thres,n_iter_min,n_iter_max,R){
    ## Assuming a network adjacency matrix of dimensions N1*N2 with rows corresponding to users and columns corresponding to products and each entry is an indicator whether the user rated the product. This framework can be generalized.
    N<-dim(net_adjacency)[1]
    
    ## initializing the arrays for parameters
    ## Initializing the cluster memberships gamma
    gamma<-array(NA_real_,dim=c(N,K,n_iter)) ## assuming there are K clusters in the single mode network
    
    ## Initializing the mixture proportions pi
    pi<-matrix(NA_real_,K,n_iter)
   
    ## Initializing the network parameters theta
    theta<-array(NA_real_,dim=c(K,K,n_iter))
    
    ## Initializing the proportional odds model intercept parameters
    delta_0<-array(NA_real_,dim=c(K,K,R-1,n_iter))
    
    gamma[,,1]<-start[[1]]
    pi[,1]<-start[[2]]
    theta[,,1]<-start[[3]]
    delta_0[,,,1]<-start[[4]]
    
    #print(theta[,,1])
    #print(omega[,1])
    #print(mu[,,,1])
    
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_curr<-10^10
    time_per_iter<-c()

    while(((error>thres)&(iter_index<n_iter_max))|(iter_index<n_iter_min)){ 
      #while((error>thres|iter_index<200)&(iter_index<250)&(iter_index==2|error<10)){ 
      #while(iter_index<21){ 
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      
      ## Calculating the logit inverse of delta_0
      logitInv_delta_0.curr<-logitInv_delta_0_cal(delta_0=delta_0[,,,iter_index-1], K=K, R=R)
      
      ## Transforming logit inverse of delta_0 to phi
      #phi.curr<-phi_cal(logitInv_delta_0=logitInv_delta_0.curr, K1=K1, K2=K2, R=R)
      
      ## Calculating the logit inverse of phi
      #logitInv_phi.curr<-logitInv_phi_cal(phi=phi.curr, K1=K1, K2=K2, R=R)
      
      ## Updating gamma_U, the mixed memberships for users
      gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr=gamma[,,iter_index-1], pi.curr=pi[,iter_index-1], theta.curr=theta[,,iter_index-1], logitInv_delta_0.curr=logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating=net_rating, N=N,K=K,R=R)
      ## Updating the mixture proportions pi_U
      pi[,iter_index]<-pi.update(gamma.curr=gamma[,,iter_index])
      
      ## Updating the network parameters theta
      theta[,,iter_index]<-theta.update(theta.curr=theta[,,iter_index-1],gamma=gamma[,,iter_index], net_adjacency=net_adjacency, N=N, K=K)
      #print('theta')
      #print(theta[,,iter_index])
      
      ## Updating the omega matrix
      delta_0[,,,iter_index]<-delta_0.update(gamma=gamma[,,iter_index], delta_0=delta_0[,,,iter_index-1],logitInv_delta_0 = logitInv_delta_0.curr,net_adjacency=net_adjacency, net_rating=net_rating,N=N, K=K, R=R)
      #print('delta_0')
      #print(delta_0[,,,iter_index])
      #Sys.sleep(1.5)
      
      ELBO_prev<-ELBO_curr
      ELBO_curr<-ELBO_conv_bipartite_stat_undir(gamma=gamma[,,iter_index], pi=pi[,iter_index], theta = theta[,,iter_index], logitInv_delta_0=logitInv_delta_0.curr, net_adjacency=net_adjacency, net_rating=net_rating,N=N,K=K)
      
      #error<-abs(((ELBO_prev-ELBO_curr)/ELBO_prev))
      
      if(!is.nan(abs(((ELBO_prev-ELBO_curr)/ELBO_prev)))){
        error<-abs(((ELBO_prev-ELBO_curr)/ELBO_prev))
      }
      
      #print('error')
      #print(error)
      
      #print("probs")
      #print(Gvals_cube_0)
      # print(apply(X = Gvals_cube_0,MARGIN = c(1,2),FUN = function(x){
      #   
      #   return(x/sum(x))
      # }))
      
      if((iter_index%%1)==0){
        #print(theta[,,iter_index])
        print('iter_index')
        print(iter_index)
        print('error')
        print(error)
        print(proc.time()-ptm)
      }
      ptm_diff<-proc.time()-ptm
      time_per_iter<-c(time_per_iter,ptm_diff[3])
      iter_index<-iter_index+1
    }
    names(time_per_iter)<-NULL
    total_iterations<-iter_index-1
    return(list(gamma,pi,theta,delta_0,time_per_iter,total_iterations))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining Model Selection functions based on converged paramters
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI_cal<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  ########################################################################################################
  ## Defining the function to evaluate clustering performance using normalized mutual information.
  # NMI_cal<-function(cluster_ids_true_U,cluster_ids_true_P,K_true_U,K_true_P,cluster_ids_est_U,cluster_ids_est_P,K1,K2,N1,N2){
  #   probab_vec_true_U<-probab_clust_cal(cluster_ids = cluster_ids_true_U,N = N1+N2,K = K_true_U)
  #   probab_vec_true_P<-probab_clust_cal(cluster_ids = cluster_ids_true_P,N = N1+N2,K = K_true_P)
  #   probab_vec_true<-c(probab_vec_true_U,probab_vec_true_P)
  #   total_entropy_true<-0
  #   for(k in 1:length(probab_vec_true)){
  #     if(probab_vec_true[k]!=0){
  #       total_entropy_true<-total_entropy_true-(probab_vec_true[k]*log(probab_vec_true[k],base = 2))
  #     }
  #   }
  #   
  #   probab_vec_est_U<-probab_clust_cal(cluster_ids = cluster_ids_est_U,N = N1+N2,K = K1)
  #   probab_vec_est_P<-probab_clust_cal(cluster_ids = cluster_ids_est_P,N = N1+N2,K = K2)
  #   probab_vec_est<-c(probab_vec_est_U,probab_vec_est_P)
  #   total_entropy_est<-0
  #   for(k in 1:length(probab_vec_est)){
  #     if(probab_vec_est[k]!=0){
  #       total_entropy_est<-total_entropy_est-(probab_vec_est[k]*log(probab_vec_est[k],base = 2))
  #     }
  #   }
  #   
  #   clust_ids_true_total<-c(cluster_ids_true_U,(cluster_ids_true_P+K_true_U))
  #   clust_ids_est_total<-c(cluster_ids_est_U,(cluster_ids_est_P+K1))
  #   
  #   conditional_entropy<-rep(NA_real_,(K1+K2))
  #   for(k in 1:(K1+K2)){
  #     probab_conditional_vec<-rep(0,(K_true_U+K_true_P))
  #     for(m in 1:(K_true_U+K_true_P)){
  #       if(sum(clust_ids_est_total==k)!=0){
  #         probab_conditional_vec[m]<-sum((clust_ids_true_total==m)&(clust_ids_est_total==k))/sum(clust_ids_est_total==k)
  #       }
  #     }
  #     plogp_conditional_vec<-rep(0,length(probab_conditional_vec))
  #     for(m in 1:length(probab_conditional_vec)){
  #       if(probab_conditional_vec[m]!=0){
  #         plogp_conditional_vec[m]<-probab_conditional_vec[m]*log(x = probab_conditional_vec[m],base = 2)
  #       }
  #     }
  #     
  #     conditional_entropy[k]<--probab_vec_est[k]*sum(plogp_conditional_vec)
  #   }
  #   
  #   mutual_info<-total_entropy_true-sum(conditional_entropy)
  #   NMI_val<-(2*mutual_info)/(total_entropy_true+total_entropy_est)
  #   return(NMI_val)
  # }
  ## Defining the functions to evaluate clustering performance using normalized mutual information.
  probab_clust_cal<-function(cluster_ids,N,K){
    probab_vec<-rep(NA_real_,K)
    for(k in 1:K){
      probab_vec[k]<-sum(cluster_ids==k)/N
    }
    return(probab_vec)
  }
  NMI_cal<-function(cluster_ids_true_U,cluster_ids_true_P,K_true_U,K_true_P,cluster_ids_est_U,cluster_ids_est_P,K1,K2,N1,N2){
    probab_vec_true_U<-probab_clust_cal(cluster_ids = cluster_ids_true_U,N = N1+N2,K = K_true_U)
    probab_vec_true_P<-probab_clust_cal(cluster_ids = cluster_ids_true_P,N = N1+N2,K = K_true_P)
    probab_vec_true<-c(probab_vec_true_U,probab_vec_true_P)
    total_entropy_true<-0
    for(k in 1:length(probab_vec_true)){
      if(probab_vec_true[k]!=0){
        total_entropy_true<-total_entropy_true-(probab_vec_true[k]*log(probab_vec_true[k],base = 2))
      }
    }
    
    probab_vec_est_U<-probab_clust_cal(cluster_ids = cluster_ids_est_U,N = N1+N2,K = K1)
    probab_vec_est_P<-probab_clust_cal(cluster_ids = cluster_ids_est_P,N = N1+N2,K = K2)
    probab_vec_est<-c(probab_vec_est_U,probab_vec_est_P)
    total_entropy_est<-0
    for(k in 1:length(probab_vec_est)){
      if(probab_vec_est[k]!=0){
        total_entropy_est<-total_entropy_est-(probab_vec_est[k]*log(probab_vec_est[k],base = 2))
      }
    }
    
    clust_ids_true_total<-c(cluster_ids_true_U,(cluster_ids_true_P+K_true_U))
    clust_ids_est_total<-c(cluster_ids_est_U,(cluster_ids_est_P+K1))
    
    conditional_entropy<-rep(NA_real_,(K1+K2))
    for(k in 1:(K1+K2)){
      probab_conditional_vec<-rep(0,(K_true_U+K_true_P))
      for(m in 1:(K_true_U+K_true_P)){
        if(sum(clust_ids_est_total==k)!=0){
          probab_conditional_vec[m]<-sum((clust_ids_true_total==m)&(clust_ids_est_total==k))/sum(clust_ids_est_total==k)
        }
      }
      plogp_conditional_vec<-rep(0,length(probab_conditional_vec))
      for(m in 1:length(probab_conditional_vec)){
        if(probab_conditional_vec[m]!=0){
          plogp_conditional_vec[m]<-probab_conditional_vec[m]*log(x = probab_conditional_vec[m],base = 2)
        }
      }
      
      conditional_entropy[k]<--probab_vec_est[k]*sum(plogp_conditional_vec)
    }
    
    mutual_info<-total_entropy_true-sum(conditional_entropy)
    NMI_val<-(2*mutual_info)/(total_entropy_true+total_entropy_est)
    return(NMI_val)
  }
  
  #################################################
  ## RASE functions
  RASE_cal<-function(param_est, param_true){
    if(is.vector(param_est)!=1){
      RASE_val<-sqrt(sum((param_est-param_true)^2)/prod(dim(param_est)))
    }else{RASE_val<-sqrt(sum((param_est-param_true)^2)/length(param_est))}
    return(RASE_val)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N<-dim(net_adjacency)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  K<-nclust ## Defining the number of clusters
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  start[[1]]<-matrix(rep(1/K,N*K),N,K)
  start[[2]]<-rep(1/K,K)
  start[[3]]<-theta_init
  start[[4]]<-delta_0_init
  if(K!=1){ param<-iterator(start=start, net_adjacency = net_adjacency, net_rating=net_rating, K=K, n_iter=1000, thres=thres, n_iter_min=n_iter_min, n_iter_max=n_iter_max, R=R)}
  
  #################################################
  ## extracting the coverged parameter values and calculating ICL
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    temp<-is.na(param[[1]][1,1,n_iter])
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  #print(param)
  
  if(K!=1){
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,,n_last]
    param_converge[[4]]<-param[[4]][,,,n_last]
    output_list<-param_converge
    cluster_ids_est<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[1]][x,])
      return(cluster_id)
    }))
    #ICL_vec<-ICL_cal(cluster_ids_est_U = cluster_ids_est_U, cluster_ids_est_P = cluster_ids_est_P, pi_U=param_converge[[3]], pi_P=param_converge[[4]] ,theta = param_converge[[5]], beta_u = as.matrix(param_converge[[6]]), beta_p = as.matrix(param_converge[[7]]), delta_0 = param_converge[[8]], net_adjacency=net_adjacency, net_rating=net_rating, cov_u_outer = cov_u_outer, cov_p_outer = cov_p_outer, N1 = N1, N2 = N2, K1 = K1, K2 = K2, R = R)
    ICL_vec<-NA
    if(sim_indicator==1){
      RI_val_U<-RI_cal(cluster_ids_est = cluster_ids_est[1:N1],cluster_ids_true = cluster_ids_true_U)
      RI_val_P<-RI_cal(cluster_ids_est = cluster_ids_est[(N1+1):(N1+N2)],cluster_ids_true = cluster_ids_true_P)
      NMI_val<-NMI_cal(cluster_ids_true_U = cluster_ids_true_U, cluster_ids_true_P = cluster_ids_true_P, K_true_U = K_true, K_true_P = K_true, cluster_ids_est_U = cluster_ids_est[1:N1], cluster_ids_est_P = cluster_ids_est[(N1+1):(N1+N2)], K1 = K1, K2 = K2, N1 = N1, N2 = N2)
      if(K==K_true){
        K_permute_mat<-do.call(rbind,permn(1:K))
        RASE_theta_mat<-matrix(NA_real_,nrow(K_permute_mat),nrow(K_permute_mat))
        for (k in 1:nrow(K_permute_mat)){
          for (l in 1:nrow(K_permute_mat)){
            theta_true_rotated<-as(as.integer(K_permute_mat[k,]), "pMatrix")%*%theta_true%*%t(as(as.integer(K_permute_mat[l,]), "pMatrix"))
            RASE_theta_mat[k,l]<-RASE_cal(param_est = param_converge[[3]], param_true = theta_true_rotated)
          }
        }
        
        RASE_theta_val<-min(RASE_theta_mat)
        rotate_indices<-which(RASE_theta_mat == min(RASE_theta_mat), arr.ind = TRUE)
        
        delta_0_true_rotated<-array(NA_real_,dim=c(K,K,(R-1)))
        for (r in 1:(R-1)){
          delta_0_true_rotated[,,r]<-as(as.integer(K_permute_mat[rotate_indices[1],]), "pMatrix")%*%delta_0_true[,,r]%*%t(as(as.integer(K_permute_mat[rotate_indices[2],]), "pMatrix"))
          
        }
        RASE_delta_0_val<-RASE_cal(param_est = param_converge[[4]], param_true = delta_0_true_rotated)
        
        output_list<-list("parameters_converged"=param_converge,"Estimated_cluster_IDs"=cluster_ids_est,"ICL"=ICL_vec,"Set1_Rand_Index"=RI_val_U,"Set2_Rand_Index"=RI_val_P, "Normalized_Mutual_Info"=NMI_val, "RASE_theta"=RASE_theta_val, "RASE_delta_0"=RASE_delta_0_val,"time_per_iter_vec"=param[[5]], "total_iterations"=param[[6]])
      }else{
        output_list<-list("parameters_converged"=param_converge,"Estimated_cluster_IDs"=cluster_ids_est,"ICL"=ICL_vec,"Set1_Rand_Index"=RI_val_U,"Set2_Rand_Index"=RI_val_P, "Normalized_Mutual_Info"=NMI_val, "time_per_iter_vec"=param[[5]], "total_iterations"=param[[6]])
      }
    }else{output_list<-list("parameters_converged"=param_converge,"Estimated_cluster_IDs"=cluster_ids_est,"ICL"=ICL_vec,"time_per_iter_vec"=param[[5]], "total_iterations"=param[[6]])}
    #}
  }
  return(output_list)
}

simulate_network_bipartite_model11<-function(N,pi,theta,delta_0,K,R){
  library(Rcpp)
  sourceCpp(file = paste0(file_path,"model_11.cpp"))
  #library(Matrix)
  net_adjacency<-matrix(0,N,N)
  net_rating<-matrix(0,N,N)
  #net_adjacency<-Matrix(0,N1,N2,sparse = T)
  #net_rating<-Matrix(0,N1,N2,sparse = T)
  if(K==1){
    clust_id_U_sampling<-rmultinom(n = N1, size = 1, prob = pi_U)
    clust_id_P_sampling<-rmultinom(n = N2, size = 1, prob = pi_P)
    clust_ids_U<-apply(X = clust_id_U_sampling,MARGIN = 2,FUN = function(x){
      y<-which(x==1)
      return(y)
    })
    clust_ids_P<-apply(X = clust_id_P_sampling,MARGIN = 2,FUN = function(x){
      y<-which(x==1)
      return(y)
    })
    for (i in 1:N1){
      for (j in 1:N2){
        exp_val<-exp(theta)
        probab<-exp_val/(1+exp_val)
        edge_sample<-rbinom(n = 1, size = 1, prob = probab)
        if(edge_sample==0){
          network[i,j]<-0
        }else if(edge_sample==1){
          G_vec<-rep(NA_real_,R)
          for (r in 1:R){
            G_vec[r]<-exp(mu[r]*omega-(mu[r])^2/2)
          }
          probab_discrete<-G_vec/sum(G_vec)
          network[i,j] <- which(as.vector(rmultinom(n = 1, size = 1, prob = probab_discrete))==1)
        }
      }
    }
  }
  else{
    clust_id_sampling<-rmultinom(n = N, size = 1, prob = pi)
    clust_ids<-apply(X = clust_id_sampling,MARGIN = 2,FUN = function(x){
      y<-which(x==1)
      return(y)
    })
    
    logitInv_delta_0<-logitInv_delta_0_cal(delta_0=delta_0, K=K, R=R)
    probab_block_ratings<-array(NA_real_,dim=c(K,K,R))
    for(k in 1:K){
      for(l in 1:K){
        for(r in 1:R){
          if(r==1){
            probab_block_ratings[k,l,r]<-logitInv_delta_0[k,l,r]
          }else if(r==R){
            probab_block_ratings[k,l,r]<-1-logitInv_delta_0[k,l,(R-1)]
          }else{
            probab_block_ratings[k,l,r]<-logitInv_delta_0[k,l,r]-logitInv_delta_0[k,l,(r-1)]
          }
        }
      }
    }
    
    for (i in 1:(N-1)){
      for (j in (i+1):N){
        exp_val<-exp(theta[clust_ids[i],clust_ids[j]])
        probab<-exp_val/(1+exp_val)
        edge_sample<-rbinom(n = 1, size = 1, prob = probab)
        if(edge_sample==1){
          net_adjacency[i,j]<-1
          net_rating[i,j] <- which(as.vector(rmultinom(n = 1, size = 1, prob = probab_block_ratings[clust_ids[i],clust_ids[j],]))==1)
        }
      }
    }
  }
  return(list(net_adjacency,net_rating,clust_ids))
}

########################################################################################################

