
#############################################
## Longitudinal IRT models 
#############################################

library(flirt)


long <- data(longitudinal)
male <- long[,1]
bin3 <- long[,2:22]

###############################
## Andersen's model 
## a theta for time point 
###############################



## no measurement invariance 
## all items are freely estimated 
model1 <- flirt(data=bin3, mul=list(on=T, dim_info=list(dim1=1:7, dim2=8:14, dim3=15:21)) , 
            control=list(nq=5)  )


# correlation 

est_cov <- model1@pars[1:6,1]  # vector of cov estimates 
cov_matrix <- matrix(c(est_cov[1]^2,est_cov[4],est_cov[5], est_cov[4], est_cov[2]^2,est_cov[6],est_cov[5],est_cov[6],est_cov[3]^2),3,3, byrow=F) # covariance matrix 
dim_info <- list(dim1=1:7, dim2=8:14, dim3=15:21) 

test1 <- std_cov(dim_info=dim_info, cov_matrix= cov_matrix)     


## measurement invariance 

I <- 7
m <- diag(I) 
comm <- rbind(m,m,m)  
colnames(comm) <- paste("i",1:I)
w2<- c(rep(0,I), rep(1,I), rep(0,I))
w3<- c(rep(0,I), rep(0,I), rep(1,I))

item_design <- cbind(comm, w2, w3)

## same item parameters across time 
model2 <- flirt(data=bin3, mul=list(on=T, dim_info=list(dim1=1:7, dim2=8:14, dim3=15:21)) , 
            item_cov=list(on=T, item_matrix_beta=comm), 
            control=list(nq=5)  )

# model 1 (non measurement invariance) fits better 
> anova(model2, model1)

Chi-squared difference: 
X2 = 22.089, df = 14, p = 0.0768
AIC difference =  -5.911 
BIC difference =  -52.087 


# plus time specific effects 
model3 <- flirt(data=bin3, mul=list(on=T, dim_info=list(dim1=1:7, dim2=8:14, dim3=15:21)) , 
            item_cov=list(on=T, item_matrix_beta=item_design), 
            control=list(nq=5)  )


## person covariate 

male <- person -1 
model4 <- flirt(data=bin3, mul=list(on=T, dim_info=list(dim1=1:7, dim2=8:14, dim3=15:21), cov_info=list(dim1=1, dim2=1, dim3=1)) , 
            item_cov=list(on=T, item_matrix_beta=item_design), 
            person_cov=list(on=T, person_matrix=male), 
            control=list(nq=5)  )
           
            
################################            
### Embretson's change model 
################################            

# using a simplex structure 

model5 <- flirt(data=bin3, mul=list(on=T, dim_info=list(dim1=1:21, dim2=8:21, dim3=15:21)) , 
            item_cov=list(on=T, item_matrix_beta=comm), 
            control=list(nq=5)  )
            

# correlation 

est_cov <- model5@pars[1:6,1]  # vector of cov estimates 
cov_matrix <- matrix(c(est_cov[1]^2,est_cov[4],est_cov[5], est_cov[4], est_cov[2]^2,est_cov[6],est_cov[5],est_cov[6],est_cov[3]^2),3,3, byrow=F) # covariance matrix 
dim_info <- list(dim1=1:21, dim2=8:21, dim3=15:21) 

test2 <- std_cov(dim_info=dim_info, cov_matrix= cov_matrix)     



###############################################             
### Bifactor model  / second-order model 
###############################################        

I <- 7
m <- diag(I) 
comm <- rbind(m,m,m)  
colnames(comm) <- paste("i",1:I)


# to obtain the general ability measure 
model6 <- flirt(data=bin3, bifac=list(on=T, dim_info=list(dim1=1:7, dim2=8:14, dim3=15:21)) , 
            item_cov=list(on=T, item_matrix_beta=comm, item_matrix_alpha=comm), 
            control=list(conv=0.1, nq=5), post=T )
            

model7 <- flirt(data=bin3, second=list(on=T, dim_info=list(dim1=1:7, dim2=8:14, dim3=15:21)) , 
            item_cov=list(on=T, item_matrix_beta=comm), 
            control=list(conv=0.1, nq=5), post=T )            
            
            
