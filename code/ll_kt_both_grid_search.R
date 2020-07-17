#altogether, as a function to source
#for testing
#i=1
#j=1

#EQUATION 1
LL_kt_grid <- function(k, theta, Assignments, Sampled_reefs, Distances, Reef_sizes, Adult_sample_proportions){

#Equation 1: Dispersal Kernel estimate
#k <- par[1]
#theta <- par[2]
d_pred <- function(par, Distances){ 
    #k <- par[1]
    #theta <- par[2]
    d <- Distances #distance matrix between patches
    disp <- exp(k)*theta*exp(-(exp(k)*d)^theta)/gamma(1/theta)
    return(disp)
}


#Grab some parameters from the inputs
NumReefs = nrow(Centroids)
NumSampledReefs = length(Sampled_reefs)
#For each reef, create the expected proportions between every reef pair
Proportions = d_pred(par, Distances)

# Inflate the values by the total area of the source reef
# (larger reefs send out more larvae, all else being equal)

Settlers = Proportions*(matrix(Reef_sizes, nrow = nrow(Distances), ncol=nrow(Distances))) #see repmat: https://www.mathworks.com/help/matlab/ref/repmat.html 
 
 
#Sum the total settlers on each reef
AllSettlers = colSums(Settlers)
AssignedSettlers = matrix(0, nrow=NumSampledReefs+1,ncol=NumSampledReefs)


#Equation 2: Account for postsettlement mortality, the proportion of settlers that come from adults on patch i **THINK MORE ABOUT THIS
### SOME THESE PROPORTIONS ARE GREATER THAN ONE, CHECK THIS. Likely because of problems before this in the code

for(i in 1:NumSampledReefs){
   This_SS_A = Adult_sample_proportions[i]
   for(j in 1:NumSampledReefs){
      SettlersFromAssignedReefs = Settlers[Sampled_reefs[i],Sampled_reefs[j]]
    #Not all settlers from assigned reefs will be assigned, because not all adults were sampled
    AssignedSettlers[i,j] = SettlersFromAssignedReefs*(This_SS_A^2 + 2*This_SS_A*(1 - This_SS_A))
    AssignedSettlers[NumSampledReefs+1,j] = AssignedSettlers[NumSampledReefs+1,j] + SettlersFromAssignedReefs*(1-This_SS_A)^2 #The three dots '...' tell matlab that the code on a given line continues on the next line.
   }
}
Unsampled = as.matrix(setdiff(1:NumReefs,Sampled_reefs))


for(j in 1:length(Sampled_reefs)){
    AssignedSettlers[NumSampledReefs+1,j] =  sum(Settlers[Unsampled,Sampled_reefs[j]]) + AssignedSettlers[NumSampledReefs+1,j] 
}
   

#Normalise values into multinomial probabilities
PredictedProportions = AssignedSettlers/(matrix(rep(t(colSums(AssignedSettlers)), NumSampledReefs+1 ), ncol = ncol(AssignedSettlers) , byrow = TRUE ))

colSums(PredictedProportions) #everything should be to the -e power because these are probabilities and the columns should sum to 1

#Loglikelihoods can't handle zeros. Make them very small instead.
PredictedProportions[PredictedProportions==0] = 1e-12
PredictedProportions[PredictedProportions<=1e-12] = 1e-12

#Re-normalise into probabilities
PredictedProportions = PredictedProportions/(matrix(rep(t(colSums(PredictedProportions)), nrow(PredictedProportions) ), ncol = ncol(PredictedProportions) , byrow = TRUE )) #this feels overengineered...could have just made a matrix of ones. whatever.



#initialize
log_like =0

for(j in 1:NumSampledReefs){ 
   # Go through reefs one-by-one: What is the LL of each reef's observations?

   ObsVector = Assignments[,j] # explanation of : A single : in a subscript position is shorthand notation for 1:end and is often used to select entire rows or columns:
   ProbVector = PredictedProportions[,j]
   
   #Calculate the log likelihood for this particular reef, given the sample that's been observed
    log_like =  log_like + sum(ObsVector*log(ProbVector))

    }
        neg_ll = -log_like

    return(neg_ll)
}