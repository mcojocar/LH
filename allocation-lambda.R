### 1. Load necessary packages ###



### 2. Establish data ####

demand <- c(0.2,0.2,0.2,0.2,0.2)

# Note everything is scaled by 10^(-3)
group_sizes <- 10^(-3)*c(152930,148584,145840,146319,146037)
price <-c(200,300)
prob_outcome <- c(0.385,0.391,0.396,0.402,0.408)
budget <- 68000000*10^(-3)
overhead_costs <- 0 #10^(-3)*budget*0.3 # vaccine programs require $$ for administration and cold chain storage (need to cite)
#budget <- budget - overhead_costs

e1=  c(0.6,0.62,0.62,0.6,0.56) # zostavax efficiency
e2= c(0.9,0.9,0.89,0.89,0.88) # shingrix efficiency
r.names = c("G1","G2","G3","G4","G5") 
c.names = c("Zostavax","Shingrix")
m.names = c("efficacy")

arr = array (c (e1,e2), dim=c (5,2,1), dimnames=list (r.names, c.names, m.names)) # make efficacy matrix
print(arr)


##################### Initialize and establish Groups #########################

G1 <- c(group_sizes[1],demand[1]*group_sizes[1])
G2 <- c(group_sizes[2],demand[2]*group_sizes[2])
G3 <- c(group_sizes[3],demand[3]*group_sizes[3])
G4 <- c(group_sizes[4],demand[4]*group_sizes[4])
G5 <- c(group_sizes[5],demand[5]*group_sizes[5])

G.data <- data.frame(G1,G2,G3,G4,G5)
#rownames(G.data)[1] <- "demand"
rownames(G.data)[1] <- "population"
rownames(G.data)[2] <- "estimated-demand"

G.data.shuffled <- G.data[,sample(ncol(G.data))] # shuffling columns into random order
#rownames(G.data.shuffled)[1] <- "demand"
rownames(G.data.shuffled)[1] <- "population"
rownames(G.data.shuffled)[2] <- "estimated-demand"

shuffled <- colnames(G.data.shuffled) # store group ordering as vector

# # transpose dataframes for easier accessing
# G.data.trans <- as.data.frame(t(G.data))
# G.data.shuffled.trans <- as.data.frame(t(G.data.shuffled))

cat("The new group ordering is: ", shuffled)


##################### 3. Allocation functions #########################

zostavax.allocation.year1 <- function() {

  #shuffled <- colnames(G.data.shuffled)
  dose_total = (budget - overhead_costs)/price[1]
  dose_initial <- dose_total
  
  dose_demand_total = 0 # initialize a counter for how many doses we estimate will be demanded by 65-70s
  max_allocations <- vector(mode="double", length=5) # create a vector to store the max number of allocated vaccines in each group
  first_allocation <- vector(mode="double", length=5)
  dd_allocation <- vector(mode="double", length=5)
  
  for (i in 1:5) {
    dose_demand_total = dose_demand_total + demand[i]*group_sizes[i] # sum
  }
  
  cat("Our budget will allow for", floor(dose_total*10^3), "doses of Zostavax. We estimate needing", dose_demand_total*10^3, ".\n")
  
        if (dose_total < dose_demand_total) {
             cat("We will have a shortage of", abs(dose_total - dose_demand_total)*10^3, "doses of Zostavax.\n")
        } else if (dose_total > dose_demand_total) {
             cat("We will have an excess of", abs(dose_total - dose_demand_total)*10^3, "doses of Zostavax.\n")
        } else if (dose_total == dose_demand_total) {
             cat("We have exactly enough doses of Zostavax.\n")
        } else {
             cat("There was a comparison error.\n")
        }
  

  lambda <- vector(mode="double", length=5)
  arrival_rate <- vector(mode="double", length=5)
  
  
    for (i in 1:5) {
      cat("\n", "----> Currently allocating in ", shuffled[i])
      
      w <- c(min(G.data.shuffled[[i]][2],dose_total)) # take whichever is smaller, the pop size or the number of doses left
      decimal_doses <- floor((w - floor(w))*10^3)
      
      cat("\n The number of Zostavax doses we can distribute to this group is w = : ", w*10^3 , "\n")
      cat("\n By our process, we first distribute floor(w) = : ", floor(w)*10^3 , "doses. \n")
      cat("\n After, we will have to account for", decimal_doses, "doses from the truncated decimal. \n")
      w <- floor(w)
      
  # compute Poisson arrival likelihood
      lambda[i] <- (1-0.2*runif(1,0,1)*(G.data.shuffled[2, i]))
      arrival_rate[i] <- 0.368 #(lambda[i]^(G.data.shuffled[2, i])*exp(-lambda[i]))/factorial(G.data.shuffled[2, i])
      
      cat("The arrival likelihood in", shuffled[i], "is", arrival_rate[i], ".\n")

      
      objective_function <- function(i,j) {
            rev <- (e1[i] * xt[j]) * prob_outcome[i] * arrival_rate[i]
            return(rev)
      }
      
      objective_function2 <- function(i,j) {
        rev <- (e1[i] * xt_dd[j]) * prob_outcome[i] * arrival_rate[i]
        return(rev)
      }
    
      
  #For the first floor(w) doses :
      xt <- vector(mode="double", length=w) # xt has has many elements as the value of w
        for (k in 1:w){
          xt[k] <- k   # fill xt
        }
      
  #For the decimal doses due to floor()  :  
      xt_dd <- vector(mode="double", length=decimal_doses) # xt_dd has has many elements as the value of decimal_doses
        for (k in 1:decimal_doses){
          xt_dd[k] <- k   # fill xt_dd
        }


  for (j in 1:w) {
          xt[j] <- objective_function(i,j) # calculate rev for this (i,j) via objective_function and store it in column j of xt
  }
  
  for (j in 1:decimal_doses) {
        xt_dd[j] <- objective_function2(i,j)
  }
      
  # The index of the max of the objective function corresponds to maxnumber of doses that are allocatable
      first_allocation[i] <- match(max(xt), xt)
        #first_allocation[i]
      dd_allocation[i] <- match(max(xt_dd), xt_dd)
        #dd_allocation[i]
      
  # Combine both sets and we get the total amount of doses we will allocate
      max_allocations[i] <- (floor(first_allocation[i]*10^3) + dd_allocation[i])*10^(-3)
        #max_allocations[i]

      
      cat("Number of doses allocated in ", shuffled[i], " is ", floor(max_allocations[i]*10^3), "\n")

      #### Remove the doses we just allocated and update new total:

      if (dose_total - max_allocations[i] >= 0) {
        dose_total <- dose_total - max_allocations[i]
        cat("Remaining number of doses to be allocated is ", dose_total, ". There are ", 5-i, "allocations left \n")
      } else {
          print("We ran out of doses in", G.data.shuffled[i], ". \n")
          break
      }
}

  ### end of for loop
  

    cat("Final vaccine allocations per group: \n", max_allocations, "\n")

    if (dose_total > 0) {
      cat("There are", floor(dose_total*10^3), "vaccine doses left. \n We estimated needing", dose_demand_total*10^3, "doses this year. \n We distributed", floor((dose_initial-dose_total)*10^3), "doses of Zostavax this year =", 100*((dose_initial-dose_total)/dose_initial), "% of stock //", 100*((dose_initial-dose_total)/dose_demand_total),"% of estimated demand. \n") 

    } else {
        cat("There are no vaccine doses left. \n")
    }
    
    
######## Coverage estimates
    
    cat("Based on the Zostavax vaccines administered, we estimate the following vaccine protection levels for 65-70 yos against shingles in Ontario: \n")
}

zostavax.allocation.year1()

# zostavax.allocation.year2 <- function() {
#   
#   z_doses = (budget - overhead_costs)/price[1]
#   
#   # DON'T FORGET TO CLEAR (SOME) OF THE PREVIOUS FUNCTIONS
#   
# }
      
# ###### Different scenarios - compare coverage and demand!
      # 1) more uptake from top of age window 
      # 2) more uptake from bottom of age window
      # 3) revised budget

#--end of file--