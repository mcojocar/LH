### 1. Load necessary packages ###
require("ggplot2")
require("data.table")

### 2. Establish data ####

demand <- c(0.2,0.2,0.2,0.2,0.2)

# Note everything is scaled by 10^(-3)
group_sizes_year1 <- 10^(-3)*c(152930,148584,145840,146319,146037) # from 2016 ON data
group_sizes_year2 <- 10^(-3)*c(156367,151966,147573,144607,144905) # from 2017 ON data
group_sizes_year3 <- 10^(-3)*c(162441,155275,150753,146197,143127) # from 2018 ON data

price <-c(200,300)
prob_outcome <- c(0.385,0.391,0.396,0.402,0.408)
budget <- 68000000*10^(-3)
overhead_costs <- 0 #10^(-3)*budget*0.3 # vaccine programs require $$ for administration and cold chain storage (need to cite)
#budget <- budget - overhead_costs
budget_doses = (budget - overhead_costs)/price[1]

e1=  c(0.6,0.62,0.62,0.6,0.56) # zostavax efficiency
e2= c(0.9,0.9,0.89,0.89,0.88) # shingrix efficiency
r.names = c("G1","G2","G3","G4","G5") 
c.names = c("Zostavax","Shingrix")
m.names = c("efficacy")

arr = array (c (e1,e2), dim=c (5,2,1), dimnames=list (r.names, c.names, m.names)) # make efficacy matrix
print(arr)

##################### Initialize and establish Groups #########################

G1 <- c(group_sizes_year1[1],demand[1]*group_sizes_year1[1])
G2 <- c(group_sizes_year1[2],demand[2]*group_sizes_year1[2])
G3 <- c(group_sizes_year1[3],demand[3]*group_sizes_year1[3])
G4 <- c(group_sizes_year1[4],demand[4]*group_sizes_year1[4])
G5 <- c(group_sizes_year1[5],demand[5]*group_sizes_year1[5])

G.data <- data.frame(G1,G2,G3,G4,G5)
#rownames(G.data)[1] <- "demand"
rownames(G.data)[1] <- "population"
rownames(G.data)[2] <- "estimated-demand"

G.data.shuffled <- G.data[,sample(ncol(G.data))] # shuffling columns into random order
#rownames(G.data.shuffled)[1] <- "demand"
rownames(G.data.shuffled)[1] <- "population"
rownames(G.data.shuffled)[2] <- "estimated-demand"

shuffled <- colnames(G.data.shuffled) # store group ordering as vector

cat("The new group ordering is: ", shuffled)


##################### 3. Allocation functions #########################

#################################################################################################################
#################################################################################################################
##################################              YEAR 1            ###############################################
#################################################################################################################
#################################################################################################################

zostavax.allocation.year1 <- function() {

  #shuffled <- colnames(G.data.shuffled)
  dose_total = (budget - overhead_costs)/price[1]
  dose_initial <- dose_total
  
  dose_demand_total = 0 # initialize a counter for how many doses we estimate will be demanded by 65-70s
  max_allocations <- vector(mode="double", length=5) # create a vector to store the max number of allocated vaccines in each group
  first_allocation <- vector(mode="double", length=5)
  dd_allocation <- vector(mode="double", length=5)
  
  for (i in 1:5) {
    dose_demand_total = dose_demand_total + demand[i]*group_sizes_year1[i] # sum
  }
  
  cat("Our budget will allow for", floor(dose_total*10^3), "doses of Zostavax in Year 1. We estimate needing", dose_demand_total*10^3, ".\n")
  
        if (dose_total < dose_demand_total) {
             cat("We will have a shortage of", abs(dose_total - dose_demand_total)*10^3, "doses of Zostavax in Year 1.\n")
        } else if (dose_total > dose_demand_total) {
             cat("We will have an excess of", abs(dose_total - dose_demand_total)*10^3, "doses of Zostavax in Year 1.\n")
        } else if (dose_total == dose_demand_total) {
             cat("We have exactly enough doses of Zostavax.\n")
        } else {
             cat("There was a comparison error.\n")
        }
  

  lambda <- vector(mode="double", length=5)
  arrival_rate <- vector(mode="double", length=5)
  #arrival_rate <- c(0.368,0.4,0.45,0.56,0.57)
  
  
    for (i in 1:5) {
      cat("\n", "----> Currently allocating in ", shuffled[i])
      
      w <- c(min(G.data.shuffled[[i]][2],dose_total)) # take whichever is smaller, the pop size or the number of doses left
      decimal_doses <- floor((w - floor(w))*10^3)
      
      cat("\n The number of Zostavax doses we can distribute to this group is w = : ", w*10^3 , "\n")
      cat("\n By our process, we first distribute floor(w) = : ", floor(w)*10^3 , "doses. \n")
      cat("\n After, we will have to account for", decimal_doses, "doses from the truncated decimal. \n")
      w <- floor(w)
      
  # compute Poisson arrival likelihood in x amount of time
      #lambda[i] <- (1-0.2*runif(1,0,1)*(G.data.shuffled[2, i]))
      lambda[i] <- (G.data.shuffled[2, i]*10^3)/52
      arrival_rate[i] <- (lambda[i]^(G.data.shuffled[2, i])*exp(-lambda[i]))/factorial(G.data.shuffled[2, i])
      
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

      
      cat("\n *** The TOTAL number of doses allocated in ", shuffled[i], " is ", floor(max_allocations[i]*10^3), ".***\n\n")

      #### Remove the doses we just allocated and update new total:

      if (dose_total - max_allocations[i] >= 0) {
        dose_total <- dose_total - max_allocations[i]
        cat("Remaining number of doses allocatable: ", dose_total*10^3, ". There are ", 5-i, "allocations left \n")
      } else {
          print("We ran out of doses here.\n")
          break
      }
}

  ### end of for loop
  

    cat("\n #### Final vaccine allocations per group in Year 1: \n", max_allocations, "####\n")

    if (dose_total > 0) {
      cat("\n There are", floor(dose_total*10^3), "vaccine doses left at the end of Year 1. \n We estimated needing", dose_demand_total*10^3, "doses this year. \n We distributed", floor((dose_initial-dose_total)*10^3), "doses of Zostavax this year =", 100*((dose_initial-dose_total)/dose_initial), "% of stock //", 100*((dose_initial-dose_total)/dose_demand_total),"% of estimated demand. \n") 

    } else {
        cat("\n There are no vaccine doses left. \n")
    }
    
    
######## Coverage estimates

# Combine max_allocations data with other group data
    G.data.shuffled1 <- rbind(G.data.shuffled, max_allocations)
    rownames(G.data.shuffled1)[3] <- "max_allocations1"
    
# unshuffle and save Group data
    G.data <- G.data.shuffled1[,order(colnames(G.data.shuffled1))]
    
# return multiple items to store in global envrironment for the next function
    year1_list <- list("G.data.shuffled1" = G.data.shuffled1, "max_allocations1" = max_allocations, "G.data1" = G.data, "leftover_doses1" = dose_total)
    list2env(year1_list ,.GlobalEnv)
    #return()
    
    cat("\n Based on the Zostavax vaccines administered, we estimate the following vaccine protection levels for 65-70 yos against shingles in Ontario: INCOMPLETE \n\n")
    
##################### 4. Plotting #########################
    
    # year_efficacy1 <- vector(mode="double", length=5)
    # 
    # for (i in 1:5) {
    #   year_efficacy1[i] <- (G.data1[[i]][3]*10^3*e1[i]) / (G.data1[[i]][3]*10^3)
    #   
    # }
    
    # df <- data.frame(Age = c("G1 \n (65-66)","G2 \n (66-67)","G3 \n (67-68)","G4 \n (68-69)","G5 \n (69-70)"), Efficacy = year_efficacy1)
    df2 <- data.frame(Age = c("G1 \n (65-66)","G2 \n (66-67)","G3 \n (67-68)","G4 \n (68-69)","G5 \n (69-70)"), Doses_administered = max_allocations1*10^3)
    
    plot1 <- ggplot(df2, aes(Age,Doses_administered)) +
      geom_col(fill="black") +
      labs(color=" ", title="Zostavax doses administered in Year 1 (2016)", x = "\n Age Group", y = "Doses \n") +
      ylim(0,31000)+
      theme_light()
    
    # plot2 <- ggplot(df, aes(Age,Efficacy)) +
    #   geom_col(fill="blue") +
    #   labs(color=" ", title="Population Coverage with Zostavax in Year 1 (2016)", x = "\n Age Group", y = "Coverage Level \n") +
    #   scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
    #   theme_light()
    
    return(list(plot1))
}

## We pass the following into the next year:
#  - subtract previously vaccinated individuals (NOT covered individuals--even if vaccine isn't effective, they still took a vaccine) from age groups
#  - leftover doses (if any)
#  - budget is irrelevant, we spent it to acquire the max number of doses in year 1
#  
# # New numbers incoming
#   - new 65-66 year olds
#   - new group sizes



#################################################################################################################
#################################################################################################################
##################################              YEAR 2            ###############################################
#################################################################################################################
#################################################################################################################


zostavax.allocation.year2 <- function() {
  
  # subtract previously vaccinated individuals (NOT covered individuals--even if vaccine isn't effective, they still took a vaccine) from age groups
  # REMEMBER WE HAD TO UNSHUFFLE FIRST TO SUBTRACT FROM THE RIGHT GROUPS --> USE G.data1 NOT G.data.shuffled

  for (i in 1:4) {
    group_sizes_year2[i+1] <- (group_sizes_year2[i]*10^3 - G.data1[[i]][3]*10^3) * 10^(-3)
  }
  
  dose_total <- leftover_doses1
  dose_initial <- dose_total
  
  ##################### Initialize and establish Groups #########################
  
  G1 <- c(group_sizes_year2[1],demand[1]*group_sizes_year2[1])
  G2 <- c(group_sizes_year2[2],demand[2]*group_sizes_year2[2])
  G3 <- c(group_sizes_year2[3],demand[3]*group_sizes_year2[3])
  G4 <- c(group_sizes_year2[4],demand[4]*group_sizes_year2[4])
  G5 <- c(group_sizes_year2[5],demand[5]*group_sizes_year2[5])
  
  G.data <- data.frame(G1,G2,G3,G4,G5)
  #rownames(G.data)[1] <- "demand"
  rownames(G.data)[1] <- "population"
  rownames(G.data)[2] <- "estimated-demand"
  
  G.data.shuffled <- G.data[,sample(ncol(G.data))] # shuffling columns into random order
  #rownames(G.data.shuffled)[1] <- "demand"
  rownames(G.data.shuffled)[1] <- "population"
  rownames(G.data.shuffled)[2] <- "estimated-demand"
  
  shuffled <- colnames(G.data.shuffled) # store group ordering as vector
  
  cat("The new group ordering is: ", shuffled)

  # DON'T FORGET TO CLEAR (SOME) OF THE PREVIOUS FUNCTIONS
  
  dose_demand_total = 0 # initialize a counter for how many doses we estimate will be demanded by 65-70s
  max_allocations <- vector(mode="double", length=5) # create a vector to store the max number of allocated vaccines in each group
  first_allocation <- vector(mode="double", length=5)
  dd_allocation <- vector(mode="double", length=5)
  
  for (i in 1:5) {
    dose_demand_total = dose_demand_total + demand[i]*group_sizes_year2[i] # sum
  }
  
  cat("\n Our budget will allow for", floor(dose_total*10^3), "doses of Zostavax in Year 2. We estimate needing", dose_demand_total*10^3, ".\n")
  
  if (dose_total < dose_demand_total) {
    cat("We will have a shortage of", abs(dose_total - dose_demand_total)*10^3, "doses of Zostavax in Year 2.\n")
  } else if (dose_total > dose_demand_total) {
    cat("We will have an excess of", abs(dose_total - dose_demand_total)*10^3, "doses of Zostavax in Year 2.\n")
  } else if (dose_total == dose_demand_total) {
    cat("We have exactly enough doses of Zostavax.\n")
  } else {
    cat("There was a comparison error.\n")
  }
  
  
  lambda <- vector(mode="double", length=5)
  arrival_rate <- vector(mode="double", length=5)
  #arrival_rate <- c(0.368,0.4,0.45,0.56,0.57)
  
  
  for (i in 1:5) {
    cat("\n", "----> Currently allocating in ", shuffled[i])
    
    w <- c(min(G.data.shuffled[[i]][2],dose_total)) # take whichever is smaller, the pop size or the number of doses left
    decimal_doses <- floor((w - floor(w))*10^3)
    
    cat("\n The number of Zostavax doses we can distribute to this group is w = : ", w*10^3 , "\n")
    cat("\n By our process, we first distribute floor(w) = : ", floor(w)*10^3 , "doses. \n")
    cat("\n After, we will have to account for", decimal_doses, "doses from the truncated decimal. \n")
    w <- floor(w)
    
    # compute Poisson arrival likelihood
    lambda[i] <- (G.data.shuffled[2, i]*10^3)/52
    arrival_rate[i] <- (lambda[i]^(G.data.shuffled[2, i])*exp(-lambda[i]))/factorial(G.data.shuffled[2, i])
    
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
    
    
    cat("\n *** The TOTAL number of doses allocated in ", shuffled[i], " is ", floor(max_allocations[i]*10^3), ".***\n\n")
    
    #### Remove the doses we just allocated and update new total:
    
    if (dose_total - max_allocations[i] >= 0) {
      dose_total <- dose_total - max_allocations[i]
      cat("Remaining number of doses allocatable: ", dose_total*10^3, ". There are ", 5-i, "allocations left \n")
    } else {
      print("We ran out of doses here.\n")
      break
    }
  }
  
  ### end of for loop
  
  
  cat("\n ##### Final vaccine allocations per group in Year 2: \n", max_allocations, "#####\n")
  
  if (dose_total > 0) {
    cat("\n There are", floor(dose_total*10^3), "vaccine doses left at the end of Year 2. \n We estimated needing", dose_demand_total*10^3, "doses this year. \n We distributed", floor((dose_initial-dose_total)*10^3), "doses of Zostavax this year =", 100*((dose_initial-dose_total)/dose_initial), "% of remaining stock //", 100*((dose_initial-dose_total)/dose_demand_total),"% of estimated demand. \n") 
    
  } else {
    cat("\n There are no vaccine doses left. \n")
  }
  
  
  ######## Coverage estimates
  
  # Combine max_allocations data with other group data
  G.data.shuffled2 <- rbind(G.data.shuffled, max_allocations)
  rownames(G.data.shuffled2)[3] <- "max_allocations2"
  
  # unshuffle and save Group data
  G.data <- G.data.shuffled2[,order(colnames(G.data.shuffled2))]
  
  # return multiple items to store in global envrironment for the next function
  year2_list <- list("G.data.shuffled2" = G.data.shuffled2, "max_allocations2" = max_allocations, "G.data2" = G.data, "leftover_doses2" = dose_total)
  list2env(year2_list ,.GlobalEnv)
  #return()
  
  cat("\n Based on the Zostavax vaccines administered, we estimate the following vaccine protection levels for 65-70 yos against shingles in Ontario: \n\n")
  
  ##################### 4. Plotting #########################
  
  # year_efficacy2 <- vector(mode="double", length=5)
  # 
  # for (i in 1:5) {
  #   year_efficacy2[i] <- (G.data2[[i]][3]*10^3*e1[i]) / (G.data2[[i]][3]*10^3)
  #   
  # }
  
  # df <- data.frame(Age = c("G1 \n (65-66)","G2 \n (66-67)","G3 \n (67-68)","G4 \n (68-69)","G5 \n (69-70)"), Efficacy = year_efficacy2)
  df2 <- data.frame(Age = c("G1 \n (65-66)","G2 \n (66-67)","G3 \n (67-68)","G4 \n (68-69)","G5 \n (69-70)"), Doses_administered = max_allocations2*10^3)
  
  plot1 <- ggplot(df2, aes(Age,Doses_administered)) +
    geom_col(fill="black") +
    labs(color=" ", title="Zostavax doses administered in Year 2 (2017)", x = "\n Age Group", y = "Doses \n") +
    ylim(0,31000)+
    theme_light()
  
  # plot2 <- ggplot(df, aes(Age,Efficacy)) +
  #   geom_col(fill="blue") +
  #   labs(color=" ", title="Population Coverage with Zostavax in Year 2 (2017)", x = "\n Age Group", y = "Coverage Level \n") +
  #   scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
  #   theme_light()
  
  return(list(plot1))

}

#################################################################################################################
#################################################################################################################
##################################              YEAR 3            ###############################################
#################################################################################################################
#################################################################################################################

zostavax.allocation.year3 <- function() {
  
  # subtract previously vaccinated individuals (NOT covered individuals--even if vaccine isn't effective, they still took a vaccine) from age groups
  # REMEMBER WE HAD TO UNSHUFFLE FIRST TO SUBTRACT FROM THE RIGHT GROUPS --> USE G.data2 NOT G.data.shuffled
  
  for (i in 1:4) {
    group_sizes_year3[i+1] <- (group_sizes_year3[i]*10^3 - G.data2[[i]][3]*10^3) * 10^(-3)
  }
  
  dose_total <- leftover_doses2
  dose_initial <- dose_total
  
  ##################### Initialize and establish Groups #########################
  
  G1 <- c(group_sizes_year3[1],demand[1]*group_sizes_year3[1])
  G2 <- c(group_sizes_year3[2],demand[2]*group_sizes_year3[2])
  G3 <- c(group_sizes_year3[3],demand[3]*group_sizes_year3[3])
  G4 <- c(group_sizes_year3[4],demand[4]*group_sizes_year3[4])
  G5 <- c(group_sizes_year3[5],demand[5]*group_sizes_year3[5])
  
  G.data <- data.frame(G1,G2,G3,G4,G5)
  #rownames(G.data)[1] <- "demand"
  rownames(G.data)[1] <- "population"
  rownames(G.data)[2] <- "estimated-demand"
  
  G.data.shuffled <- G.data[,sample(ncol(G.data))] # shuffling columns into random order
  #rownames(G.data.shuffled)[1] <- "demand"
  rownames(G.data.shuffled)[1] <- "population"
  rownames(G.data.shuffled)[2] <- "estimated-demand"
  
  shuffled <- colnames(G.data.shuffled) # store group ordering as vector
  
  cat("The new group ordering is: ", shuffled)
  
  # DON'T FORGET TO CLEAR (SOME) OF THE PREVIOUS FUNCTIONS
  
  dose_demand_total = 0 # initialize a counter for how many doses we estimate will be demanded by 65-70s
  max_allocations <- vector(mode="double", length=5) # create a vector to store the max number of allocated vaccines in each group
  first_allocation <- vector(mode="double", length=5)
  dd_allocation <- vector(mode="double", length=5)
  
  for (i in 1:5) {
    dose_demand_total = dose_demand_total + demand[i]*group_sizes_year3[i] # sum
  }
  
  cat("\n Our budget will allow for", floor(dose_total*10^3), "doses of Zostavax in Year 3. We estimate needing", dose_demand_total*10^3, ".\n")
  
  if (dose_total < dose_demand_total) {
    cat("We will have a shortage of", abs(dose_total - dose_demand_total)*10^3, "doses of Zostavax in Year 3.\n")
  } else if (dose_total > dose_demand_total) {
    cat("We will have an excess of", abs(dose_total - dose_demand_total)*10^3, "doses of Zostavax in Year 3.\n")
  } else if (dose_total == dose_demand_total) {
    cat("We have exactly enough doses of Zostavax.\n")
  } else {
    cat("There was a comparison error.\n")
  }
  
  
  lambda <- vector(mode="double", length=5)
  arrival_rate <- vector(mode="double", length=5)
  #arrival_rate <- c(0.368,0.4,0.45,0.56,0.57)
  
  
  for (i in 1:5) {
    cat("\n", "----> Currently allocating in ", shuffled[i])
    
    w <- c(min(G.data.shuffled[[i]][2],dose_total)) # take whichever is smaller, the pop size or the number of doses left
    decimal_doses <- floor((w - floor(w))*10^3)
    
    cat("\n The number of Zostavax doses we can distribute to this group is w = : ", w*10^3 , "\n")
    cat("\n By our process, we first distribute floor(w) = : ", floor(w)*10^3 , "doses. \n")
    cat("\n After, we will have to account for", decimal_doses, "doses from the truncated decimal. \n")
    w <- floor(w)
    
    # compute Poisson arrival likelihood
    lambda[i] <- (G.data.shuffled[2, i]*10^3)/52
    arrival_rate[i] <- (lambda[i]^(G.data.shuffled[2, i])*exp(-lambda[i]))/factorial(G.data.shuffled[2, i])
    
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
    
    
    cat("\n *** The TOTAL number of doses allocated in ", shuffled[i], " is ", floor(max_allocations[i]*10^3), ".***\n\n")
    
    #### Remove the doses we just allocated and update new total:
    
    if (dose_total - max_allocations[i] >= 0) {
      dose_total <- dose_total - max_allocations[i]
      cat("Remaining number of doses allocatable: ", dose_total*10^3, ". There are ", 5-i, "allocations left \n")
    } else {
      print("We ran out of doses here.\n")
      break
    }
  }
  
  ### end of for loop
  
  
  cat("\n ##### Final vaccine allocations per group in Year 3: \n", max_allocations, "#####\n")
  
  if (dose_total > 0) {
    cat("\n There are", floor(dose_total*10^3), "vaccine doses left at the end of Year 3. \n We estimated needing", dose_demand_total*10^3, "doses this year. \n We distributed", floor((dose_initial-dose_total)*10^3), "doses of Zostavax this year =", 100*((dose_initial-dose_total)/dose_initial), "% of remaining stock //", 100*((dose_initial-dose_total)/dose_demand_total),"% of estimated demand. \n") 
    
  } else {
    cat("\n There are no vaccine doses left. \n")
  }
  
  
  ######## Coverage estimates
  
  # Combine max_allocations data with other group data
  G.data.shuffled3 <- rbind(G.data.shuffled, max_allocations)
  rownames(G.data.shuffled3)[3] <- "max_allocations3"
  
  # unshuffle and save Group data
  G.data <- G.data.shuffled3[,order(colnames(G.data.shuffled3))]
  
  # return multiple items to store in global envrironment for the next function
  year3_list <- list("G.data.shuffled3" = G.data.shuffled3, "max_allocations3" = max_allocations, "G.data3" = G.data, "leftover_doses3" = dose_total)
  list2env(year3_list ,.GlobalEnv)
  #return()
 
  cat("\n Based on the Zostavax vaccines administered, we estimate the following vaccine protection levels for 65-70 yos against shingles in Ontario: \n\n")
  
  ##################### 4. Plotting #########################
  
  # year_efficacy3 <- vector(mode="double", length=5)
  # 
  # for (i in 1:5) {
  #   year_efficacy3[i] <- (G.data3[[i]][3]*10^3*e1[i]) / (G.data3[[i]][3]*10^3)
  #   
  # }
  
  # df <- data.frame(Age = c("G1 \n (65-66)","G2 \n (66-67)","G3 \n (67-68)","G4 \n (68-69)","G5 \n (69-70)"), Efficacy = year_efficacy3)
  df2 <- data.frame(Age = c("G1 \n (65-66)","G2 \n (66-67)","G3 \n (67-68)","G4 \n (68-69)","G5 \n (69-70)"), Doses_administered = max_allocations3*10^3)
  
  plot1 <- ggplot(df2, aes(Age,Doses_administered)) +
    geom_col(fill="black") +
    labs(color=" ", title="Zostavax doses administered in Year 3 (2018)", x = "\n Age Group", y = "Doses \n") +
    ylim(0,31000)+
    theme_light()
  
  # plot2 <- ggplot(df, aes(Age,Efficacy)) +
  #   geom_col(fill="blue") +
  #   labs(color=" ", title="Population Coverage with Zostavax in Year 3 (2018)", x = "\n Age Group", y = "Coverage Level \n") +
  #   scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
  #   theme_light()
  
  return(list(plot1))
}



# ###### Different scenarios - compare coverage and demand!
      # 1) more uptake from top of age window 
      # 2) more uptake from bottom of age window
      # 3) revised budget




zostavax.recap <- function() {
  
  year1_dose_total <- 0
  year2_dose_total <- 0
  year3_dose_total <- 0
  program_dose_total <- 0
  

  for (i in 1:5) {
    year1_dose_total <- year1_dose_total + G.data1[[i]][3]
    year2_dose_total <- year2_dose_total + G.data2[[i]][3]
    year3_dose_total <- year3_dose_total + G.data3[[i]][3]
  }
  
program_dose_total <- year1_dose_total + year2_dose_total + year3_dose_total
  
cat("\n\n\n SUMMARY \n\n We started with a budget of $", budget*10^3, "resulting in the purchase of", budget_doses*10^3, "doses of Zostavax at time 0.\n")
cat(" In Year 1, we administered", year1_dose_total*10^3, "doses of Zostavax.\n")
cat(" In Year 2, we administered", year2_dose_total*10^3, "doses of Zostavax.\n")
cat(" In Year 3, we administered", year3_dose_total*10^3, "doses of Zostavax.\n")

  if ((budget_doses-program_dose_total) > 0) {
    cat("By the end of the program, we have", (budget_doses-program_dose_total)*10^3, "doses left. \n")
  } else if ((budget_doses-program_dose_total) <= 0) {
    cat("We have no doses left. \n")
  } else {
    cat("SOMETHING WENT WRONG HERE! \n")
  }

cat("Total program coverage is: TBA\n")
cat("End of program.")

}

  zostavax.allocation.year1()
  zostavax.allocation.year2()
  zostavax.allocation.year3()
  zostavax.recap()

cat("\n FYI: Right now, lambda and arrival rate are set as fixed numbers because otherwise they result in NaN. NEED TO FIX!!!\n")


#--end of file--