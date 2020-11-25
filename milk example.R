# Probabilistic Dynamic programming exercise - milk distribution
  # Current stage costs are UNCERTAIN, but next period's state is CERTAIN (inventory is certain, demand and therefore revenue is uncertain)


library(readxl)
require("openxlsx")

store_data <- read_excel("C:/Users/Lia/OneDrive/GRAD SCHOOL/F20/Shingles/Shingles in R/store_data.xlsx")

###########  Initialization  ########## 

  total_gallons <- 6
  price_per_gallon <- 2
  rev_per_unsold_gallon <- 0.5
  max_demand <- 3
  

#------------ FUNCTIONS -----------##
  # expected revenue function
  
    expected_revenue <- function(gt, prob) {
      
      for (i in 2:4) {
        
        x <- 0
        for (j in 1:3) {
          
          x <- max()
          
        }
      }
      
      return(rev)
    }
#----------------------------------##    


r1 <- c(0,0,0)
r2 <- c(0,0,0)
r3 <- c(0,0,0)

store_list <- list(r1,r2,r3)

############ Step 1. Compute expected revenue earned from g_t in {1, ..., 3} ########################################################################




# r3[1] = (store3_data$gallons[[1]]*2.00)*(store3_data$probability[[1]] + store3_data$probability[[2]] + store3_data$probability[[3]] + store3_data$probability[[4]])
#        #[price (gt x price/gallon) * (probabilities of obtaining that gt or higher)] +  [price (gt x price/gallon) * (probabilities of obtaining less than gt)]
# r3[2] <- (store3_data$gallons[[2]]*2.00) * (store3_data$probability[[2]] + store3_data$probability[[3]] + store3_data$probability[[4]])
# r3[3] <- (store3_data$gallons[[3]]*2.00) * (store3_data$probability[[3]] + store3_data$probability[[4]]) + ((store3_data$gallons[[2]]*2.00) * (store3_data$probability[[2]])) + store3_data$probability[[2]]
# r3[4] <- (store3_data$gallons[[4]]*2.00) * (store3_data$probability[[4]]) +

  
  
for (i in 1:length(store_list)) {
  for (j in 1:max_demand) {
    
    
    
  }
  
  r1[i] <- something
  r2[i] <- something
  r3[i] <- something
}


############ Step 2. Compute optimal allocation of milk to stores ########################################################################
# We need to find each gt(x) of milk to store t that attains ft(x) (the max expected revenue from x gallons left)

# Start with the last store (t = max distribution = 3)
for (i in length(store_list):1)
  for (j in 1:total_gallons) { # x can be 0-6 gallons
    
    
    
  }

# 



