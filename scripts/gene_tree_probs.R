#Calculate probabilities of gene trees under different sets of coalescent branch lengths

#----MSC----#
#species tree
p_s_msc <- function(x){
  p <- 1- (2/3) * exp(-x)
  return(p)
}

#discordant gene trees = prob(discordant topology 1) + prob(discordant topology 2)
p_d_msc <- function(x){
  p <- (2/3) * exp(-x)
  return(p)
}

# one discordant topology = prob(discordant topology 1) = prob(discordant topology 2)
p_d1_msc <- function(x){
  p <- (1/3) * exp(-x)
  return(p)
}
#-----------#

#----NMSC----#
#species tree
p_s_nmsc <- function(x,t1,gamma){
  p <- (1 - gamma) * (1 - ((2/3) * exp(-t1))) + gamma * ((1/3) * exp(-x))
  return(p)
} 
#major discordance topology
p_major_nmsc <- function(x,t1,gamma){
  p <- (1 - gamma) * ((1/3) * exp(-t1)) + gamma * ((1/3) * exp(-x))
  return(p)
} 
#minor discordant topology
p_minor_nmsc <- function(x,t1,gamma){
  p <- (1 - gamma) * ((1/3) * exp(-t1)) + gamma * (1 - ((2/3) * exp(-x)))
  return(p)
} 
#------------#