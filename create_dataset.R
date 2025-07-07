library(dplyr)
library(mvtnorm)
library(spatstat)
library(ggplot2)

### Simulations

source("params.R") # Parameters
source("utilis.R") # Utilis functions


# Create a realization of the spatio temporal point pattern
set.seed(17)
sims <- list()
for (i in 0:(max.t-1)){
  point <- h_t*i+h_t/2
  sims[[(i+1)]] <- rpoispp(lambda_e,win=square(1))
}


# Convert the data from the simulation to a data frame structure
df <- data.frame()
for (i in 1:max.t){
  d_temp <- as.data.frame(sims[[i]])
  d_temp$time <- (h_t*(i-1)+h_t/2)
  d_temp$time.index <- i
  d_temp$n <- sims[[i]]$n
  df <- bind_rows(df,d_temp)
}

df$timeL <- paste0(df$time," (n=",df$n,")") # Number of points in each time stamp

# Plot the results
pplot1 <- ggplot() + geom_sf() + coord_sf(datum = NA) +
  geom_point(
    data = df, aes(x = x, y = y),
    size = 0.5
  ) +
  labs(x = "", y = "") +
  facet_wrap(~timeL) +
  theme_bw()
pplot1

ggsave(pplot1,
       file=paste0("./images/simulation.pdf"),
       width = 25, height = 20, units = "cm", dpi=320)


# Save the results for modeling
write.csv(df[,1:4], file = "./data/obs.csv", row.names = FALSE)