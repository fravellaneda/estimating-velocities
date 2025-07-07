library(dplyr)
library(mvtnorm)
library(spatstat)
library(sf)
library(sp)
library(INLA)
library(inlabru)
library(INLAspacetime)

### Simulations

source("params.R")
source("utilis.R")

# read dataset with point pattern
df <- read.csv("./data/obs.csv")

# Get the location in the proper format for modeling.
df_sf <- st_as_sf(df,coords = c("x","y"))
df$geometry <- st_geometry(df_sf)

# Create the spatial region for the point pattern
box <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0,
                           xmax = 1, ymax = 1)))
map <- st_as_sf(box, crs = st_crs(4326))

# Create the location of the prediction points
bb <- st_bbox(map)
h.t <- (1/(espatial.res-2))
x <- round(seq(bb$xmin-h.t, bb$xmax+h.t ,length.out = (espatial.res+1)),4)
dp <- as.matrix(expand.grid(x, x))
write.csv(dp, file = "./data/dp.csv", row.names = FALSE)

dpdf_sf <- st_as_sf(data.frame(x=dp[,1],y=dp[,2]),coords = c("x","y"))

bnd <- inla.nonconvex.hull(dp)

# Create the spatial mesh for the SPDE
smesh <- inla.mesh.2d(
  boundary = bnd,
  max.edge = c(0.1),
  offset = c(0.1, 0.1)
)

# Create the temporal mesh for the SPDE
tmesh <- inla.mesh.1d(
  loc = unique(df$time.index))

#dataset clean with the location and the time
df.new <- df[,c("geometry","time.index")]
names(df.new) <- c("geometry","time")

# Dataset for predictions
dpdf <- data.frame(time = unlist(lapply(1:max.t, function(x) rep(x, nrow(dp)))))
dpdf$geometry <- rep(st_geometry(dpdf_sf),max.t)

# Defining the model based on the priors and the model in params.R
stmodel <- stModel.define(
  smesh, tmesh, model,
  control.priors = list(
    prs = c(1, 0.05),
    prt = c(10, 0.05),
    psigma = c(2, 0.05)),
  constr = TRUE)

# Formula for the log of tje lcpg
formula.M <- geometry + time ~ Intercept(1) +
  field(list(space = geometry, 
             time = time),
        model = stmodel)

# Estimation of the model
print("Modeling")
t1 <-Sys.time()
fit <- lgcp(formula.M,
            data = df.new,
            domain = list(
              geometry = smesh,
              time = tmesh
            ),
            options = bru_options_set(
              verbose = TRUE,
              bru_verbose = TRUE
            )
)
t2 <- Sys.time()
print("Time modeling")
print(t2-t1)

# Create the prediction dataset for predicting
ppxl_all <- fm_cprod(dpdf_sf, data.frame(time = unique(df$time.index)))

print("Predicting")
t1 <-Sys.time()
lambda1 <- predict(
  fit,
  ppxl_all,
  ~ data.frame(time = time, lambda = exp(field + Intercept))
)
t2 <- Sys.time()
print("Time modeling")
print(t2-t1)

# Save the results of the modeling
save(list = c("lambda1"), file = "./data/lambda1.Rdata")
