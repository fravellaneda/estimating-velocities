library(mvtnorm)
library(sf)
library(INLA)
library(inlabru)
library(INLAspacetime)
library(grid)
library(ggplot2)
library(viridis)
library(patchwork)

### Simulations

source("params.R")
source("utilis.R")


#Function to compute velocities
velocities <- function(point1,point2,h.t,lambda1,per,time){
  # select the columns that will keep 
  keep <- c(time,"mean")
  
  # converting the result of the modeling in a dataset
  lambda.df <- st_drop_geometry(lambda1[,keep])
  lambda.df$x <- st_coordinates(lambda1)[,1]
  lambda.df$y <- st_coordinates(lambda1)[,2]
  
  # Keep the points neccesarly for the temporal derivatives
  lambda.df <- lambda.df[lambda.df[,time] %in% c(point1,point2),]
  rownames(lambda.df) <- 1:nrow(lambda.df)
  head(lambda.df)
  
  # Pivot the information
  pivoted <- reshape(lambda.df,
                     timevar = time,
                     idvar = c("x","y"),
                     direction = "wide"
  )
  names(pivoted) <- c("x","y","mean1","mean2")
  
  # Order the dataset and stablish and id to compute finite differences
  pivoted <- pivoted[order(pivoted$x,-pivoted$y),]
  rownames(pivoted) = 1:(nrow(pivoted))
  pivoted$id <- 1:(nrow(pivoted))
  comp.vel <- pivoted
  
  # Compute derivatives in time
  comp.vel$mean <- comp.vel$mean2
  comp.vel$mean.t = (comp.vel$mean2-comp.vel$mean1)/(h.t)
  
  # organize indices to compute derivatives in space
  comp.vel$id_up <- comp.vel$id-1
  comp.vel$id_do <- comp.vel$id+1
  comp.vel$id_r <- comp.vel$id+espatial.res+1
  comp.vel$id_l <- comp.vel$id-espatial.res-1
  comp.vel[seq(1,max(comp.vel$id),by=(espatial.res+1)),"id_up"] <- NaN
  comp.vel[seq((espatial.res+1),max(comp.vel$id),by=(espatial.res+1)),"id_do"] <- NaN
  comp.vel[which(comp.vel$id >=(espatial.res+1)*(espatial.res+1)),"id_r"] <- NaN
  comp.vel[which(comp.vel$id <=(espatial.res+1)),"id_l"] <- NaN
  
  # Compute the differences that are in space
  h.x <- unique(comp.vel$x)[2]-unique(comp.vel$x)[1]
  h.y <- unique(comp.vel$y)[2]-unique(comp.vel$y)[1]
  
  # Get the values that are up in space
  comp.vel.prov = comp.vel[,c("id","mean")]
  names(comp.vel.prov) <- c("id_up","mean_up")
  comp.vel <- merge(comp.vel, comp.vel.prov, by="id_up", all.x = TRUE)
  
  # Get the values that are down in space
  comp.vel.prov = comp.vel[,c("id","mean")]
  names(comp.vel.prov) <- c("id_do","mean_do")
  comp.vel <- merge(comp.vel, comp.vel.prov, by="id_do", all.x = TRUE)
  
  # Get the values that are on the rigth in space
  comp.vel.prov = comp.vel[,c("id","mean")]
  names(comp.vel.prov) <- c("id_r","mean_r")
  comp.vel <- merge(comp.vel, comp.vel.prov, by="id_r", all.x = TRUE)
  
  # Get the values that are on the left in space
  comp.vel.prov = comp.vel[,c("id","mean")]
  names(comp.vel.prov) <- c("id_l","mean_l")
  comp.vel <- merge(comp.vel, comp.vel.prov, by="id_l", all.x = TRUE)
  
  # Order the data set again 
  comp.vel <- comp.vel[order(comp.vel$id), ]
  row.names(comp.vel) <- 1:nrow(comp.vel)
  head(comp.vel)
  
  # Compute partial derivatives in space using forwards and backwards differentiation
  comp.vel["difx1"] <- (comp.vel[,"mean_r"]-comp.vel[,"mean"])/h.x
  comp.vel["difx2"] <- (comp.vel[,"mean"]-comp.vel[,"mean_l"])/h.x
  comp.vel["dify1"] <- (comp.vel[,"mean"]-comp.vel[,"mean_do"])/h.y
  comp.vel["dify2"] <- (comp.vel[,"mean_up"]-comp.vel[,"mean"])/h.y
  
  Mean2 <- function(x,y) {mean(c(x,y))}
  Mean4 <- function(x,y,z,w) {mean(c(x,y,z,w),na.rm=FALSE)}
  
  # Approximate derivatives in direction x and y
  comp.vel["difx_pre"] = mapply(Mean2,comp.vel$difx1,comp.vel$difx2)
  comp.vel["dify_pre"] = mapply(Mean2,comp.vel$dify1,comp.vel$dify2)
  
  # Computing approximation of the gradient
  comp.vel["difxn"] = (comp.vel["difx_pre"]/
                         (comp.vel["difx_pre"]**2+comp.vel["dify_pre"]**2)**(1/2))*h.x
  comp.vel["difyn"] = (comp.vel["dify_pre"]/
                         (comp.vel["difx_pre"]**2+comp.vel["dify_pre"]**2)**(1/2))*h.y
  
  # Computing approximation of the norm of the gradient
  comp.vel["norm1"] = mapply(Norma,comp.vel$difx1,comp.vel$dify1)
  comp.vel["norm2"] = mapply(Norma,comp.vel$difx1,comp.vel$dify2)
  comp.vel["norm3"] = mapply(Norma,comp.vel$difx2,comp.vel$dify1)
  comp.vel["norm4"] = mapply(Norma,comp.vel$difx2,comp.vel$dify2)
  comp.vel["norm"] = mapply(Mean4,comp.vel$norm1,comp.vel$norm2,
                            comp.vel$norm3,comp.vel$norm4)
  
  # Computing velocities and its direction
  comp.vel["vel"] = (abs(comp.vel["mean.t"])/comp.vel["norm"])
  comp.vel["difxn2"] = sign(comp.vel$mean.t)*comp.vel["difxn"]
  comp.vel["difyn2"] = sign(comp.vel$mean.t)*comp.vel["difyn"]
  
  # Truncating velocities for ploting purposes 
  comp.vel[,"velt"] <- comp.vel$vel
  sup = round(quantile(comp.vel$vel,per,names=FALSE,na.rm=TRUE),1)
  comp.vel[!is.na(comp.vel$vel) & comp.vel$vel>sup,"velt"] <- sup*0.999
  
  return(comp.vel)
}

# Loading datasets
df <- read.csv("./data/obs.csv")
colnames(df)[1:4] <- c("x","y","time.val","time")
load(file = "./data/lambda1.Rdata")
dp <- read.csv("./data/dp.csv")

# Ploting intensity function and the observation points
pl1 <- ggplot() +
  gg(lambda1, geom = "tile", aes(fill = mean)) +
  scale_fill_viridis_c() +
  geom_point(
    data = df, aes(x = x, y = y),
    size = 0.5
  ) +
  facet_wrap(~time) +
  coord_sf()
pl1
ggsave(pl1,
       file=paste0("./images/intensity.pdf"),
       width = 25, height = 20, units = "cm", dpi=320)


# end point is the point that you want to compute the approximation of the
# velocity, ini point is the reference for the temporal derivative
end.point <- 5
ini.point <- 1
h.t <- min(df[df$time==end.point,"time.val"])-min(
  df[df$time==ini.point,"time.val"])

# Compute true velocities and velocity estimation
true.vel <- true.values(dp,h_t*end.point,0.95)
comp.vel <- velocities(ini.point,end.point,h.t,lambda1,0.99,"time")

# Defining same scale and in this case 10 breakpoints
vlim <- range(c(true.vel$velt,comp.vel$velt), na.rm = TRUE)
breaks <- pretty(vlim, n = 10)

# True velocities plot
p1 <- ggplot(true.vel, aes(x = x, y = y)) +
  geom_contour_filled(aes(z = velt), breaks = breaks) +
  geom_segment(aes(xend = x + difxn2/2, yend = y + difyn2/2),
               arrow = arrow(length = unit(0.1, "cm")),
               size = 0.25, colour = "black") +
  scale_fill_viridis_d(drop = FALSE)

# Approximation of the velocities plot
p2 <- ggplot(comp.vel, aes(x = x, y = y)) +
  geom_contour_filled(aes(z = velt), breaks = breaks) +
  geom_segment(aes(xend = x + difxn2/2, yend = y + difyn2/2),
               arrow = arrow(length = unit(0.1, "cm")),
               size = 0.25, colour = "black") +
  scale_fill_viridis_d(drop = FALSE)  # same scale

p1 + p2
ggsave(p1+p2,
       file=paste0("./images/velocity.pdf"),
       width = 30, height = 15, units = "cm", dpi=320)

