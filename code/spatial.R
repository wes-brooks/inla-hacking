# source: https://ourcodingclub.github.io/tutorials/inla/

library( "readr" )
library( ggregplot) 
library( stringr)
library( magrittr )

Hosts = read_csv( url("https://raw.githubusercontent.com/ourcodingclub/CC-INLA/master/HostCaptures.csv"))


substr(names(Hosts), 1, 1) <- toupper(substr(names(Hosts), 1, 1)) # Giving the host names capital letters

phen <- c("Grid", "ID", "Easting", "Northing") # Base columns with spatial information we'll need

resp <- "Parasite.count" # Response variable

covar <- c("Month", # Julian month of sampling
           "Sex", # Sex
           "Smi", # Body condition
           "Supp.corrected", # Nutrition supplementation
           "Treated") # Treatment

TestHosts <- na.omit(Hosts[, c(phen, resp, covar)]) # Getting rid of NA's, picking adults
# We are using the [] to subset and only extract specific columns

# Turning variables into factors
TestHosts$Month <- as.factor(TestHosts$Month)
TestHosts$Grid <- as.factor(TestHosts$Grid)

TestHosts$Parasite.count <- round(TestHosts$Parasite.count) # Parasite counts should be integers

table(table(TestHosts$ID)) # Enough repeat samples for a mixed model?



# Setting up a custom theme
THEME <- theme(axis.text.x = element_text(size = 12,colour = "black"),
               axis.text.y = element_text(size = 12, colour = "black"),
               axis.title.x = element_text(vjust = -0.35),
               axis.title.y = element_text(vjust = 1.2)) + theme_bw()

(samp_locations <- ggplot(TestHosts, aes(Easting, Northing)) + 
    geom_jitter(aes(colour = factor(Grid))) + coord_fixed() + 
    THEME + 
    labs(colour = "Grid"))



# model selection ---------------------------------------------------------
# First without random effects ####

# Specify the formula
f0.1 <- as.formula(paste0(resp, " ~ ", # Response first
                          paste(covar, collapse = " + ") # Collapse the vector of covariates
))

# Run the model
IM0.1  <- inla(Parasite.count ~ Month + Sex + Smi + Supp.corrected + Treated, 
               family = "nbinomial", # Specify the family. Can be a wide range (see r-inla.org).
               data = TestHosts) # Specify the data

# Run the model # (This is the same thing)
IM0.1  <- inla(f0.1, 
               family = "nbinomial", # Specify the family. Can be a wide range (see r-inla.org).
               data = TestHosts) # Specify the data

# Then with an ID random effect ####

f0.2 <- as.formula(paste0(resp, " ~ ", 
                          paste(covar, collapse = " + "), 
                          " +  f(ID, model = 'iid')")) # This is how you include  a typical random effect.

IM0.2  <- inla(f0.2, 
               family = "nbinomial",
               data = TestHosts) 

summary(IM0.1)
summary(IM0.2)
Efxplot(list(IM0.1, IM0.2))

# Let's try it on our data ####

HostModelSel <- INLAModelSel(resp, covar, "ID", "iid", "nbinomial", TestHosts)

Finalcovar <- HostModelSel$Removed[[length(HostModelSel$Removed)]]

# Setting up a mesh -------------------------------------------------------
Locations = cbind(TestHosts$Easting, TestHosts$Northing) # using the sampling locations 

MeshA <- inla.mesh.2d(jitter(Locations), max.edge = c(20, 40))
MeshB <- inla.mesh.2d(Locations, max.edge = c(20, 40))
MeshC <- inla.mesh.2d(Locations, max.edge = c(10, 20))

Mesh <- MeshB

plot(MeshA)

plot(MeshB)

plot(MeshC)

points(Locations, col = "red", pch = 2)


##
# Making the A matrix

HostsA <- inla.spde.make.A(Mesh, loc = Locations) # Making A matrix
Hosts.spde = inla.spde2.pcmatern(mesh = Mesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
w.Host <- inla.spde.make.index('w', n.spde = Hosts.spde$n.spde) # making the w




# Make the stack ----------------------------------------------------------

# Making the model matrix #### 

X0 <- model.matrix(as.formula(paste0(" ~ -1 + ", paste(Finalcovar, collapse = " + "))), data = TestHosts) # make the model matrix using the final model selection formula without a response variable.

X <- as.data.frame(X0[,-which(colnames(X0)%in%c("Month7"))]) # convert to a data frame. Eliminate the base level of the first categorical variable if applicable (you will manually specify an intercept below) 

head(X)

# Making the stack ####

N <- nrow(TestHosts)

StackHost <- inla.stack(
  data = list(y = TestHosts[,resp]), # specify the response variable
  
  A = list(1, 1, 1, HostsA), # Vector of Multiplication factors for random and fixed effects              
  
  effects = list(
    
    Intercept = rep(1, N), # specify the manual intercept!
    
    X = X, # attach the model matrix
    
    ID = TestHosts$ID, # insert vectors of any random effects
    
    w = w.Host)) # attach the w 



N <- nrow(TestHosts)

GOODSTACK <- inla.stack(
  data = list(y = TestHosts[,resp]), # specify the response variable
  
  A = list(1, 1, 1, 1, HostsA), # Vector of Multiplication factors for random and fixed effects              
  
  effects = list(
    
    Intercept = rep(1, N), # specify the manual intercept!
    
    X = X, # attach the model matrix
    
    ID = TestHosts$ID, # insert vectors of any random effects
    Grid = TestHosts$Grid,
    
    w = w.Host)) # Leave


