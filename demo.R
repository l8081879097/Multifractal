source("~/Downloads/Multifractal.R") # Load multifractal function

pdsi <- read.csv("~/Downloads/PDSI_DebrisFlow.csv") # Load example data
pdsi_n <- (pdsi-min(pdsi))/(max(pdsi)-min(pdsi)) # Normalization

plot(pdsi$avg,type = "n",ylab = "PDSI")
lines(pdsi$avg) # Show the PDSI time series

mf <- multi_fractal(pdsi_n$avg,
                    q = 5,
                    step = 1,
                    nbox = 4)

plot(mf,type = "b")