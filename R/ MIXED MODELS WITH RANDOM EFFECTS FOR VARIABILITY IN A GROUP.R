data.fn <- function(nsite = 5, nyear = 40, alpha = 4.18456, 
                    beta1 = 1.90672, beta2 = 0.10852, beta3 = - 1.17121, sd.site = 0.5, 
                    sd.year = 0.2) {
  log.expected.count <- array(NA, dim = c(nyear, nsite))
  C <- log.expected.count 
  year <- 1:nyear 
  yr <- (year - 20)/20 #Standardize 
  site <- 1:nsite
  alpha.site <- rnorm(n = nsite, mean = alpha, sd = sd.site)
  eps.year <- rnorm(n = nyear, mean = 0, sd = sd.year)
  for (j in 1:nsite){
    log.expected.count[,j] <- alpha.site[j] + beta1*yr + beta2*yr^2 + beta3*yr^3 +beta3*yr^3 + eps.year
    expected.count <- exp(log.expected.count[,j])
      C[,j] <- rpois(n = nyear, lambda = expected.count) }
  
  matplot(year, C, type = "l", lty = 1, lwd = 2, main = "", las = 1, 
          ylab = "Population size", xlab = "Year") 
  return(list(nsite = nsite, nyear = nyear, alpha.site = alpha.site, 
              beta1 = beta1, beta2 = beta2, beta3 = beta3, year = year, 
              sd.site = sd.site, sd.year = sd.year, expected.count = expected.count, C = C)) }

data.1 <- data.fn(nsite = 100, nyear = 40, sd.site = 0.3, sd.year = 0.2)
win.data <- list(C = data.1$C, nsite = ncol(data.1$C), nyear = nrow(data.1$C), year = (data.1$year-20) / 20) # Note year standardized


