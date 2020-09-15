
# CLEAR AND PATHS ---------------------------------------------------------

rm(list = ls())


# FUNCTIONS ---------------------------------------------------------------

sigmoid <- function(x, a, b, c, d)
{
    y <- a + (b -a) / (1 + exp((c - x) / d))
    return(y)
}

# CREATE VARIABLES --------------------------------------------------------

x <- 1:6
a <- 0.5
b <- 1
c <- 3
d <- 1

y <- sigmoid(x, a, b, c, d)

# PLOT --------------------------------------------------------------------

par(font.axis=2, font.lab=2, lwd=3)

plot(x, y, type='b', ylim=c(0, 1), xlab='N frames', ylab='Proportion Correct',
     main='Sigmoid')
lines(x, rep(0.50, length(x)), lty=3)
