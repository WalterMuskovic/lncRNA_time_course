library(ggplot2)
library(patchwork)
library(forecast)

# Define function to generate expression data using an impulse model and adding some noise
impulse_plus_noise <- function(h0=0, h1=10, h2=0, lambda=0.1, t=seq(0,400,10), t1) { 
  t2=t1+120 # set offset time to 120 min after onset
  (1/h1)*((h0+(h1-h0)*(1/(1+exp(-lambda*(t-t1)))))*(h2+(h1-h2)*(1/(1+exp(lambda*(t-t2)))))) + rnorm(length(t), mean = 1, sd = 0.5)
}

plot_lags <- function(shift){
  # Create some gene expression data
  test_data <- data.frame(time = seq(0,400,10),
                          gene_1_expression = impulse_plus_noise(t1 = 100),
                          gene_2_expression = impulse_plus_noise(t1 = 100 + shift))
  
  # Plot expression of the two genes + ccf plot
  g1 <- ggplot(test_data, aes(x=time, y=gene_1_expression)) + geom_point()
  g2 <- ggplot(test_data, aes(x=time, y=gene_2_expression)) + geom_point()
  g3 <- ggCcf(x = test_data$gene_1_expression, y = test_data$gene_2_expression, lag.max = 20, type="correlation", plot = TRUE) + ggtitle("lagged expression")
  
  return(g1 | g2 | g3)
}

# Shift gene 2's expression forward and backwards in time relative to gene 1
set.seed(1234)
ggsave(filename = "figures_track/lag_example.png",
       plot = plot_lags(-60) / plot_lags(0) / plot_lags(60),
       device = "png", width = 25, height = 20, units = "cm")
