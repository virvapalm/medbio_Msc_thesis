# Set seed for reproducibility
set.seed(123)

# Convert the values to the log2 scale
log2_values <- log2(runif(52 * 500) + 1e-10) # Adding a small value to avoid log2(0)

# Create the data frame
df <- data.frame(matrix(log2_values, nrow = 52, ncol = 500))

# Display the first few rows of the data frame
colnames(df)<-paste('Protein', sample(35:3825, 500, replace = TRUE))
rownames(df)<-paste('Patient', sample(1:1000, 52, replace = TRUE))

#write.csv(df, 'test_data.csv')

#Meta data
df$Sex<-sample(c('female','male'), 52, replace = TRUE)
df$Age<-sample(26:95, 52, replace = TRUE)
df$Co_infection<-sample(c(TRUE, FALSE),52, replace = TRUE)
df$Alive90<-sample(c(TRUE, FALSE),52, replace = TRUE)

df%>%
  colnames()
