## Applying 'limits'

# Setting basic limits
one2ten <- 1:10

# Creating a dataframe
ggdat <- data.frame(first=one2ten, second=one2ten)

# Base plot
plot(one2ten, one2ten)

# Base ggplot
require(ggplot2)
print(qplot(first, second, data=ggdat))

# Setting typical limits
# Base
plot(one2ten, one2ten, xlim=c(-2,10))

# ggplot
print(qplot(first, second, data=ggdat) + xlim(-2, 10))

# Reverse direction
# ggplot
print(qplot(first, second, data=ggdat) + xlim(10, 1))

## 'Expansion' : To plot all data inside the plot area
plot(one2ten, one2ten, xlim=c(0,10))
par("usr")

# Specifying axes [base]
plot(one2ten, one2ten, xlim=c(0,10), xaxs="i")
plot(one2ten, one2ten, xlim=c(0,10), xaxs="i", yaxs="i")

# Specifying axes [ggplot]
print(qplot(first, second, data=ggdat) +
        scale_x_continuous(expand=c(0,0)))

