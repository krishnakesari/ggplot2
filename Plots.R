---
  title: "Plots"
output: html_document
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
library(dplyr)
library(ggplot2)
library(quantreg)
library(hexbin)
library(Hmisc)
```

```{r}
data(mtcars)
head(mtcars,3)
```

```{r}
## Convering cyl to a factor (grouping of variables in a data) varaible 
df <- mtcars[, c("mpg", "cyl", "wt")]
df$cyl <- as.factor(df$cyl)
head(df)
```

```{r}
# Basic plots
qplot (mpg, wt, data = mtcars, geom = c("point", "smooth"))

## Applying colour by a continuous numerical variable
qplot(mpg, wt, data = mtcars, color = cyl)

## Applying colour and shape by group (factor variable - df)
mtcars$cyl <- factor(mtcars$cyl)
qplot(mpg, wt, data = mtcars, colour = cyl, shape = cyl)

## Changing size of the points based on the value of continous variable
qplot(mpg, wt, data = mtcars, size = mpg)
```

```{r}
# Box plot, histogram and density plots

set.seed(1234)
wdata = data.frame(
  sex = factor(rep(c("F", "M"), each = 200)), 
  weight = c(rnorm(200, 55), rnorm(200,55)))
head(wdata, 3)

## Boxplots (fill by sex)
qplot(sex, weight, data = wdata, geom = "boxplot")

## Histogram 
qplot(weight, data = wdata, geom = "histogram")

## Denisty plot with Axes (labellling)
qplot(weight, data = wdata, geom = "density",
      xlab = "Weight (Kgs)", ylab = "Density", 
      main = "Density plot")
```

```{r}
# Custom plots (ggplot)
## Point plot
ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point(size = 1.5, shape = 18)

## Geometry function - Density plot
ggplot(wdata, aes(x= weight)) + geom_density()

## Statistics function - Density plot
ggplot(wdata, aes(x=weight)) + stat_density()

## Points and lines
ggplot(data = mtcars, aes(x=wt, y=mpg)) +
  geom_point() + 
  geom_line()

## Colouring specific data with line
ggplot(data = mtcars, aes(x=wt, y=mpg)) +
  geom_point()+
  geom_line(data = head(mtcars), color = "red")
```

```{r}
#Advanced graphs
## Log transformation in aes()
ggplot(data = mtcars, aes(x= log2(wt), y=log2(mpg))) +
  geom_point()
```

```{r}
# Using ggpoints
## Building helper function for ggpoints
ggpoints <- function(data, xName, yName) {
  p <- ggplot(data = data, aes_string(xName, yName)) +
    geom_point(color = "red") +
    geom_smooth()
  
  return(p)
}

## Scatter Plot
ggpoints(mtcars, xName = "wt", yName = "mpg")
```

```{r}
## Saving plots in pdf
pdf("myplot.pdf")
myplot <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
print(myplot)
dev.off()

## saving in png file
png("myplot.png")
print(myplot)
dev.off()

## preferred way
ggsave("myplot.png")
```

```{r}
# Plotting Continuous variables (stats)
set.seed(1234)
wdata <- data.frame(
  sex = factor(rep(c("F", "M"), each = 200)),
  weight = c(rnorm(200, 55), rnorm(200, 58)))
head(wdata, 4)

library(dplyr)
mu <- wdata %>%
  group_by(sex) %>%
  summarise(grp.mean = mean(weight))

head(mu)

## Plotting
### For one Continous Variable
a <- ggplot(wdata, aes (x=weight))
a+ geom_density() ## Density plot
a+ geom_dotplot() ## Dot plot
a+ geom_freqpoly() ## frequency polygon
a+ geom_histogram() ## Histogram plot
a+ stat_ecdf () ## Emperical Cumulative Density Function 
## a+ stat_qq() ## quantile - quantile plot

### For one Discrete Variable
a+geom_bar() ## For Bar Chart

```

```{r}
# Area Plots

a + geom_area(stat = "bin",
              color = "black", fill = "#00AFBB")

## y-value corresponds to count of weight values. 

### To have densities on y-axis
a + geom_area(aes(y= ..density..), stat = "bin")
```

```{r}
# Desnity plots with mean line

a + geom_density(color = "black", fill = "gray") +
  geom_vline(aes(xintercept = mean(weight)),
             color = "#FC4E07", linetype = "dashed", size = 1)

## Changing colours by group
a + geom_density(aes(color = sex))

## Changing colours by group with transparency
a + geom_density(aes(fill = sex), alpha = 0.4)

## Adding mean lines and colour by sex
a + geom_density(aes(color = sex), alpha = 0.4) +
  geom_vline(data = mu, aes(xintercept = grp.mean, color = sex), 
             linetype = "dashed")

## Changing colour manually
a2 <- a + geom_density(aes(color = sex)) +
  geom_vline(data = mu, aes(xintercept = grp.mean, color = sex), 
             linetype = "dashed") + theme_minimal()
a2 + scale_color_manual(values = c("#999999", "#E69F00"))

## Using Brewer color palettes
a2 + scale_color_brewer(palette = "paired")

## Using grey scale
a2 + scale_color_grey()

## With colour fill

a3 <- a + geom_density(aes(fill = sex), alpha = 0.4) + theme_minimal()
a3 + scale_fill_manual(values = c("#999999", "#E69F00"))

a3 + scale_fill_brewer(palette = "Dark2") + theme_minimal()

a3 + scale_fill_grey() + theme_minimal()
```

```{r}
# Histogram Plots
## Distribution of continuous variables by dividing into binds and counting the number of observations in each bin   using "geom_histogram()"

### Basic Plot
a + geom_histogram()

### Changing the number of bins
a + geom_histogram(bins = 50)

### Changing line colour and fill colour, add mean line
a + geom_histogram(color = "Blue", fill = "gray") +
  geom_vline(aes(xintercept = mean(weight)),
             color = "#FC4E07", linetype = "dashed", size = 1)

### Notes:  In general Y-axis is count, to change the Y-axis

a + geom_histogram(aes(y= ..density..)) + geom_histogram(color = "Blue", fill = "gray") + geom_vline(aes(xintercept = mean(weight)),
                                                                                                     color = "#FC4E07", linetype = "dashed", size = 1)

## Changing colors by group (Sex)
a + geom_histogram(aes(color = sex), fill = "white")

### Position adjustments: "identity" (overlaid)
a + geom_histogram(aes(color = sex), fill = "white", alpha = 0.6,
                   position = "identity")

### Position adjustments: "dodge" (interleaved)
a + geom_histogram(aes(color=sex), fill = "white",
                   position = "dodge") +
  geom_vline(data = mu, aes(xintercept = grp.mean, color = sex),
             linetype = "dashed")

## Changing fill colors 
### Change outline color manually
a + geom_histogram(aes(color=sex), fill = "white",
                   alpha = 0.4, position = "identity") +
  scale_color_manual(values = c("#00AFBB", "#E7B800"))

### Change and fill outline color manually 
a + geom_histogram(aes(color = sex, fill = sex),
                   aplha = 0.4, position = "identity") +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
```

```{r}
# Combining Histogram and Density plots
## Basic plot
a + geom_histogram(aes(y=..density..), color = "black", fill = "white") +
  geom_density(alpha = 0.2, fill = "#FF6666")

## Colors with groups
a + geom_histogram(aes(y=..density.., color = sex, fill = sex), 
                   alpha = 0.5, position = "identity") +
  geom_density(aes(color = sex), size = 1)
```

```{r}
# Frequency Polygon
## Histograms uses 'bars' amd frequency polygon uses "lines" (function: geom_freqpoly)

### Basic plot
a + geom_freqpoly(bins=30)+
  theme_minimal()

### Using color and linetype by sex
a + geom_freqpoly(aes(color = sex, linetype = sex)) +
  scale_color_manual(values = c("#999999", "#E69F00")) +
  theme_minimal()
### With density
a + geom_freqpoly(aes(y=..density..)) + geom_freqpoly(aes(color = sex, linetype = sex))     + scale_color_manual(values = c("#999999", "#E69F00")) +
  theme_minimal()
```

```{r}
# Dotplots for one variable
## width of a dot corresponds to the bin width
a + geom_dotplot(aes(fill=sex))
```

```{r}
# ECDF (Empirical Cumulative Density Function) plots 
## Reports for any given number the percent of individuals that are below that threshold. (function - stat_ecdf())
a + stat_ecdf(geom = "point")
a + stat_ecdf(geom = "step")
```

```{r}
# QQ Plots (Quantile - Quantile plots) (function - stat_qq() or qplot())
data(mtcars)
## Covert cyl column from a numeric to a factor variable
mtcars$cyl <- as.factor(mtcars$cyl)
head(mtcars [, c("mpg", "cyl")])

p <- ggplot(mtcars, aes(sample = mpg))

### Basic plot
p + stat_qq()

### Change point shapes by groups
#### Use custom color palettes
p + stat_qq(aes(shape = cyl, color = cyl)) +
  scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"))
```

```{r}
# Bar plots of counts
## To visualize discrete variable 
data(mpg)
ggplot(mpg, aes(fl)) +
  geom_bar(fill = "steelblue") +
  theme_minimal()
```

```{r}
# Plot two variables - X & Y : Both Continuous or Discrete
## Scatter plots: Continuous X and Y
data(mtcars)
mtcars$cyl <- as.factor(mtcars$cyl)
head(mtcars[, c("wt", "mpg", "cyl")], 3)
b <- ggplot(mtcars, aes(x=wt, y=mpg))

b + geom_point() # scatter plot
b + geom_smooth() # added smoothed line such as regression line
b + geom_quantile() # added quantile lines
b + geom_rug() # added a marginal rug
b + geom_jitter() # for avoiding overplotting
b + geom_text(label = rownames(mtcars)) # text plot

## Changing colors
b + geom_point(color = "#00AFBB", size = 2, shape = 23)

## Control point size by continous variable values
b + geom_point(aes(size=qsec), color = "#00AFBB")

## adding point labels
b + geom_point() +
  geom_text(label = rownames(mtcars), nudge_x = 0.5)

# Scatter Plot with multiple groups
## Change point shapes by level of cyl
b + geom_point(aes(shape = cyl))
## Change point shapes and colors
b + geom_point(aes(shape = cyl, color = cyl))

# Change the point color/shape/size manually
## Change the point sizes manually

b + geom_point(aes(color = cyl, shape = cyl, size = cyl)) +
  scale_size_manual(values = c(2,3,4))

## Change point shapes and colours manually

b + geom_point(aes(color = cyl, shape = cyl)) +
  scale_size_manual(values = c(3,16,17)) +
  scale_color_manual(values = c('#999999', '#E69F00', '#56B4E9'))

## using brewer color palettes

b + geom_point(aes(color = cyl, shape = cyl)) +
  scale_color_brewer(palette="Dark2") + theme_minimal()

## Using grey scale

b + geom_point(aes(color = cyl, shape = cyl)) +
  scale_color_grey() + theme_minimal()

## Adding regression line or smoothed conditional mean

### Notes: possible methods (lm, glm, gam, loess, rlm)

b + geom_point() + geom_smooth(method = lm)

## Point + regression line and removing confidence interval
b + geom_point() + geom_smooth(method = lm, se = FALSE)

## loess method: local regression fitting
b + geom_point() + geom_smooth()

### Changing point color and shapes by groups(cyl)
b + geom_point(aes(color = cyl, shape = cyl)) + 
  geom_smooth(aes(color = cyl, fill = cyl), method = lm)

### Remove confidence intervals
### Extent the regression lines: fullrange
b + geom_point(aes(color = cyl, shape = cyl)) + 
  geom_smooth(aes(color = cyl), method = lm, se = FALSE,
              fullrange = TRUE)

### Add quantile lines from a quantile regression
ggplot(mpg, aes(cty, hwy)) +
  geom_point() + geom_quantile() +
  theme_minimal()

### Add marginal rugs to a scatter plot
b + geom_point() + geom_rug()

#### Adding colors
b + geom_point(aes(color = cyl)) + 
  geom_rug(aes(color = cyl))

#### Adding marginal rugs using faithful data
data(faithful)
ggplot(faithful, aes(x = eruptions, y = waiting)) +
  geom_point() + geom_rug()

#### Jitter points to reduce overplotting
p <- ggplot(mpg, aes(displ, hwy))
##### Default scatter plot
p + geom_point()
##### Use jitter to reduce over plotting
p + geom_jitter(position = position_jitter(width = 0.5, height = 0.5))
##### Textual annotations
b + geom_text(aes(label = rownames(mtcars)), 
              size = 3)
```

```{r}
# Continuous bivariate distribution
data(diamonds)
head(diamonds[, c("carat", "price")])

c <- ggplot(diamonds, aes(carat, price))

## possible layers
c + geom_bin2d()   #adding a heatmap of 2d bin counts [rectangular bining]
c + geom_hex()   # need 'hexbin' package 
c + geom_density2d() #for adding contours from a 2d density estimate

## Add heatmap of 2d bin counts
c + geom_bin2d(bins = 15) # with number of bins = 15
c + geom_bin2d(binwidth = c(1, 1000)) # specific the width of bins

### notes: Above can be applied to hex too

## Adding colors to the density map and points
data("faithful")
# scatter plot
sp <- ggplot(faithful, aes(x = eruptions, y = waiting))
# Default plot
sp + geom_density_2d(color = "#E7B800")
## Adding points
sp + geom_point(color = "#00AFBB") +
  geom_density2d(color = "#E7B800")
## using stat_density_2d with geom = "polygon"
sp + geom_point() +
  stat_density2d(aes(fill=..level..), geom = "polygon") 

### Adding gradient colors
sp + geom_point() +
  stat_density2d(aes(fill=..level..), geom = "polygon") +
  scale_fill_gradient(low = "#00AFBB", high = "#FC4E07")

### Alternative function 
sp + stat_density2d()

## Continuous function 
data(economics)
d <- ggplot(economics, aes(x = date, y = unemploy))
d + geom_area(fill = "#00AFBB", color = "white") # Area plot
d + geom_line(color = "#E7B800") # Line plot 

### Connecitng observations by stairs - a subset of economic data is used 

set.seed(1234)
ss <- economics[sample(1:nrow(economics), 15), ]
ggplot(ss, aes(x=date, y=unemploy)) +
  geom_step(color = "#FC4E07")
```

```{r}
# Two variables: Discrete X, Discrete Y
data("diamonds")
ggplot(diamonds, aes(cut, color)) +
  geom_jitter(aes(color = cut), size = 0.5)
```

```{r}
# Plot two variables - X & Y: Discrete X, Continuous Y
data("ToothGrowth")
ToothGrowth$dose <- as.factor(ToothGrowth$dose)
head(ToothGrowth)

e <- ggplot(ToothGrowth, aes(x = dose, y = len))
e + geom_boxplot()
e + geom_violin()
e + geom_dotplot()
e + geom_line()
# e + geom_bar()
```

```{r}
# Box Plots 
## geom_boxplot(outlier.colour = "black", outlier.shape = 16, outlier.size = 2, notch = FALSE)

e + geom_boxplot() # basic box plot
e + geom_boxplot() + coord_flip() # flipping plot
e + geom_boxplot(notch = TRUE) # notched box plot
e + geom_boxplot() + 
  stat_summary(fun.y = mean, geom = "point",
               shape = 18, size = 4, color = "blue")

## Choose which items to display: group "0.5" and "2"
e + geom_boxplot() + 
  scale_x_discrete(limits = c("0.5", "2"))

## Change the default order of items
e + geom_boxplot() +
  scale_x_discrete(limits = c("2", "0.5", "1"))

e + stat_boxplot(coeff = 1.5)

### Change colors by groups
e + geom_boxplot(color = "black", fill = "steelblue")
e + geom_boxplot(aes(color = dose))
e + geom_boxplot(aes(fill = dose))

### Manual box plot outline/fill colors functions
# scale_color_manual(), scale_fill_manual() : to use custom colors
e2 <- e + geom_boxplot(aes(color = dose)) + theme_minimal()
e2 + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

# scale_color_brewer(), scale_fill_brewer() : to use color palettes
e2 + scale_color_brewer(palette = "Dark2")

# scale_color_grey(), scale_fill_grey() : to use grey color palettes
e2 + scale_color_grey()

### Change manually fill colors
# Using custom color palettes
e3 <- e + geom_boxplot(aes(fill = dose)) + theme_minimal()
e3 + scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))

#Using brewer color palettes
e3 + scale_fill_brewer(palette = "Dark2")

# Using grey scale
e3 + scale_fill_grey()

### Box plot with multiple groups
# change box plot colors by groups
e + geom_boxplot(aes(fill = supp))

# Change the position
e + geom_boxplot(aes(fill = supp), position = position_dodge(1))

# Change fill colours
e + geom_boxplot(aes(fill = supp), position = position_dodge(1)) +
  scale_fill_manual(values = c("#999999", "#E69F00"))
```

```{r}
# Violing plots

e + geom_violin()
e + geom_violin() + coord_flip() # Rotate the violin plot
e + geom_violin(trim = FALSE, fill = "steelblue")

## Adding Summary statistics
### Adding mean or median points: use fun.y=mean or fun.y=median

e + geom_violin(trim= FALSE) +
  stat_summary(fun.y = mean, geom = "point",
               shape = 23, size = 2, color = "blue")

### Add mean points +/- SD
#### use geom = "pointrange" or geom = "crossbar"

e + geom_violin(trim = FALSE) + 
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult=1),
               geom= "pointrange", color = "red")

### Combine with box plot to add median and quartiles

e + geom_violin(trim= FALSE) + 
  geom_boxplot(width = 0.2)

#### Change colors by groups
##### Outline color by dose (groups)
e + geom_violin(aes(color = dose), trim = FALSE)

##### Fill color by dose (groups)
e + geom_violin(aes(fill=dose), trim = FALSE)

##### Change manually outline colors
e2 <- e + geom_violin(aes(color = dose), trim = FALSE) + theme_minimal()
e2 + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

##### Change manually fill colors
e3 <- e + geom_violin(aes(fill = dose), trim = FALSE) + theme_minimal()
e3 + scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))

#### Violin plots with multiple groups
##### Change colors by groups
e + geom_violin(aes(fill = supp), trim = FALSE)

##### Change fill colors by groups
e + geom_violin(aes(fill = supp), trim = FALSE) +
  scale_fill_manual(values = c("#999999", "#E69F00"))
```

```{r}
# Dot plots

e + geom_dotplot(binaxis = "y", stackdir = "center")

## Changing dotsize and stack ratio
e + geom_dotplot(binaxis = "y", stackdir = "center",
                 stackratio = 1.5, dotsize = 1.1)

## add summary statistics
### Add mean or median: use fun.y = mean or fun.y = median

e + geom_dotplot(binaxis = "y", stackdir = "center") +
  stat_summary(fun.y = mean, geom = "point",
               shape = 18, size = 3, color = "red")

# Add mean points +/_ SD
## use geom = "pointrange" or geom = "crossbar"

e + geom_dotplot(binaxis = "y", stackdir = "center") +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "pointrange", color = "red")

## Combine with box plot

e + geom_boxplot() + 
  geom_dotplot(binaxis = "y", stackdir = "center")

## Combine with violin plot

e + geom_violin(trim = FALSE) +
  geom_dotplot(binaxis = 'y', stackdir = "center")

## Combine dot plot + violin plot + stat summary

e + geom_violin(trim = FALSE) +
  geom_dotplot(binaxis = 'y', stackdir = "center") +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "pointrange", color = "red")

### Changing colours by groups

# single color
e + geom_dotplot(binaxis = 'y', stackdir = 'center', color = "black",
                 fill = "lightgray") + theme_minimal()
# Change outline color by dose(groups)
e + geom_dotplot(aes(color = dose),  binaxis = 'y', stackdir = 'center',
                 fill = "white") + theme_minimal()

#Change fill color by dose (groups)
e + geom_dotplot(aes(fill = dose), binaxis = 'y', stackdir = 'center') +
  theme_minimal()

### Change manually outline colors
# use custom color pallettes

e2 <- e + geom_dotplot(aes(color = dose), binaxis = 'y',
                       stackdir = 'center', fill = "white") + theme_minimal()
e2 + scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9"))

# Use custom fill pallettes

e3 <- e + geom_dotplot(aes(fill = dose), binaxis = 'y',
                       stackdir = "center") + theme_minimal()
e3 + scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))

# Dot plot with multiple groups

e + geom_dotplot(aes(fill = supp), binaxis = 'y', stackdir = 'center',
                 position = position_dodge(0.8))

# change colors

e + geom_dotplot(aes(fill = supp), binaxis = 'y', stackdir = 'center',
                 position = position_dodge(0.8)) +
  scale_fill_manual(values = c("#999999", "#E69F00"))

# Add box plots

e + geom_boxplot(fill = "white") +
  geom_dotplot(aes(fill = supp), binaxis = 'y', stackdir = 'center')

# Change the position

e + geom_boxplot(aes(fill = supp), position = position_dodge(0.8)) +
  geom_dotplot(aes(fill = supp), binaxis = 'y', stackdir = 'center',
               position = position_dodge(0.8))


```

```{r}
# Strip Charts (For small sample sizes)
## One dimensional scatter plots

e + geom_jitter(position = position_jitter(0.2), shape = 17, size = 1.2)

## Adding summary statistics  :: stat_summary()
### Adding mean and median point : use fun.y = mean or fun.y = median

e + geom_jitter(position = position_jitter(0.2)) + 
  stat_summary(fun.y = mean, geom = "point",
               shape = 18, size = 3, color = "red")

## Adding mean points +/- Standard Deviation
### Use geom = "pointrange" or geom = "crossbar"

e + geom_jitter(position = position_jitter(0.2)) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult=1),
               geom = "pointrange", color = "red")

### Combine with boxplot

e + geom_boxplot() +
  geom_jitter(position = position_jitter(0.2))

### Combine with violin plot

e + geom_violin() +
  geom_jitter(position = position_jitter(0.2))

### Strip chart + violin plot + stat summary

e + geom_violin(trim = FALSE) +
  geom_jitter(position = position_jitter(0.2)) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "pointrange", color = "red")

### Change point shapes by groups

e + geom_jitter(aes(shape = dose), position = position_jitter(0.2))

### Change point shapes manually

e + geom_jitter(aes(shape = dose), position = position_jitter(0.2)) +
  scale_shape_manual(values = c(1,17,19))

### Change colors using single color

e + geom_jitter(position = position_jitter(0.2), color = "steelblue") + theme_minimal()

### Change colors using group

e + geom_jitter(aes(color = dose), position = position_jitter(0.2)) +
  theme_minimal()

### Changing colors manually

e3 <- e + geom_jitter(aes(color = dose), position = position_jitter(0.2)) + theme_minimal()
e3 + scale_fill_manual(values = c("#999999","#E69F00","#56B4E9"))

### Using brewer color palettes

e3 + scale_color_brewer(palette = "Dark2")

### Strip charts with multiple groups
#### Changing colors and shapes by groups

e + geom_jitter(aes(color = supp, shape = supp),
                position = position_jitter(0.2))

#### Changing the position:: interval between dotplot of the same group

e + geom_jitter(aes(color = supp, shape = supp),
                position = position_jitter(0.2))

### Changing colors + adding box plots + changing position

e + geom_jitter(aes(color = supp, shape = supp),
                position = position_jitter(0.2)) +
  scale_color_manual(values = c("#999999", "#E69F00"))

e + geom_boxplot(color = "black")+ 
  geom_jitter(aes(color = supp, shape = supp),
              position = position_jitter(0.2))

e + geom_boxplot(aes(color = supp), position = position_dodge(0.8)) +
  geom_jitter(aes(color = supp, shape = supp),
              position = position_dodge(0.8))


```

```{r}
# Line plots [ geom_line, geom_step, geom_path]

df <- data.frame(dose = c("D0.5", "D1", "D2"),
                 len = c(4.2,10,29.5))
head(df)


df2 <- data.frame(supp = rep(c("VC", "OJ"), each = 3),
                  dose = rep(c("D0.5", "D1", "D2"),2),
                  len = c(6.8, 15, 33, 4.2, 10, 29.5))
head(df2)

## Basic line plots

p <- ggplot(data = df, aes(x=dose, y= len, group = 1))
p + geom_line() + geom_point()

## Changing line type and color

p + geom_line(linetype = "dashed", color = "steelblue") +
  geom_point(color = "steelblue")

## Use geom_step()

p + geom_step() + geom_point()

## Line plot with multiple groups

p <- ggplot(df2, aes(x= dose, y = len, group = supp))

p + geom_line(aes(linetype = supp)) +
  geom_point(aes(shape = supp))

p + geom_line(aes(linetype = supp, color = supp)) +
  geom_point(aes(shape = supp, color = supp))

p + geom_line(aes(linetype = supp, color = supp)) +
  geom_point(aes(shape = supp, color = supp)) +
  scale_color_manual(values = c("#999999", "#E69F00"))


## Line with numerical X- axis

df3 <- data.frame(supp = rep(c("VC","OJ"), each = 3),
                  dose = rep(c("0.5", "1", "2"), 2),
                  len = c(6.8, 15, 33, 4.2, 10, 29.5))
df3

### X- axis treated as continuous variable

df3$dose <- as.numeric(as.vector(df3$dose))
ggplot(data = df3, aes(x = dose, y = len, group = supp, color = supp)) + geom_line() + geom_point()

### X- axis treated as discrete variable

df2$dose <- as.factor(df3$dose)
ggplot(data = df2, aes(x = dose, y = len, group = supp, color = supp)) + geom_line() + geom_point()

## Line plot with dates on X- axis: Time Series
head(economics)

### Basic line plot 
ggplot(data = economics, aes(x = date, y = pop)) +
  geom_line()

### Plotting a subset of data

ss <- subset(economics, date > as.Date("2006-1-1"))
ggplot(data = ss, aes(x = date, y = pop)) + geom_line()

### Changing line size

ggplot(data = economics, aes(x = date, y = pop, size = unemploy/pop)) + geom_line()

## Using Multiple time series data

ggplot(economics, aes(x = date)) +
  geom_line(aes(y = psavert), color = "darkred") +
  geom_line(aes(y = uempmed), color = "steelblue", linetype = "twodash") +
  theme_minimal()

## Melt chart
library(reshape2)
df <- melt(economics [, c("date", "psavert", "uempmed")], id="date")
ggplot(df, aes(x=date, y=value)) +
  geom_line(aes(color = variable), size = 1) +
  scale_color_manual(values = c("#999999", "#E69F00")) +
  theme_minimal()

## Area plot
ggplot(economics, aes(x = date)) +
  geom_area(aes(y = psavert), fill = "#999999",
            color = "#999999", alpha = 0.5) +
  geom_area(aes(y = uempmed), fill = "#E69F00",
            color = "#E69F00", alpha = 0.5) +
  theme_minimal()
```

```{r}
# Bar plots

df <- data.frame(dose = c("D0.5", "D1", "D2"),
                 len = c(4.2, 10, 29.5))

head(df)

df2 <- data.frame(supp = rep(c("VC", "OJ"), each = 3),
                  dose = rep(c("D0.5", "D1", "D2"), 2),
                  len = c(6.8, 15, 33, 4.2, 10, 29.5))
head(df2)

## Basic plot
f <- ggplot(df, aes(x = dose, y = len))
f + geom_bar(stat = "identity")

## Change fill color and addd labels to the top (vjust = -0.3)

f + geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = len), vjust = -0.3, size = 3.5) +
  theme_minimal()

## Label inside the bar ['vjust' +ve value inside and -ve value outside plot]

f + geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = len), vjust = 1.6, color = "white", size = 3.5) + theme_minimal()

## Changing colours by groups
### Plot line colours
f + geom_bar(aes(color = dose), stat = "identity", fill = "white")

### Plot fill colours
f + geom_bar(aes(fill = dose), stat = "identity")

### Change outline color with custom color (manually)

f + geom_bar(aes(color = dose), stat = "identity", fill = "white") +
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9"))

### Change fill color with custom color (manually)

f + geom_bar(aes(fill = dose), stat = "identity") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))

## Barplots with multiple groups

g <- ggplot(data = df2, aes(x = dose, y = len, fill = supp))

### Stacked bar plot
g + geom_bar(stat = "identity")

### Use position = position_dodge()

g + geom_bar(stat = "identity", position = position_dodge())

### Labelling dodged bar plot

ggplot(data = df2, aes(x = dose, y = len, fill = supp)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = len), vjust = 1.6, color = "white",
            position = position_dodge(0.9), size = 3.5)

## Adding labels to stacked plots [3 steps]
library(plyr)
### 1. Sorting data 
df_sorted <- arrange(df2, dose, supp)
df_sorted

### 2. Calculating cumulative sum of len for each dose

df_cumsum <- ddply(df_sorted, "dose", transform,
                   label_ypos = cumsum(len))
head(df_cumsum)

### 3. Create the bar plot

ggplot(data = df_cumsum, aes(x = dose, y = len, fill = supp)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = label_ypos, label = len), vjust = 1.6, 
            color = "white", size = 3.5)

## To keep the numbers in the middle of the bar

df_cumsum <- ddply(df_sorted, "dose", transform,
                   label_ypos = cumsum(len) - 0.5*len)

### Create the bar plot
ggplot(data = df_cumsum, aes(x = dose, y = len, fill = supp)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = label_ypos, label = len), vjust = 1.6,
            color = "white", size = 3.5)
```

```{r}
# Visualizing Error [Need to Re-do]

## Data format

df <- ToothGrowth
df$dose <- as.factor(df$dose)
head(df,3)

## Computing mean and standard deviation
library("dplyr")
df2 <- df %>%
  group_by(dose) %>%
  summarise(
    sd = sd(len),
    mean = mean(len)
  )
df2
```

```{r}
# Pie Charts

df <- data.frame(
  group = c("Male", "Female", "Child"),
  value = c(25, 25, 50))
head(df)

## Default plot
p <- ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0)
p

## Use custom fill color palettes
p + scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))

## Customized pie charts

blank_theme <- theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 14, face = "bold")
  )

## Apply blank theme
require(scales)
p + scale_fill_brewer("Blues") + blank_theme +
  geom_text(aes(y = value/3 + c(0, cumsum(value) [-length(value)]),
                label = percent(value/1000)), size = 5)

```

```{r}
# Graphical Primitives
require(maps)
france = map_data('world', region = "France")
ggplot(france, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = 'white', color = "Black")

## Producing path, ribbons and rectangles

h <- ggplot(economics, aes(date, unemploy))

### path
h + geom_path(size = 0.8, color = "#E46726") +
  theme_minimal()

### Combining path, ribbon and rectangle

h + geom_rect(aes(xmin = as.Date('1980-01-01'), ymin = -Inf,
                  xmax = as.Date('1985-01-01'), ymax = Inf),
              fill = "#A29B32", color = "#D8DA9E", size = 1.5) +
  geom_ribbon(aes(ymin = unemploy - 900, ymax = unemploy + 900),
              gill = "#F3BF94") +
  geom_path(size = 0.8, color = "#E46726") +
  theme_minimal()

### Add line segements and curves between points (x1, y1) and (x2, y2): 

# Create scatter plot

i <- ggplot(mtcars, aes(wt, mpg)) +
  geom_point()

## Add segment
i + geom_segment(aes(x = 2, y = 15, xend = 3, yend = 15))

## Add arrow
require(grid)
i + geom_segment(aes(x = 5, y = 30, xend = 3.5, yend = 25),
                 arrow = arrow(length = unit(0.5, "cm")))

## Add curves
i + geom_curve(aes(x = 2, y = 15, xend = 3, yend = 15))

```

```{r}
# Main titles, axis labels and legend title

## Convert the variable dose from numeric to factor variables

ToothGrowth$dose <- as.factor (ToothGrowth$dose)
p <- ggplot(ToothGrowth, aes(x = dose, y = len, fill = dose)) +
  geom_boxplot()

p + ggtitle("Main title")
p + xlab("X axis label")
p + ylab("Y axis label")
p + labs(title = "main title", x = "X axis label", y = "Y axis label") 

## Change the appearance of labels
### Values for face are one of "plain", "italic", "bold", "bold.italic"

p + theme(
  plot.title = element_text(color = "red", size = 12, face = "bold.italic"),
  axis.title.x = element_text(color = "blue", size = 12, face = "bold"),
  axis.title.y = element_text(color = "#999999", size = 12, face = "bold"))

### Hide labels

p + theme(plot.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())

## Change legend titles [labs() function]

p + labs(fill = "Dose (mg)")

```

```{r}
# Legend Position and Appearance

## Create box plot
ToothGrowth$dose <- as.factor(ToothGrowth$dose)
p <- ggplot(ToothGrowth, aes (x = dose, y = len, fill = dose)) +
  geom_boxplot()

### Changing legend position and appearance
#### Chagne legend position: "left", "top", "right", "bottom", "none"

p + theme(legend.position = "top")

#### Legend position as numeric vector c(x, y)

p + theme(legend.position = c(0.8, 0.2))

#### Remove legends
p + theme (legend.position = "none")

#### Change the appearnace of legend title and labels
p + theme (legend.title = element_text(color = "blue"),
           legend.text = element_text(color = "red"))

#### Change legend box background color
p + theme (legend.background = element_rect(fill = "lightblue"))

#### Change the order of legend items
p + scale_x_discrete(limits = c("2", "0.5", "1"))

#### set legend title and labels
p + scale_fill_discrete(name = "Dose", labels = c("A", "B", "C"))

## Using Guides : set or remove the legend for a specific aesthetic

### Prepare data
#### 1. Convert cyl and gear to factor variables
mtcars$cyl <- as.factor(mtcars$cyl)
mtcars$gear <- as.factor(mtcars$gear)

#### 2. Create a scatter plot with multiple aesthetics (guides)
##### Plot with multiple aesthetics
p <- ggplot(data = mtcars, 
            aes(x = mpg, y = wt, color = cyl, size = qsec, shape = gear)) +
  geom_point()
p
##### Chagne the order of gudies using guide_legend()

p + guides(color = guide_legend(order = 1),
           size = guide_legend(order = 2),
           shape = guide_legend(order = 3))
##### Remove a legend for a particular aesthetic (color and size)
p + guides(color = FALSE, size = FALSE)

##### In case of continous color; use guide_colourbar()
qplot(data = mpg, x = displ, y = cty, size = hwy,
      color = cyl, shape = drv) +
  guides(color = guide_colorbar(order = 1),
         alpha = guide_legend(order = 2),
         size = guide_legend(order = 3))

##### Removing legend for the point shape

p + scale_shape(guide = FALSE)

##### Remove legend for size

p + scale_size(guide = FALSE)

##### Remove legend for colors

p + scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9"),
                       guide = FALSE)
```

```{r}
# Colors

## Convert dose and cyl columns from numeric to factor variables
ToothGrowth$dose <- as.factor(ToothGrowth$dose)
mtcars$cyl <- as.factor(mtcars$cyl)

## Box plot
bp <- ggplot(ToothGrowth, aes(x = dose, y = len))
bp + geom_boxplot(fill = 'steelblue', color = "red")

## scatter plot
sp <- ggplot(mtcars, aes(x = wt, y = mpg))
sp + geom_point(color = 'darkblue')

## Changing color by group
bp + geom_boxplot(aes(fill = dose))
sp + geom_point(aes(color = cyl))

## Manual Hues
bp + geom_boxplot(aes(fill = dose))+ scale_fill_hue(l = 40, c = 35)
sp + geom_point(aes(color = cyl)) + scale_color_hue(l = 40, c = 35)

## Changing colors manually
bp + geom_boxplot(aes(fill = dose)) +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))

sp + geom_point(aes(color = cyl)) + 
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9"))

## R Color Brewer Pallette
bp + geom_boxplot(aes(fill = dose)) +
  scale_fill_brewer(palette = "Dark2")

sp + geom_point(aes(color = cyl)) +
  scale_color_brewer(palette = "Dark2")

## Using Wes Anderson Color Palettes
library(wesanderson)

bp + geom_boxplot(aes(fill = dose)) +
  scale_fill_manual(values = wes_palette(n = 3, name = "GrandBudapest"))

sp + geom_point(aes(color = cyl)) +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest"))

## Using grey scale

bp + geom_boxplot(aes(fill = dose)) +
  scale_fill_grey(start = 0.8, end = 0.2) + theme_minimal()

sp + geom_point(aes(color = cyl)) +
  scale_color_grey(start = 0.8, end = 0.2) + theme_minimal()

## Gradient and continuous colors

### Color by qsec values

sp2 <- ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point(aes(color = qsec))
sp2

### Change the low and high colors
#### Sequential color scheme

sp2 + scale_color_gradient(low = "blue", high = "red")

#### Diverge color scheme

mid <- mean(mtcars$qsec)
sp2 + scale_color_gradient2(midpoint = mid, low = "blue", mid = "white", high = "red", space = "Lab")

## Gradient between n colors

sp3 <- ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point(aes(color = mpg))

sp3 + scale_colour_gradientn(colors = rainbow(5))

```

```{r}
# Point shapes, colors and size

## Convert cyl as factor variable
mtcars$cyl <- as.factor(mtcars$cyl)

## Changing shape, color and size
ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point(shape = 18, color = "steelblue", size =4)

## Chaing shape, color, fill and size
ggplot(mtcars, aes(x= wt, y = mpg)) +
  geom_point(shape = 23, fill = "blue",
             color = "darkred", size = 3)

## Change point shapes and colors by groups
ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point(aes(shape = cyl, color = cyl))

## Controlling the size by groups
ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point(aes(size = cyl))

## Change colors and shapes manually
ggplot(mtcars, aes(x = wt, y = mpg, group = cyl)) +
  geom_point(aes(shape = cyl, color = cyl), size = 2) +
  scale_shape_manual(values = c(3, 16, 17)) +
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(legend.position = "top")

## Change colors, shapes and point size manually
ggplot(mtcars, aes(x = wt, y = mpg, group = cyl)) +
  geom_point(aes(shape = cyl, color = cyl, size = cyl)) +
  scale_shape_manual(values = c(3, 16, 17)) +
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  scale_size_manual(values = c(1.5, 2, 3)) +
  theme(legend.position = "top")
```

```{r}
# Line types

## Basic dash plot
df <- data.frame(time = c("breakfast", "Lunch", "Dinner"),
                 bill = c(10, 30, 15))
head(df)

ggplot(data = df, aes(x = time, y = bill, group = 1)) +
  geom_line(linetype = "dashed") +
  geom_point()

## Line point with multiple groups
df2 <- data.frame(sex = rep(c("Female", "Male"), each = 3),
                  time = c("breakfast", "Lunch", "Dinner"),
                  bill = c(10, 30, 15, 13, 40, 17))
head(df2)

ggplot(df2, aes(x = time, y = bill, group = sex)) +
  geom_line(aes(linetype = sex, color = sex)) +
  geom_point(aes(color = sex))

## Change line types, colors and sizes

ggplot(df2, aes(x = time, y = bill, group = sex)) +
  geom_line(aes(linetype = sex, color = sex, size = sex)) +
  geom_point() +
  scale_linetype_manual(values = c("twodash", "dotted")) +
  scale_color_manual(values = c("#999999", "#E69F00")) +
  scale_size_manual(values = c(1,1.5))
```

```{r}
# Axis limits: Minimum and Maximum Values

data(cars)
p <- ggplot(cars, aes(x = speed, y = dist)) + geom_point()
print(p)

## Change axis limits using coord_cartesian()
p + coord_cartesian(xlim = c(5,20), ylim = c(0,50))

## Set the intercept of x and y axis at (0,0)
p + expand_limits(x = 0, y = 0)

## Using scale_x_x() function:
## Change x and y axis labels and limits
p + scale_x_continuous(name = "speed of cars", limits = c(0,30)) +
  scale_y_continuous(name = "stopping distance", limits = c(0,150))
```

```{r}
# Axis transformations: log and sqrt
data(cars)
p <- ggplot(cars, aes(x = speed, y = dist)) +
  geom_point()
print(p)

## Log transformation using scale_xx()
## Possible values for trans: 'log2', 'log10', 'sqrt'
p + scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2')

### Format axis tick mark labels to show exponents
require(scales)
p + scale_y_continuous(trans = log2_trans(),
                       breaks = trans_breaks("log2", function(x) 2^x),
                       labels = trans_format("log2", math_format(2^.x)))

### Reverse coordinates
p + scale_y_reverse()

### Percent
p + scale_y_continuous(labels = percent)

### Dollar
p + scale_y_continuous(labels = dollar)

### Scientific
p + scale_y_continuous(labels = scientific)

## Log tick marks: annotation_logticks()
require(MASS)
require (scales)
data ("Animals")

### x and y axis are transformed and formatted
p2 <- ggplot(Animals, aes(x = body, y = brain)) +
  geom_point() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()

### log-log plot without log tick marks
p2

### show log tick marks
p2 + annotation_logticks()

### Log ticks on left and right
p2 + annotation_logticks(sides = "lr")

### All sides
p2 + annotation_logticks(sides="trbl")
```

```{r}
# Date Axes
set.seed(1234)
df <- data.frame(
  date = seq(Sys.Date(), len = 100, by="1 day") [sample(100,50)],
  price = runif(50)
)
df <- df[order(df$date),]
head(df)

p <- ggplot(data = df, aes(x = date, y = price)) +
  geom_line()

## Formate axis tick mark labels
### Format : month/day
require(scales)
p + scale_x_date(labels = date_format("%m/%d")) +
  theme(axis.text.x = element_text(angle = 45))

### Format: week
p + scale_x_date(labels = date_format("%w"))

### Format: Months only
p + scale_x_date(breaks = date_breaks("months"),
                 labels = date_format("%b")) +
  theme(axis.text.x = element_text(angle = 45))

## Date axis limits
data("economics")
p <- ggplot(data = economics, aes(x = date, y = psavert)) +
  geom_line(color = "steelblue")
p

### Axis limits c(min, max)
min <- as.Date("2002-1-1")
max <- max(economics$date)
p + scale_x_date(limits = c(min, max))
```

```{r}
# Axis Ticks: Customize tick marks and labels, reorder and select items - element_text(face, color, size, angle):change text style - element_blank()

data("ToothGrowth")
ToothGrowth$dose <- as.factor(ToothGrowth$dose)
p <- ggplot(ToothGrowth, aes(x = dose, y = len)) +
  geom_boxplot()

## Changing style of axis tick labels
### face can be "plain", "italic", "bold", "bold.italic"
p + theme(axis.text.x = element_text(face = "bold", color = "#993333", size = 12, angle = 45),
          axis.text.y = element_text(face = "bold", color = "blue",
                                     size = 12, angle = 45))

### Remove axis ticks and tick mark labels
p + theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank())

### Change the line type and color of axis lines
p + theme(axis.line = element_line(color = "darkblue",
                                   size = 1, linetype = "solid"))

## Discrete axes
p + scale_x_discrete(name = "Dose (mg)",
                     limits = c("2", "1", "0.5"))
### Change tick mark labels
p + scale_x_discrete(breaks = c("0.5","1","2"),
                     labels = c("D0.5","D1","D2"))
### Choose which items to display
p + scale_x_discrete(limits=c("0.5", "2"))

#### Alternative # p + xlim("0.5", "2")

## Continuous axes  
sp <- ggplot(cars, aes(x = speed, y = dist)) + geom_point()
sp
### Change x and y axis labels, and limits
sp + scale_x_continuous(name = "speed of cars", limits = c(0,30)) +
  scale_y_continuous(name = "stopping distance", limits =c(0,150))

### Set marks on y-axis: a tick mark is shown on every 50
sp + scale_y_continuous(breaks = seq(0, 150, 50))

### Tick marks can be spaced randomly
sp + scale_y_continuous(breaks = c(0, 50, 65, 75, 150))

### Removing tick mark labels and gridlines
sp + scale_y_continuous(breaks = NULL)

## Format the labels
require(scales)
### possible values for labels = percent, scientific, .....
sp + scale_y_continuous(labels = percent)
```

```{r}
# Themes and Background Colors
data ("ToothGrowth")

ToothGrowth$dose <- as.factor(ToothGrowth$dose)

p <- ggplot(ToothGrowth, aes(x = dose, y = len)) +
  geom_boxplot()
p + theme_gray()
p + theme_gray(base_size = 14)

p + theme_bw()
p + theme_linedraw()
p + theme_light()
p + theme_minimal()
p + theme_classic()
p + theme_void()
p + theme_dark()

## Customize the plot background

p + theme(panel.background = element_rect(fill = "#BFD5E3", color = "#6D9EC1", size = 2, linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "white"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid', color = "white"))

### Change the plot background color
p + theme(plot.background = element_rect(fill = "#BFD5E3"))

### Remove panel borders and grid lines

p + theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(size = 0.5, linetype = "solid",
                                   color = "black"))

## Using ggthemes
library(ggthemes)

p <- ggplot(ToothGrowth, aes(x = dose, y = len)) +
  geom_boxplot(aes(fill = dose))

### theme economist
p + theme_economist() + scale_fill_economist()

## theme stata

p + theme_stata() + scale_fill_stata()

##Creating a Theme
### 1. Set theme for the current session
theme_set(theme_grey(base_size = 20))

## theme_gray()   ----> needed to reset

## function(base_size = 11, base_family = "") { half_line <- base_size/2 theme( axis.text = element_text(size = rel(0.8), color = "grey30"),)}
```

```{r}
# Text Annotations

set.seed(1234)
df <- mtcars[sample (1: nrow(mtcars), 10),]
df$cyl <- as.factor(df$cyl)

## Simple scatter plot
sp <- ggplot(df, aes(wt, mpg, label = rownames(df))) +
  geom_point()

sp + geom_text(aes(label = rownames(df), color = cyl),
               size = 3, vjust = -1)

sp + geom_text(x = 3, y = 30, label = "Scatter Plot",
               color = "red", fontface = 2)

sp + geom_label()

## Annotation_custom : Add a static text annotation

sp2 <- ggplot(mtcars, aes(x = wt, y = mpg, label = rownames(mtcars))) +
  geom_point()

library(grid)
grob <- grobTree(textGrob("scatter plot", x = 0.1, y = 0.95, hjust = 0, gp = gpar(col = "red", fontsize = 13, fontface = "italic")))

sp2 + annotation_custom(grob) +
  facet_wrap(~ cyl, scales = "free")

## ggrepel : Avoid overlapping of text labels

library(ggrepel)
set.seed(1234)
ss <- sample(1:32, 15)
df <- mtcars[ss,]

p <- ggplot(df, aes(wt, mpg)) +
  geom_point(color = 'red') +
  theme_minimal(base_size = 10)

p + geom_text(aes(label = rownames(df)),
              size = 3.5)

## Using geom_text_repel
require(ggrepel)
set.seed(42)
p + geom_text_repel(aes(label = rownames(df)),
                    size = 3.5)

## Using ggrepel:: geom_lable_repel and change color by group

set.seed(42)
p + geom_label_repel(aes(label = rownames(df),
                         fill = factor(cyl)), color = 'white',
                     size = 3.5) +
  theme(legend.position = "bottom")

```

```{r}
# Add Straight Lines to a plot: Horizontal, Vertical and Regression lines

sp <- ggplot(data = mtcars, aes(x= wt, y = mpg)) + geom_point()

## Horizontal Segment of Line (y = 20)

sp + geom_hline(yintercept = 20, linetype = "dashed", color = "red")

## Vertical Segment of Line (x = 3)

sp + geom_vline(xintercept = 3, color = "blue", size = 1.5)

## Add regression line

sp + geom_abline(intercept = 37, slope = -5, color = "blue") +
  ggtitle("y = -5x + 37")

## Add Vertical line segment

sp + geom_segment(aes(x = 4, y = 15, xend = 4, yend = 27))

## Add Horizontal line segment

sp + geom_segment(aes(x = 2, y = 15, xend = 3, yend = 15))

## Add arrow at the end of the segment
require(grid)
sp + geom_segment(aes(x = 5, y = 30, xend = 3.5, yend = 25),
                  arrow = arrow(length = unit(0.5, "cm")))

```

```{r}
# Rotate a Plot : Flip and Reverse

hp <- qplot(x = rnorm(200), geom = "histogram")
hp

## Horizontal histogram
hp + coord_flip()

## Y- axis reversed

hp + scale_y_reverse()

```

```{r}
# Facets: Split a Plot into a Matrix of Panels

data("ToothGrowth")
ToothGrowth$dose <- as.factor(ToothGrowth$dose)

p <- ggplot(ToothGrowth, aes(x = dose, y = len, group = dose)) +
  geom_boxplot(aes(fill = dose))
p

p + facet_grid(supp~.) ## Facet in vertical direction
p + facet_grid(.~supp) ## Facet in horizontal direction
p + facet_grid(dose ~ supp) ## Facet in horizontal and vertical direction

p + facet_grid(dose ~ supp, margins = TRUE)

## Facet Scales

p + facet_grid(dose ~ supp, scales = 'free')

p + facet_grid(dose ~ supp, labeller = label_both)

## Customization:

p + facet_grid(dose ~ supp) +
  theme(strip.text.x = element_text(size = 12, color = "red",
                                    face = "bold.italic"),
        strip.text.y = element_text(size = 12, color = "red",
                                    face = "bold.italic"))

p + facet_grid(dose ~ supp) +
  theme(strip.background = element_rect(color = "black", fill = "white", size = 1.5, linetype = "solid"))


p + facet_wrap(~dose) ## Facet Wrapping

p + facet_wrap(~dose, ncol = 2)

```

```{r}
# Position Adjustements

p <- ggplot(mpg, aes(fl, fill = drv))

p + geom_bar(position = "dodge")

p + geom_bar(position = "fill")

p + geom_bar(position = "stack")

ggplot(mpg, aes(cty, hwy)) +
  geom_point(position = "jitter")

p + geom_bar(position = position_dodge(width = 1))
```

```{r}
# Coordinate Systems

p <- ggplot(mpg, aes(fl)) + geom_bar()

p + coord_cartesian(xlim = NULL, ylim = c(0, 200))

p + coord_fixed(ratio = 1/50, xlim = NULL, ylim = NULL)

p + coord_flip()

p + coord_polar(theta = "x", direction = 1)

p + coord_trans(y = "sqrt")

```

```{r}
# Extensions to ggplot2
## Arrange multiple graphs on the same page
library(gridExtra)
library(cowplot)

data("ToothGrowth")
ToothGrowth$dose <- as.factor(ToothGrowth$dose)

data("economics")
data("diamonds")

## Define a custom color set
my3cols <- c("#E7B800", "#2E9FDF", "#FC4E07")

require(cowplot)
p <- ggplot(ToothGrowth, aes(x = dose, y = len))

### Box Plot
bxp <- p + geom_boxplot(aes(color = dose)) +
  scale_color_manual(values = my3cols)
bxp
#### Adding grid lines
bxp + background_grid(major = "xy", minor = "none")
#### Adding grey theme
bxp + theme_gray()


### Dot Plot
dp <- p + geom_dotplot(aes(color = dose, fill = dose),
                       binaxis = 'y', stackdir = 'center') +
  scale_color_manual(values = my3cols) +
  scale_fill_manual(values = my3cols)
dp

## Line Plot
lp <- ggplot(economics, aes(x = date, y = psavert)) +
  geom_line(color = "#E46726")
lp

## Combine Multiple Plots
plot_grid(bxp, dp, lp, labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)

## Drawing
ggdraw() +
  draw_plot(bxp, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(dp, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(lp, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C"),
                  x = c(0, 0.5, 0), y = c(1,1,0.5), size = 15)

## Saving multi-figure plots
save_plot("mpg.pdf", bxp,
          base_aspect_ratio = 1.3)

## Saving using plot_grid
plot2by2 <- plot_grid(bxp, dp, lp, labels = c("A", "B", "C"),
                      ncol = 2, nrow = 2)
save_plot("plot2by2.pdf", plot2by2,
          ncol = 2,
          nrow = 2,
          base_aspect_ratio = 1.3)

## gridExtra Package - grid.arrange() : Arrange Multiple plots on a page
my5cols <- c("#6D9EC1", "#646567", "#A29B32", "#E46726", "#F3BF94")
data("diamonds")
brp <- ggplot(diamonds, aes(x = clarity)) +
  geom_bar(aes(fill = cut)) + scale_fill_manual(values = my5cols)

require(gridExtra)
grid.arrange(bxp, dp, lp, brp, ncol = 2, nrow = 2)

## grid.arrange() and arrangeGrob() : Change column/row span of a plot

grid.arrange(bxp, arrangeGrob(dp, brp), ncol =2)

grid.arrange(brp, bxp, dp, ncol = 2, nrow = 2,
             layout_matrix = rbind(c(1,1), c(2,3)))

## Use common legend for multiple graphs
## 1. Create the plots: p1, p2,....

get_legend <- function(myggplot){
  require(gridExtra)
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend())
}
## 2. Save the legend of the plot p1 as an external graphic element ("grob")
bxp <- p + geom_boxplot(aes(color = dose)) +
  scale_color_manual(values = my3cols)
dp <- p + geom_dotplot(aes(color = dose, fill = dose),
                       binaxis = 'y', stackdir = 'center') +
  scale_color_manual(values = my3cols) +
  scale_fill_manual(values = my3cols)

## legend <- get_legend(dp)

## 3. Remove the legends from all plots
### bxp2 <- bxp + theme(legend.position = "none")
### dp2 <- dp + theme(legend.position = "none")

## 4. Draw all the plots with only one legend in the right panel

### grid.arrange(bxp2, dp2, legend, ncol = 3, widths = c(2.3, 2.3, 0.8))

## Scatter plot with marginal density plots

my3cols <- c("#6D9EC1", "#646567", "#A29B32")

set.seed(1234)
x <- c(rnorm(350, mean = -1), rnorm(350, mean = 1.5),
       rnorm(350, mean = 4))
y <- c(rnorm(350, mean = -0.5), rnorm(350, mean = 1.7),
       rnorm(350, mean = 2.5))
group <- as.factor(rep(c(1,2,3), each = 350))

df2 <- data.frame(x, y, group)
head(df2)
```






