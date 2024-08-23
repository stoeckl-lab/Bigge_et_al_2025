rm(list = ls())
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2020 03 01
#     MODIFIED:	James Foster              DATE: 2021 03 01
#     MODIFIED: Anna Stoeckl               DATE: 2023 10 16
#
#  DESCRIPTION: Loads a text files, fits a mixed effects logistic model and
#               estimates the p-value for mean ≤ reference mean.
#               
#       INPUTS: A ".xlsx" table with columns for experiment type ("Test"),
#               individual ID ("Moth number"), proportion of correct choices ("%")
#               total number of choices ("total Choice number"), correct choices
#               ("choices cross"), incorrect choices ("choices circle").
#               User should specify h0 reference correct choice rate (line 50).
#               
#      OUTPUTS: Plots of confidence intervals (.pdf).
#
#	   CHANGES: 
#             
#             
#
#   REFERENCES: Bates D, Maechler M, Bolker B, Walker S (2015).
#               Fitting Linear Mixed-Effects Models Using lme4.
#               Journal of Statistical Software, 67(1), 1-48.
#               doi:10.18637/jss.v067.i01.
#
#    EXAMPLES:  Fill out user input (line 50), then press ctrl+shift+s to run
#
# 
#TODO   ---------------------------------------------

# Useful functions --------------------------------------------------------
# . Load package ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
require(lme4)#package for fitting GLMMs (generalised linear mixed-effects models)
require(readxl)
require(lmerTest)
require(emmeans)
require(ggplot2)
require(beeswarm)
require(multcomp)
require(broom)
install.packages("ARTool")
require(ARTool)

 
# Input Variables ----------------------------------------------------------

#  .  User input -----------------------------------------------------------
#further below code checks if data contains threshold for choices itself and uses that if so
h1 = 'greater'#alternative hypothesis, mean "greater" than reference or "less"
csv_sep = ','#Is the csv comma separated or semicolon separated? For tab sep, use "\t"


#Check the operating system and assign a logical flag (T or F)
sys_win <- Sys.info()[['sysname']] == 'Windows'
#On computers set up by JMU Würzburg, use user profile instead of home directory
if(sys_win){
  #get rid of all the backslashes
  ltp <- gsub('\\\\', '/', Sys.getenv('USERPROFILE'))#Why does windows have to make this so difficult
}else{#Root directory should be the "HOME" directory on a Mac (or Linux?)
  ltp <- Sys.getenv('HOME')#Life was easier on Mac
}

# . Select files ---------------------------------------------------------


# manually set path file
wd="F:/Anna Backup/AG Stoeckl/Anna&Ronja/Dorsal_Ventral_Imaging/ResultsNew/"; 

# NATURAL DATASETS
# add="/Natural/contrast/"
# filename <- 'contrast_data_all_long_lateral'
# 
# add="/Natural/of/"
# filename <- 'of_data_all_long_lateral'
# 
# reorderCons = c("open", "semi", "closed")
# reorderQuads = c("ventral", "lateral", "dorsal")

# TUNNEL DATASETS

# add="/Tunnel/contrast/"
# filename <- 'contrast_data_all_long'
# plotMax=.00025

# add="/Tunnel/of/"
# filename <- 'of_data_all_long'
# plotMax=2
# 
# reorderCons = c("blank", "down", "left", "right", "up", "StripeDown", "StripeUp")
# reorderQuads = c("ventral", "left", "right", "dorsal")

# 
add="/Tunnel/contrast/"
filename <- 'contrast_data_all_long'
plotMax=.00025

reorderCons = c("down", "left", "right", "up", "StripeUp")
reorderQuads = c("ventral", "left", "right", "dorsal")

# add="/Tunnel/of/"
# filename <- 'of_data_all_long'
# plotMax=2
# 
# # selection for summary data with lateral
# reorderCons = c("down", "left", "right", "Up", "StripeUp")
# reorderQuads = c("ventral", "left", "right", "dorsal")


# Read in file ------------------------------------------------------------
#read in xlsx file with multiple sheets
path_file <- paste0(wd, add,filename,'.xls')
mdata = read_excel(path_file)


# Plot the data -------------------------------------------------------------------

mdata$condition <- factor(mdata$condition , levels=reorderCons) # reorder factors in open, semi, closed
mdata$quadrant <- factor(mdata$quadrant , levels=reorderQuads) # reorder factors in open, semi, closed

setEPS()
postscript(paste0(wd, add,'stackbar_allConds.eps'))
par(mar = c(2, 8, 1, 3))  # Set margins (bottom, left, top, right)
ggplot(mdata, aes(fill=quadrant, y=magnitude, x=condition)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity") +
  ylim(0, plotMax)
dev.off()

setEPS()
postscript(paste0(wd, add,'boxplot_allData.eps'))
par(mar = c(2, 6, 1, 3))  # Set margins (bottom, left, top, right)
with(mdata, boxplot(magnitude ~  condition*quadrant,ylim=c(0,1.05*max(mdata$magnitude))))
with(mdata, beeswarm(magnitude ~ condition*quadrant, add = TRUE,
                        col = 1,   # Color
                        pch = 19,  # Symbol
                        cex = 1))
# abline(h=c(0,1,0.5),lty=c(1,1,2))
dev.off()

setEPS()
postscript(paste0(wd, add,'boxplot_conds.eps'))
with(mdata, boxplot(magnitude ~  quadrant,ylim=c(0,1.05*max(mdata$magnitude)))) # ,ylim=c(0,1)
with(mdata, beeswarm(magnitude ~ quadrant, add = TRUE,
                     col = 1,   # Color
                     pch = 19,  # Symbol
                     cex = 1))
dev.off()

# Model -------------------------------------------------------------------

#make a model with all factors
ctrl_opt = lmerControl(optimizer = 'bobyqa')

lmm.all <- lmer((magnitude) ~ quadrant * condition + (1 | scene), data = mdata,control = ctrl_opt) 
# lmm.all <- glmer(ORresponse ~ Frequency * Condition + (1 | Animal), data = mdata, family = Gamma(link = 'inverse') )


#check normality of residuals of selected model (only relevant for Gaussian models)
qqnorm(residuals(lmm.all))
qqline(residuals(lmm.all))

plot(density(residuals(lmm.all)), col = 'darkblue', main = 'Residuals', xlab = 'Residual size', ylab = 'Frequency Density')
# plot(density(summary(lmm.allbg)$residuals), col = 'darkblue', main = 'Residuals', xlab = 'Residual size', ylab = 'Frequency Density', xlim = c(-10, 10))
lines(rep(0,2), c(0,1), lwd = .25)


#also check that this model explains the data at all
lmm.null <- lmer((magnitude) ~ 1 | scene, data = mdata)#model with only random effects, how well is the data explained by just splitting it into individuals and ignoring all experimental manipulation?
lmm.null <- lmer((magnitude) ~ quadrant + condition + (1 | scene), data = mdata)#model with only random effects, how well is the data explained by just splitting it into individuals and ignoring all experimental manipulation?

AIC(lmm.all)
AIC(lmm.null)


# # for the sake of publications, you can test for change in "deviance"
print(" ")
print("Model Comparison with null model")
print(anova(lmm.all, lmm.null,test="Chisq"))
#this does a chi-squared test on whether deviance between the model and the data is significantly lower for lmm.all
#the degrees of freedom are the difference in number of factors (1)
# the statistic is a deviance and it gives a p-value that reviewers like to see
#Now you can do some post-hoc tests

#for completion put here again
resultEMS <- emmeans(lmm.all, ~ quadrant * condition ,type="response")
summary(resultEMS)


pwpp(resultEMS, by = "condition", type = "response")
coef(pairs(resultEMS))

pairResultEMS <-emmeans(lmm.all, list(pairwise ~ quadrant * condition))
print(summary(pairResultEMS))

# print(summary(emmeans(lmm.all, list(pairwise ~ quadrant),type="response")))
# print(summary(emmeans(lmm.all, ~ condition,type="response")))    

#write results to csv
result = as.data.frame(summary(resultEMS))
write.csv(result, "emmeans_all.csv" )
result_all = as.data.frame(summary(pairResultEMS[[2]]))
write.csv(result_all, file = paste0(wd, add,'emmeans_pairwise.csv') )

