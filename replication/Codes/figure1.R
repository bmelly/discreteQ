#Load the dataset for the first application
data <- read.dta(paste0(code.dir, "/Data/data_oregon.dta"))
#outcome variable
doc_num_mod_12m <- data$doc_num_mod_12m
rm(data)

#Load the dataset for the second application
data  <- read.dta(paste0(code.dir, "/Data/cpp_selected.dta"))
#drop if interviewer interviewd only blacks (8 months outcome)
for (i in setdiff(unique(data$interviewer_8months[data$black == 1]),
                  unique(data$interviewer_8months[data$black == 0]))) {
  data <- data[data$interviewer_8months != i,]
}
#outcome variable
stand_mental_score_8months <- data$stand_mental_score_8months
rm(data)

#Load the dataset for the second application
data  <- read.dta(paste0(code.dir, "/Data/cpp_selected.dta"))
#drop if interviewer interviewd only blacks (7 years outcome)
for (i in setdiff(unique(data$interviewer_7years[data$black == 1]),
                  unique(data$interviewer_7years[data$black == 0]))) {
  data <- data[data$interviewer_7years != i, ]
}
#outcome variable
stand_fullscale_iq_7years <- data$stand_fullscale_iq_7years
rm(data)

#Figure 1
pdf(
  paste0(code.dir, "/Results/figure1.pdf"),
  pointsize = 15,
  width = 6.0,
  height = 12.0
)
par(mfrow = c(3, 1))
#Top panel
barplot(
  prop.table(table(doc_num_mod_12m)),
  xlab = "Number of outpatient visits (last 6 months)",
  ylab = "Fraction",
  ylim = c(0, 0.4),
  main = "Panel A: Health care utilization"
)
#Middle panel
prop <- prop.table(table(stand_mental_score_8months))[41:72]
names.prop <- as.numeric(names(prop))
names <- rep(NA, length(prop))
for (i in-3:3){
  names[which.min(abs(names.prop - i))] <- i
}
barplot(
  prop,
  xlab = "Test score",
  ylab = "Fraction",
  names.arg = names,
  ylim = c(0, 0.12),
  main = "Panel B: Test scores for nine-month-olds"
)
#Bottom panel
prop <- prop.table(table(stand_fullscale_iq_7years))[35:110]
names.prop <- as.numeric(names(prop))
names <- rep(NA, length(prop))
for (i in-3:3){
  names[which.min(abs(names.prop - i))] <- i
}
barplot(
  prop,
  xlab = "Test score",
  ylab = "Fraction",
  names.arg = names,
  ylim = c(0, 0.05),
  main = "Panel C: Test scores for seven-year-olds"
)

dev.off()
