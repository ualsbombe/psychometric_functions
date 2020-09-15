
# CLEAN UP AND SET DIRS -------------------------------------------------------------------------------------------

rm(list = ls())

results.dir <- '/home/lau/Dropbox/paradigmer_i_aarhus/expectations/data_files' ## in Documents now
figures.dir <- '/home/lau/Dropbox/vanderbilt_experiment/article/figures/'
script.dir <- '/home/lau/Dropbox/vanderbilt_experiment/article/r_files/'


setwd(results.dir)

library(matrixStats)
library(lme4)
library(LMERConvenienceFunctions)
library(boot)
library(multcomp)
library(gplots)
library(Hmisc)
library(MCMCglmm)

r2.corr.mer <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
}

omega2 <- function(m) {
    omega <- 1 - var(residuals(m)) / (var(model.response(model.frame(m))))
    return(omega)
} ## https://stats.stackexchange.com/questions/95054/how-to-get-an-overall-p-value-and-effect-size-for-a-categorical-factor-in-a-mi


# SUBJECTS --------------------------------------------------------------------------------------------------------

subject.names <- c(
    
                    '001.v4_2016_jun_20_1556.csv',
                    '002_v4_2016_jun_22_0901.csv',
                    '003_v4_2016_jun_22_1013.csv',
                    '004_v4_2016_jun_22_1216.csv',
                    '005_v4_2016_jun_22_1531.csv',
                    '006_v4_2016_jun_23_0910.csv',
                    '007_v4_2016_jun_23_1036.csv',
                    '008_v4_2016_jun_23_1418.csv',
                    '009_v4_2016_jun_23_1558.csv',
                    '010_v4_2016_jun_28_0923.csv',
                    '011_v4_2016_jun_28_1041.csv',
                    '012_v4_2016_jun_28_1157.csv',
                    '013_v4_2016_jun_28_1333.csv',
                    '014_v4_2016_jun_29_1023.csv',
                    '015_v4_2016_jun_28_1506.csv',
                    '016_rigtig_v4_2016_jun_30_1155.csv',
                    '017_v4_2016_jun_29_1332.csv',
                    '018_v4_2016_jul_01_0840.csv',
                    '019_v4_2016_okt_12_0831.csv',
                    '020_v4_2016_okt_12_1006.csv',
                    '021_v4_2016_okt_12_1122.csv',
                    '022_v4_2016_okt_12_1258.csv',
                    '023v _v4_2016_okt_12_1548.csv',
                    '024 _v4_2016_okt_13_0956.csv',
                    '025 _v4_2016_okt_13_1122.csv',
                    '026 _v4_2016_okt_13_1304.csv',
                    '027_v4_2016_okt_13_1419.csv',
                    '028_v4_2016_okt_13_1552.csv',
                    '029_v4_2016_dec_13_1806.csv'
)


# LOAD DATA -------------------------------------------------------------------------------------------------------

for(subject.name in subject.names)
{
    subject.index <- which(subject.name == subject.names)
    dt.subject <- read.csv(subject.name, header = TRUE)
    dt.subject <- within(dt.subject,
                 {
                     task <- factor(task, levels = c('singles', 'pairs', 'quadruplet'))
                     acc <- 0
                     acc[(target.type == 'even' & obj.resp == 'e') | (target.type == 'odd' & obj.resp == 'o')] <- 1
                     acc[(target.type == 'even' & obj.resp == 'o') | (target.type == 'odd' & obj.resp == 'e')] <- 0
                     shown.digit <- 0
                     shown.digit[target.type == 'even'] <- even.digit[target.type == 'even']
                     shown.digit[target.type == 'odd'] <- odd.digit[target.type == 'odd']
                     obj.resp <- factor(obj.resp)
                     even.digit <- factor(even.digit)
                     odd.digit <- factor(odd.digit)
                     pas <- factor(pas)
                 })
    
    if(subject.index > 1) dt <- rbind(dt, dt.subject) else dt <- dt.subject

}

# GENERAL PARAMETERS ----------------------------------------------------------------------------------------------

general.par <- par(font.lab = 2, font.axis = 2, lwd = 3)
n.tasks <- length(unique(dt$task))
subjects <- unique(dt$subject)
n.subjects <- length(subjects)
n.pas <- length(unique(dt$pas))

# REACTION TIMES, (OBJECTIVE) ---------------------------------------------------------------------

par(mfrow = c(1, 5))

# single subjects
for(pas in 0:n.pas)
{

    plot(x = 1:n.tasks, y = 1:n.tasks, type = 'n', main = 'Reaction Times per Task',
         ylab = 'Reaction Time (s)', ylim = c(0, 3), xlab = 'Cue Type', xaxt = 'n')
    axis(side = 1, at = 1:n.tasks, levels(dt$task))
    summary.rts <- matrix(nrow = n.subjects, ncol = n.tasks)
    
    for(subject in subjects)
    {
        subject.index = which(subject == subjects)
        if(pas == 0) qualifier <- dt$subject == subject else qualifier <- dt$subject == subject & dt$pas == pas
        dt.subject <- dt[qualifier & dt$trial.type == 'experiment', ]
        summary.rt <- with(dt.subject, tapply(rt.obj, task, median))
        lines(1:n.tasks, summary.rt, type = 'b', col = rgb(0, 0, 0, 0.5))
        summary.rts[subject.index, ] <- summary.rt
    }
    
    # add group mean
    group.rt <- colMeans(summary.rts, na.rm = TRUE)
    lines(1:n.tasks, group.rt, type = 'b', col = 'red')
    
}

# ACCURACY, (OBJECTIVE) ---------------------------------------------------------------------

par(mfrow = c(1, 5))

# single subjects
for(pas in 0:n.pas)
{
    
    plot(x = 1:n.tasks, y = 1:n.tasks, type = 'n', main = paste('Accuracy per Task; PAS:', pas),
            ylab = 'Proportion Corrtect', ylim = c(0.0, 1), xlab = 'Cue Type', xaxt = 'n')
    axis(side = 1, at = 1:n.tasks, levels(dt$task))
    summary.accs <- matrix(nrow = n.subjects, ncol = n.tasks)
    
    for(subject in subjects)
    {
        subject.index = which(subject == subjects)
        if(pas == 0) qualifier <- dt$subject == subject else qualifier <- dt$subject == subject & dt$pas == pas
        dt.subject <- dt[qualifier & dt$trial.type == 'experiment', ]
        summary.acc <- with(dt.subject, tapply(acc, task, mean, na.rm = TRUE))
        lines(1:n.tasks, summary.acc, type = 'b', col = rgb(0, 0, 0, 0.5))
        summary.accs[subject.index, ] <- summary.acc
    }
    
    # add group mean
    group.acc <- colMeans(summary.accs, na.rm = TRUE)
    lines(1:n.tasks, group.acc, type = 'b', col = 'red')
    lines(0:4, rep(0.5, 5), lty = 3)
    
}


# NO. OF RATINGS --------------------------------------------------------------------------------------------------

dt.count <- data.frame(count = numeric(), 
                       pas = numeric(),
                       task = numeric(), 
                       subject = numeric(),
                       acc = numeric())

row.index <- 1
for(subject in subjects)
{
    for(task in levels(dt$task))
    {
        for(pas in 1:n.pas)
        {
            ## count
            count <- sum(dt$subject == subject & dt$task == task & dt$pas == pas & dt$trial.type == 'experiment')
            dt.count[row.index, ] <- NA
            dt.count$pas[row.index] <- pas
            dt.count$task[row.index] <- task
            dt.count$subject[row.index] <- subject
            dt.count$count[row.index] <- count
            ## mean accuracy
            dt.count$acc[row.index] <- sum(dt$acc[dt$subject == subject & dt$task == task &
                                                      dt$pas == pas & dt$trial.type == 'experiment']) / count
            row.index <- row.index + 1
        }
    }
}

dt.count <- within(dt.count,
                   {
                       pas <- factor(pas, levels = c(1, 2, 3, 4))
                       task <- factor(task, levels = c('singles', 'pairs', 'quadruplet'))
                       subject <- factor(subject)
                   })


# MODEL DEPENDENT VARIABLES ------------------------------------------------------------------------------------------

acc.model.path <- paste(script.dir, 'acc_models_experiment_2.RData', sep='')
do.models.exist <- file.exists(acc.model.path)

if(!do.models.exist)
{

    model.acc.full <- glmer(acc ~ task * pas + (1 | subject), data = dt, family = 'binomial', verbose = TRUE, subset = dt$trial.type == 'experiment')
    model.acc.additive <- glmer(acc ~ task + pas + (1 | subject), data = dt, family = 'binomial', verbose = TRUE, subset = dt$trial.type == 'experiment')
    model.acc.pas <- glmer(acc ~ pas + (1 | subject), data = dt, family = 'binomial', verbose = TRUE, subset = dt$trial.type == 'experiment')
    model.acc.task <- glmer(acc ~ task + (1 | subject), data = dt, family = 'binomial', verbose = TRUE, subset = dt$trial.type == 'experiment')
    model.acc.intercept <- glmer(acc ~ 1 + (1 | subject), data = dt, family = 'binomial', verbose = TRUE, subset = dt$trial.type == 'experiment')
    
    acc.models <- list(model.acc.full, model.acc.additive, model.acc.pas, model.acc.task, model.acc.intercept)
    save(acc.models, file = acc.model.path)
} else load(acc.model.path)


model.rt.full <- lmer(log(rt.obj) ~ task * pas + (1 | subject), data = dt, verbose = TRUE, subset = dt$trial.type == 'experiment')
model.rt.additive <- lmer(log(rt.obj) ~ task + pas + (1 | subject), data = dt, verbose = TRUE, subset = dt$trial.type == 'experiment')

model.count.full <- glmer(count ~ task * pas + (1 | subject), data = dt.count, family = 'poisson', verbose = TRUE)
model.count.additive <- glmer(count ~ task + pas + (1 | subject), data = dt.count, family = 'poisson', verbose = TRUE)


# TRY OUT ORDINAL LOGISTIC FOR COUNT DATA ------------------------------------------------------------------------------

# dt.multi <- subset(dt, dt$trial.type == 'experiment')
# 
# mn <- MCMCglmm(fixed = pas ~ trait * task, random = ~ us(trait):subject, rcov = ~ us(trait):units, 
#                data = dt.multi, family = 'categorical', verbose = TRUE)




# PLOT MODELS -----------------------------------------------------------------------------------------------------

dt$interaction <- with(dt, interaction(pas, task))
dt.count$interaction <- with(dt.count, interaction(pas, task))

model.acc.plot <- glmer(acc ~ interaction + 0 + (1 | subject), data = dt, family = 'binomial', verbose = TRUE, subset = dt$trial.type == 'experiment')
model.rt.plot <- lmer(log(rt.obj) ~ interaction + 0 + (1 | subject), data = dt, verbose = TRUE, subset = dt$trial.type == 'experiment')
model.count.plot <- glmer(count ~ interaction + 0 + (1 | subject), data = dt.count, family = 'poisson', verbose = TRUE)


models <- list(model.acc.plot, model.rt.plot, model.count.plot)

# INFORMATIVE RATINGS ---------------------------------------------------------------------------------------------

count.model <- models[[3]]
counts <- exp(fixef(count.model))
informative.indices <- list(c(2, 3, 4), c(6, 7, 8), c(10, 11, 12))

for(indices in informative.indices)
{
    print(sum(counts[indices[1:2]]) / sum(counts[indices]))
}


# ACTUALLy PLOT THEM -----------------------------------------------------------------------------------------------


model.index <- 1
y.lims <- list(c(0.4, 1.0), c(0.0, 1.2), c(0, 65))
y.labels <- c('Proportion Correct', 'Reaction Time (s)', 'Number of Ratings')
mains <- c('Accuracy', 'Reaction Time', 'Counts')
names <- c('acc', 'rt', 'freq')
names.args <- c('2', '4', '8')

for(model in models)
{
    figure.path <- paste(figures.dir, names[model.index], '_exp_2.jpeg', sep = '')
    jpeg(filename = figure.path, width = 8.5, height = 9, units = 'cm', res = 600)
    par(mfrow = c(1, 1), cex = 1, font.lab = 2, font.axis = 2)
    effects <- fixef(model)
    if(model.index == 1) matrix.effects <- matrix(inv.logit(effects), nrow = 4, byrow = FALSE)
    if(model.index > 1) matrix.effects <- matrix(exp(effects), nrow = 4, byrow = FALSE)
    
    se <- summary(model)$coefficients[, 2]
    if(model.index == 1)
    {
        conf.up <- inv.logit(effects + qnorm(0.975) * se)
        conf.down <- inv.logit(effects - qnorm(0.975) * se)
    }
    if(model.index > 1)
    {
        conf.up <- exp(effects + qnorm(0.975) * se)
        conf.down <- exp(effects - qnorm(0.975) * se)
    }
    matrix.conf.up <- matrix(conf.up, nrow = 4)
    matrix.conf.down <- matrix(conf.down, nrow = 4)
    h <- barplot2(matrix.effects, beside = TRUE,
             plot.ci = TRUE, ci.l = matrix.conf.down, ci.u = matrix.conf.up, ylim = y.lims[[model.index]],
             xpd = FALSE, names.arg = names.args, 
             main = mains[model.index], ci.lwd = 3,
             ylab = y.labels[model.index], xlab = 'No. of Possible Targets')
    
    if(model.index == 1)
    {
        lines(c(h[1] - 1, h[length(h)] + 1), rep(0.50, 2), lty =3, lwd = 3)

    }
    if(model.index == 3) legend('topleft', paste('PAS', 1:4), cex = 0.5, col = heat.colors(4), lty = 1, text.font = 2)
    model.index <- model.index + 1
    dev.off()
}


# POST HOC COMPARISONS --------------------------------------------------------------------------------------------
n.tests <- 12
K <- matrix(0, nrow = n.tests, ncol = length(effects))
comparisons <- list(c(1, 5), c(1, 9), c(5, 9),
                    c(2, 6), c(2, 10), c(6, 10),
                    c(3, 7), c(3, 11), c(7, 11),
                    c(4, 8), c(4, 12), c(8, 12))
comparison.index <- 1

for(comparison in comparisons)
{
    K[comparison.index, comparison[1]] <- 1
    K[comparison.index, comparison[2]] <- -1
    comparison.index <- comparison.index + 1
}

for(model in models)
{
    gh <- glht(model, linfct = K)
    s.gh <- summary(gh, test = adjusted('none'))
    print(s.gh)
}

## accuracy model
n.tests <- 3
model.test <-  glmer(acc ~ pas + 0 + (1 | subject), data = dt, family = 'binomial', verbose = TRUE, subset = dt$trial.type == 'experiment')
n.effects <- length(fixef(model.test))
gh <- glht(model.test, linfct = K)
s.gh <- summary(gh, test = adjusted('none'))
print(s.gh)

# SUMMARY STATISTICS ----------------------------------------------------------------------------------------------

K <- matrix(0, nrow = n.tests, ncol = n.effects)
K[1, 1] <- -1; K[1, 2] <- 1
K[2, 2] <- -1; K[2, 3] <- 1
K[3, 3] <- -1; K[3, 4] <- 1

summary.model.1 <- lmer(acc ~ pas + task + count + pas:task + pas:count + task:count + pas:task:count + (1 | subject), data = dt.count, verbose = TRUE)
summary.model.2 <- lmer(acc ~ pas + task + count + pas:task + pas:count + task:count + (1 | subject), data = dt.count, verbose = TRUE)
summary.model.3 <- lmer(acc ~ pas + task + count + pas:task + pas:count + (1 | subject), data = dt.count, verbose = TRUE)
summary.model.4 <- lmer(acc ~ pas + task + count + pas:task + (1 | subject), data = dt.count, verbose = TRUE)
summary.model.5 <- lmer(acc ~ pas + task + count + (1 | subject), data = dt.count, verbose = TRUE)
summary.model.6 <- lmer(acc ~ pas + task + (1 | subject), data = dt.count, verbose = TRUE)
summary.model.7 <- lmer(acc ~ pas + (1 | subject), data = dt.count, verbose = TRUE) ## winner
summary.model.8 <- lmer(acc ~ 1   + (1 | subject), data = dt.count, verbose = TRUE)

