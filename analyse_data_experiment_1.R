# CLEAN UP --------------------------------------------------------------------------------------------------------

rm(list = ls())
results.dir <- '/home/lau/Dropbox/vanderbilt_experiment/data_files'
script.dir <- '/home/lau/Dropbox/vanderbilt_experiment/article/r_files/'
figures.dir <- '/home/lau/Dropbox/vanderbilt_experiment/article/figures/'
# results.dir <- '/Users/Lau/Dropbox/VanderbiltNashville/experiment/data_files'
# script.dir <- '/Users/Lau/Dropbox/VanderbiltNashville/experiment/article/r_files/'
# figures.dir <- '/Users/Lau/Dropbox/VanderbiltNashville/experiment/article/figures/'
setwd(results.dir)


library(matrixStats)
library(lme4)
library(LMERConvenienceFunctions)
library(boot)
library(multcomp)
library(gplots)
library(Hmisc)

r2.corr.mer <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
}

omega2 <- function(m) {
    omega <- 1-var(residuals(m))/(var(model.response(model.frame(m))))
    return(omega)
}## https://stats.stackexchange.com/questions/95054/how-to-get-an-overall-p-value-and-effect-size-for-a-categorical-factor-in-a-mi

# READ IN DATA ----------------------------------------------------------------------------------------------------

subject.names <- c(
                'devin_2015_jan_29_1035.csv',
                'youngeun_2015_jan_29_1849.csv',
                'lisbeth_2015_feb_01_1624.csv',
                'co_2015_feb_17_1653.csv',
                'gn_2015_feb_18_1237.csv',
                'js_2015_feb_18_1417.csv',
                'ms_2015_feb_18_1537.csv',
                'LC_2015_feb_18_1642.csv',
                'na_2015_feb_19_0922.csv',
                'md_2015_feb_19_1136.csv',
                'cs_2015_feb_19_1247.csv',
                'ar_2015_feb_19_1407.csv',
                'vr_2015_feb_19_1633.csv',
                'jj_2015_feb_23_0807.csv',
                'ht_2015_feb_23_0925.csv',
                'mh_2015_feb_23_1153.csv',
                'th_2015_feb_23_1430.csv',
                'mz_2015_feb_23_1555.csv',
                'mg_2015_mar_02_0815.csv',
                'ra_pasted_together.csv',
                'sa_2015_mar_02_1314.csv',
                'mh2_2015_mar_02_1423.csv',
                'ao_2015_mar_02_1537.csv',
                'bg_2015_feb_18_0906.csv',
                'ab_2015_feb_23_1307.csv',
                'jr_2015_feb_23_1707.csv',
                'mr_2015_mar_02_0928.csv',
                'pp_pasted_together.csv',
                'jmt_pasted_together.csv'
)

# LOAD DATA -------------------------------------------------------------------------------------------------------

for(subject.name in subject.names)
{
    subject.index <- which(subject.name == subject.names)
    dt.subject <- read.csv(subject.name, header = TRUE)
    dt.subject <- within(dt.subject,
                         {
                             if(subject.index < 1) trial.type <- 'experiment'
                             if(subject.index < 4) seed <- 'None'
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
n.frames <- 6

# NO. OF RATINGS --------------------------------------------------------------------------------------------------

dt.count <- data.frame(count = numeric(), 
                       pas = numeric(),
                       task = numeric(), 
                       subject = numeric(),
                       acc = numeric(),
                       target.frames = numeric())

row.index <- 1
for(subject in subjects)
{
    for(task in levels(dt$task))
    {
        for(pas in 1:n.pas)
        {
            for(frame in 1:n.frames)
             {
                ## count
                count <- sum(dt$subject == subject & dt$task == task & dt$pas == pas & dt$target.frames == frame &
                                 dt$trial.type == 'experiment')
                dt.count[row.index, ] <- NA
                dt.count$pas[row.index] <- pas
                dt.count$task[row.index] <- task
                dt.count$subject[row.index] <- subject
                dt.count$target.frames[row.index] <- frame
                dt.count$count[row.index] <- count
                ## mean accuracy
                dt.count$acc[row.index] <- sum(dt$acc[dt$subject == subject & dt$task == task & dt$target.frames == frame &
                                                          dt$pas == pas & dt$trial.type == 'experiment'], na.rm = TRUE) / count
                
                row.index <- row.index + 1
                
            }

        }
    }
}

dt.count <- within(dt.count,
                   {
                       pas <- factor(pas, levels = c(1, 2, 3, 4))
                       task <- factor(task, levels = c('singles', 'pairs', 'quadruplet'))
                       target.frames <- factor(target.frames)
                       subject <- factor(subject)
                   })
# FIT PSYCHOMETRIC FUNCTION FOR ACCURACY FIXED PARAMETERS  -----------------------------------------------------

par(mfrow = c(1, 1))
yhat <- function(a, b, c, d, x)
{
    y <- a + (b - a) / (1 + exp((c - x) / d))
    return(y)
}

a <- 0.5 ## fix at chance

## initial fit to fix b and d
min.RSS <- function(data, par)
{
    with(data, sum(((a + (par[1] - a) / (1 + exp((par[2] - x) / par[3]))) - y)^2))
}

init.matrix <- matrix(NA, nrow = n.subjects, ncol = 3)
est.pars <- init.matrix
for(subject in subjects)
{
    n.subject <- which(subject == subjects)
    indices <-dt$trial.type == 'experiment' & dt$subject == subject
    data <- data.frame(x = dt$target.frames[indices], y = dt$acc[indices])
    result <- optim(par = c(1, 1, 1), min.RSS, data = data, method = 'L-BFGS-B',
                    lower = c(0.0, -Inf, -Inf), upper = c(1, Inf, Inf))
    par <- result$par
    est.pars[n.subject, 1] <- par[1]
    est.pars[n.subject, 2] <- par[2]
    est.pars[n.subject, 3] <- par[3]
}

b <- colMeans(est.pars)[1]
c <- colMeans(est.pars)[2]
d <- colMeans(est.pars)[3]

min.RSS.fix.d <- function(data, par)
{
    with(data, sum(((a + (b - a) / (1 + exp((par[1] - x) / d))) - y)^2))
}

min.RSS.fix.c <- function(data, par)
{
    with(data, sum(((a + (b - a) / (1 + exp((c - x) / par[1]))) - y)^2))
}

init.matrix <- matrix(NA, nrow = n.subjects, ncol = 1)
est.pars <- list(init.matrix, init.matrix, init.matrix)
means <- list(matrix(NA, nrow = n.subjects, ncol = 6), matrix(NA, nrow = n.subjects, ncol = 6), matrix(NA, nrow = n.subjects, ncol = 6))

settings <- c('singles', 'pairs', 'quadruplet')
for(subject in subjects)
{
    n.subject <- which(subject == subjects)
    
    for(i in 1:3)
    {
        indices <- dt$task == settings[i] & dt$trial.type == 'experiment' & dt$subject == subject
        data <- data.frame(x = dt$target.frames[indices], y = dt$acc[indices])
        result <- optim(par = c(1), min.RSS.fix.d, data = data, method = 'L-BFGS-B',
                        lower = c(-Inf), upper = c(Inf))
        par <- result$par
        est.pars[[i]][n.subject, 1] <- par[1]
        yest <- yhat(a, b, c, par[1], c(data$x))
        means[[i]][n.subject, ] <- tapply(data$y, data$x, mean) ## get points
        
    }
    
}

figure.path <- paste(figures.dir, 'psychometric_acc_fixed_exp_1.jpeg', sep = '')
# jpeg(filename = figure.path, width = 9, height = 9, units = 'cm', res = 600, pointsize = 6)

frame <- 1000/85
labels <- frame * 1:6

par(font.axis = 2, font.lab = 2, cex = 1.2)
plot(0, ylim = c(0.4, 1), xlim = c(1, 6), type = 'n', xlab = 'Target Duration (ms)', ylab = 'Proportion Correct',
     main = '', xaxt='n')
axis(1, at=1:6, labels=round(labels, 1))

for(i in 1:3)
{
    pars <- est.pars[[i]]
    mean.pars <- colMeans(pars)
    y.est.mean <- yhat(a, b, c, mean.pars[1], 1:6)
    xy <- spline(1:6, y.est.mean, method = 'hyman')
    lines(xy, lty = i, lwd = 1)
    points(1:6, colMeans(means[[i]]), pch=i) ## plot point
}

legend('topleft', c('No. of Possible Targets', '\t\t2', '\t\t4', '\t\t8'), lty = 0:3, text.font = 2, bty = 'n', lwd = 1,
       pch=0:3, pt.cex=c(0, rep(1, 3)))

# dev.off()


# FIT PSYCHOMETRIC FUNCTION FOR ACCURACY FREE PARAMETERS  -----------------------------------------------------
par(mfrow = c(1, 1))
yhat <- function(a, b, c, d, x)
{
    y <- a + (b - a) / (1 + exp((c - x) / d))
    return(y)
}


min.RSS <- function(data, par)
{
    with(data, sum(((par[1] + (par[2] - par[1]) / (1 + exp((par[3] - x) / par[4]))) - y)^2))
}

init.matrix <- matrix(NA, nrow = n.subjects, ncol = 4)
est.pars <- list(init.matrix, init.matrix, init.matrix)
means <- list(matrix(NA, nrow = n.subjects, ncol = 6), matrix(NA, nrow = n.subjects, ncol = 6), matrix(NA, nrow = n.subjects, ncol = 6))

settings <- c('singles', 'pairs', 'quadruplet')
for(subject in subjects)
{
    n.subject <- which(subject == subjects)

    for(i in 1:3)
    {
        indices <- dt$task == settings[i] & dt$trial.type == 'experiment' & dt$subject == subject
        data <- data.frame(x = dt$target.frames[indices], y = dt$acc[indices])
        result <- optim(par = c(0.5, 1, 1, 1), min.RSS, data = data, method = 'L-BFGS-B',
                        lower = c(0.0, 0.0, -Inf, -Inf), upper = c(1, 1, Inf, Inf))
        par <- result$par
        est.pars[[i]][n.subject, 1] <- par[1]
        est.pars[[i]][n.subject, 2] <- par[2]
        est.pars[[i]][n.subject, 3] <- par[3]
        est.pars[[i]][n.subject, 4] <- par[4]
        yest <- yhat(par[1], par[2], par[3], par[4], c(data$x))
        means[[i]][n.subject, ] <- tapply(data$y, data$x, mean) ## get points
        
        
    }
    
}

figure.path <- paste(figures.dir, 'psychometric_acc_exp_1.jpeg', sep = '')
jpeg(filename = figure.path, width = 9, height = 9, units = 'cm', res = 600, pointsize = 6)

frame <- 1000/85
labels <- frame * 1:6

par(font.axis = 2, font.lab = 2, cex = 1.2)
plot(0, ylim = c(0.4, 1), xlim = c(1, 6), type = 'n', xlab = 'Target Duration (ms)', ylab = 'Proportion Correct',
     main = '', xaxt='n')
axis(1, at=1:6, labels=round(labels, 1))

for(i in 1:3)
{
    pars <- est.pars[[i]]
    mean.pars <- colMeans(pars)
    y.est.mean <- yhat(mean.pars[1], mean.pars[2], mean.pars[3], mean.pars[4], 1:6)
    xy <- spline(1:6, y.est.mean, method = 'hyman')
    lines(xy, lty = i, lwd = 1)
    points(1:6, colMeans(means[[i]]), pch=i) ## plot point
}

legend('topleft', c('No. of Possible Targets', '\t\t2', '\t\t4', '\t\t8'), lty = 0:3, text.font = 2, bty = 'n', lwd = 1,
     pch=0:3, pt.cex=c(0, rep(1, 3)))

dev.off()


# TEST PARAMETERS (ANOVA) -----------------------------------------------------------------------------------------

n.tests <- 4
n.tests <- 1
comparisons <- list()
for(test in 1:n.tests)
{
    parameters <- c(est.pars[[1]][, test], est.pars[[2]][, test], est.pars[[3]][, test])
    noa <- rep(c('2', '4', '8'), each = n.subjects)
    psychometric.dt <- data.frame('par' = parameters, 'noa' = noa, 'sub' = factor(1:29))
    model.noa <-  lmer(parameters ~ noa + (1 | sub), data = psychometric.dt)
    model.null <- lmer(parameters ~ 1 + (1 | sub), data = psychometric.dt)
    comparisons[[test]] <- anova(model.noa, model.null)
    print(summary(model.null))
}

# MODEL DEPENDENT VARIABLES ------------------------------------------------------------------------------------------

acc.model.path <- paste(script.dir, 'acc_models_experiment_1.RData', sep='')
do.models.exist <- file.exists(acc.model.path)

if(!do.models.exist)
{

    model.acc.full <- glmer(acc ~ task * pas * target.frames + (1 | subject), data = dt, family = 'binomial',
                            verbose = TRUE, subset = dt$trial.type == 'experiment') ## can be killed
    
    model.acc.minus.3.way <-  glmer(acc ~ task + pas + target.frames + task:pas + task:target.frames + pas:target.frames + 
                                        (1 | subject), data = dt, family = 'binomial',
                               verbose = TRUE, subset = dt$trial.type == 'experiment') ## can be killed
    
    model.acc.minus.3.way.minus.pt <-  glmer(acc ~ task + pas + target.frames + task:pas + task:target.frames + 
                                                 (1 | subject), data = dt, family = 'binomial',
                                             verbose = TRUE, subset = dt$trial.type == 'experiment') ## CANNOT be killed
    
    model.acc.minus.3.way.minus.tt <-  glmer(acc ~ task + pas + target.frames + task:pas + pas:target.frames +
                                                 (1 | subject), data = dt, family = 'binomial',
                                              verbose = TRUE, subset = dt$trial.type == 'experiment')
    
    model.acc.minus.3.way.minus.tt.minus.tp <-  glmer(acc ~ task + pas + target.frames + pas:target.frames +
                                                 (1 | subject), data = dt, family = 'binomial',
                                             verbose = TRUE, subset = dt$trial.type == 'experiment') ## can be killed
    
    model.acc.additive.plus.pt.minus.task <-  glmer(acc ~ pas + target.frames + pas:target.frames +
                                                          (1 | subject), data = dt, family = 'binomial',
                                                      verbose = TRUE, subset = dt$trial.type == 'experiment')
    
    acc.models <- list(model.acc.full, model.acc.minus.3.way, model.acc.minus.3.way.minus.pt,
                       model.acc.minus.3.way.minus.tt, model.acc.minus.3.way.minus.tt.minus.tp,
                       model.acc.additive.plus.pt.minus.task)
    
    save(acc.models, file = acc.model.path)

} else load(acc.model.path)

freq.model.path <- paste(script.dir, 'freq_models_experiment_1.RData', sep='')
do.models.exist <- file.exists(freq.model.path)

if(!do.models.exist)
{
    model.count.full <- glmer(count ~ task * pas * target.frames + (1 | subject), data = dt.count, family = 'poisson',
                              verbose = TRUE)
    
    model.count.minus.3.way <- glmer(count ~ task + pas + target.frames + task:pas + task:target.frames + 
                                         pas:target.frames + (1 | subject), data = dt.count, family = 'poisson',
                                     verbose = TRUE)
    
    model.count.minus.3.way.minus.pt <-  glmer(count ~ task + pas + target.frames + task:pas + task:target.frames + 
                                                 (1 | subject), data = dt.count, family = 'poisson',
                                             verbose = TRUE) ## CANNOT be killed
    
    model.count.minus.3.way.minus.tt <-  glmer(count ~ task + pas + target.frames + task:pas + pas:target.frames + 
                                                   (1 | subject), data = dt.count, family = 'poisson',
                                               verbose = TRUE)
    
    model.count.minus.3.way.minus.tt.minus.tp <-  glmer(count ~ task + pas + target.frames + pas:target.frames + 
                                                   (1 | subject), data = dt.count, family = 'poisson',
                                               verbose = TRUE) ## CANNOT be killed
    
    freq.models <- list(model.count.full, model.count.minus.3.way, model.count.minus.3.way.minus.pt,
                        model.count.minus.3.way.minus.tt, model.count.minus.3.way.minus.tt.minus.tp)
    
    save(freq.models, file = freq.model.path)
} else load(freq.model.path)

# PLOT MODELS -----------------------------------------------------------------------------------------------------

full.model.path <- paste(script.dir, 'full_models_experiment_1.RData', sep='')
do.models.exist <- file.exists(full.model.path)

if(!do.models.exist)
{

    dt$interaction <- with(dt, interaction(target.frames, pas, task))
    dt.count$interaction <- with(dt.count, interaction(target.frames, pas, task))
    
    model.acc.plot <- glmer(acc ~ interaction + 0 + (1 | subject), data = dt, family = 'binomial', verbose = TRUE, subset = dt$trial.type == 'experiment')
    model.rt.plot <- lmer(log(rt.obj) ~ interaction + 0 + (1 | subject), data = dt, verbose = TRUE, subset = dt$trial.type == 'experiment')
    model.count.plot <- glmer(count ~ interaction + 0 + (1 | subject), data = dt.count, family = 'poisson', verbose = TRUE)

    
    models <- list(model.acc.plot, model.rt.plot, model.count.plot)    
    save(models, file = full.model.path)
} else load(full.model.path)


# GET PROPORTION OF GRADED RESPONSES AMONG THE INFORMATIVE RESPONSES ---------------------------------------------

effects <- fixef(models[[3]])
inv.effects <- exp(effects)

inv.effects.singles <- inv.effects[c(3, 9, 15, 21)]
inv.effects.pairs <- inv.effects[c(3, 9, 15, 21) + 24]
inv.effects.quadruplets <- inv.effects[c(3, 9, 15, 21) + 48]

print(prop.gradeds.singles <- sum(inv.effects.singles[2:3] / sum(inv.effects.singles[2:4])))
print(prop.gradeds.pairs <- sum(inv.effects.pairs[2:3] / sum(inv.effects.pairs[2:4])))
print(prop.gradeds.quadruplets <- sum(inv.effects.quadruplets[2:3] / sum(inv.effects.quadruplets[2:4])))

# ACTUALLy PLOT THEM -----------------------------------------------------------------------------------------------


model.index <- 1
y.lims <- list(c(0.30, 1.0), c(0.0, 1.0), c(0, 40))
mains <- c('Accuracy', 'Reaction Time', 'Counts')
y.labels <- matrix(c('Proportion Correct', '', '', 'Proportion Correct', '', '',
                     'Reaction Time (s)',  '', '', 'Reaction Time (s)',  '', '',
                     'Counts',             '', '', 'Counts',             '', ''), nrow = 3, byrow = TRUE)
names <- c('acc', 'rt', 'freq')
names.args <- list(c('', '', ''), c('', '', ''), c('', '', ''), c('', '', ''),c('2', '4', '8'), c('', '', ''))
x.labels <- c('', '', '', '', 'No. of Possible Targets', '')
durations <- c(11.8, 23.5, 35.3, 47.1, 58.8, 70.6)
for(model in models)
{
    
    figure.path <- paste(figures.dir, names[model.index], '_exp_1.jpeg', sep = '')
    # jpeg(filename = figure.path, width = 19, height = 13, units = 'cm', res = 600)
    par(mfrow = c(2, 3), cex = 1, font.lab = 2, font.axis = 2, mar = c(5, 4, 1, 1))
    
    for(frame in 1:6)
    {
        effects <- fixef(model)
        effects <- effects[seq(frame, length(effects), 6)]
        if(model.index == 1) matrix.effects <- matrix(inv.logit(effects), nrow = 4, byrow = FALSE)
        if(model.index > 1) matrix.effects <- matrix(exp(effects), nrow = 4, byrow = FALSE)
        
        se <- summary(model)$coefficients[, 2]
        se <- se[seq(frame, length(se), 6)]
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
                 xpd = FALSE, yaxt = 'n', ylab = y.labels[model.index, frame], xlab = x.labels[frame],
                 names.arg = names.args[[frame]],
                 main = paste('Duration', durations[frame], 'ms'),
                 ci.lwd = 3)
        if(frame == 1 | frame == 4)
        {
            if(model.index == 1) axis(side = 2, at = seq(0.4, 1.0, 0.1))
            if(model.index == 2) axis(side = 2, at = seq(0.0, 1.0, 0.5))
            if(model.index == 3) axis(side = 2, at = seq(0, 40, 10))
        }
        
        if(model.index == 1) lines(c(h[1] - 1, h[length(h)] + 1), rep(0.50, 2), lty = 3, lwd = 3)
        ## add legends
        par(cex.axis = 1.0)
        if((model.index == 1 | model.index == 3) & frame == 1)
        {
            a <- axis(side = 1, at = seq(h[1, 2], h[4, 2], 1), labels = c('1', '2', '3', '4'))
            if(model.index == 1) text(x = h[4, 1] - 1, y = 0.14, 'PAS', xpd=NA, font = 2)
            if(model.index == 3) text(x = h[4, 1] - 1, y = -9, 'PAS', xpd=NA, font = 2)
        }
    }
    
        model.index <- model.index + 1
        
        # dev.off()
}


# DIRECT COMPARISONS FREQ ----------------------------------------------------------------------------------------------

n.frames <- 6
model <- models[[3]]
n.effects <- length(fixef(model))
summaries <- list()

comparisons <- list(c(1, 25), c(1, 49),   c(25, 49),
                    c(7, 31), c(7, 55),   c(31, 55),
                    c(13, 37), c(13, 61), c(37, 61),
                    c(19, 43), c(19, 67), c(43, 67))
n.comparisons <- length(comparisons)
ps <- c()
for(frame in 1:n.frames)
{
    K <- matrix(0, nrow = n.comparisons, ncol = n.effects)
    comparison.index <- 1
    
    for(comparison in comparisons)
    {
        K[comparison.index, comparison[1] + frame - 1]  <- -1
        K[comparison.index, comparison[2] + frame - 1] <- 1
        comparison.index <- comparison.index + 1
    }
    gh <- glht(model, linfct = K)
    summaries[[frame]] <- summary(gh, test = adjusted('none'))
    ps <- c(ps, summaries[[frame]]$test$pvalues)
}

# MIDDLE RANGE STATISTICS -----------------------------------------------------------------------------------------

mid.acc.model.path <- paste(script.dir, 'mid_acc_models_experiment_1.RData', sep='')
do.models.exist <- file.exists(mid.acc.model.path)

if(!do.models.exist)
{

    mid.range.indices <- dt$trial.type == 'experiment' & (dt$target.frames == 2 | dt$target.frames == 3 | dt$target.frames == 4)
    
    model.mid.acc.full <- glmer(acc ~ task * pas * target.frames + (1 | subject), data = dt, subset = mid.range.indices,
                                verbose = TRUE, family = 'binomial') ## can be killed
    model.mid.acc.minus.3.way <- glmer(acc ~ task + pas + target.frames + task:pas + task:target.frames + pas:target.frames +
                                           (1 | subject), data = dt, subset = mid.range.indices,
                                       verbose = TRUE, family = 'binomial')
    model.mid.acc.minus.3.way.minus.pt <- glmer(acc ~ task + pas + target.frames + task:pas + task:target.frames +
                                                    (1 | subject), data = dt, subset = mid.range.indices,
                                                verbose = TRUE, family = 'binomial') ## CANNOT be killed
    
    model.mid.acc.minus.3.way.minus.tt <- glmer(acc ~ task + pas + target.frames + task:pas + pas:target.frames +
                                                    (1 | subject), data = dt, subset = mid.range.indices,
                                                verbose = TRUE, family = 'binomial')## CHOSEN model
    
    model.mid.acc.minus.3.way.minus.tt.minus.tp <- glmer(acc ~ task + pas + target.frames + pas:target.frames +
                                                    (1 | subject), data = dt, subset = mid.range.indices,
                                                verbose = TRUE, family = 'binomial')## can be killed
    
    model.mid.acc.additive.plus.pt.minus.task <- glmer(acc ~ pas + target.frames + pas:target.frames +
                                                      (1 | subject), data = dt, subset = mid.range.indices,
                                                  verbose = TRUE, family = 'binomial')## can be killed # CHOSEN MODEL
    
    mid.acc.models <- list(model.mid.acc.full, model.mid.acc.minus.3.way, model.mid.acc.minus.3.way.minus.pt,
                           model.mid.acc.minus.3.way.minus.tt, model.mid.acc.minus.3.way.minus.tt.minus.tp,
                           model.mid.acc.additive.plus.pt.minus.task)
    
    save(mid.acc.models, file = mid.acc.model.path)
} else load(mid.acc.model.path)


mid.freq.model.path <- paste(script.dir, 'mid_freq_models_experiment_1.RData', sep='')
do.models.exist <- file.exists(mid.freq.model.path)

if(!do.models.exist)
{
    
    mid.range.indices <- dt.count$target.frames == 2 | dt.count$target.frames == 3 | dt.count$target.frames == 4
    
    model.mid.freq.full <- glmer(count ~ task * pas * target.frames + (1 | subject), data = dt.count, subset = mid.range.indices,
                                verbose = TRUE, family = 'poisson') ## can be killed
    model.mid.freq.minus.3.way <- glmer(count ~ task + pas + target.frames + task:pas + task:target.frames + pas:target.frames +
                                           (1 | subject), data = dt.count, subset = mid.range.indices,
                                       verbose = TRUE, family = 'poisson')
    model.mid.freq.minus.3.way.minus.pt <- glmer(count ~ task + pas + target.frames + task:pas + task:target.frames +
                                                    (1 | subject), data = dt.count, subset = mid.range.indices,
                                                verbose = TRUE, family = 'poisson') ## CANNOT be killed
    
    model.mid.freq.minus.3.way.minus.tt <- glmer(count ~ task + pas + target.frames + task:pas + pas:target.frames +
                                                    (1 | subject), data = dt.count, subset = mid.range.indices,
                                                verbose = TRUE, family = 'poisson')## can be killed
    
    model.mid.freq.minus.3.way.minus.tt.minus.tp <- glmer(count ~ task + pas + target.frames + pas:target.frames +
                                                             (1 | subject), data = dt.count, subset = mid.range.indices,
                                                         verbose = TRUE, family = 'poisson')## CANNOT be killed
    
    model.mid.freq.additive.plus.pt.minus.task <- glmer(count ~ pas + target.frames + pas:target.frames +
                                                           (1 | subject), data = dt.count, subset = mid.range.indices,
                                                       verbose = TRUE, family = 'poisson')## can be killed # CHOSEN MODEL
    
    mid.freq.models <- list(model.mid.freq.full, model.mid.freq.minus.3.way, model.mid.freq.minus.3.way.minus.pt,
                           model.mid.freq.minus.3.way.minus.tt, model.mid.freq.minus.3.way.minus.tt.minus.tp,
                           model.mid.freq.additive.plus.pt.minus.task)
    
    save(mid.freq.models, file = mid.freq.model.path)
} else load(mid.freq.model.path)



# PLOT ACCURACY INCREASE DEPENDENT ON OBJECTIVE EVIDENCE (NO. OF FRAMES) ---------------------------------------

graphics.off()
par(font.axis=2, font.lab=2)
plot(0, 0, xlim=c(1, 6), ylim=c(0.35, 1), type='n', xlab='No. of Frames', ylab='Proportion Correct')


effects <- fixef(model.dep.on.obj.evidence)
effects <- inv.logit(effects)
          
se <- summary(model.dep.on.obj.evidence)$coefficients[, 2]
conf.up <- inv.logit(effects + qnorm(0.975) * se)
conf.down <- inv.logit(effects - qnorm(0.975) * se)

indices <- c(1, 7, 13, 19)
colours <- c('red', 'blue', 'purple', 'black')       
counter <- 1

for(index in indices)
{
    lines(1:6, effects[index:(index+5)], col=colours[counter], lwd=3, type='b')
    counter <- counter + 1
}

lines(1:6, rep(1.00, 6), lty=3, lwd=3)
lines(1:6, rep(0.50, 6), lty=3, lwd=3)

# # FIT PSYCHOMETRIC FUNCTION FOR ACCURACY PER PAS -----------------------------------------------------
par(mfrow = c(1, 1))
yhat <- function(a, b, c, d, x)
{
    y <- a + (b - a) / (1 + exp((c - x) / d))
    return(y)
}

a <- 0
b <- 1

min.RSS <- function(data, par)
{
    with(data, sum(((par[1] + (par[2] - par[1]) / (1 + exp((par[3] - x) / par[4]))) - y)^2))
}

init.matrix <- matrix(NA, nrow = n.subjects, ncol = 4)
est.pars <- list(init.matrix, init.matrix, init.matrix, init.matrix)
means <- list(matrix(NA, nrow = n.subjects, ncol = 6), matrix(NA, nrow = n.subjects, ncol = 6), 
              matrix(NA, nrow = n.subjects, ncol = 6), matrix(NA, nrow = n.subjects, ncol = 6))

settings <- 1:4
for(subject in subjects)
{
    n.subject <- which(subject == subjects)

    for(i in 1:4)
    {
        indices <- dt$pas == settings[i] & dt$trial.type == 'experiment' & dt$subject == subject
        data <- data.frame(x = dt$target.frames[indices], y = dt$acc[indices])
        result <- optim(par = c(0.5, 1, 1, 1), min.RSS, data = data, method = 'L-BFGS-B',
                        lower = c(0.0, 0.0, -Inf, -Inf), upper = c(1, 1, Inf, Inf))
        par <- result$par
        est.pars[[i]][n.subject, 1] <- par[1]
        est.pars[[i]][n.subject, 2] <- par[2]
        est.pars[[i]][n.subject, 3] <- par[3]
        est.pars[[i]][n.subject, 4] <- par[4]
        yest <- yhat(par[1], par[2], par[3], par[4], c(data$x))
        mean.vector <- rep(NA, 6)
        temp <- tapply(data$y, data$x, mean) ## get points
        for(name.index in 1:6)
        {
            if(any(names(temp) == name.index)) mean.vector[name.index] <- temp[names(temp) == name.index]
        }
        
        means[[i]][n.subject, ] <- mean.vector 

        
    }
    
}

figure.path <- paste(figures.dir, 'psychometric_acc_exp_1_follow_up.jpeg', sep = '')
jpeg(filename = figure.path, width = 9, height = 9, units = 'cm', res = 600, pointsize = 6)

par(font.axis = 2, font.lab = 2, cex = 1.2)
plot(0, ylim = c(0.0, 1.1), xlim = c(1, 6), type = 'n', xlab = 'Target Duration (ms)', ylab = 'Proportion Correct',
     main = '', xaxt='n')
axis(1, at=1:6, labels=round(labels, 1))
colours <- rainbow(4)

for(i in 1:4)
{
    pars <- est.pars[[i]]
    mean.pars <- colMeans(pars)
    y.est.mean <- yhat(mean.pars[1], mean.pars[2], mean.pars[3], mean.pars[4], 1:6)
    xy <- spline(1:6, y.est.mean, method = 'hyman')
    lines(xy, lty = 1, lwd = 3, col=colours[i])
    points(1:6, colMeans(means[[i]], na.rm = TRUE), pch=i, col=colours[i], cex=2) ## plot point

}

legend('topleft', c('Perceptual Awaraness Scale', '\t1', '\t2', '\t3', '\t4'),
         lty = 1, text.font = 2, bty = 'n', lwd = 3, col=c('white', colours), pch=0:3, pt.cex=c(0, rep(1, 4)))

dev.off()

# TEST PARAMETERS (ANOVA) -----------------------------------------------------------------------------------------

n.tests <- 4
n.tests <- 1
comparisons <- list()
for(test in 1:n.tests)
{
    parameters <- c(est.pars[[1]][, test], est.pars[[2]][, test], est.pars[[3]][, test], est.pars[[4]][, test])
    pas <- rep(c('1', '2', '3', '4'), each = n.subjects)
    psychometric.dt <- data.frame('par' = parameters, 'pas' = pas, 'sub' = factor(1:29))
    model.noa <-  lmer(parameters ~ pas + (1 | sub) + 0, data = psychometric.dt, subset=psychometric.dt$pas != '1')
    model.null <- lmer(parameters ~ 1 + (1 | sub) + 0 , data = psychometric.dt, subset=psychometric.dt$pas != '1')
    comparisons[[test]] <- anova(model.noa, model.null)
    # print(summary(model.null))
}

K <- matrix(c(-1,  1, 0,
              -1,  0, 1,
               0, -1, 1), nrow=3, byrow=TRUE)

gh <- glht(model.noa, linfct = K)
print(summary(gh, test = adjusted('none')))
