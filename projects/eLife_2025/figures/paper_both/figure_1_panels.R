## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(patchwork)


## -----------------------------------------------------------------------------
Nch= 20

fs = 50e3
nsamples= 100
noise=1e-4 *fs/nsamples


## -----------------------------------------------------------------------------
# The two vectors are zipped BY POSITION (see the read loop below): row i of `filenames`
# is labelled with row i of `algorithms`. A wrong filename errors at read.csv; a wrong
# ORDER silently mislabels a whole column. Edit them together, always.
filenames=c(
  "../data/figure_1_likelihood_diagnostic_LSE",
  "../data/figure_1_likelihood_diagnostic_NR",
  "../data/figure_1_likelihood_diagnostic_MNR",
  "../data/figure_1_likelihood_diagnostic_R",
  "../data/figure_1_likelihood_diagnostic_MR",
  "../data/figure_1_likelihood_diagnostic_VR",
  "../data/figure_1_likelihood_diagnostic_IR",
  "../data/figure_1_simulation"

)

algorithms=c("LSE"
             ,"NR"
             ,"NMR"
             ,"R"
             ,"MR"
             ,"VR"
             ,"IR"
             ,"SIM")



## -----------------------------------------------------------------------------

d=data.frame()
for (i in 1:length(filenames)){
  f <- paste0(filenames[i], ".csv")
  if (!file.exists(f)) { message("figure_1_panels: skipping missing dump ", f); next }  # LSE dump may not exist yet
  d0=read.csv(f, check.names = FALSE)
  d0$algo=algorithms[i]
  d=rbind(d,d0)
}



## -----------------------------------------------------------------------------
d_s=d%>%filter(algo=="SIM", scope=="simulation_sub")%>%select(sub_index,sample_index,step_start,step_middle,step_end,component_path,value)%>%pivot_wider(names_from = component_path, values_from = value)

d_s_mean<-d_s%>%group_by(sample_index)%>%summarize(step_start=min(step_start), step_end=max(step_end), patch_current= mean(patch_current))

d_IR=d%>%filter(algo=="IR", scope=="evolution", value_row==1| is.na(value_row),  is.na(value_col), !is.na(step_start))%>%mutate(component= component_path %>%
           str_replace("^Algo_State_Dynamic\\.", "") %>%
           str_replace("^Patch_State\\.", ""))%>%
  select(sample_index,step_start,step_middle,step_end,component,patch_current,value)%>%
  pivot_wider(names_from = component, values_from = value)


## -----------------------------------------------------------------------------
d_NR=d%>%filter(algo=="NR", scope=="evolution", value_row==1| is.na(value_row),  is.na(value_col), !is.na(step_start))%>%mutate(component= component_path %>%
           str_replace("^Algo_State_Dynamic\\.", "") %>%
           str_replace("^Patch_State\\.", ""))%>%
  select(sample_index,step_start,step_middle,step_end,component,patch_current,value)%>%
  pivot_wider(names_from = component, values_from = value)

# LSE (least squares on the mean): non-recursive like NR. Its predictive variance is a single CONSTANT
# sigma_hat^2 = SSE/n (the marginalised noise scale that makes the logL well-defined), so its row-B
# band is FLAT. Dump = figure_1.macroir's LSE block (family = 2). Columns used are VERIFIED present in
# figure_1_likelihood_diagnostic_LSE.csv: P_mean_t20_y1, y_mean, y_var, patch_current, logL.
# Empty until that dump exists (tolerant read).
d_LSE=d%>%filter(algo=="LSE", scope=="evolution", value_row==1| is.na(value_row),  is.na(value_col), !is.na(step_start))%>%mutate(component= component_path %>%
           str_replace("^Algo_State_Dynamic\\.", "") %>%
           str_replace("^Patch_State\\.", ""))%>%
  select(sample_index,step_start,step_middle,step_end,component,patch_current,value)%>%
  pivot_wider(names_from = component, values_from = value)

d_NMR=d%>%filter(algo=="NMR", scope=="evolution", value_row==1| is.na(value_row),  is.na(value_col), !is.na(step_start))%>%mutate(component= component_path %>%
           str_replace("^Algo_State_Dynamic\\.", "") %>%
           str_replace("^Patch_State\\.", ""))%>%
  select(sample_index,step_start,step_middle,step_end,component,patch_current,value)%>%
  pivot_wider(names_from = component, values_from = value)

d_R=d%>%filter(algo=="R", scope=="evolution", value_row==1| is.na(value_row),  is.na(value_col), !is.na(step_start))%>%mutate(component= component_path %>%
           str_replace("^Algo_State_Dynamic\\.", "") %>%
           str_replace("^Patch_State\\.", ""))%>%
  select(sample_index,step_start,step_middle,step_end,component,patch_current,value)%>%
  pivot_wider(names_from = component, values_from = value)

d_MR=d%>%filter(algo=="MR", scope=="evolution", value_row==1| is.na(value_row),  is.na(value_col), !is.na(step_start))%>%mutate(component= component_path %>%
           str_replace("^Algo_State_Dynamic\\.", "") %>%
           str_replace("^Patch_State\\.", ""))%>%
  select(sample_index,step_start,step_middle,step_end,component,patch_current,value)%>%
  pivot_wider(names_from = component, values_from = value)

# VR is built from MR's formulas with the RESIDUAL predictive variance instead of the
# total one, so its panels below are MR's CODE on d_VR — but not MR's numbers. y_var
# divides the innovation (chi = dy / y_var) and the covariance down-date ((N/y_var)*XTX(gS)),
# so the whole recursion diverges. Measured on this recording: y_var 26% below MR at the
# median; sample 0 identical (no agonist) and sample 1 identical in mean and covariance
# only because that prior is deterministic (gS = 0, so the update is null whatever the
# variance); posterior apart from sample 2, propagated mean from sample 3.
d_VR=d%>%filter(algo=="VR", scope=="evolution", value_row==1| is.na(value_row),  is.na(value_col), !is.na(step_start))%>%mutate(component= component_path %>%
           str_replace("^Algo_State_Dynamic\\.", "") %>%
           str_replace("^Patch_State\\.", ""))%>%
  select(sample_index,step_start,step_middle,step_end,component,patch_current,value)%>%
  pivot_wider(names_from = component, values_from = value)






## -----------------------------------------------------------------------------
common_theme = theme_classic()+ theme( legend.position = "bottom")

# Okabe-Ito colour-blind-safe semantic palette, shared across panels so patchwork guides="collect"
# merges everything into ONE legend. Colour = FILTER PHASE: predict (orange) groups the Markov step
# with the prior it produces; update (purple) groups the Bayes step with the posterior it produces
# (this operation->state link holds for MacroR/MR/VR/IR alike). SHAPE separates operation (arrow) from
# state (point/segment); ALPHA is the carried-over channel: an arrow's ORIGIN (the state carried in
# from the previous step) is drawn faded (ORIG_ALPHA) while the state the panel PRODUCES stays at full
# alpha. One mechanism for points and segments (the prior produced by the Markov panel reappears faded
# as the origin of the Bayes panel's arrow, and the posterior the other way). Hue and alpha never collide.
SEM <- c("channel current"          = "#7FCDBB",   # light green — the true resolved current (background reality)
         "Markov & prior"            = "#E69F00",   # orange      — PREDICT: Markov step, prior, AND predicted current
         "observation & innovation"  = "#0072B2",   # blue        — OBSERVE: observed current AND the innovation (obs - pred)
         "Bayes & posterior"         = "#CC79A7",   # purple      — UPDATE : Bayes step AND posterior
         "logLikelihood"             = "#000000")   # black       — the step OUTPUT/score; achromatic + colorblind-safe. Plot as a BOLD line so it does not blend with the (black) axes
sem_scale <- ggplot2::scale_colour_manual(values = SEM, limits = names(SEM), name = NULL, drop = FALSE,
                                          breaks = setdiff(names(SEM), "logLikelihood"),   # logL row is self-labelled -> drop from legend
                                          labels = function(x) ifelse(x == "Markov & prior", "Markov, prior & prediction", x))  # orange spans prior (state) + prediction (measurement)

# common x (time) range over the FULL recording (anchored at 0); each figure applies its OWN crop
# later via build_figure(), so the arrows are always the full-graph arrows, just zoomed.
XLIM <- range(c(0, d_s$step_middle, d_IR$step_start, d_IR$step_end), na.rm = TRUE)
xcommon <- coord_cartesian(xlim = XLIM)

# --- common y-scales, recomputed FROM THE DATA so they adapt if the data changes (no hardcoding) ---
.dP <- function(x) unlist(x[, grepl("^P_mean", names(x)), drop = FALSE])          # every P(open) column
.dI <- function(x) c(-x$patch_current+1, -x$y_mean+1, -x$y_mean+1 - sqrt(x$y_var), -x$y_mean+1 + sqrt(x$y_var))  # current (with +1 baseline compensation)
.algos <- list(d_NR, d_NMR, d_R, d_MR, d_VR, d_IR)
YP  <- range(unlist(lapply(.algos, .dP)), na.rm = TRUE)                           # prior & posterior share this
YI  <- range(c(-d_s$patch_current, unlist(lapply(.algos, .dI))), na.rm = TRUE)    # observation row

compact_theme <- common_theme +theme(
  axis.title.x = element_blank(),
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank(),
  legend.position = "right",
  plot.margin = margin(1, 5, 1, 5)
)
x_only_ticks=theme(axis.title.x = element_blank(), axis.text.x = element_blank())
y_only_ticks=theme(axis.title.y = element_blank(), axis.text.y = element_blank())
no_y_axis= theme( axis.text.y = element_blank(), # Removes y-axis text
         axis.title.y = element_blank(), # Removes y-axis title
          axis.ticks.y = element_blank(),
        axis.line.y = element_blank())
no_x_axis= theme( axis.text.x = element_blank(), # Removes y-axis text
         axis.title.x = element_blank(), # Removes y-axis title
          axis.ticks.x = element_blank(),
        axis.line.x = element_blank())
y_no_title=theme(axis.title.y = element_blank())


## -----------------------------------------------------------------------------
markov_arrow = arrow(length = unit(0.07, "in"), type="closed")
bayes_arrow= arrow(length = unit(0.07, "in"), type="closed")
pred_arrow = arrow(length = unit(0.07, "in"), type="closed")

# de-emphasis for an arrow's ORIGIN: the state carried in from the previous step is drawn faded, so it
# recedes to context while the full-colour mark is the state the panel produces. One value for points
# and segments alike; kept high enough that the hue still reads (not grey) after PDF downsample.
ORIG_ALPHA = 0.45




## -----------------------------------------------------------------------------
fA<-ggplot(d_IR)+geom_line(data=d_s, aes(x=step_middle, y= -patch_current, color="channel current"),
                linewidth = 0.5) +
  geom_segment(aes(x=step_start,xend=step_end,y=-patch_current+1,color="observation & innovation"), linewidth = 1)+
  geom_line(data=data.frame(role=names(SEM), x=NA_real_, y=NA_real_), aes(x=x, y=y, colour=role), na.rm=TRUE)+
  sem_scale + ylab("current (pA)") + common_theme+x_only_ticks + xcommon
fA


## ----fB_IR--------------------------------------------------------------------
fB_IR<-ggplot(d_IR)+
  geom_segment(aes(x=lag(step_start), xend = lag(step_end), y=lag(P_mean_t10_y1), yend = lag(P_mean_t20_y1), color="Bayes & posterior"), linewidth = 1, alpha=ORIG_ALPHA)+
  geom_segment(aes(x=lead(step_start), xend = lead(step_end), y=P_mean_t20_y1, yend =lead( P_mean_t11_y0), color="Markov & prior"), linewidth = 1, alpha=1)+
geom_curve(aes(x=lag(step_middle), xend = step_middle, y=lag(P_mean_t10_y1)/2+ lag(P_mean_t20_y1)/2,
   yend = lag(P_mean_t20_y1)/2+ P_mean_t11_y0/2, color="Markov & prior"),curvature = 0.6,  angle = 90, ncp=10,linetype = 1, linewidth = 0.5, arrow=markov_arrow)+
  sem_scale + ylab("P(open)") + common_theme+x_only_ticks + guides(colour="none") + xcommon

fB_IR


## -----------------------------------------------------------------------------
fC_IR<-ggplot(d_IR)+
  geom_rect(aes(xmin=step_start, xmax=step_end, ymin=-y_mean+1-sqrt(y_var), ymax=-y_mean+1+sqrt(y_var), fill="Markov & prior"), alpha=0.3)+
  geom_segment(aes(x=step_start, xend= step_end, y=-patch_current+1, color="observation & innovation"), linewidth=1, alpha=1)+
  geom_segment(aes(x=step_start, xend= step_end, y=-y_mean+1, color="Markov & prior"), linewidth=1)+
  geom_segment(aes(x=(step_start+step_end)/2, xend=(step_start+step_end)/2, y=-y_mean+1, yend=-patch_current+1, color="observation & innovation"),
               arrow=pred_arrow, linewidth=0.5, na.rm=TRUE)+
  scale_fill_manual(values = SEM, guide = "none") +
  sem_scale + ylab("predicted current (pA)") + common_theme+x_only_ticks + guides(colour="none") + xcommon
fC_IR



## -----------------------------------------------------------------------------

# Pre-calculate to avoid NA/length errors in the plot call
plot_data_clean_IR <- d_IR %>%
  mutate(
    # Handle the lag NA by replacing or filtering if necessary
    y_start = lag(P_mean_t20_y1),
    y_end = P_mean_t11_y0,
    y_mid_start = (y_start / 2) + (y_end / 2),
    y_mid_end = (P_mean_t10_y1 / 2) + (P_mean_t20_y1 / 2)
  ) %>%
  filter(!is.na(y_start)) # Remove the NA introduced by lag
fD_IR <- ggplot(plot_data_clean_IR) +
  geom_segment(aes(x = step_start, xend = step_end, y = y_start, yend = y_end, color = "Markov & prior"),
               linewidth = 1, alpha = ORIG_ALPHA) +
  geom_segment(aes(x = step_start, xend = step_end, y = P_mean_t10_y1, yend = P_mean_t20_y1, color = "Bayes & posterior"),
               linewidth = 1) +
  geom_curve(aes(x = step_middle, y = y_mid_start,
                 xend = step_middle, yend = y_mid_end, color = "Bayes & posterior"),
             linewidth = 0.5, curvature = -0.9,  angle = 180, ncp=10,
             arrow = bayes_arrow) +
  sem_scale + common_theme +
  xlab("time (ms)") + ylab("P(open)") + guides(colour="none") + xcommon+x_only_ticks
fD_IR





## -----------------------------------------------------------------------------
fE_IR<-ggplot(d_IR)+
 # geom_segment(aes(x=step_start, xend=step_end, y=logL))+
  geom_line(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  geom_point(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  

  sem_scale + ylab("logL") + common_theme+guides(colour="none") + xcommon+xlab("time (ms)")
fE_IR



## -----------------------------------------------------------------------------
fB_NR<-ggplot(d_NR)+
  geom_point(aes(x=step_middle, y=P_mean_t2_y0, color="Markov & prior"), linewidth = 1, alpha=1)+
geom_curve(aes(x=lag(step_middle), xend = step_middle, y=lag(P_mean_t2_y0),
   yend = P_mean_t2_y0, color="Markov & prior"),curvature = 0.6,  angle = 90, ncp=10,linetype = 1, linewidth = 0.5, arrow=markov_arrow)+
  sem_scale + ylab("P(open)") + common_theme+x_only_ticks + guides(colour="none") + xcommon

fB_NR


## -----------------------------------------------------------------------------
fC_NR<-ggplot(d_NR)+
  geom_errorbar(aes(x=step_middle, ymin=-y_mean+1-sqrt(y_var), ymax=-y_mean+1+sqrt(y_var),width=0.0005, color="Markov & prior"), alpha=1)+
  geom_segment(aes(x=step_start, xend= step_end, y=-patch_current+1, color="observation & innovation"), linewidth=1, alpha=1)+
  geom_point(aes(x=step_middle,  y=-y_mean+1, color="Markov & prior"), linewidth=1)+
  geom_segment(aes(x=(step_start+step_end)/2, xend=(step_start+step_end)/2, y=-y_mean+1, yend=-patch_current+1, color="observation & innovation"),
               arrow=pred_arrow, linewidth=0.5, na.rm=TRUE)+
  scale_fill_manual(values = SEM, guide = "none") +
  sem_scale + ylab("predicted current (pA)") + common_theme+x_only_ticks + guides(colour="none") + xcommon
fC_NR


## ----LSE panels: like NR (non-recursive), with a CONSTANT least-squares variance (Luciano 2026-07-23)
## LSE is least squares on the mean, noise scale marginalised to a single plug-in sigma_hat^2 = SSE/n
## (mean(r_std^2)=1 confirms y_var IS that). Drawn in NR's grammar — POINT mean + ERRORBAR — NOT
## macroir's band + interval-spanning segment. Its errorbars are all the SAME width (the flat
## sigma_hat^2), against NR/R's width that varies with the per-interval gating variance.
## fB_LSE (row A, P(open)) uses P_mean_t20_y1 — VERIFIED present in the LSE dump (its ONLY P_mean key;
## LSE is built averaging=1 so it follows IR's naming, NOT NR's P_mean_t2_y0). y_mean, y_var,
## patch_current, logL also verified present in figure_1_likelihood_diagnostic_LSE.csv.
fB_LSE<-ggplot(d_LSE)+
  geom_point(aes(x=step_middle, y=P_mean_t20_y1, color="Markov & prior"), linewidth = 1, alpha=1)+
geom_curve(aes(x=lag(step_middle), xend = step_middle, y=lag(P_mean_t20_y1),
   yend = P_mean_t20_y1, color="Markov & prior"),curvature = 0.6,  angle = 90, ncp=10,linetype = 1, linewidth = 0.5, arrow=markov_arrow)+
  sem_scale + ylab("P(open)") + common_theme+x_only_ticks + guides(colour="none") + xcommon

# fC_LSE = NR's row-B ERRORBAR style on the CONSTANT sigma_hat^2, so every errorbar is the SAME width
# (the flat least-squares noise scale). POINT mean at step_middle (NR-style), NOT an interval-spanning
# segment; innovation arrow to the observation. Same code as fC_NR — the only difference from NR is
# that y_var here is one constant, so the widths do not breathe.
fC_LSE<-ggplot(d_LSE)+
  geom_errorbar(aes(x=step_middle, ymin=-y_mean+1-sqrt(y_var), ymax=-y_mean+1+sqrt(y_var),width=0.0005, color="Markov & prior"), alpha=1)+
  geom_segment(aes(x=step_start, xend= step_end, y=-patch_current+1, color="observation & innovation"), linewidth=1, alpha=1)+
  geom_point(aes(x=step_middle,  y=-y_mean+1, color="Markov & prior"), linewidth=1)+
  geom_segment(aes(x=(step_start+step_end)/2, xend=(step_start+step_end)/2, y=-y_mean+1, yend=-patch_current+1, color="observation & innovation"),
               arrow=pred_arrow, linewidth=0.5, na.rm=TRUE)+
  scale_fill_manual(values = SEM, guide = "none") +
  sem_scale + ylab("predicted current (pA)") + common_theme+x_only_ticks + guides(colour="none") + xcommon



## -----------------------------------------------------------------------------
fE_NR<-ggplot(d_NR)+
 # geom_segment(aes(x=step_start, xend=step_end, y=logL))+
  geom_line(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  geom_point(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  

  sem_scale + ylab("logL") + common_theme+guides(colour="none") + xcommon+xlab("time (ms)")
fE_NR



## -----------------------------------------------------------------------------
fB_NMR<-ggplot(d_NMR)+
  geom_point(aes(x=step_start, y=d_NMR$P_mean_t2_y0, color="Markov & prior"), linewidth = 1, alpha=1)+
geom_curve(aes(x=lag(step_start), xend = step_start, y=lag(P_mean_t2_y0),
   yend = P_mean_t2_y0, color="Markov & prior"),curvature = 0.6,  angle = 90, ncp=10,linetype = 1, linewidth = 0.5, arrow=markov_arrow)+
  sem_scale + ylab("P(open)") + common_theme+x_only_ticks + guides(colour="none") + xcommon

fB_NMR


## -----------------------------------------------------------------------------
fC_NMR<-ggplot(d_NMR)+
  geom_rect(aes(xmin=step_start, xmax=step_end, ymin=-y_mean+1-sqrt(y_var), ymax=-y_mean+1+sqrt(y_var), fill="Markov & prior"), alpha=0.3)+
  geom_segment(aes(x=step_start, xend= step_end, y=-patch_current+1, color="observation & innovation"), linewidth=1, alpha=1)+
  geom_segment(aes(x=step_start, xend= step_end, y=-y_mean+1, color="Markov & prior"), linewidth=1)+
  geom_segment(aes(x=(step_start+step_end)/2, xend=(step_start+step_end)/2, y=-y_mean+1, yend=-patch_current+1, color="observation & innovation"),
               arrow=pred_arrow, linewidth=0.5, na.rm=TRUE)+
  scale_fill_manual(values = SEM, guide = "none") +
  sem_scale + ylab("predicted current (pA)") + common_theme+x_only_ticks + guides(colour="none") + xcommon
fC_NMR



## -----------------------------------------------------------------------------
fE_NMR<-ggplot(d_NMR)+
 # geom_segment(aes(x=step_start, xend=step_end, y=logL))+
  geom_line(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  geom_point(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  

  sem_scale + ylab("logL") + common_theme+guides(colour="none") + xcommon+xlab("time (ms)")
fE_NR



## ----fB_R---------------------------------------------------------------------
fB_R<-ggplot(d_R)+
  geom_point(aes(x=lag(step_middle), y=lag(P_mean_t15_y1), color="Bayes & posterior"), alpha=ORIG_ALPHA)+
  geom_point(aes(x=lead(step_middle), y=lead(P_mean_t15_y0), color="Markov & prior"),  alpha=1)+
geom_curve(aes(x=lag(step_middle), xend = step_middle, y=lag(P_mean_t15_y1),
   yend = P_mean_t15_y0, color="Markov & prior"),curvature = 0.6,  angle = 90, ncp=10,linetype = 1, linewidth = 0.5, arrow=markov_arrow)+
  sem_scale + ylab("P(open)") + common_theme+x_only_ticks + guides(colour="none") + xcommon

fB_R


## -----------------------------------------------------------------------------
fC_R<-ggplot(d_R)+
  geom_errorbar(aes(x=step_middle, ymin=-y_mean+1-sqrt(y_var), ymax=-y_mean+1+sqrt(y_var), width=0.0005,color="Markov & prior"), alpha=1)+
  geom_segment(aes(x=step_start, xend= step_end, y=-patch_current+1, color="observation & innovation"), linewidth=1, alpha=1)+
  geom_point(aes(x=step_middle,  y=-y_mean+1, color="Markov & prior"), linewidth=1)+
  geom_segment(aes(x=(step_start+step_end)/2, xend=(step_start+step_end)/2, y=-y_mean+1, yend=-patch_current+1, color="observation & innovation"),
               arrow=pred_arrow, linewidth=0.5, na.rm=TRUE)+
  scale_fill_manual(values = SEM, guide = "none") +
  sem_scale + ylab("predicted current (pA)") + common_theme+x_only_ticks + guides(colour="none") + xcommon
fC_R



## -----------------------------------------------------------------------------

# Pre-calculate to avoid NA/length errors in the plot call
fD_R <- d_R %>% filter(step_start > 0) %>% ggplot() +
  geom_point(aes(x = step_middle, y = P_mean_t15_y0, color = "Markov & prior"), alpha = ORIG_ALPHA) +
  geom_point(aes(x = step_middle, y = P_mean_t15_y1, color = "Bayes & posterior")) +
  geom_curve(aes(x = step_middle, xend = step_middle, y = P_mean_t15_y0, yend = P_mean_t15_y1,
                 color = "Bayes & posterior"),
             linewidth = 0.5, curvature = -0.9, angle = 180, ncp = 10, arrow = bayes_arrow) +
  sem_scale + common_theme +
  ylab("P(open)") + guides(colour = "none") + xcommon + x_only_ticks
fD_R





## -----------------------------------------------------------------------------
fE_R<-ggplot(d_R)+
 # geom_segment(aes(x=step_start, xend=step_end, y=logL))+
  geom_line(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  geom_point(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  

  sem_scale + ylab("logL") + common_theme+guides(colour="none") + xcommon+xlab("time (ms)")
fE_NR



## ----fB_MR--------------------------------------------------------------------
fB_MR<-d_MR%>%filter(sample_index>0)%>%ggplot()+
  geom_point(aes(x=step_start, y=P_mean_t1_y1, color="Bayes & posterior"), alpha=ORIG_ALPHA)+
  geom_point(aes(x=step_end, y=P_mean_t2_y1, color="Markov & prior"),  alpha=1)+
geom_curve(aes(x=step_start, xend = step_end, y=P_mean_t1_y1,
   yend = P_mean_t2_y1, color="Markov & prior"),curvature = 0.6,  angle = 90, ncp=10,linetype = 1, linewidth = 0.5, arrow=markov_arrow)+
  sem_scale + ylab("P(open)") + common_theme + x_only_ticks + guides(colour="none") + xcommon

fB_MR


## -----------------------------------------------------------------------------
fC_MR<-ggplot(d_MR)+
  geom_rect(aes(xmin=step_start, xmax=step_end, ymin=-y_mean+1-sqrt(y_var), ymax=-y_mean+1+sqrt(y_var), fill="Markov & prior"), alpha=0.3)+
  geom_segment(aes(x=step_start, xend= step_end, y=-patch_current+1, color="observation & innovation"), linewidth=1, alpha=1)+
  geom_segment(aes(x=step_start, xend= step_end, y=-y_mean+1, color="Markov & prior"), linewidth=1)+
  geom_segment(aes(x=(step_start+step_end)/2, xend=(step_start+step_end)/2, y=-y_mean+1, yend=-patch_current+1, color="observation & innovation"),
               arrow=pred_arrow, linewidth=0.5, na.rm=TRUE)+
  scale_fill_manual(values = SEM, guide = "none") +
  sem_scale + ylab("predicted current (pA)") + common_theme+x_only_ticks + guides(colour="none") + xcommon
fC_MR



## -----------------------------------------------------------------------------

# Pre-calculate to avoid NA/length errors in the plot call
fD_MR <- d_MR %>% filter(step_start > 0) %>% ggplot() +
  geom_point(aes(x=step_start, y=P_mean_t1_y1, color="Bayes & posterior"), alpha=1)+
  geom_point(aes(x=lead(step_end), y=lead(P_mean_t2_y1), color="Markov & prior"),  alpha=ORIG_ALPHA)+
geom_curve(aes(x=step_end, xend = step_end, y=P_mean_t2_y1,
   yend = lead(P_mean_t1_y1), 
                 color = "Bayes & posterior"),
             linewidth = 0.5, curvature = -0.9, angle = 180, ncp = 10, arrow = bayes_arrow) +
  sem_scale + common_theme +
  ylab("P(open)") + guides(colour = "none") + xcommon + x_only_ticks
fD_MR





## ----fB_VR--------------------------------------------------------------------
# VR's three panels are MR's code, verbatim, on d_VR. The visible difference to look for
# is the predictive band in fC_VR (narrower: it is drawn from the residual y_var), and
# then, downstream of it, a posterior and a mean that drift away from MR's because the
# variance divides the gain.
fB_VR<-d_VR%>%filter(sample_index>0)%>%ggplot()+
  geom_point(aes(x=step_start, y=P_mean_t1_y1, color="Bayes & posterior"), alpha=ORIG_ALPHA)+
  geom_point(aes(x=step_end, y=P_mean_t2_y1, color="Markov & prior"),  alpha=1)+
geom_curve(aes(x=step_start, xend = step_end, y=P_mean_t1_y1,
   yend = P_mean_t2_y1, color="Markov & prior"),curvature = 0.6,  angle = 90, ncp=10,linetype = 1, linewidth = 0.5, arrow=markov_arrow)+
  sem_scale + ylab("P(open)") + common_theme + x_only_ticks + guides(colour="none") + xcommon

fB_VR


## ----fC_VR--------------------------------------------------------------------
fC_VR<-ggplot(d_VR)+
  geom_rect(aes(xmin=step_start, xmax=step_end, ymin=-y_mean+1-sqrt(y_var), ymax=-y_mean+1+sqrt(y_var), fill="Markov & prior"), alpha=0.3)+
  geom_segment(aes(x=step_start, xend= step_end, y=-patch_current+1, color="observation & innovation"), linewidth=1, alpha=1)+
  geom_segment(aes(x=step_start, xend= step_end, y=-y_mean+1, color="Markov & prior"), linewidth=1)+
  geom_segment(aes(x=(step_start+step_end)/2, xend=(step_start+step_end)/2, y=-y_mean+1, yend=-patch_current+1, color="observation & innovation"),
               arrow=pred_arrow, linewidth=0.5, na.rm=TRUE)+
  scale_fill_manual(values = SEM, guide = "none") +
  sem_scale + ylab("predicted current (pA)") + common_theme+x_only_ticks + guides(colour="none") + xcommon
fC_VR


## ----fD_VR--------------------------------------------------------------------
fD_VR <- d_VR %>% filter(step_start > 0) %>% ggplot() +
  geom_point(aes(x=step_start, y=P_mean_t1_y1, color="Bayes & posterior"), alpha=1)+
  geom_point(aes(x=lead(step_end), y=lead(P_mean_t2_y1), color="Markov & prior"),  alpha=ORIG_ALPHA)+
geom_curve(aes(x=step_end, xend = step_end, y=P_mean_t2_y1,
   yend = lead(P_mean_t1_y1),
                 color = "Bayes & posterior"),
             linewidth = 0.5, curvature = -0.9, angle = 180, ncp = 10, arrow = bayes_arrow) +
  sem_scale + common_theme +
  ylab("P(open)") + guides(colour = "none") + xcommon + x_only_ticks
fD_VR


## -----------------------------------------------------------------------------
fE_MR<-ggplot(d_MR)+
 # geom_segment(aes(x=step_start, xend=step_end, y=logL))+
  geom_line(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  geom_point(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  

  sem_scale + ylab("logL") + common_theme+guides(colour="none") + xcommon+xlab("time (ms)")
fE_MR



## ----fig.width = 7.0, fig.height=9--------------------------------------------
# channel (truth) added to the observation row; ONE panel carries the shared legend
chan <- geom_line(data = d_s, aes(x = step_middle, y = -patch_current, color = "channel current"), linewidth = 0.5)
# realised open fraction = channel current / N_ch: the ground truth the P(open) belief is tracking
chan_P <- geom_line(data = d_s, aes(x = step_middle, y = -patch_current / Nch, color = "channel current"), linewidth = 0.5)
# legend keys (NA coords -> invisible, legend only; logLikelihood excluded via the scale breaks):
# the three PHASES are operations -> arrow keys; channel current is the truth -> a plain line key.
.ops <- c("Markov & prior", "observation & innovation", "Bayes & posterior")
legend_arrow <- geom_segment(data = data.frame(role = .ops, x = NA_real_, xend = NA_real_, y = NA_real_, yend = NA_real_),
                             aes(x = x, xend = xend, y = y, yend = yend, colour = role), na.rm = TRUE,
                             arrow = arrow(length = unit(0.05, "in"), type = "closed"))
legend_line  <- geom_line(data = data.frame(role = "channel current", x = NA_real_, y = NA_real_),
                          aes(x = x, y = y, colour = role), na.rm = TRUE)
# algorithm names spelled out (defined at first use), as a subtitle under each abbreviation
# VR's spelled-out name follows the roster's own logic: the prefix letter says WHAT the
# conditioning is about (Mean, Variance, Interval), so V = Variance. It is one string,
# changed here and nowhere else if the paper settles on another wording.
FULL <- c(LSE = "Least Squares", NR = "Non-Recursive", MNR = "Mean Non-Recursive", R = "Recursive",
          MR = "Mean Recursive", VR = "Variance Recursive", IR = "Interval Recursive")
# column headers live in their OWN strip row (colhead) so the top-left panel is free to carry the
# row's filter-phase title, facet_grid style (algorithms on top, phases on the left of each row).
colhead <- function(ab) ggplot() + theme_void() + ggtitle(ab, subtitle = FULL[[ab]]) +
  theme(plot.title    = element_text(hjust = 0.5, size = 9),
        plot.subtitle = element_text(hjust = 0.5, size = 7))
# (the legend is attached inside build_figure, to whichever column comes first)
scP  <- scale_y_continuous(limits = YP)     # prior & posterior (full-data range, window-independent)
scI  <- scale_y_continuous(limits = YI)     # observation current (full-data range, window-independent)

# CAPTION NOTE: in the observation row the predictive spread is a BAND for gmean_i / averaged-
# conductance algorithms (NMR, MR, VR, IR) and an ERRORBAR for plain-g algorithms (NR, R) — the
# band-vs-errorbar distinction encodes the averaging (av) axis, on purpose.

# CAPTION NOTE: in the Markov (prior) and Bayes (posterior) panels of the recursive algorithms
# (R, MR, VR, IR) each arrow starts from a faded mark, the state carried in from the previous step, and
# points to the full-colour state that the panel computes; the non-recursive algorithms (NR, NMR) have
# no such carried-in state.

# --- build ONE figure for a given time window and save it. -----------------------------
# COLUMNS ARE A PARAMETER (`cols`): this file is the single source of the figure, and the two
# notebooks that use it differ only in the roster they pass.
#   figure_1.Rmd      -> c("R","MR","VR","IR")            the paper figure, the recursive ladder
#   figure_1_all.Rmd  -> c("NR","MNR","R","MR","VR","IR")  every algorithm, for exploration
# Nothing below hardcodes a column name: the FIRST column of `cols` carries the y quantity and
# the legend, the middle one carries the x title, and the "no update (open loop)" placeholder is
# used for whichever of NR/MNR are present.
#
# sel = NULL is the whole recording; sel = c(i, j) crops the TIME axis to those windows (full data
# kept, so every lag/lead arrow is exactly the full-graph arrow, just zoomed). logL is the exception:
# it re-accumulates from 0 at the window start, and when zooming all columns share ONE logL axis
# (unified). The full recording keeps two logL blocks (the non-recursive pair crashes to ~-95) when
# any non-recursive column is present. y-NUMBERS live only on the first column (+ stage-label title);
# every other column keeps tick MARKS only. x-NUMBERS live on the bottom row (D). ---

# panel registries, so a roster is just a vector of names
.DAT <- list(LSE = d_LSE, NR = d_NR, NMR = d_NMR, R = d_R, MR = d_MR, VR = d_VR, IR = d_IR)
.FB  <- list(LSE = fB_LSE, NR = fB_NR, NMR = fB_NMR, R = fB_R, MR = fB_MR, VR = fB_VR, IR = fB_IR)
.FC  <- list(LSE = fC_LSE, NR = fC_NR, NMR = fC_NMR, R = fC_R, MR = fC_MR, VR = fC_VR, IR = fC_IR)
.FD  <- list(R = fD_R, MR = fD_MR, VR = fD_VR, IR = fD_IR)     # non-recursive members (LSE, NR, NMR) have no posterior update
.NAIVE <- c("LSE", "NR", "NMR")
# A roster is written in DATA keys. The mean-non-recursive algorithm is a pre-existing wart: its data
# key is NMR, its file token and its printed header are MNR (see the filenames/algorithms vectors and
# the FULL map). Keep the two apart here rather than propagating the wart into the roster.
.DISP <- c(LSE = "LSE", NR = "NR", NMR = "MNR", R = "R", MR = "MR", VR = "VR", IR = "IR")

build_figure <- function(sel, outfile, cols) {
  stopifnot(length(cols) >= 2, all(cols %in% names(.DAT)))
  unify <- !is.null(sel)
  first <- cols[1]                                  # carries the y title and the legend
  xcol  <- cols[ceiling(length(cols) / 2)]          # carries the single x title
  win <- function(df) if (is.null(sel)) df else dplyr::filter(df, sample_index %in% sel)
  dW  <- lapply(.DAT[cols], win)

  # time crop
  if (is.null(sel)) {
    xlim <- range(c(0, d_s$step_middle, d_IR$step_start, d_IR$step_end), na.rm = TRUE)
  } else {
    xlim <- range(c(win(d_s)$step_start, win(d_s)$step_end), na.rm = TRUE)
  }

  # AUTO y-ranges from the (windowed) data, SHARED across the roster per row. The full recording
  # keeps the data range (which reaches 0); a zoom uses the window's own range (non-zero) so the
  # arrows are large, padded a touch so arrowheads at the extremes are not clipped.
  pad <- function(r, f = 0.06) r + c(-1, 1) * diff(r) * f
  yP <- range(unlist(lapply(dW, .dP)), na.rm = TRUE)
  yI <- range(c(-win(d_s)$patch_current, unlist(lapply(dW, .dI))), na.rm = TRUE)
  if (!is.null(sel)) { yP <- pad(yP); yI <- pad(yI) }

  # windowed cumulative logL, re-anchored to 0 at the window start
  lw <- function(df) {
    d <- dplyr::arrange(win(df), step_end)
    tibble::tibble(t = c(min(d$step_start), d$step_end), cL = c(0, cumsum(d$logL)))
  }
  fE <- function(df) ggplot(lw(df), aes(t, cL)) +
    geom_line(aes(colour = "logLikelihood")) + geom_point(aes(colour = "logLikelihood")) +
    sem_scale + common_theme + guides(colour = "none")
  naive <- intersect(cols, .NAIVE)
  recg  <- setdiff(cols, .NAIVE)
  if (unify || !length(naive)) {          # one logL axis: zoomed, or no non-recursive column present
    yLn <- yLr <- range(unlist(lapply(.DAT[cols], function(x) lw(x)$cL)), na.rm = TRUE)
  } else {                                # full recording WITH the naive pair: two blocks
    yLn <- range(unlist(lapply(.DAT[naive], function(x) lw(x)$cL)), na.rm = TRUE)
    yLr <- range(unlist(lapply(.DAT[recg],  function(x) lw(x)$cL)), na.rm = TRUE)
  }

  noUpd <- function() ggplot() +          # empty open-loop panel, text centred in THIS window
    annotate("text", x = mean(xlim), y = mean(yP), label = "no update\n(open loop)", size = 2.5, colour = "grey60") +
    common_theme

  # per-cell placement. ONE coord per cell crops BOTH axes (x = window; y = the row's shared auto
  # range; coord CLIPS instead of dropping, so arrows at the edge survive).
  cell <- function(p, ylim, col, bottom = FALSE, tag = NULL, yl = NULL, title = NULL) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
    if (col == first) p <- p + ylab(yl) + theme(axis.title.y = element_text(size = 8))
    else if (length(naive) && col == recg[1] && bottom && !unify) p <- p + y_no_title  # 2nd logL axis
    else p <- p + y_only_ticks
    if (bottom) {
      if (col == xcol) p <- p + xlab("time (ms)") else p <- p + theme(axis.title.x = element_blank())
    } else p <- p + x_only_ticks
    if (!is.null(tag)) p <- p + labs(tag = tag)
    if (!is.null(title)) p <- p + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0, size = 8.5, face = "plain"))
    p
  }

  # the legend rides on the first column's B panel, whichever algorithm that is
  legendise <- function(p) p + legend_arrow + legend_line +
    guides(colour = guide_legend(override.aes = list(linetype = 1, shape = NA, linewidth = 1.0)))

  rowA <- lapply(cols, function(cl) {
    p <- .FB[[cl]] + chan_P
    if (cl == first) cell(legendise(p), yP, cl, tag = "A", yl = expression(P[open]), title = "Markov & prior")
    else cell(p, yP, cl)
  })
  rowB <- lapply(cols, function(cl) {
    p <- .FC[[cl]] + chan
    if (cl == first) cell(p, yI, cl, tag = "B", yl = expression("current (pA)"), title = "prediction & innovation")
    else cell(p, yI, cl)
  })
  rowC <- lapply(cols, function(cl) {
    p <- if (cl %in% .NAIVE) noUpd() else .FD[[cl]] + chan_P
    if (cl == first) cell(p, yP, cl, tag = "C", yl = expression(P[open]), title = "Bayes & posterior")
    else cell(p, yP, cl)
  })
  rowD <- lapply(cols, function(cl) {
    yl_ <- if (cl %in% naive) yLn else yLr
    if (cl == first) cell(fE(.DAT[[cl]]), yl_, cl, bottom = TRUE, tag = "D",
                          yl = expression(Sigma~"logL"), title = "logLikelihood")
    else cell(fE(.DAT[[cl]]), yl_, cl, bottom = TRUE)
  })

  g <- do.call(wrap_plots, c(
    lapply(cols, function(cl) colhead(.DISP[[cl]])), rowA, rowB, rowC, rowD,
    list(ncol = length(cols), heights = c(0.01, 1, 1, 1, 1), guides = "collect"))) &
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8),
          legend.key.width = unit(1, "lines"),
          legend.key.spacing.x = unit(3, "pt"),
          plot.margin = margin(1.5, 2, 1.5, 2)) &
    scale_x_continuous(labels = function(x) x * 1000)   # seconds -> ms

  ggsave(outfile, g, width = 7.0, height = 7.5)
  g
}
