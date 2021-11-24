library(data.table)
library(ggplot2)
library(haven)
library(survminer)
library(survival)
library(survey)
library(grid)
library(gridExtra)
library(flextable)
library(officer)
gLegend<-function(a.plot){
  if ("ggplot" %in% class(a.plot)) {
    tmp <- ggplot_gtable(ggplot_build(a.plot))
  } else if ("grob" %in% class(a.plot)) {
    tmp <- .gplot
  }
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

## Set path to local repo.
repo <- 'C:/Users/ncgra/Dropbox/Penn/repos/dem_res_2021/'

## Load clean input data.
cj_surv <- readRDS(paste0(repo,'/clean_psid_data.RDS'))

## KM function
km_curve_plot <- function(i, d=cj_surv, plot_ci=F, return='gg', ref_cat='race_eth=White, as.factor(get(edu_var))=4', custom_ylimits=NULL, drop_hispanic=FALSE) {
  v <- arg_table[i,v]
  target_sex <- arg_table[i,target_sex]
  target_edu <- arg_table[i,target_edu]
  target_weight <- arg_table[i,target_weight]
  if(target_weight=='Longitudinal weight') {
    weight_var <- 'svywt_orig_cds'
    cj_surv_sub <- d[race_eth %in% c('White','Black','Hispanic') & !is.na(get(weight_var)) & id_tas %in% cds$id,]
  }
  if(target_weight=='Back fill weight') {
    weight_var <- 'full_long_weight'
    cj_surv_sub <- d[race_eth %in% c('White','Black','Hispanic') & !is.na(get(weight_var)) & id_tas %in% cds$id,]
  }
  if(target_weight=='CDS1997 weight') {
    weight_var <- 'cds1997_weight'
    cj_surv_sub <- d[race_eth %in% c('White','Black','Hispanic') & !is.na(get(weight_var)) & id_tas %in% cds$id,]
  }
  message(target_edu)
  edu_var <- 'none'
  if(target_sex!='Both') cj_surv_sub <- cj_surv_sub[ta_dem_male==ifelse(target_sex=='Male',1,0) & !is.na(get(weight_var)),]
  if(target_edu=='Respondent education') edu_var <- 'ta_edu_agg'
  if(target_edu=='Maternal education') edu_var <- 'mom_edu_agg'
  if(target_edu=='Maternal edu detailed') edu_var <- 'mom_edu'
  if(target_edu=='Highest education') edu_var <- 'highest_edu_agg'
  if(target_edu=='Highest edu detailed') edu_var <- 'highest_edu'
  if(target_edu!='AllEdu') {
    ## Weighted KM
    message(edu_var)
    if(drop_hispanic) cj_surv_sub <- cj_surv_sub[race_eth!='Hispanic']
    mort <- survfit(Surv(get(paste0('ta_cj_f_',v,'_age')),get(v)) ~ race_eth + as.factor(get(edu_var)),
                    data=cj_surv_sub,
                    weights=get(weight_var))
    ## Unweighted KM (for sample size table)
    mort_unweighted <- survfit(Surv(get(paste0('ta_cj_f_',v,'_age')),get(v)) ~ race_eth + as.factor(get(edu_var)),
                               data=cj_surv_sub)
  }
  if(target_edu=='AllEdu') {
    ## Weighted KM
    message(edu_var)
    mort <- survfit(Surv(get(paste0('ta_cj_f_',v,'_age')),get(v)) ~ race_eth,
                    data=cj_surv_sub,
                    weights=get(weight_var))
    ## Unweighted KM (for sample size table)
    mort_unweighted <- survfit(Surv(get(paste0('ta_cj_f_',v,'_age')),get(v)) ~ race_eth,
                               data=cj_surv_sub)
  }
  ## Assign color palette
  color_cb <- color_cb[strata %in% names(mort$strata),]
  color_cb[, strata := factor(strata, levels=names(mort$strata))]
  color_cb <- color_cb[order(strata)]
  plot_colors <- color_cb[, color]
  color_labels <- color_cb[, label]
  ## Delete Hispanic CIs
  if(plot_ci){
    if(target_edu=='AllEdu') {
      mort$lower[(mort$strata[['race_eth=Black']]+1):(mort$strata[['race_eth=Black']]+mort$strata[['race_eth=Hispanic']])] <- NA
      mort$upper[(mort$strata[['race_eth=Black']]+1):(mort$strata[['race_eth=Black']]+mort$strata[['race_eth=Hispanic']])] <- NA
    }
  }
  ## Make plots of curves
  if(plot_ci) ylimits <- c(0,1-min(mort$lower[mort$time==26]))
  if(!plot_ci) ylimits <- c(0,1-min(mort$surv[mort$time==28]))
  if(!is.null(custom_ylimits)) ylimits <- custom_ylimits
  km_plot <- function(km_fit,edu_var) {
    ggsurvplot(km_fit,
               fun='event',
               data = cj_surv_sub,
               size = 2,
               xlim = c(12,26),
               ylim = ylimits,
               palette = plot_colors,
               conf.int = plot_ci,
               pval = FALSE,
               risk.table = TRUE,
               risk.table.y.text.col = T,
               risk.table.y.text = FALSE,
               legend.labs = color_labels,
               risk.table.height = 0.25, 
               break.time.by = 2,
               ggtheme = theme_bw())
  }
  gg <- km_plot(mort,edu_var)
  gg <- gg + labs(x='Age',y='Cumulative risk')
  gg_unweighted <- km_plot(mort_unweighted,edu_var)
  ## Grab unweighted sample size table.
  gg$table <- gg_unweighted$table
  ## Add title.
  if(v=='arrest') var <- 'Arrest'
  if(v=='prob') var <- 'Probation'
  if(v=='incar') var <- 'Incarceration'
  gg$plot <- gg$plot + ggtitle(paste0(var,': ',target_sex,', ',target_edu,', ',target_weight))
  ## Get p-values from log-rank test (survey weighted) 
  if(target_edu=='AllEdu') {
    logrank_test <- function(compare_race,var) {
      message(paste0('ta_cj_f_',var,'_age'))
      des <- svydesign(id=~id_tas, weight=as.formula(paste0('~',weight_var)), data=cj_surv_sub[race_eth %in% compare_race,])
      if(v=='arrest') p <- summary(svycoxph(Surv(ta_cj_f_arrest_age,arrest) ~ as.factor(race_eth), design=des))$coefficients[,6]
      if(v=='prob') p <- summary(svycoxph(Surv(ta_cj_f_prob_age,prob) ~ as.factor(race_eth), design=des))$coefficients[,6]
      if(v=='incar') p <- summary(svycoxph(Surv(ta_cj_f_incar_age,incar) ~ as.factor(race_eth), design=des))$coefficients[,6]
      ps <- data.table(race_eth=paste(compare_race,collapse='-'), p=p)
      return(ps)
    }
    if(!drop_hispanic) pvals <- rbindlist(lapply(list(c('White','Black'),c('White','Hispanic'),c('Black','Hispanic')), logrank_test, var=v))
    if(drop_hispanic) pvals <- logrank_test(c('White','Black'), var=v)
  }
  if(target_edu!='AllEdu') {
    logrank_test <- function(compare_edu, compare_race) {
      des <- svydesign(id=~id_tas, weight=~cds1997_weight, data=cj_surv_sub[race_eth_edu %in% c(ref_cat, paste0('race_eth=',compare_race,', as.factor(get(edu_var))=',compare_edu))])
      if(v=='arrest') p <- summary(svycoxph(Surv(ta_cj_f_arrest_age,arrest) ~ factor(race_eth), design=des))$coefficients[,6]
      if(v=='prob') p <- summary(svycoxph(Surv(ta_cj_f_prob_age,prob) ~ factor(race_eth), design=des))$coefficients[,6]
      if(v=='incar') p <- summary(svycoxph(Surv(ta_cj_f_incar_age,incar) ~ factor(race_eth), design=des))$coefficients[,6]
      return(data.table(race_eth=compare_race, edu=compare_edu, p=p))
    }
    if(grepl('White',ref_cat)) {
      compare_races <- c('Black','Hispanic')
      if(drop_hispanic) compare_races <- c('Black')
      cj_surv_sub[, race_eth := factor(race_eth, levels=c('White','Black','Hispanic'))]
    }
    if(grepl('Black',ref_cat)) {
      compare_races <- c('White','Hispanic')
      if(drop_hispanic) compare_races <- c('White')
      cj_surv_sub[, race_eth := factor(race_eth, levels=c('Black','White','Hispanic'))]
    }
    if(grepl('Hispanic',ref_cat)) {
      compare_races <- c('White','Black')
      cj_surv_sub[, race_eth := factor(race_eth, levels=c('Hispanic','Black','White'))]
    }
    pvals <- rbindlist(lapply(2:4, logrank_test, compare_race=compare_races[1]))
    if(!drop_hispanic) {
      pvals_black <- rbindlist(lapply(2:4, logrank_test, compare_race=compare_races[1]))
      pvals_hisp <- rbindlist(lapply(2, logrank_test, compare_race=compare_races[2]))
      pvals <- rbind(pvals_black, pvals_hisp)
    }
  }
  if(return=='pvalue') return(pvals)
  if(return=='gg') return(gg)
  if(return=='mort') return(mort)
}
## Estimate KM survival curves for each indicator and strata.
arg_table <- expand.grid(list(c('arrest','prob','incar'),
                              c('Both','Male','Female'),
                              c('AllEdu','Respondent education','Maternal education','Maternal edu detailed','Highest education','Highest edu detailed'),
                              c('CDS1997 weight','Longitudinal weight')))
color_cb <- data.table(strata=c('race_eth=Black','race_eth=Hispanic','race_eth=White',
                                paste0('race_eth=Black, as.factor(get(edu_var))=',1:4),
                                paste0('race_eth=Hispanic, as.factor(get(edu_var))=',1:4),
                                paste0('race_eth=White, as.factor(get(edu_var))=',1:4),
                                paste0('race_eth=Black, weight=',c('Baseline ','Attrition')),
                                paste0('race_eth=Hispanic, weight=',c('Baseline ','Attrition')),
                                paste0('race_eth=White, weight=',c('Baseline ','Attrition')),
                                paste0('race_eth=Black, as.factor(get(edu_var))=',c('Cohort 84-90','Cohort 91-97')),
                                paste0('race_eth=Hispanic, as.factor(get(edu_var))=',c('Cohort 84-90','Cohort 91-97')),
                                paste0('race_eth=White, as.factor(get(edu_var))=',c('Cohort 84-90','Cohort 91-97'))),
                       label=c('Black','Hispanic','White',
                               paste0('Black, ',c('<HS','HS','Some college','>=College')),
                               paste0('Hispanic, ',c('<HS','HS','Some college','>=College')),
                               paste0('White, ',c('<HS','HS','Some college','>=College')),
                               paste0('Black, ',c('Baseline','Attrition')),
                               paste0('Hispanic, ',c('Baseline','Attrition')),
                               paste0('White, ',c('Baseline','Attrition')),
                               paste0('Black, ',c('Cohort 84-90','Cohort 91-97')),
                               paste0('Hispanic, ',c('Cohort 84-90','Cohort 91-97')),
                               paste0('White, ',c('Cohort 84-90','Cohort 91-97'))),
                       color=c('#cb181d','#2171b5','#238b45',
                               '#fee0d2','#fc9272','#ef3b2c','#a50f15',
                               '#deebf7','#9ecae1','#4292c6','#08519c',
                               '#e5f5e0','#a1d99b','#41ab5d','#006d2c',
                               '#a50f15','#fc9272',
                               '#08519c','#9ecae1',
                               '#006d2c','#a1d99b',
                               '#a50f15','#fc9272',
                               '#08519c','#9ecae1',
                               '#006d2c','#a1d99b'))
arg_table <- as.data.table(arg_table)
setnames(arg_table, c('v','target_sex','target_edu','target_weight'))
arg_table <- arg_table[, lapply(.SD, as.character)]
todays_date <- gsub('-','_',Sys.Date())
dir.create(paste0(repo,'/',todays_date))
cj_surv[, race_eth_edu := paste0('race_eth=',race_eth,', as.factor(get(edu_var))=',mom_edu_agg)]

## Figure 1: 
  # A. Men, arrests by race
  # B. Men, probation by race
  # C. Men, incarceration by race
  # D. Women, arrests by race
  # E. Women, probation by race
  # F. Women, incarceration by race
all_plots <- lapply(1:9, km_curve_plot, d=cj_surv, plot_ci=T, return='gg')
target_plots <- all_plots[4:9]
legend_plot <- ggplot() + geom_line(data=cj_surv, aes(x=max_age,y=max_age,color=race_eth)) +
  scale_color_manual(values=c('#cb181d','#2171b5','#238b45'),labels=c('Black','Latinx','White'), name='', guide=guide_legend(ncol=3, override.aes = list(fill=NA,size=2))) + 
  theme_bw() + 
  theme(legend.text = element_text(size=16))
target_plots[[1]] <- target_plots[[1]]$plot + ggtitle('(A) Men, arrests') + theme(legend.position = 'none') + scale_fill_manual(values=c('#cb181d',NA,'#238b45'),labels=c('Black','Latinx','White'))
target_plots[[2]] <- target_plots[[2]]$plot + ggtitle('(B) Men, probation') + theme(legend.position = 'none') + scale_fill_manual(values=c('#cb181d',NA,'#238b45'),labels=c('Black','Latinx','White'))
target_plots[[3]] <- target_plots[[3]]$plot + ggtitle('(C) Men, incarceration') + theme(legend.position = 'none') + scale_fill_manual(values=c('#cb181d',NA,'#238b45'),labels=c('Black','Latinx','White'))
target_plots[[4]] <- target_plots[[4]]$plot + ggtitle('(D) Women, arrests') + theme(legend.position = 'none') + scale_fill_manual(values=c('#cb181d',NA,'#238b45'),labels=c('Black','Latinx','White'))
target_plots[[5]] <- target_plots[[5]]$plot + ggtitle('(E) Women, probation') + theme(legend.position = 'none') + scale_fill_manual(values=c('#cb181d',NA,'#238b45'),labels=c('Black','Latinx','White'))
target_plots[[6]] <- target_plots[[6]]$plot + ggtitle('(F) Women, incarceration') + theme(legend.position = 'none')+ scale_fill_manual(values=c('#cb181d',NA,'#238b45'),labels=c('Black','Latinx','White'))
target_plots[[7]] <- gLegend(legend_plot)
png(paste0(repo,'/',todays_date,'/Figure1.png'),height=9.8,width=16.8,res=900,units = 'in')
lay <- rbind(c(1,1,1,1,2,2,2,2,3,3,3,3),
             c(1,1,1,1,2,2,2,2,3,3,3,3),
             c(1,1,1,1,2,2,2,2,3,3,3,3),
             c(4,4,4,4,5,5,5,5,6,6,6,6),
             c(4,4,4,4,5,5,5,5,6,6,6,6),
             c(4,4,4,4,5,5,5,5,6,6,6,6),
             c(7,7,7,7,7,7,7,7,7,7,7,7))
grid.arrange(grobs=target_plots,
             layout_matrix=lay)
dev.off()

## Figure 2:
# A. Men, arrests by race & highest ed
# B. Men, probation by race & highest ed
# C. Men, incarceration by race & highest ed
cj_surv[, race_eth_edu := paste0('race_eth=',race_eth,', as.factor(get(edu_var))=',highest_edu_agg)]
cj_surv[, race_eth_edu := factor(race_eth_edu, levels=c(paste0('race_eth=Black, as.factor(get(edu_var))=',2:4),
                                                        paste0('race_eth=White, as.factor(get(edu_var))=',2:4),
                                                        paste0('race_eth=Hispanic, as.factor(get(edu_var))=',2:4)))]
edu_colors <- color_cb[strata %in% c(paste0('race_eth=Black, as.factor(get(edu_var))=',2:4),
                                     paste0('race_eth=White, as.factor(get(edu_var))=',2:4))]
edu_colors <- merge(edu_colors, data.table(strata=c(c(paste0('race_eth=Black, as.factor(get(edu_var))=',2:4),paste0('race_eth=White, as.factor(get(edu_var))=',2:4))),order=1:6))
edu_colors <- edu_colors[order(order)]
edu_labels <- edu_colors[strata %in% c(paste0('race_eth=Black, as.factor(get(edu_var))=',2:4),paste0('race_eth=White, as.factor(get(edu_var))=',2:4)), label]
edu_colors <- edu_colors[,color]
names(edu_colors) <- c(paste0('race_eth=Black, as.factor(get(edu_var))=',2:4),paste0('race_eth=White, as.factor(get(edu_var))=',2:4))
target_plots <- lapply(40:42, km_curve_plot, d=cj_surv[!(race_eth_edu %in% c('race_eth=Hispanic, as.factor(get(edu_var))=4',
                                                                             'race_eth=Hispanic, as.factor(get(edu_var))=3',
                                                                             'race_eth=Hispanic, as.factor(get(edu_var))=2',
                                                                             'race_eth=Hispanic, as.factor(get(edu_var))=NA')),], plot_ci=F, return='gg', drop_hispanic=TRUE)
legend_plot <- ggplot() + geom_line(data=cj_surv[!(race_eth_edu %in% c('race_eth=Hispanic, as.factor(get(edu_var))=4',
                                                                       'race_eth=Hispanic, as.factor(get(edu_var))=3',
                                                                       'race_eth=Hispanic, as.factor(get(edu_var))=2')) & !is.na(race_eth_edu),], aes(x=max_age,y=max_age,color=race_eth_edu)) +
  scale_color_manual(values=edu_colors, labels=edu_labels, name='', guide=guide_legend(ncol=2, override.aes = list(fill=NA,size=2))) + 
  theme_bw() + 
  theme(legend.text = element_text(size=16))
target_plots[[1]] <- target_plots[[1]]$plot + ggtitle('(A) Men, arrests') + theme(legend.position = 'none')
target_plots[[2]] <- target_plots[[2]]$plot + ggtitle('(B) Men, probation') + theme(legend.position = 'none')
target_plots[[3]] <- target_plots[[3]]$plot + ggtitle('(C) Men, incarceration') + theme(legend.position = 'none')
target_plots[[4]] <- gLegend(legend_plot)
png(paste0(repo,'/',todays_date,'/Figure2.png'),height=14,width=9,res=900,units = 'in')
lay <- rbind(c(1,1,1,1),
             c(1,1,1,1),
             c(1,1,1,1),
             c(2,2,2,2),
             c(2,2,2,2),
             c(2,2,2,2),
             c(3,3,3,3),
             c(3,3,3,3),
             c(3,3,3,3),
             c(4,4,4,4))
grid.arrange(grobs=target_plots,
             layout_matrix=lay)
dev.off()

## Table 1
## Spits out many versions of the table with p-values calculated against different baseline race-eth-edu comparisons
for(r in c('Black','White','Hispanic')) {
if(r %in% c('Black','White')) edu_levels <- c('2','3','4')
if(r %in% c('Hispanic')) edu_levels <- c('2','3','4')
for(e in edu_levels) {
if(e=='2') e_name <- 'High school'
if(e=='3') e_name <- 'Some college'
if(e=='4') e_name <- 'College'
target_args <- arg_table[c(40:45),]
all_plots <- lapply(c(40:45), km_curve_plot, d=cj_surv, return='gg', ref_cat=paste0('race_eth=',r,', as.factor(get(edu_var))=',e), drop_hispanic=FALSE)
all_pvals <- lapply(c(40:45), km_curve_plot, d=cj_surv, return='pvalue', ref_cat=paste0('race_eth=',r,', as.factor(get(edu_var))=',e), drop_hispanic=FALSE)
get_ages <- function(i, reshape_ss=F, add_ss=F) {
  d <- as.data.table(all_plots[[i]]$data.survplot)
  d <- d[time %in% c(18,26),c('time','surv','strata','upper','lower')]
  d <- merge(d, color_cb[, c('strata','label')], by='strata')
  d[, surv := as.character(round((1 - surv)*100))]
  d[, strata := NULL]
  d[, sex := target_args[i, target_sex]]
  d[, cj := target_args[i, v]]
  d[, edu_type := target_args[i, target_edu]]
  ## Add p-values (comparing entire curves)
  pvals <- all_pvals[[i]]
  if('AllEdu' %in% d[, edu_type]) pvals[, label := race_eth]
  if(!('AllEdu' %in% d[, edu_type])) {
    pvals[edu==2, label := paste0(race_eth,', HS')]
    pvals[edu==3, label := paste0(race_eth,', Some college')]
    pvals[edu==4, label := paste0(race_eth,', >=College')]
  }
  for(r in pvals[,label]) {
    p <- ifelse(pvals[label==r,p]<0.05, '*', '')
    p <- ifelse(pvals[label==r,p]<0.01, '**', p)
    p <- ifelse(pvals[label==r,p]<0.001, '***', p)
    if(r=="Black-Hispanic") {
      p <- gsub('*','^',p)
      this_r <- 'Black'
    }
    this_r <- gsub('White-','',r)
    ## Continuous p-value
    #p <- round(pvals[label==r,p],3)
    p <- '' ## Do not include in final tables per Dem Research guidelines
    d[label==this_r, surv := paste0(surv,'\n(',round((1-upper)*100),'-',round((1-lower)*100),')',p)]
  }
  ## Add sample sizes
  if(add_ss) {
    ss <- as.data.table(all_plots[[i]]$table$data)
    setnames(ss, 'strata', 'label')
    if(reshape_ss) {
      ss <- dcast(ss, label ~ time, value.var = 'n.risk')
      ss[, n.risk := paste0(`0`,'-',`26`)]
      ss[, time := 26]
    }
    d <- merge(d, ss[, c('n.risk','time','label')], by=c('time','label'))
  d[, surv := paste0(surv, '\n(N=', n.risk, ')')]
  }
  return(d)
}
## Table 1
all_mort <- rbindlist(lapply(c(1:6),get_ages,reshape_ss=F,add_ss=F))
all_mort[, race := tstrsplit(label,', ',keep=1)]
all_mort[, edu := tstrsplit(label,', ',keep=2)]
all_mort[edu_type=='Respondent education', edu_type := 'Resp']
all_mort[edu_type=='Highest education', edu_type := 'Parent']
all_mort <- dcast(all_mort, time + cj + race ~ edu_type + edu + sex, value.var = 'surv')
all_mort[, cj := factor(cj, levels=c('arrest','prob','incar'))]
all_mort <- all_mort[order(cj,time,race)]
all_mort <- all_mort[time %in% c(26),]
setcolorder(all_mort, c('cj','time','race',
                        'Parent_HS_Male','Parent_Some college_Male','Parent_>=College_Male',
                        'Parent_HS_Female','Parent_Some college_Female','Parent_>=College_Female'))
all_mort[race=='White', race_order := 1]
all_mort[race=='Black', race_order := 2]
all_mort[race=='Hispanic', race_order := 3]
all_mort <- all_mort[order(cj,race_order)]
all_mort[, race_order := NULL]
for(v in c('Parent_Some college_Male','Parent_>=College_Male','Parent_Some college_Female','Parent_>=College_Female')) all_mort[race=='Hispanic', (v) := '--']
first_header <- as.list(c(rep('',3),rep(c('High school','Some college','College'),2)))
second_header <- as.list(c(rep('',3), rep('Male',3), rep('Female',3)))
names(first_header) <- names(all_mort)
names(second_header) <- names(all_mort)
ft <- flextable(all_mort, theme_fun = theme_booktabs) %>%
  delete_part(part = "header") %>%
  add_header(values=first_header) %>%
  add_header(values=second_header) %>%
  merge_h(part = "header") %>%
  hline_top(part = 'body', border = fp_border(color='black',width=2)) %>%
  hline_top(part = 'header', border = fp_border(color='black',width=2)) %>%
  hline(j=3:dim(all_mort)[2],part = 'header', border = fp_border(color='black',width=2)) %>%
  align(align = 'center', part='body') %>%
  align(align = 'center', part='header') %>%
  align(j=1, align = 'left', part='body') %>%
  add_footer_lines(c('* p < 0.05; ** p < 0.01; *** p < 0.001')) %>%
  flextable::font(fontname = 'Times New Roman', part='all') %>%
  fontsize(size = 12, part='all')
doc <- read_docx() %>%
  body_add_flextable(value = ft, split = TRUE) %>%
  body_end_section_landscape() %>% # a landscape section is ending here
  print(target = paste0(repo,'/',todays_date,'/Table1_',r,'_',e_name,'.docx'))
}
}
