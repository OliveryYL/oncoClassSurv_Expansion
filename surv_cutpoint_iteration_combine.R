#Step1:迭代拆分为多个预后分组####
#surv_cutpoint迭代
surv_cutpoint.iteration<-function(data,
                                  time = NULL,
                                  event = NULL,
                                  variables,
                                  minprop = 0.1,
                                  cut.n=2,
                                  show.survplot=TRUE,
                                  pval = TRUE,
                                  show.pval.detail=T,
                                  fun = "pct",
                                  palette= NULL,
                                  surv.median.line = "hv",
                                  legend="top",
                                  conf.int=TRUE,
                                  size=1.3,
                                  legend.labs=names(survfit.dt$strata),
                                  legend.title="Cluster",
                                  xlab="Time (months)",
                                  ylab='Overall survival (%)',
                                  risk.table=TRUE,
                                  break.time.by = 6,
                                  risk.table.title="Number at risk",
                                  risk.table.height=0.3,
                                  risk.table.y.text = TRUE,
                                  ggtheme = theme_survminer(),
                                  progressbar = TRUE){
  #安装并加载R包
  if(!require("survival",quietly = T)){
    warning('Need to install "survival" package !')
    install.packages("survival")}
  if(!require("survminer",quietly = T)){
    warning('Need to install "survminer" package !')
    install.packages("survminer")}
  #初始值
  surv_data<-data
  cut.point.stat<-data.frame()

  for (i in 1:(cut.n-1)) {
    cat(paste0("Cut.number=",i+1,"\n"))
    survcutpoint<-surv_cutpoint(data =surv_data,
                                time = time,
                                event = event,
                                variables = variables,
                                minprop = minprop,
                                progressbar=progressbar)
    cut.value.i<-survcutpoint$cutpoint$cutpoint
    cut.statistic.i<-survcutpoint$cutpoint$statistic
    cut.point.stat.i<-data.frame(cutpoint=cut.value.i,statistic=cut.statistic.i)
    #保存cut.point.stat
    cut.point.stat<-rbind(cut.point.stat,cut.point.stat.i)
    #当样本量太少时需提前终止循环
    if(is.na(cut.value.i)){warning("End for small samples!")}
    #基于所有clusters进行后续样本量筛选
    surv_data<-data
    surv_data$cut.x<-cut(surv_data[,variables],breaks = c(-Inf,cut.point.stat$cutpoint,+Inf))
    #挑选出数量最多的cluster
    #从数量较多的cluster中进行下一次循环
    subset.max.id<-table(surv_data$cut.x)%>%which.max()%>%names()
    surv_data<-subset(surv_data,cut.x==subset.max.id)
  }
  cut.point.stat$iteration=rownames(cut.point.stat)
  cut.point.stat<-cut.point.stat[order(cut.point.stat$cutpoint),]
  data[,"Cluster"]<-cut(data[,variables],breaks = c(-Inf,cut.point.stat$cutpoint,+Inf))


  if(show.survplot){
    survfit.dt<-survfit(Surv(time = data[,time],event = data[,event]) ~Cluster,data = data)

    if(show.pval.detail==T){
      #log=rank p value计算
      survfit.diff<-survdiff(Surv(time = data[,time],event = data[,event]) ~Cluster,data = data)
      surv.pvalue1<-1-pchisq(survfit.diff$chisq,df=1)
      surv.pvalue2<-format(signif(surv.pvalue1,digits = 3),scientific=T)
      surv.addstar<-ifelse(surv.pvalue1<0.001,"***",ifelse(surv.pvalue1<0.01,"**",ifelse(surv.pvalue1<0.05,"*","(NS.)")))
      surv.pvalue<-paste("Log-rank P: ",surv.pvalue2,surv.addstar,sep = "")
      #!!!#
      pval<-surv.pvalue
    }

    ggsurv.p<-ggsurvplot(fit = survfit.dt,
                         data =data ,
                         pval = pval,
                         fun = fun,
                         palette= palette,
                         surv.median.line = surv.median.line,
                         legend=legend,
                         conf.int=conf.int,
                         size=size,
                         legend.labs=legend.labs,
                         legend.title=legend.title,
                         xlab=xlab,
                         ylab=ylab,
                         risk.table=risk.table,
                         break.time.by = break.time.by,
                         risk.table.title=risk.table.title,
                         risk.table.height=risk.table.height,
                         risk.table.y.text = risk.table.y.text,
                         ggtheme = ggtheme)
    print(ggsurv.p)
    cut.point.stat<-list(cut.point.stat=cut.point.stat,ggsurv.p=ggsurv.p,data=data)
  }
  return(cut.point.stat)
}

# 示例：运行上述迭代函数
# dt<-surv_cutpoint.iteration(data =tcga.surv.dt,
#                             time = "OS",
#                             event = "status",
#                             variables = "OS",
#                             minprop = 0.1,
#                             cut.n = 17,show.survplot = T)


#Step2:自动合并彼此不显著的clusters####
#根据两两比较的pvalue，将彼此不显著的clusters前后合并
surv_cutpoint.combine<-function(data,
                                #method.check="pairwise",
                                #method.by="pairwise.pval",
                                n.iteration.max=5,
                                time = time,
                                event = event,
                                cut.pval=0.05,pval.max=0.05,
                                combine.dist.max=1){#如果距离太远，则不进行合并
  message("Combine method according to pairwise pvalue.")#使用最近邻的方法
  if(!require("survminer",quietly = T)){
    warning('Need to install "survminer" package !')
    install.packages("survminer")}
  if(!require("dplyr",quietly = T)){
    warning('Need to install "dplyr" package !')
    install.packages("dplyr")}
  if(!require("stringr",quietly = T)){
    warning('Need to install "stringr" package !')
    install.packages("stringr")}
  #初始计算pairwise.pval
  surv.f<-as.formula(paste0("Surv(",time,",",event,") ~Cluster"))
  surv.pairwise.x <- pairwise_survdiff(surv.f,data =data)
  select.pairs.condition<-!surv.pairwise.x$p.value<cut.pval
  while.condition<-which(select.pairs.condition)%>%length()>0

  if(!while.condition){warning("All clusters is significant!");break}

  needtocombine.index.0<-order(-surv.pairwise.x$p.value)[!is.na(surv.pairwise.x$p.value)]
  needtocombine.index.0<-needtocombine.index.0[1:sum(select.pairs.condition,na.rm = T)]
  needtocombine.index<-needtocombine.index.0

  n.iteration=1
  skip.n=0
  while (while.condition&n.iteration<=n.iteration.max) {
    surv.pair.locate.max<-needtocombine.index[1+skip.n]
    if(is.na(surv.pair.locate.max)){warning("Up to limmt for NA value!"); break}
    mtx.nrow<-nrow(surv.pairwise.x$p.value)
    mtx.ncol<-ncol(surv.pairwise.x$p.value)
    locate.col<-ceiling(surv.pair.locate.max/mtx.nrow)
    locate.row<-surv.pair.locate.max%%mtx.nrow
    locate.col.name<-colnames(surv.pairwise.x$p.value)[locate.col]
    if (locate.row==0) {
      locate.row.name<-tail(rownames(surv.pairwise.x$p.value),1)
    }else{
      locate.row.name<-rownames(surv.pairwise.x$p.value)[locate.row]
    }
    need.combine<-c(locate.col.name,locate.row.name)

    cat(paste0("Iteration=",n.iteration,"; "))
    cat(paste0("Location=",needtocombine.index.0[n.iteration],"; "))
    cat(paste(need.combine))
    cat(sep = "\n")

    #正则表达式提取区间名称中的数字
    ex.pattern<-"\\((\\d+\\.?\\d*),\\s*(\\d+\\.?\\d*)\\]|\\((\\-Inf\\.?\\d*),\\s*(\\d+\\.?\\d*)\\]|\\((\\d+\\.?\\d*),\\s*(\\s*Inf)\\]"
    need.combine1<-str_match(need.combine[1],ex.pattern )[-1]%>%na.omit()%>%as.numeric() %>%
      rlang::set_names(c("lower", "upper"))
    need.combine2<-str_match(need.combine[2],ex.pattern )[-1]%>%na.omit()%>%as.numeric() %>%
      rlang::set_names(c("lower", "upper"))

    combine.new.lower<-min(need.combine1,need.combine2)
    combine.new.upper<-max(need.combine1,need.combine2)

    new.interval.label<-paste0("(",combine.new.lower,",",combine.new.upper,"]")
    #替换levels(合并levels)
    which.levels<-match(need.combine,levels(data$Cluster))

    #判断是否存在交集
    interval_intersect <- function(a, b) {
      a[1] <= b[2] && b[1] <= a[2]
    }
    interval.condition<-interval_intersect(a = need.combine1,b = need.combine2)

    if(interval.condition){
      #替换levels(合并levels)
      levels(data$Cluster)[which.levels]<-new.interval.label}else{
        other.cut.nums<-setdiff(c(need.combine1,need.combine2),c(combine.new.lower,combine.new.upper))
        if(length(other.cut.nums)>1){
          check.dist<-abs(other.cut.nums[1]-other.cut.nums[2])<combine.dist.max
          check.p.max<-surv.pairwise.x$p.value[surv.pair.locate.max]>pval.max
          if(check.dist|check.p.max){
            #替换levels(合并levels)
            levels(data$Cluster)[which.levels]<-new.interval.label
          }else{skip.n=skip.n+1;message(paste0("Total.skip=",skip.n))}
        }else{
          #替换levels(合并levels)
          levels(data$Cluster)[which.levels]<-new.interval.label
        }
      }

    #更新pairwise结果
    surv.pairwise.x <- pairwise_survdiff(surv.f,data =data)
    select.pairs.condition<-!surv.pairwise.x$p.value<cut.pval
    while.condition<-which(select.pairs.condition)%>%length()>0
    needtocombine.index<-order(-surv.pairwise.x$p.value)[!is.na(surv.pairwise.x$p.value)]
    needtocombine.index<-needtocombine.index[1:sum(select.pairs.condition,na.rm = T)]

    #终止条件1：所有clusters彼此之间pairwise均显著
    if(while.condition==FALSE){
      cat(paste0("The iteration is ",n.iteration,".\n"));
      warning("Done! All clusters are significant with each other!");
      break}

    #终止条件2：达到最大迭代次数
    n.iteration=n.iteration+1
  }
  return(data)
}

#示例：运行合并clusters的代码
# dt.pair<-surv_cutpoint.combine(data =dt$data,time = "OS",event = "status",cut.pval = 0.001,n.iteration.max = 15 )

#运行后再次绘制生存曲线，注意观察是否仍然存在置信区间重叠的clusters。
# surv.pair.x<-pairwise_survdiff(Surv(OS,status) ~Cluster,data = dt.pair)
# select.pairs.condition<-!surv.pair.x$p.value<cut.pval
# select.pairs.condition
#
#
# survfit.dt.x<-survfit(Surv(OS,status) ~Cluster,data = dt.pair)
# ggsurv.p<-ggsurvplot(fit = survfit.dt.x,data =dt.pair ,pval = T,conf.int = T,
#                      fun = "pct",surv.median.line = "hv",
#                      xlab="Time (months)",
#                      ylab='Overall survival (%)',
#                      risk.table=TRUE,
#                      break.time.by = 12,
#                      risk.table.title="Number at risk",
#                      risk.table.height=.3,
#                      risk.table.y.text = T)
# ggsurv.p

#Step3:手动合并clusters（可选项）####
#如果仍然存在重叠，可以进一步选择部分距离过近的clusters进行手动合并。


#x.res.out
x.res.out<-function(x){
  ex.pattern<-"\\((\\d+\\.?\\d*),\\s*(\\d+\\.?\\d*)\\]|\\((\\-Inf\\.?\\d*),\\s*(\\d+\\.?\\d*)\\]|\\((\\d+\\.?\\d*),\\s*(\\s*Inf)\\]"
  x.res<-stringr::str_match(x,ex.pattern )[-1]%>%na.omit()%>%as.numeric() %>%
    rlang::set_names(c("lower", "upper"))
  return(x.res)
}

#Step3.1:check.cluster####
check.cluster<-function(data,vars="Cluster"){
  all.levels.list<-as.list(levels(data[,"Cluster"]))
  l.x<-lapply(all.levels.list, x.res.out)
  df.x<-do.call(rbind,l.x)
  save.intersect<-list()
  for (m in 1:(length(all.levels.list)-1)) {
    for (n in (m+1):length(all.levels.list)) {
      if(df.x[m,"upper"]>df.x[n,"lower"]){
        intersect.pair<-list(c(m,n))
        save.intersect<-c(save.intersect,intersect.pair)
      }
    }
  }

  intersect.summary<-do.call(rbind,save.intersect)

  if(!is.null(intersect.summary)){
    #绘图显示建议合并的Clusters
    plot.window(xlim = c(min(intersect.summary),max(intersect.summary)),
                ylim = c(0,1))
    par(mfrow=c(1,1),mar = c(4, 3, 2, 1))
    point.x<-unique(as.numeric(intersect.summary))
    point.y<-rep(0.5,length(point.x))
    seg.axis.x=seq(min(point.x),max(point.x),by=1)

    plot(x=point.x,y=point.y,type = "p",
         main = "Merge Cluster",
         xlab = "Cluster Number",xaxt="n",
         ylab = "",yaxt="n",
         pch=16,lty=4,col="black", cex=1)
    axis(1,at=seg.axis.x)
    arrows(intersect.summary[,1], c(0.5,0.5,0.5),
           intersect.summary[,2] , c(0.5,0.5,0.5),
           col= 1:dim(intersect.summary)[1],
           length = 0.1,angle = 30,lwd = 2)
    abline(v=point.x,col="red",lty="dashed",lwd=0.5)
  }else{warning("All clusters have no overlapping intervals.")}
  return(intersect.summary)
}
#perform check.cluster:
#check.cluster(data = tcga.dt.pair1)

#Step3.2:merge.cluster####
merge.cluster<-function(data,vars="Cluster",merge.list){
  old.labels<-list()
  newlabels<-c()
  for (i in merge.list) {
    old.label.i<-levels(data[,"Cluster"])[i]
    old.list.label<-paste0(i,collapse = "_")
    old.labels[[old.list.label]]<-old.label.i

    old.i<-lapply(as.list(old.label.i), x.res.out)%>%do.call(rbind,.)
    new.label.i<-paste0("(",min(old.i),",",max(old.i),"]")
    newlabels<-append(newlabels,new.label.i)
  }

  for (otn in 1:length(old.labels)) {
    old.match.id<-match(old.labels[[otn]],levels(data[,"Cluster"]))
    levels(data[,"Cluster"])[old.match.id]<-newlabels[otn]
  }
  return(data)
}
# 示例：perform merge.cluster:
# tcga.dt.pair.new<-merge.cluster(data =tcga.dt.pair,vars = "Cluster",
#                                 merge.list = list(c(4:5),c(6:10)) )




