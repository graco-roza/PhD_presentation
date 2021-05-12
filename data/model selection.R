#model selection function of the paper Graco-Roza 2019
# fm0 = linear
#fm1 =  piecewise
#fm2 = quadratic


mod_sel = function(y, x, ylim=range(y), xlim=range(x), y_label='response',z, scale=c('Large','Small')) {
  
 df<-data.frame(x,y) 
  require(betareg)
  require(strucchange)  
 
  if(scale == 'Large') {  
  dams=c(199.78,255.43,609.47,907.49,980.04)
  
 if (max(y) < 1 ){
   y.break = breakpoints(y ~ 1)
  fm0 <- betareg(y ~ x,type="ML") 
  fm1 <- betareg(y ~ breakfactor(y.break, breaks = length(y.break$breakpoints)), type="ML")
  fm2 <- betareg(y ~ poly(c(x), 2),type="ML")
  models <- list(fm0,fm1,fm2)
  Modnames <- paste(1:length(models))
  val = aictab(cand.set = models, modnames = Modnames, sort = TRUE)
  rsq<-round(models[[as.numeric(row.names(val)[1])]]$pseudo.r.squared,2)
} else {
  y.break = breakpoints(y ~ x)
  fm0 <- lm(y ~ x) 
  fm1 <- lm(y ~ breakfactor(y.break, breaks = length(y.break$breakpoints)))
  fm2 <- lm(y ~ poly(c(x), 2))  
  models <- list(fm0,fm1,fm2)
  Modnames <- paste(1:length(models))
  val = aictab(cand.set = models, modnames = Modnames, sort = TRUE)
  rsq<-round(summary(models[[as.numeric(row.names(val)[1])]])$adj.r.squared,2)
  }

    fig<- ggplot(data=df,aes(x=x,y=y))+
    geom_line(linetype=2, colour='gray60',alpha=.2)+
    geom_point(colour='gray60', size=2, shape=19,alpha=.2)+
    ylim(ylim)+
    theme_tufte(base_size=12,base_family='sans')+
    ylab(y_label)+
    xlab('Watercourse distance (Km)')+
      geom_vline(aes(xintercept=dams[1]),linetype=3)+
      geom_vline(aes(xintercept=dams[2]),linetype=3)+
      geom_vline(aes(xintercept=dams[3]),linetype=3)+
      geom_vline(aes(xintercept=dams[4]),linetype=3)+
      geom_vline(aes(xintercept=dams[5]),linetype=3)+
      annotate('text',dams[1]-100,max(ylim)*.9,label="PB",colour='black')+
      annotate('text',dams[2]+100,max(ylim)*.9,label="ST",colour='black')+
      annotate('text',dams[3]-100,max(ylim)*.9,label="FN",colour='black')+
      annotate('text',dams[4]-100,max(ylim)*.9,label="AN",colour='black')+
      annotate('text',dams[5]+100,max(ylim)*.9,label="IP",colour='black')+
      geom_line(aes(y=fitted(models[[as.numeric(row.names(val)[1])]]),x=x), size=1, colour='black', linetype=2)+
    annotate('text',-Inf,-Inf, hjust=-.13,vjust=0,label=paste0("\u0394","AIC= ",round(val[2,4],2)),colour='red',fontface=2)+
    annotate('text',-Inf,-Inf, hjust=-.2,vjust=0,label=paste0("R","\u00B2","= ",rsq,'\n'),colour='red',fontface=2)


} else if(scale == 'Small') {

  df<-data.frame(x,y)
  rf=factor(coord_large$Reservoir[!is.na(coord_large$Dist_dam)])
  rf = factor(rf,levels(rf)[c(5,7,4,6,2,1,3)])

  if (max(y) < 1 ){
    fm0 <- betareg(y ~ x,type="ML") 
    fm1 <- betareg(y ~ x+I(pmax(x-0.1, 0)), type="ML")
    fm2 <- betareg(y ~ poly(c(x), 2),type="ML")
    models <- list(fm0,fm1,fm2)
    Modnames <- paste(1:length(models))
    val = aictab(cand.set = models, modnames = Modnames, sort = TRUE)
    rsq<-round(models[[as.numeric(row.names(val)[1])]]$pseudo.r.squared,2)
    
    fig<- ggplot(data=df,aes(x=x,y=y))+
      geom_line(linetype=2, colour='gray60',alpha=.2)+
      geom_point(colour='gray60', size=2, shape=19,alpha=.2)+
      ylim(ylim)+
      theme_tufte(base_size=12,base_family='sans')+
      ylab(y_label)+
      xlab('Watercourse distance (Km)')+ 
      geom_line(aes(y=fitted(models[[as.numeric(row.names(val)[1])]]),x=x), size=1, colour='black', linetype=2)+
      annotate('text',-Inf,-Inf, hjust=-.13,vjust=0,label=paste0("\u0394","AIC= ",round(val[2,4],2)),colour='red',fontface=2)+
      annotate('text',-Inf,-Inf, hjust=-.2,vjust=0,label=paste0("R","\u00B2","= ",rsq,'\n'),colour='red',fontface=2)
    
  } else {
    fm0 <- lme4::lmer(y ~ x +(1|rf), REML='FALSE') 
    fm1 <- lme4::lmer(y ~ x+I(pmax(x-0.7, 0))+(1|rf), REML='FALSE')
    fm2 <- lme4::lmer(y ~ poly(c(x), 2)+(1|rf), REML='FALSE')  
    models <- list(fm0,fm1,fm2)
    Modnames <- paste(1:length(models))
    val = aictab(cand.set = models, modnames = Modnames, sort = TRUE)
    rsq<-round(r.squaredGLMM(models[[as.numeric(row.names(val)[1])]])[2],2)
    
    fig<-ggplot(data=df,aes(x=x,y=y))+
      geom_line(linetype=2, colour='gray60',alpha=.2)+
      geom_point(colour='gray60', size=2, shape=19,alpha=.2)+
      ylim(range(y))+
      theme_tufte(base_size=12,base_family='sans')+
      ylab(y_label)+
      xlab('Watercourse distance (Km)')+
      geom_vline(aes(xintercept=0),linetype=3)+
      geom_line(aes(y=predict(models[[as.numeric(row.names(val)[1])]]),x=x, group=rf), size=1, colour='black', linetype=3)+
      geom_line(aes(y=predict(models[[as.numeric(row.names(val)[1])]],re.form=NA),x=x), size=2, colour='black', linetype=1)+  
#      scale_size_manual(name="Predictions", values=c("Subjects"=1, "Population"=3)) +
      annotate('text',-Inf,-Inf, hjust=-.13,vjust=0,label=paste0("\u0394","AIC= ",round(val[2,4],2)),colour='red',fontface=2)+
      annotate('text',-Inf,-Inf, hjust=-.2,vjust=0,label=paste0("R","\u00B2","= ",rsq,'\n'),colour='red',fontface=2)
    
          }
  
}

return(fig)
  } 
