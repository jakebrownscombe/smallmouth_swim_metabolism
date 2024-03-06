
rsquared.lme=function(modlist) {
  #Iterate over each model in the list
  do.call(rbind,lapply(modlist,function(i) {
    #For models fit using lm
    if(class(i)=="lm") {
      Rsquared.mat=data.frame(Class=class(i),Marginal=summary(i)$r.squared,
                              Conditional=NA,AIC=AIC(i)) } 
    #For models fit using lme4
    else if(class(i)=="mer" | class(i)=="merMod" | class(i)=="merLmerTest") {
      #Get variance of fixed effects by multiplying coefficients by design matrix
      VarF=var(as.vector(fixef(i) %*% t(i@X))) 
      #Get variance of random effects by extracting variance components
      VarRand=colSums(do.call(rbind,lapply(VarCorr(i),function(j) j[1])))
      #Get residual variance
      VarResid=attr(VarCorr(i),"sc")^2
      #Calculate marginal R-squared (fixed effects/total variance)
      Rm=VarF/(VarF+VarRand+VarResid)
      #Calculate conditional R-squared (fixed effects+random effects/total variance)
      Rc=(VarF+VarRand)/(VarF+VarRand+VarResid)
      #Bind R^2s into a matrix and return with AIC values
      Rsquared.mat=data.frame(Class=class(i),Marginal=Rm,Conditional=Rc,
                              AIC=AIC(update(i,REML=F))) } 
    #For model fit using nlme  
    else if(class(i)=="lme") {
      #Get design matrix of fixed effects from model
      Fmat=model.matrix(eval(i$call$fixed)[-2],i$data)
      #Get variance of fixed effects by multiplying coefficients by design matrix
      VarF=var(as.vector(fixef(i) %*% t(Fmat)))
      #Get variance of random effects by extracting variance components
      VarRand=sum(suppressWarnings(as.numeric(VarCorr(i)[rownames(VarCorr(i))!=
                                                           "Residual",1])),na.rm=T)
      #Get residual variance
      VarResid=as.numeric(VarCorr(i)[rownames(VarCorr(i))=="Residual",1])
      #Calculate marginal R-squared (fixed effects/total variance)
      Rm=VarF/(VarF+VarRand+VarResid)
      #Calculate conditional R-squared (fixed effects+random effects/total variance)
      Rc=(VarF+VarRand)/(VarF+VarRand+VarResid)
      #Bind R^2s into a matrix and return with AIC values
      Rsquared.mat=data.frame(Class=class(i),Marginal=Rm,Conditional=Rc,
                              AIC=AIC(update(i,method="ML")))
    } else { print("Function requires models of class lm, lme, mer, or merMod") 
    } } ) ) }

