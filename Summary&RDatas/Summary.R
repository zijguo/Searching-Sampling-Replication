########################### Homo Settings #############################
### modified at 5/19: add Bias-aware
homo_settings <- function(VIO.str.options){
  IV.str = 0.5 
  setting.options = c("S1", "S2", "S3", "S4", "S5")
  # VIO.str.options = c(0.2, 0.4)
  n.options = c(500, 1000, 2000, 5000)
  sim.round.options = c(1,2,3,4,5)
  
  summary.Cov.mat = summary.Leng.mat = matrix(NA, nrow=length(setting.options)*
                                                length(VIO.str.options)*
                                                length(n.options), ncol=10)
  colnames(summary.Cov.mat) = c("TSLS","BA","TSHT", "CIIV", "Sear-TSHT", "Sear-CIIV",
                                "Samp-TSHT", "Samp-CIIV", "Union-1", "Union-2")
  summary.rule.mat = matrix(NA, nrow=length(setting.options)*length(VIO.str.options)*length(n.options), ncol=4)
  colnames(summary.rule.mat) = c("Sear-TSHT","Sear-CIIV","Samp-TSHT", "Samp-CIIV")
  
  for(i.setting in 1:length(setting.options)){
    setting = setting.options[i.setting]
    for(i.VIO.str in 1:length(VIO.str.options)){
      VIO.str = VIO.str.options[i.VIO.str]
      for(i.n in 1:length(n.options)){
        n = n.options[i.n]
        ind = (i.setting-1)*length(VIO.str.options)*length(n.options)+
          (i.VIO.str-1)*length(n.options) + i.n
        
        CI.mat.whole = matrix(NA, nrow=500, ncol=6)
        Cov.mat.whole = Leng.mat.whole = matrix(NA, nrow=500, ncol=9)
        cover.union05.whole = length.union05.whole = matrix(NA, nrow=500, ncol=10)
        Rule.mat.whole = matrix(NA, nrow=500, ncol=4)
        for(sim.round in sim.round.options){
          datafile = paste("RDatas/Homo/Homo-Setting",setting,"-Strength",IV.str,"-Violation",VIO.str,
                           "-n",n,"-SimRound",sim.round,".RData",sep="")
          load(datafile)
          CI.mat.whole[((sim.round-1)*100+1):(sim.round*100),] = CI.mat; colnames(CI.mat.whole) = colnames(CI.mat)
          Cov.mat.whole[((sim.round-1)*100+1):(sim.round*100),] = Cov.mat
          Leng.mat.whole[((sim.round-1)*100+1):(sim.round*100),] = Leng.mat
          cover.union05.whole[((sim.round-1)*100+1):(sim.round*100),] = cover.union05
          length.union05.whole[((sim.round-1)*100+1):(sim.round*100),] = length.union05
          Rule.mat.whole[((sim.round-1)*100+1):(sim.round*100),] = Rule.mat
        }
        
        points = (CI.mat.whole[,"TSHT-L"] + CI.mat.whole[,"TSHT-U"])/2
        bias = points - 1
        sd = sqrt(var(bias))
        temp = mean(bias) / sd
        cv_alpha = sqrt(qchisq(1-0.05,df=1, ncp=temp^2))
        CI.new = cbind(points - sd*cv_alpha, points + sd*cv_alpha)
        cov.new = mean((CI.new[,1] < 1)*(CI.new[,2] > 1))
        len.new = 2*sd*cv_alpha
        
        summary.Cov.mat[ind, 2] = cov.new
        summary.Cov.mat[ind, c(1,3,4,5,6,7,8)] = colMeans(Cov.mat.whole)[c(1,2,3,4,5,7,8)]
        summary.Cov.mat[ind, 9] = mean(na.omit(cover.union05.whole[,4]))
        summary.Cov.mat[ind, 10] = mean(na.omit(cover.union05.whole[,9]))
        summary.Leng.mat[ind, 2] = len.new
        summary.Leng.mat[ind, c(1,3,4,5,6,7,8)] = colMeans(Leng.mat.whole)[c(1,2,3,4,5,7,8)]
        summary.Leng.mat[ind,9] = mean(na.omit(length.union05.whole[,4]))
        summary.Leng.mat[ind,10] = mean(na.omit(length.union05.whole[,9]))
        summary.rule.mat[ind,] = colMeans(Rule.mat.whole)
      }
    }
  }
  return(list(summary.Cov.mat=summary.Cov.mat, summary.Leng.mat=summary.Leng.mat,
              summary.rule.mat=summary.rule.mat))
}
## tau0.2 ##
out.1 = homo_settings(c(0.2))
setting.options = c("S1", "S2", "S3", "S4", "S5")
n.options = c(500, 1000, 2000, 5000)
out.table.1 = data.frame(matrix(NA, nrow=length(setting.options)*
                                  length(n.options)*2, ncol=2+10+1))
colnames(out.table.1) = c("Set","n","TSLS","BA","TSHT","CIIV","Sear-TSHT","Sear-CIIV",
                          "Samp-TSHT", "Samp-CIIV","Test","Union-1", "Union-2")
out.table.1[,1] = rep(rep(1:5, rep(4,5)), 2)
out.table.1[,2] = rep(n.options, 10)
out.table.1[,-c(1,2,11)] = rbind(out.1$summary.Cov.mat, out.1$summary.Leng.mat)
out.table.1[,11] = c(out.1$summary.rule.mat[,1], rep(NA, 20))

library(kableExtra)
kbl(out.table.1, "latex", digits=2, align="c")%>%
  add_header_above(c(" "=2,"Oracle"=2," "=2, "Searching"=2, "Sampling"=2," ","Union"=2))%>%
  kable_styling(latex_options = c("scale_down"))%>%
  collapse_rows(columns=1:2)

## tau0.4 ##
out.2 = homo_settings(c(0.4))
setting.options = c("S1", "S2", "S3", "S4", "S5")
n.options = c(500, 1000, 2000, 5000)
out.table.2 = data.frame(matrix(NA, nrow=length(setting.options)*
                                  length(n.options)*2, ncol=2+10+1))
colnames(out.table.2) = c("Set","n","TSLS","BA","TSHT","CIIV","Sear-TSHT","Sear-CIIV",
                          "Samp-TSHT", "Samp-CIIV","Test","Union-1", "Union-2")
out.table.2[,1] = rep(rep(1:5, rep(4,5)), 2)
out.table.2[,2] = rep(n.options, 10)
out.table.2[,-c(1,2,11)] = rbind(out.2$summary.Cov.mat, out.2$summary.Leng.mat)
out.table.2[,11] = c(out.2$summary.rule.mat[,1], rep(NA, 20))

library(kableExtra)
kbl(out.table.2, "latex", digits=2, align="c")%>%
  add_header_above(c(" "=2,"Oracle"=2," "=2, "Searching"=2, "Sampling"=2," ","Union"=2))%>%
  kable_styling(latex_options = c("scale_down"))%>%
  collapse_rows(columns=1:2)

########################### Hetero Settings #############################
### modified at 5/19: add Bias-aware
hetero_settings <- function(VIO.str.options){
  IV.str = 0.5 
  setting.options = c("S1", "S2", "S3", "S4", "S5")
  #VIO.str.options = c(0.2, 0.4)
  n.options = c(500, 1000, 2000, 5000)
  sim.round.options = c(1)
  
  summary.Cov.mat = summary.Leng.mat = matrix(NA, nrow=length(setting.options)*
                                                length(VIO.str.options)*
                                                length(n.options), ncol=8)
  summary.rule.mat = matrix(NA, nrow=length(setting.options)*length(VIO.str.options)*length(n.options), ncol=4)
  colnames(summary.rule.mat) = c("Sear-TSHT","Sear-CIIV","Samp-TSHT", "Samp-CIIV")
  
  for(i.setting in 1:length(setting.options)){
    setting = setting.options[i.setting]
    for(i.VIO.str in 1:length(VIO.str.options)){
      VIO.str = VIO.str.options[i.VIO.str]
      for(i.n in 1:length(n.options)){
        n = n.options[i.n]
        ind = (i.setting-1)*length(VIO.str.options)*length(n.options)+
          (i.VIO.str-1)*length(n.options) + i.n
        
        CI.mat.whole = matrix(NA, nrow=500, ncol=12)
        Cov.mat.whole = Leng.mat.whole = matrix(NA, nrow=500, ncol=7)
        Rule.mat.whole = matrix(NA, nrow=500, ncol=4)
        for(sim.round in sim.round.options){
          datafile = paste("RDatas/Hetero/Hetero-Setting",setting,"-Strength",IV.str,"-Violation",VIO.str,
                           "-n",n,"-SimRound",sim.round,".RData",sep="")
          load(datafile)
          CI.mat.whole[((sim.round-1)*500+1):(sim.round*500),] = CI.mat; colnames(CI.mat.whole) = colnames(CI.mat)
          Cov.mat.whole[((sim.round-1)*500+1):(sim.round*500),] = Cov.mat[,c(2,4,6,7,8,10,11)]
          Leng.mat.whole[((sim.round-1)*500+1):(sim.round*500),] = Leng.mat[,c(2,4,6,7,8,10,11)]
          Rule.mat.whole[((sim.round-1)*500+1):(sim.round*500),] = Rule.mat
        }
        
        points = (CI.mat.whole[,"TSHT2-L"] + CI.mat.whole[,"TSHT2-U"])/2
        bias = points - 1
        sd = sqrt(var(bias))
        temp = mean(bias) / sd
        cv_alpha = sqrt(qchisq(1-0.05,df=1, ncp=temp^2))
        CI.new = cbind(points - sd*cv_alpha, points + sd*cv_alpha)
        cov.new = mean((CI.new[,1] < 1)*(CI.new[,2] > 1))
        len.new = 2*sd*cv_alpha
        
        summary.Cov.mat[ind, 2] = cov.new
        summary.Cov.mat[ind, c(1,3,4,5,6,7,8)] = colMeans(Cov.mat.whole)
        summary.Leng.mat[ind, 2] = len.new
        summary.Leng.mat[ind, c(1,3,4,5,6,7,8)] = colMeans(Leng.mat.whole)
        summary.rule.mat[ind,] = colMeans(Rule.mat.whole)
      }
    }
  }
  return(list(summary.Cov.mat=summary.Cov.mat, summary.Leng.mat=summary.Leng.mat,
              summary.rule.mat=summary.rule.mat))
}
## tau0.2 ##
out.1 = hetero_settings(c(0.2))
setting.options = c("S1", "S2", "S3", "S4", "S5")
n.options = c(500, 1000, 2000, 5000)
out.table.1 = data.frame(matrix(NA, nrow=length(setting.options)*
                                  length(n.options)*2, ncol=2+8+1))
colnames(out.table.1) = c("Set","n","TSLS","BA","TSHT","CIIV","Sear-TSHT","Sear-CIIV",
                          "Samp-TSHT", "Samp-CIIV","Test")
out.table.1[,1] = rep(rep(1:5, rep(4,5)), 2)
out.table.1[,2] = rep(n.options, 10)
out.table.1[,-c(1,2,11)] = rbind(out.1$summary.Cov.mat, out.1$summary.Leng.mat)
out.table.1[,11] = c(out.1$summary.rule.mat[,1], rep(NA, 20))

library(kableExtra)
kbl(out.table.1, "latex", digits=2, align="c")%>%
  add_header_above(c(" "=2,"Oracle"=2," "=2, "Searching"=2, "Sampling"=2," "))%>%
  kable_styling(latex_options = c("scale_down"))%>%
  collapse_rows(columns=1:2)

## tau0.4
out.2 = hetero_settings(c(0.4))
setting.options = c("S1", "S2", "S3", "S4", "S5")
n.options = c(500, 1000, 2000, 5000)
out.table.2 = data.frame(matrix(NA, nrow=length(setting.options)*
                                  length(n.options)*2, ncol=2+8+1))
colnames(out.table.2) = c("Set","n","TSLS","BA","TSHT","CIIV","Sear-TSHT","Sear-CIIV",
                          "Samp-TSHT", "Samp-CIIV","Test")
out.table.2[,1] = rep(rep(1:5, rep(4,5)), 2)
out.table.2[,2] = rep(n.options, 10)
out.table.2[,-c(1,2,11)] = rbind(out.2$summary.Cov.mat, out.2$summary.Leng.mat)
out.table.2[,11] = c(out.2$summary.rule.mat[,1], rep(NA, 20))

library(kableExtra)
kbl(out.table.2, "latex", digits=2, align="c")%>%
  add_header_above(c(" "=2,"Oracle"=2," "=2, "Searching"=2, "Sampling"=2," "))%>%
  kable_styling(latex_options = c("scale_down"))%>%
  collapse_rows(columns=1:2)

################### CIIV ##################
CIIV_setting <- function(setting.options){
  IV.str = 0.4
  #setting.options = c("CIIV1","CIIV2")
  VIO.str.options = c(0.2, 0.4)
  n.options = c(500, 1000, 2000, 5000)
  
  summary.Cov.mat = summary.Leng.mat = matrix(NA, nrow=length(setting.options)*
                                                length(VIO.str.options)*
                                                length(n.options), ncol=9)
  colnames(summary.Cov.mat) = c("TSLS","BA","TSHT", "CIIV", "Sear-TSHT", "Sear-CIIV",
                                "Samp-TSHT", "Samp-CIIV", "Union-1")
  summary.rule.mat = matrix(NA, nrow=length(setting.options)*length(VIO.str.options)*length(n.options), ncol=4)
  colnames(summary.rule.mat) = c("Sear-TSHT","Sear-CIIV","Samp-TSHT", "Samp-CIIV")
  
  for(i.setting in 1:length(setting.options)){
    setting = setting.options[i.setting]
    for(i.VIO.str in 1:length(VIO.str.options)){
      VIO.str = VIO.str.options[i.VIO.str]
      for(i.n in 1:length(n.options)){
        n = n.options[i.n]
        
        CI.mat.whole = matrix(NA, nrow=500, ncol=6)
        Cov.mat.whole = Leng.mat.whole = matrix(NA, nrow=500, ncol=9)
        cover.union05.whole = length.union05.whole = matrix(NA, nrow=500, ncol=10)
        Rule.mat.whole = matrix(NA, nrow=500, ncol=4)
        for(sim.round in 1:10){
          datafile = paste("RDatas/CIIV/CIIV-Setting",setting,"-Strength",IV.str,"-Violation",VIO.str,
                           "-n",n,"-SimRound",sim.round,".RData",sep="")
          load(datafile)
          CI.mat.whole[((sim.round-1)*50+1): (sim.round*50), ] = CI.mat
          Cov.mat.whole[((sim.round-1)*50+1): (sim.round*50), ] = Cov.mat
          Leng.mat.whole[((sim.round-1)*50+1): (sim.round*50), ] = Leng.mat
          cover.union05.whole[((sim.round-1)*50+1): (sim.round*50), ] = cover.union05
          length.union05.whole[((sim.round-1)*50+1): (sim.round*50), ] = length.union05
          Rule.mat.whole[((sim.round-1)*50+1): (sim.round*50), ] = Rule.mat
        }
        
        ind = (i.setting-1)*length(VIO.str.options)*length(n.options)+
          (i.VIO.str-1)*length(n.options) + i.n
        colnames(CI.mat.whole) = colnames(CI.mat)
        points = (CI.mat.whole[,"TSHT-L"] + CI.mat.whole[,"TSHT-U"])/2
        bias = points - 1
        sd = sqrt(var(bias))
        temp = mean(bias) / sd
        cv_alpha = sqrt(qchisq(1-0.05,df=1, ncp=temp^2))
        CI.new = cbind(points - sd*cv_alpha, points + sd*cv_alpha)
        cov.new = mean((CI.new[,1] < 1)*(CI.new[,2] > 1))
        len.new = 2*sd*cv_alpha
        
        summary.Cov.mat[ind,c(1,3,4,5,6,7,8)] = colMeans(Cov.mat.whole)[c(1,2,3,4,5,7,8)]
        summary.Cov.mat[ind,2] = cov.new
        summary.Cov.mat[ind,9] = mean(na.omit(cover.union05[,4]))
        summary.Leng.mat[ind,c(1,3,4,5,6,7,8)] = colMeans(Leng.mat.whole)[c(1,2,3,4,5,7,8)]
        summary.Leng.mat[ind, 2] = len.new
        summary.Leng.mat[ind,9] = mean(na.omit(length.union05[,4]))
        summary.rule.mat[ind,] = colMeans(Rule.mat)
      }
    }
  }
  return(list(summary.Cov.mat=summary.Cov.mat, summary.Leng.mat=summary.Leng.mat,
              summary.rule.mat=summary.rule.mat))
}
VIO.str.options = c(0.2, 0.4)
n.options = c(500, 1000, 2000, 5000)

## CIIV1 ##
out.1 = CIIV_setting(c("CIIV1"))
out.table.1 = data.frame(matrix(NA, nrow=length(VIO.str.options)*
                                  length(n.options)*2, ncol=2+9+1))
colnames(out.table.1) = c("tau","n","TSLS","BA","TSHT","CIIV","Sear-TSHT","Sear-CIIV",
                          "Samp-TSHT","Samp-CIIV","Check","Union")
out.table.1[,1] = rep(rep(VIO.str.options, rep(length(n.options), length(VIO.str.options))) ,2)
out.table.1[,2] = rep(rep(n.options, length(VIO.str.options)), 2)
out.table.1[,c(3,4,5,6,7,8,9,10,12)] = rbind(out.1$summary.Cov.mat, out.1$summary.Leng.mat)
out.table.1[,11] = c(out.1$summary.rule.mat[,1], rep(NA, length(VIO.str.options)*
                                                       length(n.options)))
library(kableExtra)
kbl(out.table.1, "latex", digits=2, align="c")%>%
  add_header_above(c(" "=2, "Oracle"=2, " "=2, "Searching"=2, "Sampling"=2, " "=2))%>%
  kable_styling(latex_options = c("scale_down"))%>%
  collapse_rows(columns=1)

## CIIV2 ##
out.1 = CIIV_setting(c("CIIV2"))
out.table.1 = data.frame(matrix(NA, nrow=length(VIO.str.options)*
                                  length(n.options)*2, ncol=2+9+1))
colnames(out.table.1) = c("tau","n","TSLS","BA","TSHT","CIIV","Sear-TSHT","Sear-CIIV",
                          "Samp-TSHT","Samp-CIIV","Check","Union")
out.table.1[,1] = rep(rep(VIO.str.options, rep(length(n.options), length(VIO.str.options))) ,2)
out.table.1[,2] = rep(rep(n.options, length(VIO.str.options)), 2)
out.table.1[,c(3,4,5,6,7,8,9,10,12)] = rbind(out.1$summary.Cov.mat, out.1$summary.Leng.mat)
out.table.1[,11] = c(out.1$summary.rule.mat[,1], rep(NA, length(VIO.str.options)*
                                                       length(n.options)))
library(kableExtra)
kbl(out.table.1, "latex", digits=2, align="c")%>%
  add_header_above(c(" "=2, "Oracle"=2, " "=2, "Searching"=2, "Sampling"=2, " "=2))%>%
  kable_styling(latex_options = c("scale_down"))%>%
  collapse_rows(columns=1)

################# Hetero-S2 More taus plots ####################

hetero_S2 <- function(n.options){
  setting = "S2"
  IV.str = 0.5
  VIO.str.options = c(0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5)
  #n.options = c(500)
  sim.round.options = c(1)
  nsim.perround=500/length(sim.round.options)
  
  summary.new = matrix(NA, nrow=length(VIO.str.options)*
                         length(n.options), ncol=2)
  summary.Cov.mat = summary.Leng.mat = matrix(NA, nrow=length(VIO.str.options)*
                                                length(n.options), ncol=12)
  summary.Rule.mat = matrix(NA, nrow=length(VIO.str.options)*
                              length(n.options), ncol=4)
  summary.Point.mat = matrix(NA, nrow=length(VIO.str.options)*
                               length(n.options), ncol=2)
  for(i.VIO.str in 1:length(VIO.str.options)){
    VIO.str = VIO.str.options[i.VIO.str]
    for(i.n in 1:length(n.options)){
      n = n.options[i.n]
      
      CI.mat.whole = matrix(NA, nrow=500, ncol=12)
      Cov.mat.whole = Leng.mat.whole = matrix(NA, nrow=500, ncol=12)
      Rule.mat.whole = matrix(NA, nrow=500, ncol=4)
      # colnames(Cov.mat.whole) = colnames(Leng.mat.whole) = c("oracle", "TSHT", "CIIV", "Sear-TSHT", "Sear-CIIV", "Sear-Comb",
      #                                                        "Samp-TSHT", "Samp-CIIV", "Samp-Comb")
      colnames(Rule.mat.whole) = c("Sear-TSHT","Sear-CIIV","Samp-TSHT", "Samp-CIIV")
      for(i.sim.round in 1:length(sim.round.options)){
        sim.round = sim.round.options[i.sim.round]
        datafile = paste("RDatas/hetero-S2/Hetero-Setting",setting,"-Strength",IV.str,"-Violation",VIO.str,
                         "-n",n,"-SimRound",sim.round,".RData",sep="")
        load(datafile)
        colnames(CI.mat.whole) = colnames(CI.mat)
        CI.mat.whole[((i.sim.round-1)*nsim.perround+1):(i.sim.round*nsim.perround),] = CI.mat
        Cov.mat.whole[((i.sim.round-1)*nsim.perround+1):(i.sim.round*nsim.perround),] = Cov.mat
        Leng.mat.whole[((i.sim.round-1)*nsim.perround+1):(i.sim.round*nsim.perround),] = Leng.mat
        Rule.mat.whole[((i.sim.round-1)*nsim.perround+1):(i.sim.round*nsim.perround),] = Rule.mat
      }
      
      ind = (i.VIO.str-1)*length(n.options) + i.n
      summary.Cov.mat[ind,] = colMeans(Cov.mat.whole)
      summary.Leng.mat[ind,] = colMeans(Leng.mat.whole)
      summary.Rule.mat[ind,] = colMeans(Rule.mat.whole)
      summary.Point.mat[ind,] = c(mean(CI.mat.whole[,2]-CI.mat.whole[,1]), mean(CI.mat.whole[,4]-CI.mat.whole[,3]))
      
      points = (CI.mat.whole[,"TSHT2-L"]+CI.mat.whole[,"TSHT2-U"])/2
      bias = points - 1
      sd = sqrt(var(bias))
      temp = mean(bias) / sd
      cv_alpha = sqrt(qchisq(1-0.05,df=1, ncp=temp^2))
      CI.new = cbind(points - sd*cv_alpha, points + sd*cv_alpha)
      cov.new = mean((CI.new[,1] < 1)*(CI.new[,2] > 1))
      len.new = 2*sd*cv_alpha
      summary.new[ind, ] = c(cov.new, len.new)
    }
  }
  df.cov = cbind(VIO.str.options, summary.Cov.mat[,c(4,6)], summary.new[,1], summary.Cov.mat[,c(7,10)])
  df.len = cbind(VIO.str.options, summary.Leng.mat[,c(4,6)], summary.new[,2], summary.Leng.mat[,c(7,10)])
  colnames(df.cov) = colnames(df.len) = c("VIO","TSHT","CIIV","Oracle-BA","Searching","Sampling")
  df.cov = data.frame(df.cov)
  df.len = data.frame(df.len)
  return(list(df.cov=df.cov, df.len=df.len))
}

library(reshape2)
## n=500
out.1 = hetero_S2(c(500))
df.cov.1 = melt(out.1$df.cov, id="VIO")
df.len.1 = melt(out.1$df.len, id="VIO")
## n=2000
out.2 = hetero_S2(c(2000))
df.cov.2 = melt(out.2$df.cov, id="VIO")
df.len.2 = melt(out.2$df.len, id="VIO")
## n=5000
out.3 = hetero_S2(c(5000))
df.cov.3 = melt(out.3$df.cov, id="VIO")
df.len.3 = melt(out.3$df.len, id="VIO")

## plot
library(reshape2)
library(ggplot2)
library(ggpubr)
colnames(df.cov.1)[2] = colnames(df.cov.2)[2] = colnames(df.cov.3)[2] = 
  colnames(df.len.1)[2] = colnames(df.len.2)[2] = colnames(df.len.3)[2]= "Method"
p1 = ggplot(df.cov.1, aes(x = factor(VIO), y=value, group=Method))+
  geom_line(aes(color=Method), size=0.9)+
  geom_point(aes(shape=Method, color=Method), size=2)+
  labs(x="Violation Stength", y="Coverage")
p2 = ggplot(df.len.1, aes(x = factor(VIO), y=value, group=Method))+
  geom_line(aes(color=Method), size=0.9)+
  geom_point(aes(shape=Method, color=Method), size=2)+
  labs(x="Violation Stength", y="Length")
p3 = ggplot(df.cov.2, aes(x = factor(VIO), y=value, group=Method))+
  geom_line(aes(color=Method), size=0.9)+
  geom_point(aes(shape=Method, color=Method), size=2)+
  labs(x="Violation Stength", y="Coverage")
p4 = ggplot(df.len.2, aes(x = factor(VIO), y=value, group=Method))+
  geom_line(aes(color=Method), size=0.9)+
  geom_point(aes(shape=Method, color=Method), size=2)+
  labs(x="Violation Stength", y="Length")
p5 = ggplot(df.cov.3, aes(x = factor(VIO), y=value, group=Method))+
  geom_line(aes(color=Method), size=0.9)+
  geom_point(aes(shape=Method, color=Method), size=2)+
  labs(x="Violation Stength", y="Coverage")
p6 = ggplot(df.len.3, aes(x = factor(VIO), y=value, group=Method))+
  geom_line(aes(color=Method), size=0.9)+
  geom_point(aes(shape=Method, color=Method), size=2)+
  labs(x="Violation Stength", y="Length")

p12 = ggarrange(p1,p2, ncol=2, nrow=1, legend="none")
p12 = annotate_figure(p12, top = text_grob("n=500", 
                                           color = "black", face = "bold", size = 10))
p34 = ggarrange(p3,p4, ncol=2, nrow=1, legend="none")
p34 = annotate_figure(p34, top = text_grob("n=2000", 
                                           color = "black", face = "bold", size = 10))
p56 = ggarrange(p5,p6, ncol=2, nrow=1, legend="bottom", common.legend = TRUE)
p56 = annotate_figure(p56, top = text_grob("n=5000", 
                                           color = "black", face = "bold", size = 10))
ggarrange(p12,p34,p56, nrow=3, heights = c(0.8, 0.8, 1))
######################## High d ###########################
### pz=200
n.options = c(200, 300, 1000, 2500)
VIO.str.options = c(0.5, 1, 2)
sim.round.options = seq(1,25)
nsim.per = 500/length(sim.round.options)

summary.new = matrix(NA, nrow=length(n.options)*length(VIO.str.options), ncol=2)
summary.CI.mat = matrix(NA, nrow=length(n.options)*length(VIO.str.options), ncol=4)
colnames(summary.CI.mat) = c("oracle-L","oracle-U","TSHT-L","TSHT-U")
summary.Cov.mat = summary.Leng.mat = matrix(NA, nrow=length(n.options)*length(VIO.str.options), ncol=4)
colnames(summary.Cov.mat) =  c("oracle", "TSHT", "Sear", "Samp")
summary.point.mat = matrix(NA, nrow=length(n.options)*length(VIO.str.options), ncol=2)
colnames(summary.point.mat) = c("oracle", "TSHT")

for(i.VIO.str in 1:length(VIO.str.options)){
  VIO.str = VIO.str.options[i.VIO.str]
  for(i.n in 1:length(n.options)){
    n = n.options[i.n]
    
    whole_Cov.mat = whole_Leng.mat = whole_CI.mat = matrix(NA, nrow=500, ncol=4)
    colnames(whole_Cov.mat) = colnames(whole_Leng.mat) = c("oracle", "TSHT", "Sear", "Samp")
    colnames(whole_CI.mat) = c("oracle-L","oracle-U","TSHT-L","TSHT-U")
    for(i.sim.round in 1:length(sim.round.options)){
      sim.round = sim.round.options[i.sim.round]
      datafile = paste("RDatas/pz200/Highd-n", n, "-pz200-Violation", VIO.str, "-SimRound", sim.round, ".RData", sep="")
      load(datafile)
      whole_Cov.mat[((sim.round-1)*nsim.per+1):(sim.round*nsim.per), ] = Cov.mat
      whole_Leng.mat[((sim.round-1)*nsim.per+1):(sim.round*nsim.per), ] = Leng.mat
      whole_CI.mat[((sim.round-1)*nsim.per+1):(sim.round*nsim.per), ] = CI.mat
    }
    ind = (i.VIO.str-1)*length(n.options)+i.n
    points = (whole_CI.mat[,"TSHT-L"] + whole_CI.mat[,"TSHT-U"])/2
    bias = points - 1
    sd = sqrt(var(bias))
    temp = mean(bias) / sd
    cv_alpha = sqrt(qchisq(1-0.05,df=1, ncp=temp^2))
    CI.new = cbind(points - sd*cv_alpha, points + sd*cv_alpha)
    cov.new = mean((CI.new[,1] < 1)*(CI.new[,2] > 1))
    len.new = 2*sd*cv_alpha
    summary.new[ind, ] = c(cov.new, len.new)
    
    summary.Cov.mat[ind, ] = colMeans(whole_Cov.mat)
    summary.Leng.mat[ind, ] = colMeans(whole_Leng.mat)
    summary.point.mat[ind, ] = c(mean((whole_CI.mat[,2] + whole_CI.mat[,1])/2), mean((whole_CI.mat[,4] + whole_CI.mat[,3])/2))
  }
}
out.all.2 = data.frame(matrix(NA, nrow=length(n.options)*length(VIO.str.options), ncol=2+3+4*2))
colnames(out.all.2) = c("tau", "n",rep(c("cov","len"),2), c("bias","cov","len"), rep(c("cov","len"),2))
out.all.2[,1] = rep(VIO.str.options, rep(length(n.options),length(VIO.str.options)))
out.all.2[,2] = rep(n.options, length(VIO.str.options))
out.all.2[,c(3,8,10,12)] = summary.Cov.mat[,c(1,2,3,4)] #cov
out.all.2[,c(4,9,11,13)] = summary.Leng.mat[,c(1,2,3,4)] #len
out.all.2[,7] = summary.point.mat[,2] - 1
out.all.2[,c(5,6)] = summary.new

library(kableExtra)
out.all.2 = round(out.all.2, 2)
kbl(out.all.2, "latex", align="c")%>%
  add_header_above(c(" "=2,"Oracle-TSLS"=2,"Oracle-BA"=2,"TSHT"=3,"Searching"=2,"Sampling"=2))%>%
  kable_styling(latex_options = c("scale_down"))%>%
  collapse_rows(columns=1:2)

### pz=100
n.options = c(200, 300, 1000, 2500)
VIO.str.options = c(0.5, 1, 2)
sim.round.options = seq(1,25)
nsim.per = 500/length(sim.round.options)

summary.new = matrix(NA, nrow=length(n.options)*length(VIO.str.options), ncol=2)
summary.CI.mat = matrix(NA, nrow=length(n.options)*length(VIO.str.options), ncol=4)
colnames(summary.CI.mat) = c("oracle-L","oracle-U","TSHT-L","TSHT-U")
summary.Cov.mat = summary.Leng.mat = matrix(NA, nrow=length(n.options)*length(VIO.str.options), ncol=4)
colnames(summary.Cov.mat) =  c("oracle", "TSHT", "Sear", "Samp")
summary.point.mat = matrix(NA, nrow=length(n.options)*length(VIO.str.options), ncol=2)
colnames(summary.point.mat) = c("oracle", "TSHT")

for(i.VIO.str in 1:length(VIO.str.options)){
  VIO.str = VIO.str.options[i.VIO.str]
  for(i.n in 1:length(n.options)){
    n = n.options[i.n]
    
    whole_Cov.mat = whole_Leng.mat = whole_CI.mat = matrix(NA, nrow=500, ncol=4)
    colnames(whole_Cov.mat) = colnames(whole_Leng.mat) = c("oracle", "TSHT", "Sear", "Samp")
    colnames(whole_CI.mat) = c("oracle-L","oracle-U","TSHT-L","TSHT-U")
    for(i.sim.round in 1:length(sim.round.options)){
      sim.round = sim.round.options[i.sim.round]
      datafile = paste("RDatas/pz100/Highd-n", n, "-Violation", VIO.str, "-SimRound", sim.round, ".RData", sep="")
      load(datafile)
      whole_Cov.mat[((sim.round-1)*nsim.per+1):(sim.round*nsim.per), ] = Cov.mat
      whole_Leng.mat[((sim.round-1)*nsim.per+1):(sim.round*nsim.per), ] = Leng.mat
      whole_CI.mat[((sim.round-1)*nsim.per+1):(sim.round*nsim.per), ] = CI.mat
    }
    ind = (i.VIO.str-1)*length(n.options)+i.n
    points = (whole_CI.mat[,"TSHT-L"] + whole_CI.mat[,"TSHT-U"])/2
    bias = points - 1
    sd = sqrt(var(bias))
    temp = mean(bias) / sd
    cv_alpha = sqrt(qchisq(1-0.05,df=1, ncp=temp^2))
    CI.new = cbind(points - sd*cv_alpha, points + sd*cv_alpha)
    cov.new = mean((CI.new[,1] < 1)*(CI.new[,2] > 1))
    len.new = 2*sd*cv_alpha
    summary.new[ind, ] = c(cov.new, len.new)
    
    summary.Cov.mat[ind, ] = colMeans(whole_Cov.mat)
    summary.Leng.mat[ind, ] = colMeans(whole_Leng.mat)
    summary.point.mat[ind, ] = c(mean((whole_CI.mat[,2] + whole_CI.mat[,1])/2), mean((whole_CI.mat[,4] + whole_CI.mat[,3])/2))
  }
}

out.all.1 = data.frame(matrix(NA, nrow=length(n.options)*length(VIO.str.options), ncol=2+3+4*2))
colnames(out.all.1) = c("tau", "n",rep(c("cov","len"),2), c("bias","cov","len"), rep(c("cov","len"),2))
out.all.1[,1] = rep(VIO.str.options, rep(length(n.options),length(VIO.str.options)))
out.all.1[,2] = rep(n.options, length(VIO.str.options))
out.all.1[,c(3,8,10,12)] = summary.Cov.mat[,c(1,2,3,4)] #cov
out.all.1[,c(4,9,11,13)] = summary.Leng.mat[,c(1,2,3,4)] #len
out.all.1[,7] = summary.point.mat[,2] - 1
out.all.1[,c(5,6)] = summary.new

library(kableExtra)
out.all.1 = round(out.all.1, 2)
kbl(out.all.1, "latex",digits=2, align="c")%>%
  add_header_above(c(" "=2,"Oracle-TSLS"=2,"Oracle-BA"=2,"TSHT"=3,"Searching"=2,"Sampling"=2))%>%
  kable_styling(latex_options = c("scale_down"))%>%
  collapse_rows(columns=1:2)

################# Homo-S2 More taus plots ####################
homo_taus <- function(n.options){
  setting = "S2"
  IV.str = 0.5
  VIO.str.options = c(0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5)
  #n.options = c(500)
  sim.round.options = c(1,2,3,4,5)
  nsim.perround=500/length(sim.round.options)
  
  summary.new = matrix(NA, nrow=length(VIO.str.options)*
                         length(n.options), ncol=2)
  summary.Cov.mat = summary.Leng.mat = matrix(NA, nrow=length(VIO.str.options)*
                                                length(n.options), ncol=9)
  summary.cover.union.mat = summary.length.union.mat = matrix(NA, nrow=length(VIO.str.options)*
                                                                length(n.options), ncol=2)
  summary.Rule.mat = matrix(NA, nrow=length(VIO.str.options)*
                              length(n.options), ncol=4)
  summary.Point.mat = matrix(NA, nrow=length(VIO.str.options)*
                               length(n.options), ncol=2)
  colnames(summary.Cov.mat) = c("oracle", "TSHT", "CIIV",
                                "Sear-TSHT", "Sear-CIIV", "Sear-Comb",
                                "Samp-TSHT", "Samp-CIIV", "Samp-Comb")
  method=c("AR","CLR","TSLS","SarganTSLS","SarganCLR")
  
  
  for(i.VIO.str in 1:length(VIO.str.options)){
    VIO.str = VIO.str.options[i.VIO.str]
    for(i.n in 1:length(n.options)){
      n = n.options[i.n]
      
      CI.mat.whole = matrix(NA, nrow=500, ncol=4)
      Cov.mat.whole = Leng.mat.whole = matrix(NA, nrow=500, ncol=9)
      Rule.mat.whole = matrix(NA, nrow=500, ncol=4)
      colnames(Cov.mat.whole) = colnames(Leng.mat.whole) = c("oracle", "TSHT", "CIIV", "Sear-TSHT", "Sear-CIIV", "Sear-Comb",
                                                             "Samp-TSHT", "Samp-CIIV", "Samp-Comb")
      colnames(Rule.mat.whole) = c("Sear-TSHT","Sear-CIIV","Samp-TSHT", "Samp-CIIV")
      length.union05.whole<-matrix(NA, nrow=500, ncol=10)
      cover.union05.whole<-matrix(NA, nrow=500, ncol=10)
      colnames(length.union05.whole)=c(method,method)
      colnames(cover.union05.whole)=c(method,method)
      for(i.sim.round in 1:length(sim.round.options)){
        sim.round = sim.round.options[i.sim.round]
        datafile = paste("RDatas/homo_case_S2/Homo-Setting",setting,"-Strength",IV.str,"-Violation",VIO.str,
                         "-n",n,"-SimRound",sim.round,".RData",sep="")
        load(datafile)
        colnames(CI.mat.whole) = colnames(CI.mat)
        CI.mat.whole[((i.sim.round-1)*nsim.perround+1):(i.sim.round*nsim.perround),] = CI.mat
        Cov.mat.whole[((i.sim.round-1)*nsim.perround+1):(i.sim.round*nsim.perround),] = Cov.mat
        Leng.mat.whole[((i.sim.round-1)*nsim.perround+1):(i.sim.round*nsim.perround),] = Leng.mat
        Rule.mat.whole[((i.sim.round-1)*nsim.perround+1):(i.sim.round*nsim.perround),] = Rule.mat
        length.union05.whole[((i.sim.round-1)*nsim.perround+1):(i.sim.round*nsim.perround),] = length.union05
        cover.union05.whole[((i.sim.round-1)*nsim.perround+1):(i.sim.round*nsim.perround),] = cover.union05
      }
      
      ind = (i.VIO.str-1)*length(n.options) + i.n
      summary.Cov.mat[ind,] = colMeans(Cov.mat.whole)
      summary.Leng.mat[ind,] = colMeans(Leng.mat.whole)
      summary.Rule.mat[ind,] = colMeans(Rule.mat.whole)
      summary.cover.union.mat[ind,] = colMeans(na.omit(cover.union05.whole[,c(4,9)]))
      summary.length.union.mat[ind,] = colMeans(na.omit(length.union05.whole[,c(4,9)]))
      summary.Point.mat[ind,] = c(mean(CI.mat.whole[,2]-CI.mat.whole[,1]), mean(CI.mat.whole[,4]-CI.mat.whole[,3]))
      
      points = (CI.mat.whole[,"TSHT-L"]+CI.mat.whole[,"TSHT-U"])/2
      bias = points - 1
      sd = sqrt(var(bias))
      temp = mean(bias) / sd
      cv_alpha = sqrt(qchisq(1-0.05,df=1, ncp=temp^2))
      CI.new = cbind(points - sd*cv_alpha, points + sd*cv_alpha)
      cov.new = mean((CI.new[,1] < 1)*(CI.new[,2] > 1))
      len.new = 2*sd*cv_alpha
      summary.new[ind, ] = c(cov.new, len.new)
    }
  }
  df.cov = cbind(VIO.str.options, summary.Cov.mat[,c(2,3)], summary.new[,1], summary.Cov.mat[,c(4,8)])
  df.len = cbind(VIO.str.options, summary.Leng.mat[,c(2,3)], summary.new[,2], summary.Leng.mat[,c(4,8)])
  colnames(df.cov) = colnames(df.len) = c("VIO","TSHT","CIIV","Oracle-BA","Searching","Sampling")
  df.cov = data.frame(df.cov)
  df.len = data.frame(df.len)
  return(list(df.cov=df.cov, df.len=df.len))
}
library(reshape2)
out.1 = homo_taus(c(500))
df.cov.1 = melt(out.1$df.cov, id="VIO")
df.len.1 = melt(out.1$df.len, id="VIO")

out.2 = homo_taus(c(2000))
df.cov.2 = melt(out.2$df.cov, id="VIO")
df.len.2 = melt(out.2$df.len, id="VIO")

out.3 = homo_taus(c(5000))
df.cov.3 = melt(out.3$df.cov, id="VIO")
df.len.3 = melt(out.3$df.len, id="VIO")

## plot
library(reshape2)
library(ggplot2)
library(ggpubr)
colnames(df.cov.1)[2] = colnames(df.cov.2)[2] = colnames(df.cov.3)[2] =  
  colnames(df.len.1)[2] = colnames(df.len.2)[2] = colnames(df.len.3)[2] = "Method"
p1 = ggplot(df.cov.1, aes(x = factor(VIO), y=value, group=Method))+
  geom_line(aes(color=Method), size=0.9)+
  geom_point(aes(shape=Method, color=Method), size=2)+
  scale_y_continuous(breaks = c(0.4,0.6,0.8,1.0), limits=c(0.4, 1.0))+
  labs(x="Violation Stength", y="Coverage")
p2 = ggplot(df.len.1, aes(x = factor(VIO), y=value, group=Method))+
  geom_line(aes(color=Method), size=0.9)+
  geom_point(aes(shape=Method, color=Method), size=2)+
  labs(x="Violation Stength", y="Length")
p3 = ggplot(df.cov.2, aes(x = factor(VIO), y=value, group=Method))+
  geom_line(aes(color=Method), size=0.9)+
  scale_y_continuous(breaks=c(0.4,0.6,0.8,1.0), limits=c(0.4, 1.0))+
  geom_point(aes(shape=Method, color=Method), size=2)+
  labs(x="Violation Stength", y="Coverage")
p4 = ggplot(df.len.2, aes(x = factor(VIO), y=value, group=Method))+
  geom_line(aes(color=Method), size=0.9)+
  geom_point(aes(shape=Method, color=Method), size=2)+
  labs(x="Violation Stength", y="Length")
p5 = ggplot(df.cov.3, aes(x = factor(VIO), y=value, group=Method))+
  geom_line(aes(color=Method), size=0.9)+
  scale_y_continuous(breaks=c(0.4,0.6,0.8,1.0))+
  geom_point(aes(shape=Method, color=Method), size=2)+
  labs(x="Violation Stength", y="Coverage")
p6 = ggplot(df.len.3, aes(x = factor(VIO), y=value, group=Method))+
  geom_line(aes(color=Method), size=0.9)+
  geom_point(aes(shape=Method, color=Method), size=2)+
  labs(x="Violation Stength", y="Length")

#ggarrange(p1,p2,p3,p4,p5,p6, ncol=2, nrow=3, common.legend = TRUE, legend="bottom")

p12 = ggarrange(p1,p2, ncol=2, nrow=1, legend="none")
p12 = annotate_figure(p12, top = text_grob("n=500", 
                                           color = "black", face = "bold", size = 10))
p34 = ggarrange(p3,p4, ncol=2, nrow=1, legend="none")
p34 = annotate_figure(p34, top = text_grob("n=2000", 
                                           color = "black", face = "bold", size = 10))
p56 = ggarrange(p5,p6, ncol=2, nrow=1, legend="bottom", common.legend = TRUE)
p56 = annotate_figure(p56, top = text_grob("n=5000", 
                                           color = "black", face = "bold", size = 10))
ggarrange(p12,p34,p56, nrow=3, heights = c(0.8, 0.8, 1))
################# Homo S2 table-1 ####################
setting = "S2"
IV.str = 0.5
VIO.str.options = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
n.options = c(500, 2000)
sim.round.options = c(1,2,3,4,5)
nsim.perround=500/length(sim.round.options)

summary.new = matrix(NA, nrow=length(VIO.str.options)*
                       length(n.options), ncol=2)
summary.Cov.mat = summary.Leng.mat = matrix(NA, nrow=length(VIO.str.options)*
                                              length(n.options), ncol=9)
summary.cover.union.mat = summary.length.union.mat = matrix(NA, nrow=length(VIO.str.options)*
                                                              length(n.options), ncol=2)
summary.Point.mat = matrix(NA, nrow=length(VIO.str.options)*
                             length(n.options), ncol=2)
colnames(summary.Cov.mat) = c("oracle", "TSHT", "CIIV",
                              "Sear-TSHT", "Sear-CIIV", "Sear-Comb",
                              "Samp-TSHT", "Samp-CIIV", "Samp-Comb")
method=c("AR","CLR","TSLS","SarganTSLS","SarganCLR")

for(i.VIO.str in 1:length(VIO.str.options)){
  VIO.str = VIO.str.options[i.VIO.str]
  for(i.n in 1:length(n.options)){
    n = n.options[i.n]
    
    CI.mat.whole = matrix(NA, nrow=500, ncol=4)
    Cov.mat.whole = Leng.mat.whole = matrix(NA, nrow=500, ncol=9)
    colnames(Cov.mat.whole) = colnames(Leng.mat.whole) = c("oracle", "TSHT", "CIIV", "Sear-TSHT", "Sear-CIIV", "Sear-Comb",
                                                           "Samp-TSHT", "Samp-CIIV", "Samp-Comb")
    length.union05.whole<-matrix(NA, nrow=500, ncol=10)
    cover.union05.whole<-matrix(NA, nrow=500, ncol=10)
    colnames(length.union05.whole)=c(method,method)
    colnames(cover.union05.whole)=c(method,method)
    for(i.sim.round in 1:length(sim.round.options)){
      sim.round = sim.round.options[i.sim.round]
      datafile = paste("RDatas/homo_case_S2/Homo-Setting",setting,"-Strength",IV.str,"-Violation",VIO.str,
                       "-n",n,"-SimRound",sim.round,".RData",sep="")
      load(datafile)
      CI.mat.whole[((i.sim.round-1)*nsim.perround+1):(i.sim.round*nsim.perround),] = CI.mat
      Cov.mat.whole[((i.sim.round-1)*nsim.perround+1):(i.sim.round*nsim.perround),] = Cov.mat
      Leng.mat.whole[((i.sim.round-1)*nsim.perround+1):(i.sim.round*nsim.perround),] = Leng.mat
      length.union05.whole[((i.sim.round-1)*nsim.perround+1):(i.sim.round*nsim.perround),] = length.union05
      cover.union05.whole[((i.sim.round-1)*nsim.perround+1):(i.sim.round*nsim.perround),] = cover.union05
    }
    ind = (i.VIO.str-1)*length(n.options) + i.n
    
    # bias-aware
    colnames(CI.mat.whole) = colnames(CI.mat)
    points = (CI.mat.whole[,"TSHT-L"]+CI.mat.whole[,"TSHT-U"])/2
    bias = points - 1
    sd = sqrt(var(bias))
    temp = mean(bias) / sd
    cv_alpha = sqrt(qchisq(1-0.05,df=1, ncp=temp^2))
    CI.new = cbind(points - sd*cv_alpha, points + sd*cv_alpha)
    cov.new = mean((CI.new[,1] < 1)*(CI.new[,2] > 1))
    len.new = 2*sd*cv_alpha
    summary.new[ind, ] = c(cov.new, len.new)
    summary.Cov.mat[ind,] = colMeans(Cov.mat.whole)
    summary.Leng.mat[ind,] = colMeans(Leng.mat.whole)
    summary.cover.union.mat[ind,] = colMeans(na.omit(cover.union05.whole[,c(4,9)]))
    summary.length.union.mat[ind,] = colMeans(na.omit(length.union05.whole[,c(4,9)]))
    summary.Point.mat[ind,] = c(mean(CI.mat.whole[,2]-CI.mat.whole[,1]), mean(CI.mat.whole[,4]-CI.mat.whole[,3]))
  }
}
out.table = data.frame(matrix(NA, nrow=length(VIO.str.options)*length(n.options), ncol=2+1+5*2))
colnames(out.table) = c("tao","n","T",rep(c("cov","len"),5))
out.table[,1] = rep(VIO.str.options, rep(length(n.options), length(VIO.str.options)))
out.table[,2] = rep(n.options, length(VIO.str.options))
out.table[,c(4,6,12)] = summary.Cov.mat[,c(2,3,7)] #cov:TSHT, CIIV, Samp
out.table[,c(5,7,13)] = summary.Leng.mat[,c(2,3,7)] #len:TSHT, CIIV, Samp
out.table[,c(8,9)] = summary.new[,c(1,2)] #bias-aware cov&len
out.table[,10] = summary.cover.union.mat[,1] #cov: union pz-1
out.table[,11] = summary.length.union.mat[,1] #len: union pz-1

load("RDatas/S2_thres.RData")
out.table[,3] = thres_table[,"n-maxV"]
out.table = round(out.table, 3)
library(kableExtra)
kbl(out.table, "latex", digits=2, align="c")%>%
  add_header_above(c(" "=3, "TSHT"=2,"CIIV"=2,"Bias-aware"=2,"Union"=2,"Sampling"=2))%>%
  kable_styling(latex_options = c("scale_down"))%>%
  collapse_rows(columns=1:2)

################# Homo S2 Points ###################
library(ggplot2)
library(ggpubr)
### n=500
CI.mat.whole = matrix(NA, nrow=500, ncol=4)
colnames(CI.mat.whole) = c("TSHT-L","TSHT-U","CIIV-L","CIIV-U")
load("RDatas/CIIV_CI.RData")
CI.mat.whole = (CI.mat.list$`n=500VIO=0.2`)
df = data.frame(TSHT = rowMeans(CI.mat.whole[,1:2]), CIIV=rowMeans(CI.mat.whole[,3:4]))
colMeans(df) # 1.022044 1.028433 

p<-ggplot(df,aes(x=TSHT)) + geom_histogram(binwidth=0.016,color="black", fill="white")
p1<-p+ geom_vline(aes(xintercept=1),color="black", linetype="solid", size=0.5)+
  geom_vline(aes(xintercept=mean(TSHT)),color="black", linetype="dashed", size=0.5)+
  xlim(c(0.9, 1.2))
p<-ggplot(df,aes(x=CIIV)) + geom_histogram(binwidth=0.016,color="black", fill="white")
p2<-p+ geom_vline(aes(xintercept=1),color="black", linetype="solid", size=0.5)+
  geom_vline(aes(xintercept=mean(CIIV)),color="black", linetype="dashed", size=0.5)+
  xlim(c(0.9, 1.2))

### n=2000
CI.mat.whole = matrix(NA, nrow=500, ncol=4)
colnames(CI.mat.whole) = c("TSHT-L","TSHT-U","CIIV-L","CIIV-U")
load("RDatas/CIIV_CI.RData")
CI.mat.whole = (CI.mat.list$`n=2000VIO=0.2`)
df = data.frame(TSHT = rowMeans(CI.mat.whole[,1:2]), CIIV=rowMeans(CI.mat.whole[,3:4]))
colMeans(df) # 1.029843 1.024400 

p<-ggplot(df,aes(x=TSHT)) + geom_histogram(binwidth=0.016,color="black", fill="white")
p3<-p+ geom_vline(aes(xintercept=1),color="black", linetype="solid", size=0.5)+
  geom_vline(aes(xintercept=mean(TSHT)),color="black", linetype="dashed", size=0.5)+
  xlim(c(0.9, 1.2))
p<-ggplot(df,aes(x=CIIV)) + geom_histogram(binwidth=0.016,color="black", fill="white")
p4<-p+ geom_vline(aes(xintercept=1),color="black", linetype="solid", size=0.5)+
  geom_vline(aes(xintercept=mean(CIIV)),color="black", linetype="dashed", size=0.5)+
  xlim(c(0.9, 1.2))

ggarrange(p1, p2, p3,p4, nrow = 2, ncol=2)

######################## Illustrate Sampling M ###########################
VIO.str = 0.4
n.options = c(500, 1000, 2000)
setting.options = c("S1","S2")
LUa.options = c(1,2,3)
IV.str = 0.5

summary.Cov.mat = summary.Leng.mat = summary.Diff.mat = matrix(NA, nrow=length(n.options)*length(setting.options)*length(LUa.options), ncol=4)
colnames(summary.Cov.mat) = colnames(summary.Leng.mat) = c(500,1000,1500,2000)
for(i.setting in 1:length(setting.options)){
  for(i.n in 1:length(n.options)){
    for(i.LUa in 1:length(LUa.options)){
      setting = setting.options[i.setting]
      n = n.options[i.n]
      LUa = LUa.options[i.LUa]
      ind = (i.setting-1)*length(n.options)*length(LUa.options) +
        (i.n-1)*length(LUa.options) + i.LUa
      Cov.mat.Samp.whole = Leng.mat.Samp.whole = Point.mat.Samp.whole = matrix(NA, nrow=500, ncol=4)
      for(sim.round in 1:10){
        datafile = paste("RDatas/simu_illustrate_samplingM/Illustrate_SampM-Homo-Setting",setting,"-Strength",IV.str,"-Violation",VIO.str,
                         "-n",n,"-LUa",LUa,"-SimRound",sim.round,".RData",sep="")
        load(datafile)
        Cov.mat.Samp.whole[((sim.round-1)*50+1):(sim.round*50),] = Cov.mat.Samp
        Leng.mat.Samp.whole[((sim.round-1)*50+1):(sim.round*50),] = Leng.mat.Samp
        Point.mat.Samp.whole[((sim.round-1)*50+1):(sim.round*50),] = Point.mat.Samp
      }
      summary.Cov.mat[ind, ] = colMeans(Cov.mat.Samp.whole)
      summary.Leng.mat[ind, ] = colMeans(Leng.mat.Samp.whole)
      if(i.LUa==1){
        Benchmark.CI.mat = cbind((Point.mat.Samp.whole[,2]-Leng.mat.Samp.whole[,2]/2),
                                 (Point.mat.Samp.whole[,2]-Leng.mat.Samp.whole[,2]/2))
      }
      for(i.M in 1:4){
        compare.CI.mat = cbind((Point.mat.Samp.whole[,i.M]-Leng.mat.Samp.whole[,i.M]/2),
                               (Point.mat.Samp.whole[,i.M]-Leng.mat.Samp.whole[,i.M]/2))
        upper_diff = mean(abs(compare.CI.mat[,2] - Benchmark.CI.mat[,2]))
        lower_diff = mean(abs(compare.CI.mat[,1] - Benchmark.CI.mat[,1]))
        summary.Diff.mat[ind, i.M] = upper_diff+lower_diff
      }
    }
  }
}
### With coverage
out.table = matrix(NA, nrow=length(setting.options)*length(n.options)*length(LUa.options), ncol=3+4*2)
out.table[,1] = rep(1:2, rep(length(n.options)*length(LUa.options),length(setting.options)))
out.table[,2] = rep(rep(n.options,rep(length(LUa.options), length(n.options))), length(setting.options))
out.table[,3] = rep(LUa.options, length(setting.options)*length(n.options))
out.table[,c(4,6,8,10)] = summary.Cov.mat
out.table[,c(5,7,9,11)] = summary.Leng.mat
# out.table = out.table + 0.0001
colnames(out.table) = c("set","n","LUa",rep(c("cov","len"),4))
library(kableExtra)
kbl(out.table, "latex", digits=3, align="c")%>%
  add_header_above(c(" "=3, "M=500"=2, "M=1000"=2, "M=1500"=2,"M=2000"=2))%>%
  kable_styling(latex_options = c("scale_down"))%>%
  collapse_rows(columns=1:3)

# #### Without Coverage
# out.table = matrix(NA, nrow=length(setting.options)*length(n.options)*length(LUa.options), ncol=3+4*2)
# out.table[,1] = rep(1:2, rep(length(n.options)*length(LUa.options),length(setting.options)))
# out.table[,2] = rep(rep(n.options,rep(length(LUa.options), length(n.options))), length(setting.options))
# out.table[,3] = rep(LUa.options, length(setting.options)*length(n.options))
# out.table[,c(4,6,8,10)] = summary.Leng.mat
# out.table[,c(5,7,9,11)] = summary.Diff.mat
# colnames(out.table) = c("set","n","LUa",rep(c("len","diff"),4))
# library(kableExtra)
# kbl(out.table, "latex", digits=3, align="c")%>%
#   add_header_above(c(" "=3, "M=500"=2, "M=1000"=2, "M=1500"=2,"M=2000"=2))%>%
#   kable_styling(latex_options = c("scale_down"))%>%
#   collapse_rows(columns=1:3)

######################## Illustrate LUa ###########################
LUa_illustrate <- function(setting, n){
  sim.round.options = seq(1,20)
  VIO.str = 0.4
  IV.str = 0.5 
  Cov.mat.Sear.whole = Leng.mat.Sear.whole = Point.mat.Sear.whole = 
    Cov.mat.Samp.whole = Leng.mat.Samp.whole = Point.mat.Samp.whole = 
    matrix(NA, nrow=500, ncol=7)
  for(sim.round in sim.round.options){
    datafile = paste("RDatas/simu_illustrate_LUa/Illustrate_LUa-Homo-Setting",setting,"-Strength",IV.str,"-Violation",VIO.str,
                     "-n",n,"-SimRound",sim.round,".RData",sep="")
    load(datafile)
    Cov.mat.Sear.whole[((sim.round-1)*25+1):(sim.round*25),] = Cov.mat.Sear
    Leng.mat.Sear.whole[((sim.round-1)*25+1):(sim.round*25),] = Leng.mat.Sear
    Point.mat.Sear.whole[((sim.round-1)*25+1):(sim.round*25),] = Point.mat.Sear
    Cov.mat.Samp.whole[((sim.round-1)*25+1):(sim.round*25),] = Cov.mat.Samp
    Leng.mat.Samp.whole[((sim.round-1)*25+1):(sim.round*25),] = Leng.mat.Samp
    Point.mat.Samp.whole[((sim.round-1)*25+1):(sim.round*25),] = Point.mat.Samp
  }
  
  out.table = matrix(NA, nrow=7, ncol=7)
  colnames(out.table) = c("Method","Sear-cov","Sear-len","Sear-diff",
                          "Samp-cov","Samp-len","Samp-diff")
  out.table[,1] = seq(7)
  out.table[,2] = colMeans(Cov.mat.Sear.whole)
  out.table[,3] = colMeans(Leng.mat.Sear.whole)
  for(i in 1:7){
    upper_diff = mean(abs((Point.mat.Sear.whole[,i]+Leng.mat.Sear.whole[,i]/2) - 
                            (Point.mat.Sear.whole[,1]+Leng.mat.Sear.whole[,1]/2)))
    lower_diff = mean(abs((Point.mat.Sear.whole[,i]-Leng.mat.Sear.whole[,i]/2) - 
                            (Point.mat.Sear.whole[,1]-Leng.mat.Sear.whole[,1]/2)))
    out.table[i,4] = upper_diff + lower_diff
  }
  out.table[,5] = colMeans(Cov.mat.Samp.whole)
  out.table[,6] = colMeans(Leng.mat.Samp.whole)
  for(i in 1:7){
    upper_diff = mean(abs((Point.mat.Samp.whole[,i]+Leng.mat.Samp.whole[,i]/2) - 
                            (Point.mat.Samp.whole[,1]+Leng.mat.Samp.whole[,1]/2)))
    lower_diff = mean(abs((Point.mat.Samp.whole[,i]-Leng.mat.Samp.whole[,i]/2) - 
                            (Point.mat.Samp.whole[,1]-Leng.mat.Samp.whole[,1]/2)))
    out.table[i,7] = upper_diff + lower_diff
  }
  return(out.table)
}
out.table.1 = LUa_illustrate("S1",500)
out.table.2 = LUa_illustrate("S1",2000)
out.table.3 = LUa_illustrate("S2",500)
out.table.4 = LUa_illustrate("S2",2000)

## table1: only searching
out.table = cbind(1, cbind(out.table.1[,1:4], out.table.2[,2:4]),
                  cbind(out.table.3[,2:4], out.table.4[,2:4]))

colnames(out.table)[c(1,2)] = c("LU","a")
library(kableExtra)
kbl(out.table, "latex", digits=3, align="c")%>%
  kable_styling(latex_options = c("scale_down"))%>%
  add_header_above(c(" "=2, "n=500"=3, "n=2000"=3, "n=500"=3, "n=2000"=3))%>%
  add_header_above(c(" "=2, "S1"=6, "S2"=6))

## table2: only sampling
out.table = cbind(rep(c(1,2),c(7,7)),
                  rbind(cbind(out.table.1[,c(1,5:6)], out.table.2[,5:6]),
                        cbind(out.table.3[,c(1,5:6)], out.table.4[,5:6])))
colnames(out.table)[1] ="Set"
library(kableExtra)
kbl(out.table, "latex", digits=3, align="c")%>%
  kable_styling(latex_options = c("scale_down"))%>%
  add_header_above(c(" "=2, "n=500"=2, "n=2000"=2))%>%
  collapse_rows(columns=1)

######################## Illustrate LUa additional ############################
LUa_illustrate_add <- function(setting, n){
  sim.round.options = seq(1,20)
  VIO.str = 0.4
  IV.str = 0.5 
  Cov.mat.Sear.whole = Leng.mat.Sear.whole = Point.mat.Sear.whole = Time.mat.Sear.whole =
    matrix(NA, nrow=500, ncol=12)
  for(sim.round in sim.round.options){
    datafile = paste("C:/Users/Zhenyu Wang/Dropbox/Zhenyu/Instrumental Variable/code_0506/RDatas/simu_illustrate_LUa_add/Illustrate_LUa-Additional-Homo-Setting"
                     ,setting,"-Strength",IV.str,"-Violation",VIO.str,
                     "-n",n,"-SimRound",sim.round,".RData",sep="")
    load(datafile)
    Cov.mat.Sear.whole[((sim.round-1)*25+1):(sim.round*25),] = Cov.mat.Sear
    Leng.mat.Sear.whole[((sim.round-1)*25+1):(sim.round*25),] = Leng.mat.Sear
    Point.mat.Sear.whole[((sim.round-1)*25+1):(sim.round*25),] = Point.mat.Sear
    Time.mat.Sear.whole[((sim.round-1)*25+1):(sim.round*25),] = Time.mat.Sear
  }
  
  out.table = matrix(NA, nrow=12, ncol=5)
  colnames(out.table) = c("Method","Sear-cov","Sear-len","Sear-diff",
                          "Sear-time")
  out.table[,1] = seq(12)
  out.table[,2] = colMeans(Cov.mat.Sear.whole)
  out.table[,3] = colMeans(Leng.mat.Sear.whole)
  for(i in 1:12){
    upper_diff = mean(abs((Point.mat.Sear.whole[,i]+Leng.mat.Sear.whole[,i]/2) - 
                            (Point.mat.Sear.whole[,1]+Leng.mat.Sear.whole[,1]/2)))
    lower_diff = mean(abs((Point.mat.Sear.whole[,i]-Leng.mat.Sear.whole[,i]/2) - 
                            (Point.mat.Sear.whole[,1]-Leng.mat.Sear.whole[,1]/2)))
    out.table[i,4] = upper_diff + lower_diff
  }
  out.table[,5] = colMeans(Time.mat.Sear.whole)
  return(out.table)
}
out.table.1 = LUa_illustrate_add("S1",500)
out.table.2 = LUa_illustrate_add("S1",2000)
out.table.3 = LUa_illustrate_add("S2",500)
out.table.4 = LUa_illustrate_add("S2",2000)

## table1: only searching
out.table = cbind(1, cbind(out.table.1[,1:4], out.table.2[,2:4]),
                  cbind(out.table.3[,2:4], out.table.4[,2:4]))

colnames(out.table)[c(1,2)] = c("LU","a")
library(kableExtra)
kbl(out.table, "latex", digits=3, align="c")%>%
  kable_styling(latex_options = c("scale_down"))%>%
  add_header_above(c(" "=2, "n=500"=3, "n=2000"=3, "n=500"=3, "n=2000"=3))

######################## Illustrate sampling Lambda ###########################

## tau=0.4 ##
setting.options = c("S1","S2","S3","S4","S5")
n.options = c(500, 1000, 2000)
IV.str = 0.5
VIO.str = 0.4
out.table = matrix(NA, nrow=15, ncol=5*2)
for(i.setting in 1:5){
  for(i.n in 1:3){
    setting = setting.options[i.setting]
    n = n.options[i.n]
    datafile = paste("RDatas/simu_illustrate_samplingLambda/Illustrate_SampLambda-Homo-Setting",setting,"-Strength",IV.str,"-Violation",VIO.str,
                     "-n",n,".RData",sep="")
    load(datafile)
    
    ind = (i.setting-1)*3 + i.n
    out.table[ind, c(1,3,5,7,9)] = colMeans(Cov.mat.Samp)[1:5]
    out.table[ind, c(2,4,6,8,10)] = colMeans(Leng.mat.Samp)[1:5]
  }
}
out.table = out.table+0.0001
col_setting = rep(1:5, rep(length(n.options), length(setting.options)))
col_n = rep(n.options, length(setting.options))
out.table = cbind(col_setting, col_n, out.table)
colnames(out.table) = c("Set","n",rep(c("cov","len"),5))
library(kableExtra)
kbl(out.table, "latex", digits=3, align="c")%>%
  kable_styling(latex_options = c("scale_down"))%>%
  add_header_above(c(" "=2, "prop1"=2, "prop2"=2,"prop3"=2,"prop4"=2,"prop5"=2))%>%
  collapse_rows(columns=1)

## tau=0.2 ##
setting.options = c("S1","S2","S3","S4","S5")
n.options = c(500, 1000, 2000)
IV.str = 0.5
VIO.str = 0.2
out.table = matrix(NA, nrow=15, ncol=5*2)
for(i.setting in 1:5){
  for(i.n in 1:3){
    setting = setting.options[i.setting]
    n = n.options[i.n]
    datafile = paste("RDatas/simu_illustrate_samplingLambda/Illustrate_SampLambda-Homo-Setting",setting,"-Strength",IV.str,"-Violation",VIO.str,
                     "-n",n,".RData",sep="")
    load(datafile)
    
    ind = (i.setting-1)*3 + i.n
    out.table[ind, c(1,3,5,7,9)] = colMeans(Cov.mat.Samp)[1:5]
    out.table[ind, c(2,4,6,8,10)] = colMeans(Leng.mat.Samp)[1:5]
  }
}
out.table = out.table+0.0001
col_setting = rep(1:5, rep(length(n.options), length(setting.options)))
col_n = rep(n.options, length(setting.options))
out.table = cbind(col_setting, col_n, out.table)
colnames(out.table) = c("Set","n",rep(c("cov","len"),5))
library(kableExtra)
kbl(out.table, "latex", digits=3, align="c")%>%
  kable_styling(latex_options = c("scale_down"))%>%
  add_header_above(c(" "=2, "prop1"=2, "prop2"=2,"prop3"=2,"prop4"=2,"prop5"=2))%>%
  collapse_rows(columns=1)

#################### Summary S0 ###########################
## table-1
IV.str.weak.options = c(0.05, 0.075, 0.1)
n.options = c(500,1000,2000)
VIO.str = 0.5
VIO.str = 0.5
IV.str = 0.5
setting = "S0"
Threshold = 1
summary.Cov.mat = summary.Leng.mat = matrix(NA, nrow=length(IV.str.weak.options)*
                                              length(n.options), ncol=4)
colnames(summary.Cov.mat) = colnames(summary.Leng.mat) = c("TSHT", "CIIV", "Sear-TSHT","Samp-TSHT")
summary.include.mat = matrix(NA, nrow=length(IV.str.weak.options)*
                               length(n.options), ncol=3)
colnames(summary.include.mat) = c("V1","V2","V7")
for(i.IV.weak in 1:length(IV.str.weak.options)){
  for(i.n in 1:length(n.options)){
    IV.str.weak = IV.str.weak.options[i.IV.weak]
    n = n.options[i.n]
    datafile = paste("RDatas/simu_illustrate_S0/Homo-Setting", setting, "-weak", IV.str.weak,"-Thres",Threshold, "-Strength", IV.str, "-Violation", VIO.str, "-n", n, ".RData", sep="")
    load(datafile)
    ind = (i.IV.weak-1)*length(n.options) + i.n
    summary.Cov.mat[ind, ] = colMeans(Cov.mat[,c(2,3,4,7)])
    summary.Leng.mat[ind, ] = colMeans(Leng.mat[,c(2,3,4,7)])
    summary.include.mat[ind, 1] = mean(apply(TSHT.SHat.mat, MARGIN = 1, FUN=function(X) 1%in%X))
    summary.include.mat[ind, 2] = mean(apply(TSHT.SHat.mat, MARGIN = 1, FUN=function(X) 2%in%X))
    summary.include.mat[ind, 3] = mean(apply(TSHT.SHat.mat, MARGIN = 1, FUN=function(X) 7%in%X))
  }
}
out.table = matrix(NA, nrow=length(IV.str.weak.options)*
                     length(n.options), ncol=2+2*2+3)
colnames(out.table) = c("gamma_1","n",rep(c("cov","len"),2), "V1", "V2", "V7")
out.table[,1] = rep(IV.str.weak.options, rep(length(n.options), length(IV.str.weak.options)))
out.table[,2] = rep(n.options, length(IV.str.weak.options))
out.table[,c(3,5)] = summary.Cov.mat[,c(1,2)]
out.table[,c(4,6)] = summary.Leng.mat[,c(1,2)]
out.table[,c(7,8,9)] = summary.include.mat

library(kableExtra)
kbl(out.table, "latex", digits=3, align="c")%>%
  add_header_above(c(" "=2, "TSHT"=2,"CIIV"=2,"Include in SHat"=3))%>%
  kable_styling(latex_options = c("scale_down"))%>%
  collapse_rows(columns=1)

## table-2
IV.str.weak.options = c(0.05, 0.075, 0.1)
n.options = c(500, 1000, 2000)
Threshold.options = c(1,2)
VIO.str = 0.5
IV.str = 0.5
setting= "S0"
summary.new = summary.Cov.mat.1 = summary.Leng.mat.1 = 
  summary.Cov.mat.2 = summary.Leng.mat.2 = matrix(NA, nrow=length(IV.str.weak.options)*
                                                    length(n.options), ncol=2)
colnames(summary.Cov.mat.1) = colnames(summary.Leng.mat.1) = 
  colnames(summary.Cov.mat.2) = colnames(summary.Leng.mat.2) = c("Sear-TSHT","Samp-TSHT")
for(i.IV.weak in 1:length(IV.str.weak.options)){
  for(i.n in 1:length(n.options)){
    IV.str.weak = IV.str.weak.options[i.IV.weak]
    n = n.options[i.n]
    for(i.thres in 1:length(Threshold.options)){
      Threshold = Threshold.options[i.thres]
      datafile = paste("RDatas/simu_illustrate_S0/Homo-Setting", setting, "-weak", IV.str.weak,"-Thres",Threshold, "-Strength", IV.str, "-Violation", VIO.str, "-n", n, ".RData", sep="")
      load(datafile)
      ind = (i.IV.weak-1)*length(n.options) + i.n
      if(Threshold==1){
        points = (CI.mat[,"TSHT-L"] + CI.mat[,"TSHT-U"])/2
        bias = points - 1
        sd = sqrt(var(bias))
        temp = mean(bias) / sd
        cv_alpha = sqrt(qchisq(1-0.05,df=1, ncp=temp^2))
        CI.new = cbind(points - sd*cv_alpha, points + sd*cv_alpha)
        cov.new = mean((CI.new[,1] < 1)*(CI.new[,2] > 1))
        len.new = 2*sd*cv_alpha
        summary.new[ind, ] = c(cov.new, len.new)
        summary.Cov.mat.1[ind, ] = colMeans(Cov.mat[,c(4,7)])
        summary.Leng.mat.1[ind, ] = colMeans(Leng.mat[,c(4,7)])
      }else if(Threshold==2){
        summary.Cov.mat.2[ind, ] = colMeans(Cov.mat[,c(4,7)])
        summary.Leng.mat.2[ind, ] = colMeans(Leng.mat[,c(4,7)])
      }
    }
  }
}
out.table = matrix(NA, nrow=length(IV.str.weak.options)*
                     length(n.options), ncol=2+2*2*2)
colnames(out.table) = c("gamma_1","n",rep(c("cov","len"),4))
out.table[,1] = rep(IV.str.weak.options, rep(length(n.options), length(IV.str.weak.options)))
out.table[,2] = rep(n.options, length(IV.str.weak.options))
out.table[,c(3,5)] = summary.Cov.mat.1
out.table[,c(4,6)] = summary.Leng.mat.1
out.table[,c(7,9)] = summary.Cov.mat.2
out.table[,c(8,10)] = summary.Leng.mat.2

kbl(out.table, "latex", digits=3, align="c")%>%
  add_header_above(c(" "=2, "Searching"=2, "Sampling"=2, "Searching"=2, "Sampling"=2))%>%
  add_header_above(c(" "=2, "Threshold-1"=4, "Threshold-2"=4))%>%
  kable_styling(latex_options = c("scale_down"))%>%
  collapse_rows(columns=1)

############### Bootstrap & Filtering ###################
## Bootstrap or not
n.options = c(500, 1000, 2000)
VIO.str.options = c(0.2, 0.4)
setting.options = c("S1","S2","S3","S4","S5")

Cov.mat.summary = Leng.mat.summary = matrix(NA,nrow=3*2*5, ncol=2)
for(i.setting in 1:5){
  setting = setting.options[i.setting]
  for(i.VIO.str in 1:2){
    VIO.str = VIO.str.options[i.VIO.str]
    for(i.n in 1:length(n.options)){
      n = n.options[i.n]
      CI.mat.whole = matrix(NA, nrow=500, ncol=4)
      Cov.mat.whole = Leng.mat.whole = matrix(NA, nrow=500, ncol=2)
      for(round in 1:10){
        filename = paste("RDatas/Bootstrap/Bootstrap-Homo-Setting",
                         setting,"-Strength0.5-Violation",VIO.str,
                         "-n",n,"-SimRound",round,".RData", sep="")
        load(filename)
        CI.mat.whole[((round-1)*50+1):(round*50),] = CI.mat
        Cov.mat.whole[((round-1)*50+1):(round*50),] = Cov.mat
        Leng.mat.whole[((round-1)*50+1):(round*50),] = Leng.mat
      }
      ind = (i.setting-1)*2*length(n.options) + (i.VIO.str-1)*length(n.options) + i.n
      Cov.mat.summary[ind, ] = colMeans(Cov.mat.whole)
      Leng.mat.summary[ind, ] = colMeans(Leng.mat.whole)
    }
  }
}

Cov.mat.1 = Cov.mat.summary
Leng.mat.1 = Leng.mat.summary

## Filtering or not
n.options = c(500, 1000, 2000)
VIO.str.options = c(0.2, 0.4)
setting.options = c("S1","S2","S3","S4","S5")

Cov.mat.summary = Leng.mat.summary = matrix(NA,nrow=3*2*5, ncol=2)
for(i.setting in 1:5){
  setting = setting.options[i.setting]
  for(i.VIO.str in 1:2){
    VIO.str = VIO.str.options[i.VIO.str]
    for(i.n in 1:length(n.options)){
      n = n.options[i.n]
      CI.mat.whole = matrix(NA, nrow=500, ncol=4)
      Cov.mat.whole = Leng.mat.whole = matrix(NA, nrow=500, ncol=2)
      for(round in 1:10){
        filename = paste("RDatas/Filtering/Filtering-Homo-Setting",
                         setting,"-Strength0.5-Violation",VIO.str,
                         "-n",n,"-SimRound",round,".RData", sep="")
        load(filename)
        CI.mat.whole[((round-1)*50+1):(round*50),] = CI.mat
        Cov.mat.whole[((round-1)*50+1):(round*50),] = Cov.mat
        Leng.mat.whole[((round-1)*50+1):(round*50),] = Leng.mat
      }
      ind = (i.setting-1)*2*length(n.options) + (i.VIO.str-1)*length(n.options) + i.n
      Cov.mat.summary[ind, ] = colMeans(Cov.mat.whole)
      Leng.mat.summary[ind, ] = colMeans(Leng.mat.whole)
    }
  }
}

Cov.mat.2 = Cov.mat.summary
Leng.mat.2 = Leng.mat.summary

col_setting = rep(1:5, rep(length(n.options)*length(VIO.str.options), length(setting.options)))
col_VIO.str = rep(rep(VIO.str.options, rep(length(n.options), length(VIO.str.options))), length(setting.options))
col_n = rep(n.options, length(setting.options)*length(VIO.str.options))

out.table.1 = cbind(col_setting, col_VIO.str, col_n,
                    cbind(Cov.mat.1[,1], Leng.mat.1[,1], Cov.mat.1[,2], Leng.mat.1[,2])+0.0001)
colnames(out.table.1) = c("set","tau","n","cov","len","cov","len")
library(kableExtra)
kbl(out.table.1, "latex",digits=3, align="c")%>%
  add_header_above(c(" "=3,"TRUE"=2, "FALSE"=2))%>%
  kable_styling(latex_options = c("scale_down"))%>%
  collapse_rows(columns=1:2)

out.table.2 = cbind(col_setting, col_VIO.str, col_n,
                    cbind(Cov.mat.2[,1], Leng.mat.2[,1], Cov.mat.2[,2], Leng.mat.2[,2])+0.0001)
colnames(out.table.2) = c("set","tau","n","cov","len","cov","len")
library(kableExtra)
kbl(out.table.2, "latex",digits=3, align="c")%>%
  add_header_above(c(" "=3,"TRUE"=2, "FALSE"=2))%>%
  kable_styling(latex_options = c("scale_down"))%>%
  collapse_rows(columns=1:2)