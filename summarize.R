flex=rep(0,5);S=rep(0,5);para=rep(0,5)
for(j in 1:5){
  load(sprintf('flex_cv%d.rda',j)); flex[j]=cv_error
  load(sprintf('para_cv%d.rda',j)); para[j]=cv_error
  load(sprintf('S_cv%d.rda',j)); S[j]=cv_error
}
sink('CV_result.txt')
print(c(flex=mean(flex),para=mean(para),S=mean(S)))
sink()