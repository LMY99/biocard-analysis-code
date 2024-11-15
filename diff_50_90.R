age_59 <- c(50,90)
basis_59 <- ibs(age_59, knots=knot.list[[1]], Boundary.knots = boundary.knot, 
             degree=2, intercept=TRUE)
basis_59 <- basis_59[,3:(ncol(basis_59)-2)]
incre_59 <- 
  foreach(i=1:K,.combine=rbind) %do% {
    curves <- basis_59 %*% coefs[-(1:4),i,(R/2+1):R]
    incre <- apply(curves, 2, function(x) x[2]-x[1])
    incre_est <- c(avg=mean(incre), HDInterval::hdi(incre)) 
    incre_est
  } %>% as.data.frame() %>%
  mutate_all(round, digits=2) %>%
  mutate(biomarker=c("MMSE","LogMem","DSST","ENT-THICK","HIPPO","ENT-VOL",
                     "MTL","SPARE-AD","t-tau","p-tau181","AB42/AB40 Ratio")) %>%
  remove_rownames() %>%
  column_to_rownames('biomarker')
