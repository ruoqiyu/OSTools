overlap<-function(score,z,threshold=0.5){
  pooled.sd=sqrt(var(score[z==1])/2+var(score[z==0])/2)
  min.treated.score=min(score[z==1])
  max.control.score=max(score[z==0])
  # How many treated and control subjects lack overlap by Hansen's criterion
  no.treated.lack.overlap=sum(score[z==1]>(max.control.score+threshold*pooled.sd))
  no.control.lack.overlap=sum(score[z==0]<(min.treated.score-threshold*pooled.sd))
  # If there are subjects who lack overlap, record their indices to remove
  which.remove=NULL
  if (no.treated.lack.overlap+no.control.lack.overlap>0){
    which.remove=which((score>(max.control.score+threshold*pooled.sd))|
                         (score<(min.treated.score-threshold*pooled.sd)))
  }
  list(no.treated.lack.overlap=no.treated.lack.overlap,
       no.control.lack.overlap=no.control.lack.overlap,
       index.lack.overlap=which.remove)
}
