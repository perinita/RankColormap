## Optimal Rank Decomposition
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tictoc)

######################
#### Initialize
######################

#Convex multiplier dataframe
#stepsize = 0.01
#df = expand.grid(seq(0,1,stepsize),seq(0,1,stepsize))
meshsize = 10^3
df = expand.grid(seq(0,1,length.out=sqrt(2*meshsize)),seq(0,1,length.out=sqrt(2*meshsize)))
df = df[df$Var1+df$Var2<=1,]
Lambda = data.frame(lambda1=df$Var1, lambda2=df$Var2)
Lambda$lambda3 = 1-Lambda$lambda1-Lambda$lambda2

# Transformation function for equilateral triangle
EquiTransform <- function(df) {
  #2 columns of a dataframe (x and y)
  df2 <- df
  df2[,1] <- df2[,1]-0.5*(1-df2[,2])
  df2[,2] <- 0.5*sqrt(3)*df2[,2]
  return(df2)
}

######################
#### Generate Instances
######################

#number of items:
n=5 #5, 10, 15, 20

#Input ranks
# NOTE: the values in the "rank" vector will be treated as continuous variables. 
# Therefore must interpret as rank[i]=1 if item i is ranked FIRST. rank[i] gives the "rating" of item i

# Instance 1: similar ranks
rank1 = 1:n #ranked in order of items: 1st item is FIRST, 2nd item is SECOND, etc.
rank2 = c(2:n,1) #nth item is ranked FIRST, 1st item is ranked SECOND, etc.
rank3 = c(1,3,2,4:n) #1st item is FIRST, 2nd item is THIRD, 3rd item is SECOND, rest follow position

# Instance 2: randomized ranks
set.seed(123)
rank1 = sample(n)
rank2 = sample(n)
rank3 = sample(n)

# Instance 3: randomized ratings
set.seed(1)
rank1 = runif(n)
rank2 = runif(n)
rank3 = runif(n)
# rank from ratings
rank1=order(rank1)
rank2=order(rank2)
rank3=order(rank3)

######################
#### Label Heuristic
######################


#Compute Labels for convex multipliers
if(TRUE) {
  tic()
  Lambda$Label = ""
  for(i in 1:dim(Lambda)[1]) {
    rating = Lambda$lambda1[i]*rank1 + Lambda$lambda2[i]*rank2 + Lambda$lambda3[i]*rank3
    #nonlinear rating:
    #rating = Lambda$lambda1[i]*rank1 + Lambda$lambda2[i]*rank2 + f(Lambda$lambda3[i])*rank3
    rr = rank(rating, ties.method="first") #was order
    Lambda$Label[i] = paste(rr,collapse = "")
  }
  Lambda$Label <- as.factor(Lambda$Label)
  toc()
}
summary(Lambda)
length(levels(Lambda$Label)) #number of rankings/labels found

######################
#### Plot Heuristic Labels
######################
equiLambda = EquiTransform(Lambda)
# Basic plot
g1 <- ggplot() + 
  geom_point(data=equiLambda, aes(x=lambda1,y=lambda2,color=Label)) + 
  theme_void() +
  theme(legend.position="none") +
  #add bias axes
  geom_segment(aes(x = -0.5, y = 0, xend = 0.25, yend = 0.25*sqrt(3)), color="gray") + 
  geom_segment(aes(x = 0.5, y = 0, xend = -0.25, yend = 0.25*sqrt(3)), color="gray") + 
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 0.5*sqrt(3)), color="gray") +
  xlim(c(-0.6,0.6)) + 
  ylim(c(-0.05,0.95))

# annotations
labels = data.frame(x=c(0.5,0,-0.5),
                    y=c(0,0.5*sqrt(3),0),
                    deltay=c(-0.05,0.05,-0.05),
                    lab=c("r1", "r2", "r3"))

g1 + 
  geom_point(data=labels, aes(x=x, y=y)) +
  geom_text(data=labels, aes(x=x, y=y+deltay, label=lab))

# Histogram for areas

Areadf <- Lambda %>% count(Label) #%>% mutate(count = n())
totalN = sum(Areadf$n)
Areadf$Percent = round(100*Areadf$n/totalN,1)

g2 <- ggplot() + 
  geom_bar(data=Areadf, aes(x=reorder(Label,-n),y=Percent,fill=Label),stat="identity") +
  geom_text(data=Areadf, aes(x=reorder(Label,-n),y=Percent+1,label=paste0(Percent,"%")),size=3) +
  xlab("") + ylab("") +
  theme(legend.position="none",
        axis.text.x = element_text(angle = -45, vjust = 0, hjust=0.1)) 
g2

######################
#### Compute Line Segments
######################

#describe lines by border points (bp), i.e., where they intersect with border of triangle
InLambda <- function(lambda) {
  if(sum(lambda<0)>0) return(FALSE)
  if(sum(lambda)>1) return(FALSE)
  return(TRUE)
}
AddIfUnique <- function(df,row) {
  nrow=dim(df)[1]
  ncol=dim(df)[2]
  counter=rep(0,nrow)
  for(i in 1:ncol) {
    counter=counter+(df[,i]==row[i])
  }
  if(sum(counter==ncol)==0) df=rbind(df,row)
  return(df)
}
Linedf2 = data.frame(bp1x = 0, bp1y = 0, bp2x = 0, bp2y = 0, slope = 0, intercept = 0,
                     itemA = 0,itemB = 0, label="blank")
IntersectPts = data.frame(x=0, y=0)
IRextremepts <- data.frame(lambda1=0, lambda2=0, label = "blank")

if(TRUE) {
  tic()
  #Step 1: Compute Line segments
  for(A in 1:(n-1)) {
    for(B in (A+1):n){
      AbetterthanB=0
      if(rank1[A]<rank1[B]) AbetterthanB=AbetterthanB+1
      if(rank2[A]<rank2[B]) AbetterthanB=AbetterthanB+1
      if(rank3[A]<rank3[B]) AbetterthanB=AbetterthanB+1
      if(AbetterthanB==0 || AbetterthanB==3) next
      delta1=rank1[A]-rank1[B]
      delta2=rank2[A]-rank2[B]
      delta3=rank3[A]-rank3[B]
      if(delta2-delta3==0) { #then vertical line
        val = -delta3/(delta1-delta3)
        slope=NaN
        incpt=val
        bp1x = val
        bp1y = 0
        bp2x = val
        bp2y = 1-val
        
      } else {
        slope=-(delta1-delta3)/(delta2-delta3)
        incpt=-delta3/(delta2-delta3)
        p1=c(-1,-1) #3 dummy points. the one that is not replaced will not be returned
        p2=c(-1,-1)
        p3=c(-1,-1)
        if(delta1!=delta3) {
          p1=c(-delta3/(delta1-delta3), 0)
        } 
        if(delta2!=delta3) {
          p2=c(0,-delta3/(delta2-delta3))
        }
        if(delta1!=delta2) {
          p3[1]=-delta2/(delta1-delta2)
          p3[2]=1-p3[1]
        }
        if(InLambda(p1)==FALSE) {
          bp1x = p2[1]
          bp1y = p2[2]
          bp2x = p3[1]
          bp2y = p3[2]  
        } else if(InLambda(p2)==FALSE) {
          bp1x = p1[1]
          bp1y = p1[2]
          bp2x = p3[1]
          bp2y = p3[2]  
        } else if(InLambda(p3)==FALSE) {
          bp1x = p1[1]
          bp1y = p1[2]
          bp2x = p2[1]
          bp2y = p2[2]  
        } else { #all are in Lambda, 2 should be equivalent
          if(sum(p1==p2)==2) {
            bp1x = p2[1]
            bp1y = p2[2]
            bp2x = p3[1]
            bp2y = p3[2]  
          } else if(sum(p2==p3)==2) {
            bp1x = p1[1]
            bp1y = p1[2]
            bp2x = p3[1]
            bp2y = p3[2]  
          } else if(sum(p1==p3)==2) {
            bp1x = p1[1]
            bp1y = p1[2]
            bp2x = p2[1]
            bp2y = p2[2]  
          }
        }
      }
      newrow = data.frame(bp1x = bp1x, bp1y = bp1y, bp2x = bp2x, bp2y = bp2y, 
                          slope = slope, intercept = incpt,
                          itemA = A,itemB = B, label=paste(c(A,B),collapse = "&"))
      Linedf2 = rbind(Linedf2,newrow)
    }
  }
  Linedf2 <- Linedf2[-1,]
  Linedf2$label<-as.factor(Linedf2$label)
  toc()
  #Step 2: Compute Intersection Points
  #add boundary points
  tic()
  for(l1 in 1:dim(Linedf2)[1]) {
    IntersectPts=AddIfUnique(IntersectPts,c(Linedf2$bp1x[l1],Linedf2$bp1y[l1]))
    IntersectPts=AddIfUnique(IntersectPts,c(Linedf2$bp2x[l1],Linedf2$bp2y[l1]))
  }
  #add corner points
  IntersectPts=rbind(IntersectPts,c(0,0))
  IntersectPts=rbind(IntersectPts,c(1,0))
  IntersectPts=rbind(IntersectPts,c(0,1))
  
  for(l1 in 1:(dim(Linedf2)[1]-1)) {
    for(l2 in 2:dim(Linedf2)[1]) {
      # if either endpoint are equal, then skip
      l1bp1 = c(Linedf2$bp1x[l1],Linedf2$bp1y[l1])
      l1bp2 = c(Linedf2$bp2x[l1],Linedf2$bp2y[l1])
      l2bp1 = c(Linedf2$bp1x[l2],Linedf2$bp1y[l2])
      l2bp2 = c(Linedf2$bp2x[l2],Linedf2$bp2y[l2])
      if(sum(l1bp1==l2bp1)==2) next
      if(sum(l1bp1==l2bp2)==2) next
      if(sum(l1bp2==l2bp1)==2) next
      if(sum(l1bp2==l2bp2)==2) next
      # for slope NA (vertical), then compute intersection
      if(is.na(Linedf2$slope[l1]) && is.na(Linedf2$slope[l2])) next
      if(is.na(Linedf2$slope[l1])) {
        xval = Linedf2$intercept[l1]
        yval = Linedf2$intercept[l2]+Linedf2$slope[l2]*xval
        if(xval>=0 && yval>=0 && xval+yval<=1) IntersectPts=AddIfUnique(IntersectPts,c(xval,yval))
        next
      }
      if(is.na(Linedf2$slope[l2])) {
        xval = Linedf2$intercept[l2]
        yval = Linedf2$intercept[l1]+Linedf2$slope[l1]*xval
        if(xval>=0 && yval>=0 && xval+yval<=1) IntersectPts=AddIfUnique(IntersectPts,c(xval,yval))
        next 
      }
      # if slopes equal, then skip
      if(Linedf2$slope[l1]==Linedf2$slope[l2]) next
      # otherwise, compute intersection
      xval = (Linedf2$intercept[l1]-Linedf2$intercept[l2])/(Linedf2$slope[l2]-Linedf2$slope[l1])
      yval = Linedf2$intercept[l2]+Linedf2$slope[l2]*xval
      if(xval>=0 && yval>=0 && xval+yval<=1) IntersectPts=AddIfUnique(IntersectPts,c(xval,yval))
      
    }
  }
  IntersectPts = IntersectPts[-1,]
  toc()
  #Step 3: Assign rankings to intersection points
  tic()
  for(lab in levels(Lambda$Label)) {
    #print(lab)
    num=0
    subdf <- subset(Lambda,Label==lab)
    meanx = mean(subdf$lambda1)
    meany = mean(subdf$lambda2)
    meanrating = meanx*rank1+meany*rank2+(1-meanx-meany)*rank3
    meanrank = rank(meanrating, ties.method="first") #meanrank is the true ranking for the IR
    for(i in 1:dim(IntersectPts)[1]) {
      lambda=c(IntersectPts$x[i], IntersectPts$y[i], 1-IntersectPts$x[i]-IntersectPts$y[i])
      lambdarating = lambda[1]*rank1+lambda[2]*rank2+lambda[3]*rank3
      #test if meanrank is true for lambdarating
      test=TRUE
      j=1
      while(test==TRUE && j<length(meanrank)) {
        #if(lambdarating[meanrank[j]] - lambdarating[meanrank[j+1]] > 0.001) test=FALSE
        if(lambdarating[which(meanrank==j)] - lambdarating[which(meanrank==j+1)] > 0.001) test=FALSE
        j=j+1
      }
      if(test==TRUE) {
        IRextremepts=rbind(IRextremepts,c(IntersectPts$x[i],IntersectPts$y[i], label=lab))
        num=num+1
      }
    }
    #print(paste(lab,":",num))
  }
  IRextremepts = IRextremepts[-1,]
  IRextremepts$lambda1 <- as.numeric(IRextremepts$lambda1)
  IRextremepts$lambda2 <- as.numeric(IRextremepts$lambda2)
  IRextremepts$label <- as.factor(IRextremepts$label)
  toc()
  #Step 4: Calculate convex hull for each label
  tic()
  IR_hull <- IRextremepts %>%
    group_by(label) %>%
    slice(chull(lambda1, lambda2))
  toc()
}
length(unique(IR_hull$label))

#Plot computed line segments and intersection points  
g3 <- ggplot() + 
  #labels from heuristic:
  geom_point(data=Lambda, aes(x=lambda1,y=lambda2,color=Label)) + 
  #line segments: 
  geom_segment(data=Linedf2,aes(x=bp1x,y=bp1y,xend=bp2x,yend=bp2y,group=label)) +
  #intersection points: 
  geom_point(data=IntersectPts, aes(x=x, y=y)) +
  theme(legend.position="none") 
g3 

#Plot convex hulls
g4 <- ggplot() + 
  geom_polygon(data = IR_hull, aes(x=lambda1, y=lambda2, fill=label)) +
  theme(legend.position="none") 
g4

#Final Beautiful Plot
equiLambda=EquiTransform(Lambda)
equiIR_hull=EquiTransform(IR_hull)
equiLines = cbind(EquiTransform(Linedf2[,1:2]),EquiTransform(Linedf2[,3:9]))

g5 <- ggplot() + 
  geom_polygon(data = equiIR_hull, aes(x=lambda1, y=lambda2, fill=label)) + 
  geom_segment(aes(x = -0.5, y = 0, xend = 0.25, yend = 0.25*sqrt(3)), color="gray") + 
  geom_segment(aes(x = 0.5, y = 0, xend = -0.25, yend = 0.25*sqrt(3)), color="gray") + 
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 0.5*sqrt(3)), color="gray") +
  #geom_segment(data=equiLines,aes(x=bp1x,y=bp1y,xend=bp2x,yend=bp2y,group=label), size=1) +
  geom_path(aes(x=c(-0.505,0.505,0,-0.505), y=c(0,0,0.5*sqrt(3),0)),color="black", size=2) +
  theme_void() +
  theme(legend.position="none") +
  ylim(c(0,0.5*sqrt(3)))
g5


######################
#### Compute Centroids
######################

# Compute centroids
IR_centers <- IRextremepts %>%
  group_by(label) %>%
  summarise(CentroidX = mean(lambda1), CentroidY = mean(lambda2))

#Plot IRs by convex hull and centroid
g6 <- ggplot() + 
  #geom_point(data=Lambda, aes(x=lambda1,y=lambda2,color=Label)) + 
  geom_polygon(data = IR_hull, aes(x=lambda1, y=lambda2, fill=label)) + 
  geom_point(data = IR_centers, aes(x=CentroidX, y=CentroidY), color="gray") +
  #geom_segment(data=Linedf2,aes(x=bp1x,y=bp1y,xend=bp2x,yend=bp2y,group=label)) +
  #geom_point(data=IntersectPts, aes(x=x, y=y)) +
  xlim(c(0,1)) + 
  ylim(c(0,1))
g6

######################
#### Plot Heatmap by Rank
######################

item=1
IR_hull$Position = 0
for(r in 1:dim(IR_hull)[1]) {
  lab=IR_hull$label[r]
  IR_hull$Position[r] = unlist(gregexpr(item, lab))
}
equiLambda=EquiTransform(Lambda)
equiIR_hull=EquiTransform(IR_hull)
equiLines = cbind(EquiTransform(Linedf2[,1:2]),EquiTransform(Linedf2[,3:9]))

g7 <- ggplot() + 
  geom_polygon(data = equiIR_hull, aes(x=lambda1, y=lambda2, group=label, fill=Position)) + 
  geom_segment(data=equiLines,aes(x=bp1x,y=bp1y,xend=bp2x,yend=bp2y,group=label), size=1) +
  geom_path(aes(x=c(-0.505,0.505,0,-0.505), y=c(0,0,0.5*sqrt(3),0)),color="black", size=2) +
  theme_void() +
  scale_fill_gradient(trans = 'reverse',limits=c(5,1)) +
  #theme(legend.position="none") +
  ylim(c(0,0.5*sqrt(3)))
g7

######################
#### Plot Heatmap by Distance from Centroid
######################

Lambda$CentroidX = 0
Lambda$CentroidY = 0
Lambda$CentroidDist = 0

for(i in 1:dim(Lambda)[1]) {
  point = c(Lambda$lambda1[i], Lambda$lambda2[i])
  centroid = as.numeric(IR_centers[which(as.character(IR_centers$label)==as.character(Lambda$Label[i])),
                                   c("CentroidX","CentroidY")])
  Lambda$CentroidX[i]=centroid[1]
  Lambda$CentroidY[i]=centroid[2]
  Lambda$CentroidDist[i]=sqrt(sum((point - centroid) ^ 2))
}
Lambda$CentroidDistScaled=-Lambda$CentroidDist/max(Lambda$CentroidDist)

equiLambda=EquiTransform(Lambda)
equiIR_hull=EquiTransform(IR_hull)
equiLines = cbind(EquiTransform(Linedf2[,1:2]),EquiTransform(Linedf2[,3:9]))

#Plot Heatmap for distance from centroid
g4 <- ggplot() + 
  #geom_point(data=Lambda, aes(x=lambda1,y=lambda2,color=Label,alpha=CentroidDist)) + 
  geom_polygon(data = equiIR_hull, aes(x=lambda1, y=lambda2, fill=label)) + 
  geom_point(data=equiLambda, aes(x=lambda1, y=lambda2, alpha=CentroidDist), color="white") +
  #geom_point(data = IR_centers, aes(x=CentroidX, y=CentroidY), color="gray") +
  geom_segment(data=equiLines,aes(x=bp1x,y=bp1y,xend=bp2x,yend=bp2y,group=label), size=1) +
  #geom_point(data=IntersectPts, aes(x=x, y=y)) +
  #xlim(c(0,1)) + 
  #ylim(c(0,1))
  geom_path(aes(x=c(-0.505,0.505,0,-0.505), y=c(0,0,0.5*sqrt(3),0)),color="black", size=2) +
  theme_void() +
  theme(legend.position="none") +
  ylim(c(0,0.5*sqrt(3)))
g4
