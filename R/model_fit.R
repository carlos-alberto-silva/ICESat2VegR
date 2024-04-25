randomForestRegression <- function(
    featurecollection,
    property,
    train_properties,
    nTrees = 500,
    mtry = NULL,
    nodesize = 5,
    sampsize = .632,
    seed = floor((runif(1) * 2147483647 * 2) - 2147483647)) {
  if (is.null(mtry)) {
    mtry <- max(1, floor(length(featurecollection$first()$propertyNames()$getInfo()) / 3))
    classifier <- ee$Classifier$smileRandomForest(
      numberOfTrees = nTrees,
      variablesPerSplit = mtry,
      minLeafPopulation = nodesize,
      bagFraction = sampsize,
      seed = seed
    )$
      train(
      features = featurecollection,
      inputProperties = train_properties,
      classProperty = "h_canopy"
    )$
      setOutputMode("REGRESSION")

    return(classifier)
  }

  ee$Classifier$smileRandomForest(
    numberOfTrees = nTrees,
    variablesPerSplit = mtry,
    minLeafPopulation = nodesize,
    bagFraction = sampsize,
    seed = seed
  )$
    train(
    features = featurecollection,
    classProperty = "h_canopy"
  )$
    setOutputMode("REGRESSION")
}



# modeling LOOCV
model_fit<-function(x=INV[,2:6],y=INV[,7], method="lm", LOOCV = FALSE, ...){

  maxy<-max(y)

  par(mfrow=c(1,2))

  

  if  (method=="lm"){model<- lm(y ~ ., data = x)}
  if  (method=="rf"){model<- randomForest(y=y, x=x, ...)}
  if  (method=="knn.euclidean"){model<-yai(x=x,y=y,method="euclidean", ...)}
  if  (method=="knn.mahalanobis"){model<-yai(x=x,y=y,method="mahalanobis", ...)}
  if  (method=="knn.ica"){model<-yai(x=x,y=y,method="ica", ...)}
  if  (method=="knn.msn"){model<-yai(x=x,y=y,method="msn", ...)}
  if  (method=="knn.msn2"){model<-yai(x=x,y=y,method="msn2", ...)}
  if  (method=="knn.randomForest"){model<-yai(x=x,y=y,method="randomForest", ...)}
  if  (method=="knn.raw"){model<-yai(x=x,y=y,method="raw", ...)}
  if  (method=="nnt"){
    yscales<-(y-mean(y))/sd(y)
    model<- nnet(yscales ~ ., data=x, ...)
  }
  if  (method=="svm"){model<-svm(x=x,y=y, ...)}

  if (LOOCV == TRUE) {
    pb <- txtProgressBar(min = 0, max = nrow(x), style = 3)
  if ( any(method == c("lm","rf","svm"))) {
    pred_model<-as.numeric(paste0(predict(model)))
    pred_model<-cbind(method=rep(method,nrow(x)), obs=y, pred=pred_model)
    stats_model<-cbind(method=rep(method,6), stats_model(y,as.numeric(pred_model[,3]),main=paste("model",method)))
    predloocv<-NULL
    for ( i in 1:nrow(x)){
      setTxtProgressBar(pb, i)
      if ( method=="lm") {model_i<-lm(y[-i] ~ ., data = x[-i,])}
      if ( method=="rf") {model_i<-randomForest(y=y[-i], x=x[-i,], ntree=1000)}
      if ( method=="svm") {model_i<-svm(x=x[-i,], y=y[-i])}

      pred_i<-predict(model_i, newdata = x[i,])
      predloocv<-rbind(predloocv,cbind(method=method,obs=y[i], pred=pred_i))
    }
    stats_loocv<-cbind(method=rep(method,6), stats_model(as.numeric(predloocv[,2]),as.numeric(predloocv[,3]),main=paste("LOOCV",method)))
  }


  if ( any(method == c("knn.euclidean","knn.mahalanobis","knn.ica","knn.msn","knn.msn2","knn.randomForest","knn.raw"))) {
    pred_model<-yaImpute::impute(model,vars=yaImpute::yvars(model))
    
    # yy = yaImpute::yvars(model)[1]
    for (yy in yaImpute::yvars(model)) {
      x11()
      par(mfrow=c(1,2))
      stats_model<-cbind(method=rep(method,6), stats_model(y[,yy],as.numeric(pred_model[,yy]),main=paste("model:",method, "var:", yy)))
      
      predloocv<-NULL
      # i = 1
      for ( i in 1:nrow(x)){
        setTxtProgressBar(pb, i)
        rownames(x) = NULL
        model_i<- yai(x=x[-i,],y=y[-i, yy],method=gsub("knn.","",method), k = 1)
        new.t <- newtargets(model_i,x[i,])
        pred_i<-yaImpute::impute(new.t)[,1]
        
        predloocv<-rbind(predloocv,cbind(method=method,obs=y[i,yy], pred=pred_i))
      }


      stats_loocv<-cbind(method=rep(method,6), stats_model(as.numeric(predloocv[,2]),as.numeric(predloocv[,3]),main=paste("LOOCV",method)))
    }
  }


  if  (method=="nnt"){
    yscales<-(y-mean(y))/sd(y)
    model<- nnet(yscales ~ ., data=x, size=40,linout=T, trace=FALSE,skip=FALSE,maxit = 90, MaxNWts = 2000)
    pred<-predict(model);pred<-pred*sd(y)+mean(y) # removing negative predictions
    pred_model<-cbind(method=rep(method,nrow(x)), obs=y, pred=pred)
    stats_model<-cbind(method=rep(method,6), stats_model(y,as.numeric(pred_model[,3]),main=paste("model",method)))
    predloocv<-NULL
    for ( i in 1:nrow(x)){
      setTxtProgressBar(pb, i)
      model_i<- nnet(yscales[-i] ~ ., data=x[-i,], size=40,linout=T, trace=FALSE,skip=FALSE,maxit = 90, MaxNWts = 2000)
      pred_i<-predict(model_i, newdata=x[i,]);pred_i<-pred_i*sd(y)+mean(y) # removing negative predictions
      predloocv<-rbind(predloocv,cbind(method=method,obs=y[i], pred=pred_i))
    }
    predloocv<-na.omit(predloocv)
    stats_loocv<-cbind(method=rep(method,6), stats_model(as.numeric(predloocv[,2]),as.numeric(predloocv[,3]),main=paste("LOOCV",method)))
  }
  colnames(pred_model)<-c("method","obs","pred")
  colnames(predloocv)<-colnames(pred_model)
  colnames(stats_model)<-c("method","stat","value")
  colnames(stats_loocv)<-c("method","stat","value")
  close(pb)
  return(list(method = method,model=model,predModel=as.data.frame(pred_model),
              predLOOCV=as.data.frame(predloocv),
              statModel=as.data.frame(stats_model),
              statLOOCV=as.data.frame(stats_loocv)))
  } 

  return (list(method = method,model = model))
}