

### Function returns only the price/value of a computation
fitobPrice <- function ( xmlConfigFileName, scriptName )   {
   # todo: check the correctness of the arguments, if those files exist
   ###if ( !is.string(scriptName) || !is.string(xmlConfigFileName) ){
   ###}
   price <- 0;
   price <- .Call("mainScriptPrice", xmlConfigFileName, scriptName, PACKAGE = "fitob");
   return(price);
}


### Function returns only the price/value of a computation with Monte-Carlo
fitobPriceMC <- function ( xmlConfigFileName, scriptName )   {
   # todo: check the correctness of the arguments, if those files exist
   ###if ( !is.string(scriptName) || !is.string(xmlConfigFileName) ){
   ###}
   price <- 0;
   price <- .Call("mainScriptMCPrice", xmlConfigFileName, scriptName, PACKAGE = "fitob");
   return(price);
}


### Function returns the price/value of a computation and the resulting mesh as well
fitobPriceMesh <- function ( xmlConfigFileName, scriptName, levelIn )   {

   # todo: check the correctness of the arguments, if those files exist 
   ###if ( !is.string(scriptName) || !is.string(xmlConfigFileName) ){
   ###}
   outList <- .Call("mainScriptPriceFG", xmlConfigFileName, scriptName, as.integer(levelIn), PACKAGE = "fitob");
   return(outList);
}


### Plots the resulting mesh from the FITOB pricing
fitobMeshPlot <- function ( inList )   {

   # we only plot for 1 and 2 dimensions
  
   # one dimension
   if ( inList[[2]] == 1) {
      # one dimension
      plot(c(inList[[6]]), inList[[4]], xlab=inList[[7]][1], ylab="Value", main="FITOB");
   }

   # two dimensions   
   if ( inList[[2]] == 2) {
      # create the vectors for the axis
      v1 = array(inList[[5]]);
      v2 = array(inList[[5]]);
      j1 = 1; 
      j2 = 1;
      for (i in 1:(2*inList[[5]])) {
        if (i %% 2 > 0){
           v1[j1] = inList[[6]][i];
           j1 = j1 + 1;
        }
        else{
           v2[j2] = inList[[6]][i];
           j2 = j2 + 1;
        }
      }
      # copy the results into the created matrix 
      z <- matrix(0, inList[[5]], inList[[5]]);
      for (i in 1:inList[[5]]) {
        for (j in 1:inList[[5]]) {
           z[j,i] = inList[[4]][ (i - 1) * inList[[5]] + j ];
        }
      }
      
      # two dimensional plot
      persp(v1, v2, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
            ltheta = 120, shade = 0.75, ticktype = "detailed",
            xlab = inList[[7]][1], ylab = inList[[7]][2], zlab = "Value") -> res
   }
   return(0);
}

# ---- usage ------
# require(fitob);
# liM <- fitobPriceMesh("asian1D.xml" , "asian1DScript", 5);
# lisMy <- fitobPriceMesh("Heston2D.xml" , "Heston2D.thetas",6);