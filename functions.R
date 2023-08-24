convertToLongFormat <- function (df)
{
  df.target <- df[c("SampleName", "Green_Channel_NumberOfPositiveDroplets", "Green_Channel_NumberOfNegativeDroplets", "Marker")]
  
  colnames(df.target) <- c("Sample", "Positives", "Negatives", "Marker")
  
  df.target$TargetType <- "Unknown"
  
  df.ref <- df[c("SampleName", "Red_Channel_NumberOfPositiveDroplets", "Red_Channel_NumberOfNegativeDroplets", "Marker")]
  
  colnames(df.ref) <- c("Sample", "Positives", "Negatives", "Marker")
  
  df.ref$TargetType <- "Reference"
  
  return(rbind(df.target,df.ref))
}


extractGroup <- function(df, assay)
{
  if (assay == "dual_color")
  {
    df$Group <- unlist(lapply(strsplit(df$Sample, split = "[_#-]+"), function(vector) vector[2]))  
  }
  if (assay == "dual_color_all_tag_only")
  {
    df$Group <- unlist(lapply(strsplit(df$Sample, split = "[#]+"), function(vector) vector[2]))
  }
  if (ASSAY == "triple_color" | ASSAY=="triple_color_2023_screen")
  {
    df$Group <- df$Sample
  }
  
  return(df)
}


computeNegativeRate <- function(df)
{
  return (1 - (df$Positives / (df$Positives + df$Negatives)))
}


computeNormalizedRatio <- function (df.target, df.reference)
{
  return (log(df.target$NegativeRate) / log(df.reference$NegativeRate))
}


computeConfidenceInterval <- function (df)
{
  df.target <- df[df$TargetType=="TARGET",]
  df.reference <- df[df$TargetType=="REFERENCE",]
  
  df.target$Ratio <- computeNormalizedRatio(df.target, df.reference)
  df <- df.target
  
  df$MinCI95 <- 0
  df$MaxCI95 <- 0
  
  G2 <- log(1-df.target$Positives/(df.target$Positives+df.target$Negatives))/log(1-df.reference$Positives/(df.reference$Positives+df.reference$Negatives))
 
  df$MinCI95 <- exp(log(G2)-(1.96)*sqrt(df.target$Positives/((1-df.target$Positives/(df.target$Positives+df.target$Negatives))*(df.target$Positives+df.target$Negatives)^2*(log(1-df.target$Positives/(df.target$Positives+df.target$Negatives)))^2)+df.reference$Positives/((1-df.reference$Positives/(df.reference$Positives+df.reference$Negatives))*(df.reference$Positives+df.reference$Negatives)^2*(log(1-df.reference$Positives/(df.reference$Positives+df.reference$Negatives)))^2)))
  
  df$MaxCI95 <- exp(log(G2)+(1.96)*sqrt(df.target$Positives/((1-df.target$Positives/(df.target$Positives+df.target$Negatives))*(df.target$Positives+df.target$Negatives)^2*(log(1-df.target$Positives/(df.target$Positives+df.target$Negatives)))^2)+df.reference$Positives/((1-df.reference$Positives/(df.reference$Positives+df.reference$Negatives))*(df.reference$Positives+df.reference$Negatives)^2*(log(1-df.reference$Positives/(df.reference$Positives+df.reference$Negatives)))^2)))
  
  return(df)
}


computeRelativeError <- function (df)
{
  df$RelErr <- (df$MaxCI95 - df$MinCI95) / (2 * df$Ratio)
  return(df)
}