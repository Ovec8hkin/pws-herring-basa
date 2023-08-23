#  data_reader.R
#  Created by John T. Trochta
#  This file reads in .DAT and .CTL files created for ADMB .TPL files
#  The use must ensure not to have within line comments in the files (nope, fixed that)

data.reader <- function(filename) {
    #  The user needs to make sure there is a blank line at the end of
    #  filename and that each data type (vector number or matrix) is
    #  separated by a blank line
    
    # This is kind of convoluted
    text <- readLines(filename)
    values <- grep("^\\s{0,2}[0-9]", text)
    signed.values <- grep("^\\s{0,2}[-]", text)
    read.these <- sort(c(values, signed.values))
    nlines <- length(text) 
    indices <- seq(1:nlines)
    indices <- indices[read.these]
    first.differences <- c(diff(indices),5)# This accounts for the last data
    data.types <- length(first.differences[first.differences>1])
    
    data <- vector("list", data.types)
    j <- 1
    temp <- NA
    for(i in 1:length(indices)){
        temp.1 <- scan(filename, skip=indices[i]-1, nlines=1, quiet=TRUE, flush=FALSE)
        if(first.differences[i]>1 | all(first.differences[1:(length(first.differences)-1)]==1)){
            if(all(is.na(temp))){
                data[[j]] <- temp.1
                j <- j+1
            }else{
                temp <- rbind(temp, temp.1)
                data[[j]] <- temp
                rownames(data[[j]]) <- NULL
                temp <- NA
                j <- j+1
            }
        } else{
            if(all(is.na(temp))){
                temp <- temp.1
            }else{
                temp <- rbind(temp, temp.1)
            }
        }
    }
    return(data)
}

# Function for simultaneously creating names of variables/elements within a list
listN <- function(...){
    anonList <- list(...)
    names(anonList) <- as.character(substitute(list(...)))[-1]
    return(anonList)
}
