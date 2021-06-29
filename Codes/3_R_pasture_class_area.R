####
####
getBlocksOffset <- function(colOffset, rowOffset, nBlock, nColBlock = NULL, nRowBlock = NULL) {
    if(!is.null(nColBlock)){
        print("blocos por colunas")
        #blocks by col
        blocksOffset <- data.frame(iBlock = 1, colOffset, rowOffset)
        for(iBlock in 2:nBlock){
            blockOffset_i <- c(iBlock, (colOffset + ((iBlock-1) * nColBlock)), rowOffset)
            blocksOffset <-rbind(blocksOffset, blockOffset_i)
        }
        return(blocksOffset)
    } else if(!is.null(nRowBlock)){
        print("blocos por linhas")
        #blocks by col
        blocksOffset <- data.frame(iBlock = 1, colOffset, rowOffset)
        for(iBlock in 2:nBlock){
            blockOffset_i <- c(iBlock, colOffset, (rowOffset + ((iBlock-1) * nRowBlock)))
            blocksOffset <-rbind(blocksOffset, blockOffset_i)
        }
        return(blocksOffset)
    } else {
        print("O numero de linhas ou colunas por blocks nao informado")

    }
}

getOutputExtent <- function(referenceRaster, blockOffset, blockSize) {
    extentBlock <- extent(referenceRaster, blockOffset[2], c((blockOffset[2]-1) + blockSize[2]),
                          blockOffset[1], c((blockOffset[1]-1) + blockSize[1]))
    return(extentBlock)
}

# area cvp
npixel_cvp <- function(data){
  data <- na.omit(data)
  if(nrow(data) < 1){
     c1 <- NA
     c2 <- NA
     c3 <- NA
     c4 <- NA
    }else{
     c1 <- length(data[data == 1]) 
     c2 <- length(data[data == 2]) 
     c3 <- length(data[data == 3]) 
     c4 <- length(data[data == 4])
     }
     return(c(c1, c2, c3, c4))
 }
###
###
procRasterBlocks <- function(inputFile, outputFile, nRowBlock) {
        STTot <- Sys.time()

        #preparing blocks
        referenceRaster <- raster(inputFile)
        nColBlock <- ncol(referenceRaster)
        colOffset <- 1
        rowOffset <- 1
        nRow <- nrow(referenceRaster)
        print(referenceRaster)
		nRowBlock <- nRowBlock
		print(nRowBlock)

        nBlock <- ceiling(nRow / nRowBlock)
        print(paste0("Numero de blocos = ", nBlock))

        iBlockSize <- c(nColBlock, nRowBlock)
        blocksOffset <- getBlocksOffset(colOffset, rowOffset, nBlock, nColBlock = NULL, nRowBlock)

        #process data by blocks
        df_cvp <- data.frame(bloco = 1:nBlock, c1 = NA, c2 = NA, c3 = NA, c4 = NA)
        for (bloco in 1:nBlock) {
 
            iBlockOffset <- as.numeric(blocksOffset[bloco, 2:3])

            STBloco <- Sys.time()
            outputExtent <- getOutputExtent(referenceRaster, iBlockOffset, iBlockSize)

            #read data of the block
            STRead <- Sys.time()
            blockData = data.frame(crop(referenceRaster, outputExtent)[])
            totSTRead <- paste0("time to read block ", bloco, " = ", Sys.time() - STRead)
            cvp_pixels <- npixel_cvp(blockData)
            df_cvp[bloco,] <- c(bloco, cvp_pixels)
            print(c(bloco, cvp_pixels))
        	rm(blockData)
        	gc(reset = TRUE)
        }
        write.csv(df_cvp, outputFile, row.names = FALSE)
        print(paste0("time to execute all blocks = ", Sys.time() - STTot))
}

#Run
suppressWarnings(suppressMessages(library(raster)))

###
# pa_br_pasture_quality_col5_2010_2018
inputDir <- path.files
outputDir <- inputDir

###
###
inputFileList <- Sys.glob(file.path(inputDir, "*biome*cvp*.tif"))
(outputFileList <- gsub(".tif", ".csv", inputFileList))
nRowBlock <- 500
njobs <- length(inputFileList)

for(i in 1:njobs){
    inputFile <- inputFileList[i]
    outputFile <- outputFileList[i]
    procRasterBlocks(inputFile, outputFile, nRowBlock)
    }

####
####
