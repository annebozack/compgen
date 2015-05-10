##################################################################################
##							MATRIX FACTORIZATION METHOD						    ##
##################################################################################

library(quadprog)
library(matrix)

#################### GENERATING QUADPROG MATRICES (X CONSTANT) ###################
# Generate Amat

genAmatX = function(nCellTypes, nSamp){
	AmatX = matrix(1, ncol=nCellTypes)
	Amat1 = AmatX
	for (i in 2:nSamp){
		AmatX = bdiag(AmatX, Amat1)
		}
	AmatX = as.matrix(AmatX)
	AmatX = t(AmatX)
	I = diag(nCellTypes)
	I1 = I
	for (i in 2:nSamp){
		I = bdiag(I, I1)
		}
	I = as.matrix(I)
	AmatX = cbind(AmatX, I)
	return(AmatX)
}

# Generate bvec

genBvecX = function(nCellTypes, nSamp){
	bvecX = matrix(1, nrow=nSamp)
	bvec1 = matrix(0, nrow=(nSamp*nCellTypes))
	bvecX = matrix(c(bvecX, bvec1, ncol=1))
	bvecX = bvecX[1:(nCellTypes*nSamp+nSamp),]
	bvecX = matrix(bvecX, ncol=1)
	return(bvecX)
}
	
# Generate dvec
genDvecX = function(nSamp, samples, X){
	dvecX = crossprod(X, samples[,1])
	for (i in 2:nSamp){
		dvecX = matrix(c(dvecX, crossprod(X, samples[,i])), ncol=1)
	}
	return(dvecX)
}

# Generate Dmat

genDmatX = function(nSamp, X){
	DmatX = crossprod(X, X)
	Dmat1 = DmatX
	for (i in 2:nSamp){
		DmatX = bdiag(DmatX, Dmat1)
	}
	return(DmatX)	
}

#################### GENERATING QUADPROG MATRICES (W CONSTANT) ###################
# Generate Amat
genAmatW = function(nCellTypes, nSites){
	I1 = diag(nCellTypes*nSites)
	I2 = (-1)*diag(nCellTypes*nSites)
	AmatW = rbind(I1, I2)
	AmatW = t(AmatW)
	return(AmatW)
}

# Generate bvec
genBvecW = function(nCellTypes, nSites){
	bvecW = matrix(c((matrix(0, nrow=(nCellTypes*nSites), ncol=1)), matrix(-1, nrow=(nCellTypes*nSites), ncol=1)), ncol=1)
	return(bvecW)
}

# Generate Dmat
genDmatW = function(nSamp, nSites, W){
	DmatW = crossprod(t(W[1:6,]))
	for (i in 2:nSamp){
		DmatW = DmatW + crossprod(t(W[(6*i-5):(6*i)]))
	}
	I = diag(nSites)
	DmatW = kronecker(DmatW, I)
	return(DmatW)
}

# Generate dvec
genDvecW = function(nSamp, nCellTypes, W, samples){
	y = t(crossprod(t(W[1:nCellTypes,]), t(samples[,1])))
	for (i in 2:nSamp){
		y = y + t(crossprod(t(W[(6*i-5):(6*i)]), t(samples[,i])))
	}
	dvecW = matrix(y[,1])
	for (i in 2:ncol(y)){
		dvecW = matrix(c(dvecW, y[,i]))
	}
	return(dvecW)
}

############################### SOLVING FOR W ###############################
AmatX = genAmatX(6, 100)
bvecX = genBvecX(6, 100)
solveW = function(nSamp, samples, X, AmatX, bvecX){
	dvecX = genDvecX(nSamp, samples, X)
	DmatX = genDmatX(nSamp, X)
	W = solve.QP(DmatX, dvecX, AmatX, bvec=bvecX, meq=nSamp)
	W = as.matrix(W$solution)
	return(W)
}

############################### SOLVING FOR X ###############################
solveX = function(W, nSamp, nSites, nCellTypes, samples, AmatW, bvecW){
	DmatW = genDmatW(nSamp, nSites, W)
	dvecW = genDvecW(nSamp, nCellTypes, W, samples)
	QPX = solve.QP(DmatW, dvecW, AmatW, bvec=bvecW)
	X = QPX$solution
	X = matrix((as.matrix(X)), ncol=6)
	return(X)
}

##################### ITERATION SOLVING FOR X AND W ##########################
iteration = function(X, W, nSamp, nSites, nCellTypes, samples, AmatX, bvecX, AmatW, bvecW, nIterations, verbose=TRUE){
	for (i in 1:nIterations){
		Wnew = solveW(nSamp, samples, X, AmatX, bvecX)
		normW = norm(W-Wnew)
		if(verbose) cat("W norm", normW, "\n")
		W = Wnew
		Xnew = solveX(W, nSamp, nSites, nCellTypes, samples, AmatW, bvecW)
		normX = norm(X-Xnew)
		if(verbose) cat("X norm", normX, "\n\n")
		X = Xnew
		}
	return(list(X, W))
}

######################## PROCESSING ITERATION RESULTS ########################
convertTestW = function(test2){
	testW = matrix(test2, nrow=6)
	testW = t(testW)
	return(testW)
}

######################## ITERATION WITH AVERAGING #############################
avIteration = function(nSamp, nSites, nCellTypes, samples, nIterations, nAverages, verbose=TRUE){
	AmatX = genAmatX(nCellTypes, nSamp)
	bvecX = genBvecX(nCellTypes, nSamp)
	AmatW = genAmatW(nCellTypes, nSites)
	bvecW = genBvecW(nCellTypes, nSites)
	avW = matrix(0, nrow=nSamp, ncol=nCellTypes)
	for (i in 1:nAverages){
		X = runif((nSites*nCellTypes), 0, 1)
		X = matrix(X, ncol=nCellTypes)
		W = matrix(0, nrow=(nSamp*nCellTypes), ncol=1)
		estXW = iteration(X, W, nSamp, nSites, nCellTypes, samples, AmatX, bvecX, AmatW, bvecW, nIterations, verbose=TRUE)
		X = estXW[[1]]
		W = convertTestW(estXW[[2]])
		W = cbind(W[,which.max(X[1,])], W[,which.max(X[2,])], W[,which.max(X[3,])], W[,which.max(X[4,])], W[,which.max(X[5,])], W[,which.max(X[6,])])
		avW = avW + W
		X = cbind(X[,which.max(X[1,])], X[,which.max(X[2,])], X[,which.max(X[3,])], X[,which.max(X[4,])], X[,which.max(X[5,])], X[,which.max(X[6,])])
	}
	avW = avW/nAverages
	return(list("W" = avW, "X" = X))
}

# If verbose = TRUE, prints norm(previous W = new W) and norm(previous X - new X)

######################## GENERATING SIMULATED DATA #############################
# Simulated data is included in data file.  This code is provided to demonstrate how 
# simulated data was generated.  The refernce data set, FlowSorted.Blood.450k, was
# downloaded with the FlowSorted.Blood.450k R package (source("http://bioconductor.org/
# biocLite.R"), biocLite("FlowSorted.Blood.450k")).

# # Process data to obtain betas (methylation percentages)
# source("http://bioconductor.org/biocLite.R")
# biocLite("minfi")
# require(minfi)
# Mset = preprocessIllumina(FlowSorted.Blood.450k)
# betas = getBeta(Mset)

# # CpG sites approximating criteria for unicity
# topCpg = c("cg09993145", "cg26986871", "cg08206665", "cg08864105", "cg22651103", "cg01710351", "cg24641737", "cg03429643", "cg10768932", "cg18886071", "cg08560387", "cg14317609", "cg02606840", "cg24838825", "cg07725527", "cg21761853", "cg23954655", "cg10825315", "cg21221767", "cg13203135", "cg07489413", "cg19455189", "cg11334870", "cg17343167", "cg00776174", "cg13437525", "cg26986871", "cg23942604", "cg18440048", "cg21268578", "cg19176391", "cg22632840", "cg25753024", "cg25518868", "cg08130572", "cg09804858", "cg20236089", "cg27166718", "cg03155112", "cg09121680", "cg24804707", "cg08560387", "cg24838825", "cg27471192", "cg03699843", "cg08580187", "cg00854273", "cg23551494", "cg14595786", "cg09228833")

# topCpgBetas = betas[topCpg,]
# avGran = t(t(rowMeans(topCpgBetas[,c(7, 8, 9, 31, 32, 33)])))
# avCD4 = t(t(rowMeans(topCpgBetas[,c(10, 11, 12, 34, 35, 36)])))
# avCD8 = t(t(rowMeans(topCpgBetas[,c(13, 14, 15, 37, 38, 39)])))
# avCD19 = t(t(rowMeans(topCpgBetas[,c(16, 17, 18, 40, 41, 42)])))
# avCD14 = t(t(rowMeans(topCpgBetas[,c(19, 20, 21, 43, 44, 45)])))
# avCD56 = t(t(rowMeans(topCpgBetas[,c(22, 23, 24, 46, 47, 48)])))

# meanTopCpgBetas = cbind(avGran, avCD4, avCD8, avCD19, avCD14, avCD56)

# headers = cbind("Gran", "CD4", "CD8", "CD19", "CD14", "CD56")
# colnames(meanTopCpgBetas) = headers

# # Generating a random matrix of cell type proportions for 100 samples
# randMatrix = function(nSamp){
	# x = matrix(, ncol=6)
	# r = NULL
	# for (i in 1:nSamp){
		# r = sort(runif(5, 0, 1))
		# r = c(r[1], (r[2]-r[1]), (r[3]-r[2]), (r[4]-r[3]), (r[5]-r[4]), (1-r[5]))
		# x = rbind(x, r)
	# }
	# x = x[2:(nSamp+1),]
	# return(x)
# }

# matrix100 = randMatrix(100)

# # Generating 100 samples
# genSamples = function(nBetas, nSamp, randMatrix, betas){
	# x = matrix(, nrow=nBetas, ncol=nSamp)
	# for (i in 1:nSamp){
		# x[,i] = (t(t((betas[,"Gran"])*randMatrix[i,1])) + t(t((betas[,"CD4"])*randMatrix[i,2])) + t(t((betas[,"CD8"])*randMatrix[i,3])) + t(t((betas[,"CD19"])*randMatrix[i,4])) + t(t((betas[,"CD14"])*randMatrix[i,5])) + t(t((betas[,"CD56"])*randMatrix[i,6])))
	# }
	# return(x)
# }

# samplesSim = genSamples(50, 100, matrix100, meanTopCpgBetas)

################# ESTIMATING X AND W WITH SIMULATED DATA #######################
estXW = avIteration(100, 50, 6, samples, 300, 10)

estimatedW = estXW$W
estimatedX = estXW$X

# Correlation of estimated W with actual used to generate data
correlation = cor(matrix(estimatedW, ncol=1), matrix(matrix100, ncol=1))


#################################################################################
# The following code was used to evaluated the six whole blood samples contained in 
# the reference data set and 36 cord whole blood samples obtained from Columbia 
# University’s Center for Children’s Environmental Health.

################## DATA SET INCLUDING WHOLE BLOOD FROM REFERECE #################
# WBCpg = betas[topCpg,]
# WBCpg = cbind(WBCpg[,1], WBCpg[,2], WBCpg[,3], WBCpg[,25], WBCpg[,26], WBCpg[,27])
# WBCpg = cbind(WBCpg, samplesSim[,7:100])

######## PERFORMING HOUSEMAN METHOD ON REFERENCE WHOLE BLOOD DATA ###############
# wh.WBC = which(FlowSorted.Blood.450k$CellType == "WBC")
# RGset = FlowSorted.Blood.450k[, wh.WBC]
# sampleNames(RGset) = paste(RGset$CellType, c(seq(along = wh.WBC), seq(along = wh.PBMC)), sep = "_")
# counts = estimateCellCounts(RGset, meanPlot = FALSE)
# counts = cbind(counts[,6], counts[,2], counts[,1], counts[,4], counts[,5], counts[,3])

################# ESTIMATING X AND W WITH REFERENCE WB DATA #######################
# WBestXW = avIteration(100, 50, 6, WBCpg, 300, 5)

# WBestimatedW = WBestXW$W
# WBestimatedW = WBestimatedW[1:6,]
# correlation = cor(matrix(WBestimatedW, ncol=1), matrix(counts, ncol=1))


################## DATA SET INCLUDING CHILDRENS CENTER DATA #######################
# setwd("/Users/anne/Documents/School_Work_Fall_2014/Rotation/")
# baseDir = "idat"
# targets = read.csv(file.path(baseDir, "Samples.csv"))
# targets$Basename = file.path(baseDir, targets$Sentrix.Barcode, paste0(targets$Sentrix.Barcode, "_", targets$Sample.Section))
# RGset <- read.450k.exp(base = baseDir, targets = targets)
# Mset.norm <- preprocessIllumina(RGset, bg.correct = TRUE, normalize = "controls", reference = 2)
# childCentBetas = getBeta(Mset.norm)

# childCentCpg = childCentBetas[topCpg,]
# childCentCpg = cbind(childCentCpg, samplesSim[,37:100])

############ PERFORMING HOUSEMAN METHOD ON CHILDRENS CENTER DATA ##################
# counts = estimateCellCounts(RGset, meanPlot = FALSE)
# counts = cbind(counts[,6], counts[,2], counts[,1], counts[,4], counts[,5], counts[,3])

################ ESTIMATING X AND W WITH CHILDRENS CENTER DATA #####################
# childCentXW = avIteration(100, 50, 6, childCentCpg, 300, 5)

# childCentestimatedW = childCentXW$W
# childCentestimatedW = childCentestimatedW[1:36,]
# correlation = cor(matrix(childCentestimatedW, ncol=1), matrix(counts, ncol=1))