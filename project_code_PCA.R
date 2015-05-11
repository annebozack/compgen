##################################################################################
##									PCA METHOD								    ##
##################################################################################

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

# matrix500 = randMatrix(500)

# # Generating 500 samples
# genSamples = function(nBetas, nSamp, randMatrix, betas){
	# x = matrix(, nrow=nBetas, ncol=nSamp)
	# for (i in 1:nSamp){
		# x[,i] = (t(t((betas[,"Gran"])*randMatrix[i,1])) + t(t((betas[,"CD4"])*randMatrix[i,2])) + t(t((betas[,"CD8"])*randMatrix[i,3])) + t(t((betas[,"CD19"])*randMatrix[i,4])) + t(t((betas[,"CD14"])*randMatrix[i,5])) + t(t((betas[,"CD56"])*randMatrix[i,6])))
	# }
	# return(x)
# }

# samples500 = genSamples(50, 500, matrix500, meanTopCpgBetas)

################################# DETERMINE PCs ##################################
res = prcomp(t(samples500), center = TRUE, scale = FALSE)

cumsum(res$sdev^2/sum(res$sdev^2))

pc.use = 5 # explains 100% of variance
truncPC = res$x[,1:pc.use]

#################################### SHIFT PCs ####################################
# PCs are shifted to be positive. This step is for the ease of constructing the 
# initial hyperplanes.

twos = matrix(2, nrow=500, ncol=5)
shiftTruncPC = truncPC + twos

############################ FUNCTION TO BE MINIMIZED #############################
detFun = function(t){
	A = matrix(t[1:30], nrow = 6)
	b = -(matrix(t[31:36], nrow = 6))
	A1 = A[1:5,]
	b1 = b[1:5,]
	A2 = A[2:6,]
	b2 = b[2:6,]
	A3 = rbind(A[1,], A[3:6,])
	b3 = rbind(b[1,], matrix(b[3:6,], nrow=4))
	A4 = rbind(A[1:2,], A[4:6,])
	b4 = rbind(matrix(b[1:2,], nrow=2), matrix(b[4:6,], nrow=3))
	A5 = rbind(A[1:3,], A[5:6,])
	b5 = rbind(matrix(b[1:3,], nrow=3), matrix(b[5:6,], nrow=2))
	A6 = rbind(A[1:4,], A[6,])
	b6 = rbind(matrix(b[1:4,], nrow=4), b[6,])
	vertices = cbind(mldivide(A1, b1), mldivide(A2, b2), mldivide(A3, b3), mldivide(A4, b4), mldivide(A5, b5), mldivide(A6, b6))
	X = rbind(vertices[,1], vertices[,2], vertices[,3], vertices[,4], vertices[,5], vertices[,6])
	# Test if points are contained within hyperplane
	cellProp = cart2bary(X, P) 
	# Penalty if points are outside of hyperplane (any barycentric coordinates < 0) 
	if (min(cellProp) < 0) area = 1000000 else area = abs(det(rbind(vertices, rep.int(1, 6))))
	print(area)
	return(area)
}

########################### CONSTRAINTS FOR OPTIMIZATION #############################
# ui %*% theta - ci >= 0, ak is contructed to ensure all points are on the 
# same side of the hyperplane.

a = t(shiftTruncPC)
a = rbind(a, rep.int(1, 500))
ak = kronecker(t(a), diag(6))

ci = rep.int(0, 3000)

# Global variable P is constructed to determine barycentric coordinates within detFun
P = a[1:5,]
P <<- t(P)

######################### PARAMETERS OF INITIAL HYPERPLANES ###########################
t = c(1, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 15)

################# OPTIMIZATION TO FIND HYPERPLANES OF MINIMUM VOLUME ###################
minDet = constrOptim(t, detFun, grad=NULL, outer.eps = 1e-6, control=list(maxit=10000), ui=ak, ci=ci)

########################### FINDING VERTICES OF HYPERPLANES #############################
detPoints = function(x){
	A = matrix(x[1:30], nrow = 6)
	b = -(matrix(x[31:36], nrow = 6))
	A1 = A[1:5,]
	b1 = b[1:5,]
	A2 = A[2:6,]
	b2 = b[2:6,]
	A3 = rbind(A[1,], A[3:6,])
	b3 = rbind(b[1,], matrix(b[3:6,], nrow=4))
	A4 = rbind(A[1:2,], A[4:6,])
	b4 = rbind(matrix(b[1:2,], nrow=2), matrix(b[4:6,], nrow=3))
	A5 = rbind(A[1:3,], A[5:6,])
	b5 = rbind(matrix(b[1:3,], nrow=3), matrix(b[5:6,], nrow=2))
	A6 = rbind(A[1:4,], A[6,])
	b6 = rbind(matrix(b[1:4,], nrow=4), b[6,])
	verticies = cbind(mldivide(A1, b1), mldivide(A2, b2), mldivide(A3, b3), mldivide(A4, b4), mldivide(A5, b5), mldivide(A6, b6))
	return(vertices)
}

pnts = detPoints(minDet$par)

############################ FINDING BARYCENTRIC COORDINATES ##############################
X = rbind(pnts[,1], pnts[,2], pnts[,3], pnts[,4], pnts[,5], pnts[,6])
cellProp = cart2bary(X, P)

################# DETERMINING ORDER OF COLUMNS IN BARYCENTRIC COORDINATES ##################
# The columns of the barycentric coordinates must be reordered to correspond to the appropriate
# columns of the matrix used to generate the data.

norm(matrix((cellProp[,1]-matrix500[,1])/matrix500[,1]))
norm(matrix((cellProp[,1]-matrix500[,2])/matrix500[,2]))
norm(matrix((cellProp[,1]-matrix500[,3])/matrix500[,3]))
norm(matrix((cellProp[,1]-matrix500[,4])/matrix500[,4]))
norm(matrix((cellProp[,1]-matrix500[,5])/matrix500[,5]))
norm(matrix((cellProp[,1]-matrix500[,6])/matrix500[,6]))
norm(matrix((cellProp[,2]-matrix500[,1])/matrix500[,1]))
norm(matrix((cellProp[,2]-matrix500[,2])/matrix500[,2]))
norm(matrix((cellProp[,2]-matrix500[,3])/matrix500[,3]))
norm(matrix((cellProp[,2]-matrix500[,4])/matrix500[,4]))
norm(matrix((cellProp[,2]-matrix500[,5])/matrix500[,5]))
norm(matrix((cellProp[,2]-matrix500[,6])/matrix500[,6]))
norm(matrix((cellProp[,3]-matrix500[,1])/matrix500[,1]))
norm(matrix((cellProp[,3]-matrix500[,2])/matrix500[,2]))
norm(matrix((cellProp[,3]-matrix500[,3])/matrix500[,3]))
norm(matrix((cellProp[,3]-matrix500[,4])/matrix500[,4]))
norm(matrix((cellProp[,3]-matrix500[,5])/matrix500[,5]))
norm(matrix((cellProp[,3]-matrix500[,6])/matrix500[,6]))
norm(matrix((cellProp[,4]-matrix500[,1])/matrix500[,1]))
norm(matrix((cellProp[,4]-matrix500[,2])/matrix500[,2]))
norm(matrix((cellProp[,4]-matrix500[,3])/matrix500[,3]))
norm(matrix((cellProp[,4]-matrix500[,4])/matrix500[,4]))
norm(matrix((cellProp[,4]-matrix500[,5])/matrix500[,5]))
norm(matrix((cellProp[,4]-matrix500[,6])/matrix500[,6]))
norm(matrix((cellProp[,5]-matrix500[,1])/matrix500[,1]))
norm(matrix((cellProp[,5]-matrix500[,2])/matrix500[,2]))
norm(matrix((cellProp[,5]-matrix500[,3])/matrix500[,3]))
norm(matrix((cellProp[,5]-matrix500[,4])/matrix500[,4]))
norm(matrix((cellProp[,5]-matrix500[,5])/matrix500[,5]))
norm(matrix((cellProp[,5]-matrix500[,6])/matrix500[,6]))
norm(matrix((cellProp[,6]-matrix500[,1])/matrix500[,1]))
norm(matrix((cellProp[,6]-matrix500[,2])/matrix500[,2]))
norm(matrix((cellProp[,6]-matrix500[,3])/matrix500[,3]))
norm(matrix((cellProp[,6]-matrix500[,4])/matrix500[,4]))
norm(matrix((cellProp[,6]-matrix500[,5])/matrix500[,5]))
norm(matrix((cellProp[,6]-matrix500[,6])/matrix500[,6]))

cellProp = cbind(cellProp[,6], cellProp[,2], cellProp[,1], cellProp[,4], cellProp[,5], cellProp[,3])

correlation = cor(as.vector(cellProp), as.vector(matrix500))
# correlation = 0.3811694


###################### PCA METHOD FOR THREE CELL TYPES ###########################

# ######################## GENERATING SIMULATED DATA #############################
# # Simulated data is included in data file.  This code is provided to demonstrate how 
# # simulated data was generated.  The refernce data set, FlowSorted.Blood.450k, was
# # downloaded with the FlowSorted.Blood.450k R package (source("http://bioconductor.org/
# # biocLite.R"), biocLite("FlowSorted.Blood.450k")).

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

# threeCellTypeBetas = cbind(avGran, avCD4, avCD8)

# headers = cbind("Gran", "CD4", "CD8")
# colnames(threeCellTypeBetas) = headers

# # Generating a random matrix of cell type proportions for 100 samples
# randMatrix3 = function(nSamp){
	# x = matrix(, ncol=3)
	# r = NULL
	# for (i in 1:nSamp){
		# r = sort(runif(2, 0, 1))
		# r = c(r[1], (r[2]-r[1]), (1-r[2]))
		# x = rbind(x, r)
	# }
	# x = x[2:(nSamp+1),]
	# return(x)
# }

# matrix500_3 = randMatrix(500)

# # Generating 500 samples
# genSamples3 = function(nBetas, nSamp, randMatrix, betas){
	# x = matrix(, nrow=nBetas, ncol=nSamp)
	# for (i in 1:nSamp){
		# x[,i] = (t(t((betas[,"Gran"])*randMatrix[i,1])) + t(t((betas[,"CD4"])*randMatrix[i,2])) + t(t((betas[,"CD8"])*randMatrix[i,3])))
	# }
	# return(x)
# }

# samples500_3 = genSamples3(50, 500, matrix500_3, threeCellTypeBetas)

# ################################# DETERMINE PCs ##################################
# res = prcomp(t(samples500_3), center = TRUE, scale = FALSE)

# cumsum(res$sdev^2/sum(res$sdev^2))

# pc.use = 2 # explains 100% of variance
# truncPC = res$x[,1:pc.use]

# #################################### SHIFT PCs ####################################
# # PCs are shifted to be positive. This step is for the ease of constructing the 
# # initial hyperplanes.

# ones = matrix(1, nrow=500, ncol=2)
# shiftTruncPC = truncPC + ones

# ############################ FUNCTION TO BE MINIMIZED #############################
# detFun3 = function(t){
	# A = matrix(t[1:6], nrow = 3)
	# b = -(matrix(t[7:9], nrow = 3))
	# A1 = A[1:2,]
	# b1 = b[1:2,]
	# A2 = A[2:3,]
	# b2 = b[2:3,]
	# A3 = rbind(A[1,], A[3,])
	# b3 = rbind(b[1,], b[3,])
	# vertices = cbind(mldivide(A1, b1), mldivide(A2, b2), mldivide(A3, b3))
	# X = rbind(vertices[,1], vertices[,2], vertices[,3])
	# cellProp = cart2bary(X, P) 
	# if (min(cellProp) < 0) area = 1000000 else area = abs(det(rbind(vertices, rep.int(1, 3))))
	# print(area)
	# return(area)
# }

# ########################### CONSTRAINTS FOR OPTIMIZATION #############################
# # ui %*% theta - ci >= 0, ak is contructed to ensure all points are on the 
# # same side of the hyperplane.

# a = t(shiftTruncPC)
# a = rbind(a, rep.int(1, 500))
# ak = kronecker(t(a), diag(3))

# ci = rep.int(0, 1500)

# # Global variable P is constructed to determine barycentric coordinates within detFun
# P = a[1:2,]
# P <<- t(P)

# ######################### PARAMETERS OF INITIAL HYPERPLANES ###########################
# t = c(1, 0, -1, 0, 1, -1, 0, 0, 6)

# ################# OPTIMIZATION TO FIND HYPERPLANES OF MINIMUM VOLUME ###################
# minDet = constrOptim(t, detFun3, grad=NULL, outer.eps = 1e-6, control=list(maxit=10000), ui=ak, ci=ci)

# ########################### FINDING VERTICES OF HYPERPLANES #############################
# detPoints3 = function(x){
	# A = matrix(x[1:6], nrow = 3)
	# b = -(matrix(x[7:9], nrow = 3))
	# A1 = A[1:2,]
	# b1 = b[1:2,]
	# A2 = A[2:3,]
	# b2 = b[2:3,]
	# A3 = rbind(A[1,], A[3,])
	# b3 = rbind(b[1,], b[3,])
	# vertices = cbind(mldivide(A1, b1), mldivide(A2, b2), mldivide(A3, b3))
	# return(vertices)
# }

# pnts = detPoints3(minDet$par)

# ############################ FINDING BARYCENTRIC COORDINATES ##############################
# X = rbind(pnts[,1], pnts[,2], pnts[,3])
# cellProp = cart2bary(X, P)

# ################# DETERMINING ORDER OF COLUMNS IN BARYCENTRIC COORDINATES ##################
# # The columns of the barycentric coordinates must be reordered to correspond to the appropriate
# # columns of the matrix used to generate the data.

# head(cellProp)
# head(matrix500_3)

# cellProp = cbind(cellProp[,3], cellProp[,2], cellProp[,1])

# correlation = cor(as.vector(cellProp), as.vector(matrix500_3))
# # correlation = 0.9999

# plot(cellProp[,1], matrix500_3[,1], pch=19, col=rgb(190,190,190,100,maxColorValue=255), xlim=c(0,.7), ylim=c(0,.7), main="Simulated data cell type proportions (3 cell types)", ylab="Actual proportions", xlab="Estimated proportions")
# points(cellProp[,2], matrix500_3[,2], pch=19, col=rgb(255,140,0,100,maxColorValue=255))
# points(cellProp[,3], matrix500_3[,3], pch=19, col=rgb(154,205,50,100,maxColorValue=255))
# abline(0,1, col="gray")