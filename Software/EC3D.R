#'dir.in' is the directory containing images/shapes (e.g. patient MRIs).

#'out.file' specifies where the final list of Euler Characteristic (EC) matrices should be saved. Note that the user is not specifying a directory, but the file name directly.

#'img.dir' determines the subfolders where the images are loacated. This can be a pattern.

#'stepsize' is the number of sublevel sets to be taken during the filtration of each image.

#'rotstep' is the number of directions to be taken.

#'first.only' refers to the issue of multiple unconnected shapes or regions in the image (e.g. multifocal tumors with mutliple masses). If TRUE, the function will only compute the EC for the 'first' identified shape.

#'verbose' is a boolean variable that prints out the individual image file names. If FALSE, the code will print the index of the patient and the patient identifier.

ecf <- function(in.dir,out.file,img.dir=NULL,stepsize=100,rotstep=72,first.only=TRUE,verbose=FALSE){

#This script is for .png images (e.g. slices of MRI Scans)
#Assumptions: The input image is assumed to have no background noise.
#Output #1: EC is a cell, EC{i} is the *integrated* EC curve for ith direction.
#Output #2: Shapes is the resulting structual array where:
#.Name = Name of the shape and the slice number (e.g. apple_1)
#.EC = A matrix of EC curves for that shape and slice number (#columns = #rotations)
# R version, converted from EC3D.m

# Load packages
usePackage("png")
usePackage("lattice")

    if (strsplit(in.dir,"")[[1]][length(strsplit(in.dir,"")[[1]])] != "/" ) {
        in.dir <- paste(in.dir,"/",sep='')
    }

#Load functions to compute Euler Characteristics
#source("Code/ECFunctions_mao/mao_ECscale.R")
#source("Code/ECFunctions_mao/mao_ECTmany.R")
#source("Code/ECFunctions_mao/mao_ECTone.R")

# Get the current working directory, useful for resetting
startdir <- getwd()

# Set the stepsize for calculating EC
#stepsize <- 100

# Set save directory for results
# Assume that the code is being executed somehwere containing a folder called
# "Images" in which the subject specific images are stored. The output will be 
# stored in a folder called "Output"

#out.dir = paste(startdir,"/Output",sep="")
#if ( file.exists(out.dir) != TRUE ) dir.create(out.dir)

# Set working directory to Image folder, get a vector containing all
# patient identifiers

setwd(in.dir)
TCGA_patients = dir()
n.patients <- length(TCGA_patients)
setwd("../")

# Set up desired number of rotations
#rotstep = 72
theta = seq(from= -pi,to=pi,by=2*pi/rotstep)
d1 = cos(theta);d2=sin(theta)
d= rbind(d1,d2)

MRI_list <- list(NULL)
    for ( k in 1:n.patients) {
        print(paste("K:",k))
        print(TCGA_patients[k])

        MRI_list[[k]] <- list(name = TCGA_patients[k])
#    MRI_list[[k]]$name <- TCGA_patients[k]
        MRI_list[[k]]$EC <- NULL
    # Go grab images for patient k
        subdir = paste(startdir,"/",in.dir,"/",TCGA_patients[k],"/",img.dir,sep="")
        setwd(subdir)
    # Some patients have no images, hence checking whether images are present
        if ( length(dir(pattern='.png')) > 0 ) {
            slices <- dir(pattern="*.png")
            EC = matrix(0,stepsize+1,dim(d)[2]-1)
            ec.slices <- array(0,dim=c(dim(EC),length(dir(pattern='.png'))))
            realslice = 0
            
        # Compute the EC for all slices of partient k
        # j goes through the images/slices
            for ( j in 1:length(slices) ) {
                if(verbose == TRUE ) print(slices[j])
                I = readPNG(slices[j])
            # Only use images that have non-zero regions
                if ( sum(I) > 10 ) {
                    BW = round(I)
                    components <- connected.components(BW)
                    c.unique <- sort(unique(c(components)))[-1]
                    
                    if ( first.only == TRUE ) {
                        ones = which(BW==1,arr.ind=TRUE)
                        for(q in 1:nrow(ones)){
                          if(abs(ones[q,2]-ones[q+1,2])<2&abs(ones[q,1]-ones[q+1,1])<2)
                            break
                        }
                        start.c=ones[q,]
                        Z <- boundary.trace(BW,start.c)
                        Z <- Z[dim(Z)[1]:1,]
                        E <- cbind(1:dim(Z)[1],c(2:dim(Z)[1],1))
                        complex <- list(NULL)
                        complex$V=Z
                        complex$E=E
                        complex$F=NULL
                        complex$T=NULL
                        Z=t(t(Z) - apply(Z,2,mean))
                        Z=Z/norm(Z,type='F')
                        ec <- matrix(0,stepsize+1,dim(d)[2]-1)
                        for ( i in 2:dim(d)[2] ) {
                            fun = Z%*%as.matrix(d[,i])
                            C = gEuler(complex,fun,100)
                            ec[,i-1] <- C[,2]
                            ec.slices[,i-1,j] <- C[,2]
                        }
                    } else {
                        ec <- matrix(0,stepsize+1,dim(d)[2]-1)
                        for ( comp.ind in c.unique ) {
                            ones = which(components==comp.ind,arr.ind=TRUE)
                            if(nrow(ones)>5){
                              for(q in 1:nrow(ones)){
                                if(abs(ones[q,2]-ones[q+1,2])<2&abs(ones[q,1]-ones[q+1,1])<2)
                                  break
                              }
                              start.c=ones[q,]
                              
                              Z <- boundary.trace(BW,start.c)
                              Z <- Z[dim(Z)[1]:1,]
                              E <- cbind(1:dim(Z)[1],c(2:dim(Z)[1],1))
                              complex <- list(NULL)
                              complex$V=Z
                              complex$E=E
                              complex$F=NULL
                              complex$T=NULL
                              Z=t(t(Z) - apply(Z,2,mean))
                              Z=round(Z/norm(Z,type='F'),4)
                              for ( i in 2:dim(d)[2] ) {
                                fun = Z%*%as.matrix(d[,i])
                                C = gEuler(complex,fun,100)
                                ec[,i-1] <- ec[,i-1] +  C[,2]
                              }
                            }
                        }
                    }
                # End of angle loop over i
                    EC = EC + ec #Integration step 1, over different angles
                    realslice = realslice + 1                
                }
            # End of non-zero image if            
            }
        # End of slice loop over j for patient k
            MRI_list[[k]]$EC <- EC/realslice
        }
    # End of non-zero slice if
    }
# End of patient loop over k

    setwd(startdir)
    save(file=out.file,MRI_list)
    setwd(startdir)
    
}

boundary.trace <- function(bw,start){

    current.position <- start
    pos.list <- matrix(0,1,2)
    pos.list[1,] <- current.position
    current.position <- boundary.move(bw,current.position,1)
    pos.list <- rbind(pos.list,current.position[1:2])
    while ( sum(current.position[1:2]==start) != 2 ) {
        current.position <- boundary.move(bw,current.position[1:2],current.position[3])
        pos.list <- rbind(pos.list,current.position[1:2])
    }
    return(pos.list)
}

# This will include diagonal movement
# direction of movement is ordered as
# 1 - right
# 2 - downright
# 3 - down
# 4 - downleft
# 5 - left
# 6 - upleft
# 7 - up
# 8 - upright

boundary.move <- function(bw,pos,direction){

    # check all directly adjacent cells
    # move according to priorities right, followed by clockwise
    c.state <- bw[(pos[1]-1):(pos[1]+1),(pos[2]-1):(pos[2]+1)]
#    print(c.state)
    changes <- matrix(0,8,2)
    changes[1,] <- c(0,1)
    changes[2,] <- c(1,1)
    changes[3,] <- c(1,0)
    changes[4,] <- c(1,-1)
    changes[5,] <- c(0,-1)
    changes[6,] <- c(-1,-1)
    changes[7,] <- c(-1,0)
    changes[8,] <- c(-1,1)
    seq.vec <- c(1:8,1:8,1:8)
    
    for ( j in 1:8 ) {
        if ( direction == j ) {
            for ( i in seq.vec[(9+j-3):((9+j-3)+8)] ) {
                if ( bw[pos[1]+changes[i,1],pos[2]+changes[i,2]] == 1 ) {
                    return(c(pos[1]+changes[i,1],pos[2]+changes[i,2],i))
                }
            }
        }
    }
}


connected.components <- function(BW){

    lab.mat <- matrix(0,dim(BW)[1],dim(BW)[2])
    new.lab <- 1
    equiv.list <- list(NULL)
    equiv.list[[1]] <- 1
    
    for ( i in 1:dim(BW)[1] ) {        
        for ( j in 1:dim(BW)[2] ) {
            if (BW[i,j] == 1) {

                i.range <- 0:2-1
                j.range <- 0:2-1
                if ( i==1 ) i.range <- 0:1
                if ( i==dim(BW)[1] ) i.range <- 0:1-1
                if ( j==1 ) j.range <- 0:1
                if ( j==dim(BW)[2] ) j.range <- 0:1-1
#                print(i + i.range)
#                print(j + j.range)
                current.neighbours <- BW[i+i.range,j+j.range]
                current.neighbours[2,2] <- 0
                if ( sum(current.neighbours) > 0 ) {

                    tmp.n <- which(current.neighbours == 1)
                    tmp.labs <- unique(lab.mat[i+i.range,j+j.range][tmp.n])
                    tmp.labs <- tmp.labs[which(tmp.labs>0)]
#                    print(tmp.labs)
                    if ( length(tmp.labs) > 0 ) {
                        lab.mat[i,j] <- tmp.labs[1]                                            
                        for ( k in tmp.labs ) {
                            equiv.list[[k]] <- unique(c(equiv.list[[k]],
                                                        unlist(equiv.list[tmp.labs])))
                        }
                    } else {
                        lab.mat[i,j] <- new.lab
                        equiv.list[[new.lab]] <- new.lab
                        new.lab <- new.lab+1
                    }
                    
                } else {

                    lab.mat[i,j] <- new.lab
                    equiv.list[[new.lab]] <- new.lab
                    new.lab <- new.lab+1
                }

                
                
            } else {

            }
        }
    }

    equiv.list <- lapply(equiv.list,sort)


    for ( i in 1:length(equiv.list) ) {
        
        lab.mat[which(lab.mat==i)] <- equiv.list[[i]][1]
        
    }
    
    
#    print(equiv.list)
    return(lab.mat)
}

#Input: complex: is a structure of finite dimensional simplicial complex (here is for 3 dim).
#       complex.V is n1x3 vertices
#       complex.E is n2x2 edges
#       complex.F is n3x3 faces
#       complex.T=[]; for 2D case.
#       fun: The function value on the vertices
#       stepsize: discrete Euler Characteristic curve into such dimensional
#       vector. (e.g. 100);
#Output: stepsize <- by <- 2 matrix. Euler Characteristics (second column) for scales of function values (first).

gEuler <- function(complex,fun,stepsize){

    # The thresholding rounding is in there to make sure this produces
    # the exact same results as the matlab code, which rounds to 4 
    # decimal places
    
    V = complex$V
    E = complex$E
    F = complex$F
    T = complex$T

    if ( dim(V)[1] != length(fun) ) {
        print('The size of function should be same as the number of vertices')
        return(0)
    }

    fe = sapply(1:dim(V)[1],function(ind) max(fun[E[ind,]]))
    if ( length(F) > 0 ) {
        ff = sapply(1:dim(V)[1],function(ind) max(fun[F[ind,]]))
    } else {
        ff = NULL
    }
    
    if ( length(F) > 0 ) {
        ft = sapply(1:dim(V)[1],function(ind) max(fun[T[ind,]]))
    } else {
        ft = NULL
    }
    threshold = seq(from=min(fun),to=max(fun),length.out=stepsize+1)
    kai = matrix(0,stepsize+1,2)
    kai[,1] <- threshold
    for ( i in 1:length(threshold) ) {
        v = length(which(fun<=threshold[i]))
        e = length(which(fe<=threshold[i]))
        f = length(which(ff<=threshold[i]))      
        t = length(which(ft<=threshold[i]))
        kai[i,2] <- v - e + f - t;
    }
    return(kai)
}

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
