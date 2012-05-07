# clusterR using mmap
# Original code by Robert Hijimans, mmap integration by Jonathan Greenberg
#' @export

calc_hpc <- function(x, fun, args=NULL, filename='', cl=NULL, m=2, disable_cl=FALSE,verbose=FALSE,...) {
	require("raster")
	require("snow")
	require("mmap")
	require("ff")
	
	# Do some file checks up here.
	
	if(disable_cl)
	{
		nodes <- 1
	} else
	{
		if (is.null(cl)) {
			cl <- getCluster()
			on.exit( returnCluster() )
		}
		nodes <- length(cl)
	}

	# We should test a single pixel here to see the size of the output...
	
	# We are going to pull out the first row and first two pixels to check the function...
	r_check <- crop(x, extent(x, r1=1, r2=1, c1=1,c2=2))
	
	if(!is.null(args)) {
		r_check_function <- getValues(do.call(fun, c(r_check, args)))
	} else
	{
		r_check_function <- getValues(fun(r_check)) 
	}
	
	if(class(r_check_function)=="numeric")
	{
		outbands=1
	} else
	{
		outbands=dim(r_check_function)[2]
	}
	
	
	outdata_ncells=nrow(x)*ncol(x)*outbands
	if(filename=="")
	{	
		filename <- tempfile()
	} 
	
	# How about using ff?
	out<-ff(vmode="double",length=outdata_ncells,filename=filename)
	finalizer(out) <- close
	close(out)
	
	m <- max(1, round(m))
	tr <- blockSize(x, minblocks=nodes*m )
	if (tr$n < nodes) {
		nodes <- tr$n
	}
	
	tr$row2 <- tr$row + tr$nrows - 1

	tr_out=list(row=((tr$row-1)*outbands+1))
	tr_out$row2=((tr$row2)*outbands)
	
	i=1:tr$n
	
	if(disable_cl)
	# Use only for debugging.
	{
		mapply(function(i,fun,args,x,tr,filename,outbands) 
			{
				r <- crop(x, extent(x, r1=tr$row[i], r2=tr$row2[i], c1=1, c2=ncol(x)))
				if(!is.null(args)) {
					r <- fun(r) 
				} else
				{
					r <- do.call(fun, c(r, args))
				}
				out <- mmap(filename,mode=real64())
				cellStart=((cellFromRowCol(x,row=tr$row[i],col=1))-1)*outbands+1
				cellEnd=((cellFromRowCol(x,row=tr$row2[i],col=ncol(x))))*outbands
				out[cellStart:cellEnd] <- as.vector(t(getValues(r)))
#				out[cellFromRow(x,tr$row[i]:tr$row2[i])] <- as.vector(getValues(r))
				munmap(out)
				return(NULL)
			},
			i,MoreArgs=list(fun=fun,x=x,tr=tr,args=args,filename=filename,outbands=outbands))

	} else
	{
		clusterMap(cl,function(fun,i,args,x,tr,filename,outbands) 
			{
				r <- crop(x, extent(x, r1=tr$row[i], r2=tr$row2[i], c1=1, c2=ncol(x)))
				if(!is.null(args)) {
					r <- fun(r) 
				} else
				{
					r <- do.call(fun, c(r, args))
				}
				out <- mmap(filename,mode=real64())
				cellStart=((cellFromRowCol(x,row=tr$row[i],col=1))-1)*outbands+1
				cellEnd=((cellFromRowCol(x,row=tr$row2[i],col=ncol(x))))*outbands
				# Disable transpose for BIL?
				out[cellStart:cellEnd] <- as.vector(t(getValues(r)))
#				out[cellFromRow(x,tr$row[i]:tr$row2[i])] <- as.vector(getValues(r))
				munmap(out)
				return(NULL)
			},
			i,MoreArgs=list(fun=fun,x=x,tr=tr,args=args,filename=filename,outbands=outbands))
	}
		
	# Let's see if we can trick raster into making us a proper header...
	if(outbands > 1) 
	{ 
		reference_raster=brick(raster(x,layer=1),nl=outbands) 
	} else
	{
		if(nlayers(x) > 1) { reference_raster=raster(x,layer=1) } else
		{ reference_raster=x }	
	}
	outraster_base <- writeStart(reference_raster,filename=paste(filename,".grd",sep=""),datatype="FLT8S",bandorder="BIP",...)
	suppressWarnings(outraster_base <- writeStop(outraster_base))
	file.remove(paste(filename,".gri",sep=""))
	file.rename(filename,paste(filename,".gri",sep=""))
	outraster=brick(paste(filename,".grd",sep=""))
	return(outraster)
}

