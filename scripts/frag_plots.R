#!/usr/bin/R

pdf("simPlots.pdf")

args<-commandArgs(TRUE)

if (length(args)==0){
  stop("Please provide file name\n", call.=FALSE)
} else if (length(args)>2){
  stop("This script only takes one argument (input file).\n");
}

print(paste("Reading from file: ", args[1]))

table<-read.table(file=args[1], header=T);

table<-read.table(file="sim.tsv",header=T);
print("Creating plots...");

impulsePlot <- function(x) {
	#Impulse plot with automatic axes
	plot(table[[x]],table$Fragment_length,type="h", main=x,
	xlab="Fragment length",ylab="Frequency")
		
	#Impulse plot xlim=100:5000
	plot(table[[x]],table$Fragment_length,type="h", main=x,
	xlab="Fragment length",ylab="Frequency",xlim=range(0:5000))

	#Impulse plot xlim=100:5000
	plot(table[[x]],table$Fragment_length,type="h", main=x,
	xlab="Fragment length",ylab="Frequency",xlim=range(0:1000))
	
	#Impulse plot xlim=100:5000
	plot(table[[x]],table$Fragment_length,type="h", main=x,
	xlab="Fragment length",ylab="Frequency",xlim=range(0:500))
	
}	

num<-0;
splines<-list();
preds<-list();
titles<-list();
ylen<-length(table$Sum);
for(i in names(table)){

	if (i != "Fragment_length"){
		if (i != "Sum"){
			impulsePlot(i);	
		}
		if (length(unique(table[[i]])) > 4){
			if (i == "Sum"){
				sumSpline<- loess(table[[i]]~table$Fragment_length, span=0.02)		
				sumPred <- predict(sumSpline)
			}else{
				num<-num+1;
				titles[[num]]<-i;
				splines[[num]] <-  loess(table[[i]]~table$Fragment_length, span=0.02)			
				preds[[num]] <- predict(splines[[num]])
			}
		}else{
			print(paste("too few values for", i, ", leaving out of final plot"));
		}
	}
}

colors<-rainbow(num);
t<-unlist(titles);

plot(sumPred,type="l",main="Total recovered fragments",
xlab="Fragment size", ylab="Frequency",lwd=1.5)
for (i in 1:num){
	lines(preds[[i]],col=colors[[i]],lwd=1.5)
}
legend("topright", c("Total fragments", t), col=c("black",colors),lwd=1.5)

plot(sumPred,type="l",main="Total recovered fragments",
xlab="Fragment size", ylab="Frequency",lwd=1.5, xlim=c(0,1000))
for (i in 1:num){
	lines(preds[[i]],col=colors[[i]],lwd=1.5)
}
legend("topright", c("Total fragments", t), col=c("black",colors),lwd=1.5)
