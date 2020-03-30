dpi <- 300

# Fig. 2
pdf(file = "Fig2.pdf",width = 6,height = 4)
png(file = "Fig2.png",width = 6*dpi,height = 4*dpi,res=dpi)
dates <- date_matching_cumulative_date-date_matching_i+seq(0,n_steps)
plot(dates,rep(0,n_steps+1),type="n",ylim=c(0,1),xlab="Date",ylab="Proportion",xaxt="n",xlim = as.Date(c("2020-03-16", "2020-07-1")))
rangeDates <- dates[dates >= as.Date(c("2020-03-16")) & dates <= as.Date(c("2020-07-01"))]
selectedDates <- rangeDates[format(rangeDates, "%d") %in% c("01", "10", "20")]
axis(1, selectedDates, labels = F , cex.axis = .7)
text(selectedDates, par("usr")[3]-0.07, 
     srt = 40, adj= 1, xpd = TRUE,
     labels = format(selectedDates, "%d %b '%y"), cex=0.75)
lines(dates,(state[,"R"]-cumsum(deaths))/pop_size,col="green",type="l")
lines(dates,rowSums(state[,IStages])/pop_size,col="red",type="l")
lines(dates,cumsum(deaths)/pop_size*1e1,col="black",type="l")
legend(x="topleft",
       legend = c("Recovered","Infected","Cum. mort. x 10"),
       col=c("green","red","black"),lty=c(1,1,1),bg=rgb(1,1,1,alpha = 0.7),cex=0.9)
dev.off()

# Fig. 3
for(s in 1:10){set.seed(s);source("Covid-19-Policy-Simulator.R"); simplify <- T;title(s)};rm("simplify")

pdf(file = "Fig3.pdf",width = 6.5,height = 6.5)
#png(file = "Fig3.png",width = 6.5*dpi,height = 6.5*dpi, res = dpi)
selection <- rev(c(7,10,9,8))
par(mar=c(4,4.1,0.2,2.1)) # default:  c(5.1, 4.1, 0, 2.1).
#par(mfrow=c(4,1))
layout(matrix(c(4,3,2,1),nrow = 4, ncol = 1, byrow = TRUE))
for(s in seq(selection)){
  set.seed(selection[s]);
  source("Covid-19-Policy-Simulator.R"); 
  mtext(paste0("(", letters[5-s], ")"), cex = 0.8, side = 3, adj = 0.05,line = -1.3);
  simplify <- T;
};
rm("simplify")
dev.off()

set.seed(selection[2]);
source("Covid-19-Policy-Simulator.R"); 
final_mortality

dev.off()

# Fig. 4
for(s in 1:30){set.seed(s);source("Covid-19-Policy-Simulator.R"); simplify <- T;title(s)};rm("simplify")

pdf(file = "Fig4.pdf",width = 6.5,height = 6.5)
png(file = "Fig4.png",width = 6.5*dpi,height = 6.5*dpi, res = dpi)
selection <- rev(c(27,30,4,26))
par(mar=c(4,4.1,0.2,2.1)) # default:  c(5.1, 4.1, 0, 2.1).
#par(mfrow=c(4,1))
layout(matrix(c(4,3,2,1),nrow = 4, ncol = 1, byrow = TRUE))
for(s in seq(selection)){
  set.seed(selection[s]);
  source("Covid-19-Policy-Simulator.R"); 
  mtext(paste0("(", letters[5-s], ")"), cex = 0.8, side = 3, adj = 0.05,line = -1.3);
  simplify <- T;
};
rm("simplify")
dev.off()

set.seed(selection[2]);
source("Covid-19-Policy-Simulator.R"); 
final_mortality

dev.off()

plot(final_immunization,final_mortality,cex=min(1,200/n_runs), xlab = "Population ")
vanilla.run <- data.frame(final_immunization,final_mortality)
mean(final_immunization/final_mortality < 1.1*min(final_immunization/final_mortality))
head(sort(final_mortality/final_immunization))
mean(final_mortality/final_immunization > 0.03 )

# Fig. 4
for(s in 1:40){set.seed(s);source("Covid-19-Policy-Simulator.R"); simplify <- T;title(s)};rm("simplify")
dev.off()


