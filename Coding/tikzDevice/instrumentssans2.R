require(ggplot2)
require(tikzDevice)
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}",
                               "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}",
                               "\\usepackage{amssymb}"))
tikz('instrumentssans2.tex')
cols <- c("Estimate"=rgb(128,0,128,maxColorValue=255),"Confidence"=rgb(0,0.75,1))
BCE <- ACE[ACE$Exp<6,]
BCE <- data.frame(Exp = ordered(BCE$Exp),BCE[,2:5])
p <- ggplot(BCE, aes(xmin=cl,xmax=cu,ymin=Exp-0.4,ymax=Exp+0.4))
p + theme_bw() + 
  theme(text = element_text(family = "CM Sans"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        title = element_text(size = 14)) +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = 'black'),
          panel.background = element_blank(),
          axis.line = element_line(colour = 'black')) +
#   theme(legend.position="none") +
  theme(axis.title.y=element_text(vjust=4)) +
  theme(axis.title.x=element_text(vjust=-2)) +
  theme(plot.title = element_text(vjust=2)) +
  geom_rect(aes(xmin=cl,xmax=cu,ymin=as.numeric(Exp)-0.4,ymax=as.numeric(Exp)+0.4,fill="Confidence")) +
  geom_rect(aes(xmin=lower,xmax=upper,ymin=as.numeric(Exp)-0.4,ymax=as.numeric(Exp)+0.4,fill="Estimate")) +
  scale_fill_manual(name="",values=cols,labels=c("Confidence","Estimate")) +
    #   geom_vline(xintercept = 0) +
  scale_x_continuous(breaks = seq(-0.550,0.450,0.2)) +  
  scale_y_discrete(limits=c(1,2,3,4,5), labels = c("Same gender",                                            
                                               "Male-Male",
                                               "Female-Female",
                                               "Male-Male and \n Female-Female",
                                               "All permutations")) +
  xlab("$ACE_n(D\\rightarrow Y)$") +
  ylab("Definition of instrumental variable") +
  labs(title="Enriching the support of an instrumental variable")
dev.off()