require(tikzDevice)
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}",
                               "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}",
                               "\\usepackage{amssymb}"))
tikz('allin.tex')
p <- ggplot(df1, aes(xmin=cl0,xmax=cu0,ymin=cl1,ymax=cu1))
p + theme_bw() + 
  theme(text = element_text(family = "CM Roman"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black'),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black')) +
  theme(legend.position="none") +
  geom_rect(data=df2,aes(xmin=cl0,xmax=cu0,ymin=cl1,ymax=cu1),fill="red") +
#   geom_rect(data=df2,aes(xmin=u0,xmax=cu0,ymin=cl1,ymax=cu1),fill="red") +
#   geom_rect(data=df2,aes(xmin=l0,xmax=u0,ymin=cl1,ymax=l1),fill="red") +
#   geom_rect(data=df2,aes(xmin=l0,xmax=u0,ymin=u1,ymax=cu1),fill="red") +
  geom_line(aes(x=c(0.537,0.716),y=0.476),color="blue") +
  scale_x_continuous(breaks = round(seq(0.371, 0.723, by = 0.1),1)) +
  scale_y_continuous(breaks = round(seq(0.186, 0.777, by = 0.1),1)) +
  xlab("$\\mathbb{E}_n[Y(0)]$") +
  ylab("$\\mathbb{E}_n[Y(1)]$") 
dev.off()