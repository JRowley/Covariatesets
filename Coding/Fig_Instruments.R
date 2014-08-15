p <- ggplot(ACE, aes(xmin=cl,xmax=cu,ymin=Exp-0.4,ymax=Exp+0.4))
q <- p + theme_bw() + 
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
  geom_rect(aes(xmin=cl,xmax=lower,ymin=Exp-0.4,ymax=Exp+0.4),fill="red",alpha=0.5) +
  geom_rect(aes(xmin=lower,xmax=upper,ymin=Exp-0.4,ymax=Exp+0.4),fill="blue",alpha=0.5) +
  geom_rect(aes(xmin=upper,xmax=cu,ymin=Exp-0.4,ymax=Exp+0.4),fill="red",alpha=0.5) +
  #   geom_vline(xintercept = 0) +
  scale_x_continuous(breaks = round(seq(min(ACE$lower), max(ACE$upper), by = 0.1),1)) +
  scale_y_continuous(breaks = round(seq(min(ACE$Exp), max(ACE$Exp), by = 1),1)) +
  xlab(expression(ACE[n]*(D%->%Y)))
ggsave("C:/Users/Jeffro/Documents/GitHub/Covariatesets/Diagrams/Instruments.pdf", q)
embed_fonts("C:/Users/Jeffro/Documents/GitHub/Covariatesets/Diagrams/Instruments.pdf")