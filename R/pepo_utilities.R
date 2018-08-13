calc_pepo <- function(t1, t2, t3, mu, model = 'strict') {

  pe <- pr_hemiplasy(t2, t1, t1, t1 + t2, t3, mu, mode = model)
  po <- pr_homoplasy(t2, t1, t1, t1 + t2, t3, mu, mode = model)
 #print(pe)
 print(po)
  pe / po
}


pepo_contourplot <- function(logmu, t2range, t1val, t3val, bins=21, model='strict'){
  if(length(logmu)!=2 | !(logmu[1] < logmu[2]) ) {stop('"logmu" is not c(lower, upper)')  }
  if(length(t2range)!=2 | !(t2range[1] < t2range[2]) ) {stop('"t2range" is not c(lower, upper)')  }
  muvec = seq(logmu[1],logmu[2], length.out = bins)
  t2vec = seq(t2range[1],t2range[2], length.out = bins)

  grill <- data_frame( t1 = t1val, t2= t2vec, t3 = t3val, mu = muvec)

  grillpt <- grill%>%complete(., t1, t2, t3, mu)%>%
    mutate(mu = 10^mu)%>%
    mutate(pepo = pmap_dbl(list(.$t1, .$t2, .$t3, .$mu), calc_pepo))
  grillpt
}

#tmp <- pepo_contourplot(logmu=c(-5,0), t2range=c(0.0001,4), t1val=1, t3val=5, bins=51) %>%
#  mutate(pepo = log10(pepo))%>%
#  mutate(catpepo = cut(pepo, breaks=c(-Inf, -2, -1, log10(1/2),log10(2), 1, 2, Inf)))

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

# ggplot(tmp,aes(x=t2, y=mu, z=pepo))+
#   geom_tile(aes(fill=catpepo))+
#   stat_contour(breaks=0, color='grey40')+
#   scale_y_log10(labels=scientific_10)+
#   scale_fill_carto_d(type = 'diverging', palette = "Geyser")+
#   theme_classic()+
#   labs(x=expression(italic('t')[2]), y=expression(mu))+
#   theme(legend.position = 'none',
#         axis.title.y = element_text(angle = 0, vjust=0.5))

