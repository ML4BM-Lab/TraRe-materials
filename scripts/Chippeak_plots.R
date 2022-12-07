library(dplyr)
library(xlsx)
library(ggplot2)


args = commandArgs(trailingOnly=TRUE)
method <- args[1]

stats <- read.xlsx(paste0('~/output/tables/Supplementary_data_ChIPseq_',method,'.xlsx'),sheetIndex = 1)
if (method!="TraRe"){
    stats <- stats %>% filter(GRNinChip>=5)
    if(method=="ARACNE"){
        # because regulons are more diverse, the plot looks nicer if we orther for the number of edges
         stats <- stats %>% mutate(percentage=GRNinChip/(GRNinChip+GRNNoChip))  %>% 
                       mutate(Totalregulongenes=GRNinChip+GRNNoChip)%>% 
                       arrange(desc(Totalregulongenes))
    }else{stats <- stats %>% mutate(percentage=GRNinChip/(GRNinChip+GRNNoChip)) %>% arrange(desc(percentage))} # for GRNboost it is easier to plot by percentage
    stats$NA. <- factor(stats$NA.,levels=stats$NA.)
    stats$text_position <- stats$GRNinChip-10
    stats$text_position[seq(1,nrow(stats),2)] <- stats$GRNinChip[seq(1,nrow(stats),2)]+20
    stats$text_position <- ifelse(stats$text_position<0, 0,stats$text_position)
    stats$point <- ifelse(stats$signif==1,1,NA)
    mylevels=c("GRNNoChip","GRNinChip")
}else{
    stats <- stats %>% filter(TraReinChip>=5) # filter for chip coincidence > 5
    # because regulons are more diverse, the plot looks nicer if we orther for the number of edges
    stats <- stats %>% mutate(percentage=TraReinChip/(TraReinChip+TraReNoChip))  %>% 
                       mutate(Totalregulongenes=TraReinChip+TraReNoChip)%>% 
                       arrange(desc(Totalregulongenes))
    stats$NA. <- factor(stats$NA.,levels=stats$NA.)
    stats$text_position <- stats$TraReinChip-20
    stats$text_position[seq(1,nrow(stats),2)] <- stats$TraReinChip[seq(1,nrow(stats),2)]+20
    stats$text_position <- ifelse(stats$text_position<0, 0,stats$text_position)
    stats$point <- ifelse(stats$signif==1,1,NA)
    mylevels=c("TraReNoChip","TraReinChip")
}

sig_tfs <- as.character(stats$NA.[stats$signif==1])

targets_with_Chip <- read.xlsx(paste0('~/output/tables/Supplementary_data_ChIPseq_',method,'.xlsx'),sheetIndex = 2)
targets_with_Chip <- targets_with_Chip  %>% filter(TF%in%as.character(unique(stats$NA.))) # remove those filtered TF as well
targets_with_Chip$TF <- factor(targets_with_Chip$TF,levels=stats$NA.)
targets_with_Chip$sig <- factor(ifelse(targets_with_Chip$TF%in%sig_tfs,1,0))


pdf(paste0("~/output/figures/",method,"_regulon_ChIPseq_validation.pdf"))
ggplot()+
geom_bar(data=targets_with_Chip,aes(x=TF,fill=factor(ChIPstatus_str, levels=mylevels)),width = 0.9)+
theme_bw()+
theme(axis.text.x=element_text(angle = 60, hjust=1))+
scale_fill_brewer("Set4")+
ylab("Number of regulon targets")+xlab("")+
guides(fill=guide_legend(title="ChIP status"))+
geom_text(data=stats,size=2.5, aes(label = paste0(round(percentage,3)*100,"%"),x=NA.,y=  text_position),
         stat= "identity", vjust = -.5)+
geom_point(data=stats,aes(x=NA.,y=point),shape=8)
dev.off()


