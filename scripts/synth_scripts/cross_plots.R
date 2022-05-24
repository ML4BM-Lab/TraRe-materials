library(ggplot2)
library(grid)

#crossplots
cross_plots <- function(linker_output_synth,synth_modules,pathtosave,quadgraphTH){
  

  #linker_output, synth_mods, NrModules=10, TH=0.5, sankey=FALSE, quadgraphTH=0
  eval_linker_list <- eval_linker(linker_output_synth,synth_modules,
                                  TH=0.50,quadgraphTH=quadgraphTH)
  
  #Decompose list  
  confmat <- eval_linker_list[[1]]
  tp_fp_fn <- eval_linker_list[[2]]
  quad_graphs_c <- t(eval_linker_list[[3]])
  message('Initial rows are ', nrow(quad_graphs_c))
  
  clust_adjmatrix <- eval_linker_list[[4]]
  
  #Delete True Negatives
  TN <- which(apply(quad_graphs_c,1,function(x) if (x[1]==0 & x[2]==0){return(TRUE)}else{return(FALSE)}))
  
  quad_graphs_c <- quad_graphs_c[-TN,]
  message('Final rows are ', nrow(quad_graphs_c))
  
  #Get precision and recall
  precision <- tp_fp_fn['TP']/(tp_fp_fn['TP']+tp_fp_fn['FP'])
  recall <- tp_fp_fn['TP']/(tp_fp_fn['TP']+tp_fp_fn['FN'])
  
  #add the column of confusion matrix
  quad_graphs_c <- cbind(quad_graphs_c,confmat[confmat!='TN'])
  
  #change colnames
  colnames(quad_graphs_c) <- c('X','Y','ConfMat')
  
  #fit the linear model
  fit_model <- lm(as.numeric(quad_graphs_c[,2])~as.numeric(quad_graphs_c[,1]))
  fit_model_sum <- summary(fit_model)
  
  #plot the quads
  r2adj <- round(fit_model_sum$adj.r.squared,2)
  r2adj_title <- paste0("R2 adjusted: ",r2adj)
  
  prec_recall <- paste0(" Precision: ",round(precision,2)," Recall: ",round(recall,2))
  
  pdf(paste0(getwd(),pathtosave))
  # Basic scatter plot
  scatter <- ggplot(as.data.frame(quad_graphs_c), aes(x=as.numeric(X), y=as.numeric(Y),color=ConfMat)) + 
    geom_point(size=2) + 
    geom_abline(intercept = fit_model$coefficients[1], slope = fit_model$coefficients[2])+
    labs(x='Synthetic driver weights', y = 'Inferred driver weights') +
    theme_classic() +
    ggtitle(paste0(r2adj_title,prec_recall))+
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = hues::iwanthue(length(unique(quad_graphs_c[,3])),hmax=120,cmin = 30,cmax = 80,lmin=35,lmax=80))
  print(scatter)
  dev.off()
  
  return(c(r2adj,precision,recall))
  #return(list(adjmat=clust_adjmatrix,quad_graph=quad_graphs_c, metrics = tp_fp_fn,scatter =scatter))
  
}

#bar plot associated
bar_plot <- function(values,colors=c("#4aac8d","#6a7fce","#6ba647")){
  x <- c('VBSR','VBSR','VBSR','LASSO','LASSO','LASSO','LM','LM','LM')
  y <- c('R2','Precision','Recall','R2','Precision','Recall','R2','Precision','Recall')
  
  df <- data.frame(x,y,values)
  colnames(df) <- c('Models','Metrics','Score')
  
  p <- ggplot(data=df, aes(x=Metrics, y=Score, fill=Models)) +
    geom_bar(stat="identity", color="black", position=position_dodge())+
    theme_classic()+
    scale_fill_manual(name='Models',
                      values = colors)
  
  pdf(file=paste0(getwd(),'/output/figures/crossplot_barplot.pdf'))
  print(p)
  dev.off()
}


#Synthetic modules
synth_modules <- readRDS(paste0(getwd(),'/input_data/synth_modules_0.01_to_100.rds'))

#Load inferred 
#define the metric array
metrics <- c()
for (model in c('VBSR','LASSOmin','LM')){
  
  
  #Read the model
  linker_output_synth <- readRDS(paste0(getwd(),'/input_data/models_noise_variance_1.2/',model,'_noise_variance_1.2.rds'))
  
  
  #Generate cross_plots
  cross_p_list <- cross_plots(linker_output_synth=linker_output_synth,
                              synth_modules=synth_modules$noise_1.2,
                              pathtosave=paste0('/output/figures/crossplot_',model,'.pdf'),
                              quadgraphTH=0.5)
  #add the metrics
  metrics <- c(metrics,as.numeric(cross_p_list))
  
}

#Call the barplot
bar_plot(values=round(metrics,2),colors=c("#5f92db",
                                          "#36a7aa",
                                          "#355a95"))
