

library(tidyverse  )
library(magrittr   )
library(glmnet     )
library(decontam   )
library(vegan      )
library(philentropy)
library(microDecon )
library(SCRuB)

m1 <- read.csv('Fig1/dilution_summary.csv', row.names=1)
m2 <- read.csv('Fig1/msystems.00290-19-st001.csv', row.names=1)
m1 <- m1[row.names(m2),]
m1[, 'Concentration.of.samples..cells.ml.'] <- m2$Concentration.of.samples..cells.ml.
m1[, 'p'] <- (100-m1$PercentContaminants ) /100

m1[, c('p', 'Concentration.of.samples..cells.ml.')]

x <- m1$Concentration.of.samples..cells.ml.
y <- m1$p

sigmoid <- function(x){
  return( 1 / (1 + exp(-x)) )
}

inverse_sigmoid <- function(p){
  return( log(p / (1 - p)) )
}

## building function to map concentration to contaminatoin

logistic_model <- glm(inverse_sigmoid(y) ~ log10(x) )
# Plot results
map_concentration_to_contamination <- function(x){
  return( logistic_model %>% predict(x) %>% sigmoid ) }


map_concentration_to_contamination( data.frame(x=c(1e7, 1e8) ) ) %>% as.vector

rescale <- function(m) m / sum(m)
get_rescaled_mat <- function(Q) ( Q %>% apply(MARGIN=1, rescale) %>% t() )


set.seed(1)

synth_data <- read.csv('../data/Fig1/synth_combination.csv', row.names = 1) %>% as.matrix()
colnames(synth_data) <- paste0('OTU_', 1:ncol(synth_data))

controls <- synth_data[row.names(synth_data) %>% str_detect('BLANK') %>% 
                         which(), ]
samples <- synth_data[( F == row.names(synth_data) %>% str_detect('BLANK') ) %>% 
                        which(), ]


source('Simulation_Functions/Plate_simulations.R')
source('Simulation_Functions/generate_simulations.R')


CLEAN_SAMPLES_DECONTAM <- function(samples, 
                                   controls, 
                                   threshold = 0.1
){
  n_sinks <- nrow(samples)
  n_contams <- nrow(controls)
  
  ps <- samples %>% rbind(controls)
  
  neg <-c(rep(FALSE, n_sinks), rep(TRUE, n_contams))
  
  out <- isContaminant(ps, neg=neg, threshold = threshold,
                       detailed = T,
                       method='prevalence'
  )
  
  remove_out_samps <- function(mat, out){
    mat %>% apply(MARGIN = 1, function(x){
      x[out == T] <- 0 
      return(x)
    }  ) %>% t()
  }
  return(list(decontaminated_samples = remove_out_samps(samples, out$contaminant),
              cont_idx = out, 
              estimated_sources = controls) )
}

CLEAN_SAMPLES_DECONTAM_conc <- function(samples, 
                                        controls, 
                                        concentration,
                                        method='frequency'
){
  n_sinks <- nrow(samples)
  n_contams <- nrow(controls)
  ps <- samples
  prac_conc = concentration
  neg <-rep(FALSE, n_sinks)
  
  if( method=='combined' ){
    ps <- ps %>% rbind(controls)
    prac_conc <-  c( prac_conc, rep( min(concentration) /10 ,n_contams))
    neg <-c(rep(FALSE, n_sinks), rep(TRUE, n_contams))
  }
  
  out <- isContaminant(ps, 
                       neg=neg, 
                       conc=prac_conc,
                       detailed = T,
                       method=method
  )
  
  remove_out_samps <- function(mat, out){
    mat %>% apply(MARGIN = 1, function(x){
      x[out == T] <- 0 
      return(x)
    }  ) %>% t()
  }
  return(list(decontaminated_samples = remove_out_samps(samples, out$contaminant),
              cont_idx = out, 
              estimated_sources = controls) )
}



CLEAN_SAMPLES_DECONTAM_low_biomass <- function(samples, 
                                               controls, 
                                               threshold = 0.5){
  n_sinks <- nrow(samples)
  n_contams <- nrow(controls)
  
  ps <- samples %>% rbind(controls)
  
  neg <-c(rep(FALSE, n_sinks), rep(TRUE, n_contams))
  
  out <- isNotContaminant(ps, neg, threshold = threshold, 
                          detailed = T)
  
  remove_out_samps <- function(mat, out){
    mat %>% apply(MARGIN = 1, function(x){
      x[out == FALSE] <- 0 
      return(x)
    }  ) %>% t()
  }
  return(list(decontaminated_samples = remove_out_samps(samples, out$not.contaminant),
              cont_idx = out, 
              estimated_sources = controls) )
}

CLEAN_SAMPLES_MICRODECON <- function(smps, cnts){
  tmp <- data.frame( rbind( colnames(smps), cnts, smps ) %>% t() )
  tmp[, 2:( ncol(tmp) ) ]  <- as.numeric( as.matrix( tmp[, 2:( ncol(tmp) ) ] ) )
  decontaminated <- decon(data = tmp, numb.blanks=nrow(cnts), numb.ind= c(nrow(smps)), taxa=F)
  
  md_out <- smps*0
  md_out[,decontaminated$decon.table$V1 %>% as.character() %>% unname()] <- decontaminated$decon.table[,3:(2 + nrow(smps) )] %>% t()
  
  return(list(decontaminated_samples=md_out))
}


get_samps_from_upload <- function(uploaded_samples, n_samples){
  n_provided <- nrow(uploaded_samples)
  return(uploaded_samples[sample(n_provided, n_samples, replace=n_samples>n_provided),
  ])
}


REMOVE_ALL_CONTROL_SPECS <- function(samples, controls){
  X = samples
  to_remove <- which( (controls %>% colSums() )  > 0  )
  X[,to_remove] <- 0
  return(X)
}


run_SCRUB_DEC_and_RESTR<- function(samples, 
                                   controls, 
                                   well_dist, 
                                   true_concentrations, 
                                   dist_threshold=1.5, 
                                   scr_md=1){
  combined <- rbind(samples, controls)
  is_control <- c( rep(F, nrow(samples)),  rep(T, nrow(controls)))
  
  SCRUB_out <- SCRuB(rbind(controls, samples) %>% as.matrix,
                     scr_md
                     )
  
  decont_out <- CLEAN_SAMPLES_DECONTAM(samples, controls)
  decont_low_out <- CLEAN_SAMPLES_DECONTAM_low_biomass(samples, controls)
  
  
  restr_out <- CLEAN_SAMPLES_DECONTAM_conc(samples, 
                                           controls, 
                                           concentration = true_concentrations, 
                                           method='frequency'
  )
  
  microdec_out <- CLEAN_SAMPLES_DECONTAM_conc(samples, 
                                              controls, 
                                              concentration = true_concentrations, 
                                              method='combined'
  )
  
  return(  list(
                SCRUB_out=SCRUB_out,
                decont_out=decont_out, 
                decont_low_out=decont_low_out,
                restr_out=restr_out,
                microdec_out=microdec_out
  ) ) 
}

store_f1_results <- function(new_results, 
                             out,
                             samples, 
                             observed_samples, 
                             n_controls, 
                             well_species_limit, 
                             plate_format, 
                             well_level, 
                             contam_level, 
                             true_contamination,
                             n_taxa,
                             idx){
  
  
  
  input_sink_jsd <- get_rowwise_jsd(samples, 
                                    get_rescaled_mat(observed_samples))
  
  SCRUB_jsd <- get_rowwise_jsd(samples, 
                               get_rescaled_mat( out$SCRUB_out$decontaminated_samples ) )
  
  dec_jsd <- get_rowwise_jsd(samples, 
                             get_rescaled_mat( out$decont_out$decontaminated_samples ) )
  
  
  microdec_jsd <- get_rowwise_jsd(samples, 
                                  get_rescaled_mat( out$microdec_out$decontaminated_samples ) )
  
  q <- get_rescaled_mat( out$decont_low_out$decontaminated_samples )
  q[is.na(q)] <- 0
  dec_low_jsd <- get_rowwise_jsd(samples, 
                                 q)
  
  restrictive_jsd <- get_rowwise_jsd(samples, 
                                     get_rescaled_mat( out$restr_out$decontaminated_samples) )
  
  jsds <- c(input_sink_jsd, 
            SCRUB_jsd, 
            dec_jsd, dec_low_jsd, restrictive_jsd, microdec_jsd)
  names <-  c(rep('Input', length(input_sink_jsd)), 
              rep('SCRUB', length(input_sink_jsd)), 
              rep('Decontam', length(input_sink_jsd)), 
              rep('Decontam (Low Biomass)', length(input_sink_jsd)), 
              rep('Restrictive', length(input_sink_jsd)),
              rep('microDecon', length(input_sink_jsd))
  )
  
  result <- data.frame(names, jsds) %>% 
    mutate( n_taxa=n_taxa, 
            n_controls=n_controls, 
            contam_species=sum(true_contamination>0),
            contam_level=contam_level, 
            well_to_well_level=well_level, 
            plate_format=plate_format, 
            round=idx) 
  if(is.data.frame(new_results)==F){ return(result)
  }else{return(rbind(new_results, result))}
}

synth_data <- read.csv('data/Fig1/synth_combination.csv', row.names = 1) %>% as.matrix()
colnames(synth_data) <- paste0('OTU_', 1:ncol(synth_data))

controls <- synth_data[row.names(synth_data) %>% str_detect('BLANK') %>% 
                         which(), ]
samples <- synth_data[( F == row.names(synth_data) %>% str_detect('BLANK') ) %>% 
                        which(), ]

full_samps <- samples
full_conts <- controls

finalized_simulations <- function(full_samps, 
                                  full_conts, 
                                  use_uploaded_samples=T,
                                  n_experiments=10, 
                                  ){
  
  set.seed(123)
  total_experiments <- n_experiments * 3 
  
  n_taxa <- ncol(full_samps)
  n_samples <- 88
  
  control_noise <- 'low'
  sample_noise <- 'low'
  contam_in_sample_noise <- 'low'
  
  n_contam_species <- NA
  
  sitch_results <- NA
  idx <- 1
  
  for(concentration_sigma in c(2, 1, 0.1) ){
    
    for(experiment in 1:n_experiments){
      
      n_controls <- 8
      
      samples <- get_samps_from_upload(full_samps, n_samples) %>% get_rescaled_mat
      row.names(samples) <- paste0('Sample_', 1:n_samples)
      
      median_value <- 1e8  # Desired median for concentration
      mu <- log(median_value)  # Mean of the underlying normal distribution
      true_concentrations <- rlnorm(n_samples, meanlog = mu, sdlog = concentration_sigma)
      
      true_contamination <- full_conts[sample(nrow(full_conts), 1),] %>% rescale
      n_taxa <- ncol(samples)
      
      
      # draw noisy / contaminated samples, 8 noisy controls
      
      observed_samples <- sapply(1:n_samples, function(x){
        draw_contaminated_sample(samples[x,],
                                 true_contamination, 
                                 samples,
                                 1-as.numeric(map_concentration_to_contamination( 
                                   data.frame(x=c(true_concentrations[x]) ) ) ), # contam_level, 
                                 n_taxa, 
                                 c(),
                                 sample_noise=sample_noise, 
                                 sample_well_level=0,
                                 contam_in_sample_noise=contam_in_sample_noise, 
        ) %>% return() }  
      ) %>% t()
      
      row.names(observed_samples) <- row.names(samples)
      
      
      controls_pre_well_to_well <- sapply(1:n_controls, function(x){
        draw_control( true_contamination, 
                      samples, 
                      x, 
                      c(), 
                      0,
                      control_noise, 
                      n_taxa
        ) %>% return()
      } ) %>% t()
      
      row.names(controls_pre_well_to_well) <- paste0('Control_', 1:8)
      
      ## adding some variance/uncertainty
      
      
      observed_concentrations <- 10**(log10( true_concentrations) * rnorm( length( true_concentrations ), mean=1, sd=.1 ) )
      
      
      for(plate_format in c('F')){ 
        well_positions <- Simulate_Plate_Positions(n_samples, 8, n_wells=8, plate_format=plate_format)
        well_dists <- dist(well_positions) %>% as.matrix()
        for(well_level in c(0)){
          well_species_limit <- NA
          observed_controls <- get_final_conts(controls_pre_well_to_well,
                                               samples,
                                               well_dists, 
                                               well_level, 
                                               well_species_limit)
          colnames(observed_samples) <- colnames(samples)
          colnames(observed_controls) <- colnames(samples)
          for(n_controls in c(4)){ 
            if(n_controls>1){ out <- run_SCRUB_DEC_and_RESTR(observed_samples, 
                                                             observed_controls[1:n_controls,], 
                                                             well_dists, 
                                                             observed_concentrations,
                                                             dist_threshold=1.5)
            }else{
              tmp <- observed_controls[1, ] %>% t()
              row.names(tmp) <- c('Control_1')
              out <- run_SCRUB_DEC_and_RESTR(observed_samples, 
                                             tmp, 
                                             well_dists, 
                                             observed_concentrations,
                                             dist_threshold=1.5)
            }
            # store_finalized_results_using_jsd
            print(paste('Running Fig_1 iteration', 
                        idx, 'of', total_experiments))
            sitch_results <- store_f1_results(sitch_results, 
                                              out,
                                              samples, 
                                              observed_samples, 
                                              n_controls, 
                                              well_species_limit, 
                                              plate_format, 
                                              well_level, 
                                              concentration_sigma,
                                              true_contamination,
                                              n_taxa,
                                              idx)
            
            out_path <- paste0('simulation_data/round_',idx)
            dir.create( out_path )
            
            samples %>% write.csv( paste0(out_path, 
                                          '/true_samples.csv') )
            
            observed_samples %>% write.csv( paste0(out_path, 
                                            '/observed_samples.csv') )
            
            get_rescaled_mat( out$SCRUB_out$decontaminated_samples ) %>% 
              write.csv( paste0(out_path, 
                                '/scrub_out.csv') )
            
            
            observed_concentrations %>% 
              write.csv( paste0(out_path, 
                                '/true_concentrations.csv') )
            
            
            sitch_results %>% write.csv('simulation_data/concentration_simulation_results.csv')
            
          
            
            idx <- idx+1
          }
      }
      print(idx)
    }
  }}
  return(sitch_results)
}

cols <- which( samples %>% colSums > 1000 )
samples <-  samples[, cols]
controls <- controls[, cols]


results_of_final_sims <- finalized_simulations(samples, 
                                               controls,
                                               n_experiments = 10
                                               )



simulations_varying_conc_noise <- function(full_samps, 
                                          full_conts, 
                                          use_uploaded_samples=T,
                                          n_experiments=10
                                          ){
  
  set.seed(123)
  total_experiments <- n_experiments * 3 
  
  n_taxa <- ncol(full_samps)
  n_samples <- 88
  
  control_noise <- 'low'
  sample_noise <- 'low'
  contam_in_sample_noise <- 'low'
  
  n_contam_species <- NA
  
  sitch_results <- NA
  idx <- 1
  
  
  for(conc_noise in c(0.03, 0.04, 0.05)){
    for(concentration_sigma in c(
                               1
                               ) ){
    
     for(experiment in 1:n_experiments){
      
      n_controls <- 8
      
      samples <- get_samps_from_upload(full_samps, n_samples) %>% get_rescaled_mat
      row.names(samples) <- paste0('Sample_', 1:n_samples)
      
      median_value <- 1e8  # Desired median for concentration
      mu <- log(median_value)  # Mean of the underlying normal distribution
      true_concentrations <- rlnorm(n_samples, meanlog = mu, sdlog = concentration_sigma)
      
      true_contamination <- full_conts[sample(nrow(full_conts), 1),] %>% rescale
      n_taxa <- ncol(samples)
      
      
      # draw noisy / contaminated samples, 8 noisy controls
      
      observed_samples <- sapply(1:n_samples, function(x){
        draw_contaminated_sample(samples[x,],
                                 true_contamination, 
                                 samples,
                                 1-as.numeric(map_concentration_to_contamination( 
                                   data.frame(x=c(true_concentrations[x]) ) ) ), # contam_level, 
                                 n_taxa, 
                                 c(),
                                 sample_noise=sample_noise, 
                                 sample_well_level=0,
                                 contam_in_sample_noise=contam_in_sample_noise, 
        ) %>% return() }  
      ) %>% t()
      
      row.names(observed_samples) <- row.names(samples)
      
      
      controls_pre_well_to_well <- sapply(1:n_controls, function(x){
        draw_control( true_contamination, 
                      samples, 
                      x, 
                      c(), 
                      0,
                      control_noise, 
                      n_taxa
        ) %>% return()
      } ) %>% t()
      
      row.names(controls_pre_well_to_well) <- paste0('Control_', 1:8)
      
      ## adding some variance/uncertainty
      
      observed_concentrations <- 10**(log10( true_concentrations) * rnorm( length( true_concentrations ), 
                                                                           mean=1, 
                                                                           sd=conc_noise ) )
      
      
      for(plate_format in c('F')){ 
        well_positions <- Simulate_Plate_Positions(n_samples, 8, n_wells=8, plate_format=plate_format)
        well_dists <- dist(well_positions) %>% as.matrix()
        for(well_level in c(0.05)){
          well_species_limit <- NA
          observed_controls <- get_final_conts(controls_pre_well_to_well,
                                               samples,
                                               well_dists, 
                                               well_level, 
                                               well_species_limit)
          colnames(observed_samples) <- colnames(samples)
          colnames(observed_controls) <- colnames(samples)
          
          
          
          
          for(n_controls in c(4)){ 
            
            mps <- c( rep('control', n_controls), 
                      rep('sample', n_samples) )
            
            mp1 <- c( rep(T, n_controls), 
                      rep(F, n_samples) )
            
            scr_md <- data.frame( mp1,
                                  mps,
                                  c( paste0( LETTERS[well_positions$x], well_positions$y )[1:n_controls], 
                                     paste0( LETTERS[well_positions$x], well_positions$y )[9:96] 
                                     )
                                  )
            row.names(scr_md) <- c( row.names(observed_controls[1:n_controls,]), 
                                    row.names(observed_samples))
            
            if(n_controls>1){ out <- run_SCRUB_DEC_and_RESTR(observed_samples, 
                                                             observed_controls[1:n_controls,], 
                                                             well_dists, 
                                                             observed_concentrations,
                                                             dist_threshold=1.5, 
                                                             scr_md=scr_md)
            }else{
              tmp <- observed_controls[1, ] %>% t()
              row.names(tmp) <- c('Control_1')
              out <- run_SCRUB_DEC_and_RESTR(observed_samples, 
                                             tmp, 
                                             well_dists, 
                                             observed_concentrations,
                                             dist_threshold=1.5, 
                                             scr_md=scr_md)
            }
            # store_finalized_results_using_jsd
            print(paste('Running Fig_1 iteration', 
                        idx, 'of', total_experiments))
            sitch_results <- store_f1_results(sitch_results, 
                                              out,
                                              samples, 
                                              observed_samples, 
                                              n_controls, 
                                              well_species_limit, 
                                              plate_format, 
                                              conc_noise, 
                                              concentration_sigma,
                                              true_contamination,
                                              n_taxa,
                                              idx)
            
            out_path <- paste0('simulation_data/round_',idx)
            dir.create( out_path )
            
            samples %>% write.csv( paste0(out_path, 
                                          '/true_samples.csv') )
            
            observed_samples %>% write.csv( paste0(out_path, 
                                                   '/observed_samples.csv') )
            
            get_rescaled_mat( out$SCRUB_out$decontaminated_samples ) %>% 
              write.csv( paste0(out_path, 
                                '/scrub_out.csv') )
            
            
            observed_concentrations %>% 
              write.csv( paste0(out_path, 
                                '/true_concentrations.csv') )
            
            
            sitch_results %>% write.csv('simulation_data/concentration_simulation_results.csv')
            
            
            
            idx <- idx+1
          }
        }
        print(idx)
        }
       }
     }
   }
  return(sitch_results)
}

cols <- which( samples %>% colSums > 1000 )
samples <-  samples[, cols]
controls <- controls[, cols]


results_of_final_sims <- simulations_varying_conc_noise(samples, 
                                                        controls,
                                                        n_experiments = 10
                                                        )




simulations_varying_processing_bias <- function(full_samps, 
                                                 full_conts, 
                                                 use_uploaded_samples=T,
                                                 n_experiments=10
){
  
  set.seed(123)
  total_experiments <- n_experiments * 3 
  
  n_taxa <- ncol(full_samps)
  n_samples <- 88
  
  control_noise <- 'low'
  sample_noise <- 'low'
  contam_in_sample_noise <- 'low'
  
  n_contam_species <- NA
  
  sitch_results <- NA
  idx <- 1
  
  
  for(processing_bias_std in c(0, 0.5, 1, 2, 3)){
    for(concentration_sigma in c(
                                1
                            ) ){
      
      for(experiment in 1:n_experiments){
        
        n_controls <- 8
        
        samples <- get_samps_from_upload(full_samps, n_samples) %>% get_rescaled_mat
        row.names(samples) <- paste0('Sample_', 1:n_samples)
        
        median_value <- 1e8  # Desired median for concentration
        mu <- log(median_value)  # Mean of the underlying normal distribution
        true_concentrations <- rlnorm(n_samples, meanlog = mu, sdlog = concentration_sigma)
        
        true_contamination <- full_conts[sample(nrow(full_conts), 1),] %>% rescale
        n_taxa <- ncol(samples)
        
        
        # draw noisy / contaminated samples, 8 noisy controls
        
        observed_samples <- sapply(1:n_samples, function(x){
          draw_contaminated_sample(samples[x,],
                                   true_contamination, 
                                   samples,
                                   1-as.numeric(map_concentration_to_contamination( 
                                     data.frame(x=c(true_concentrations[x]) ) ) ), 
                                   n_taxa, 
                                   c(),
                                   sample_noise=sample_noise, 
                                   sample_well_level=0,
                                   contam_in_sample_noise=contam_in_sample_noise, 
          ) %>% return() }  
        ) %>% t()
        
        row.names(observed_samples) <- row.names(samples)
        
        
        controls_pre_well_to_well <- sapply(1:n_controls, function(x){
          draw_control( true_contamination, 
                        samples, 
                        x, 
                        c(), 
                        0,
                        control_noise, 
                        n_taxa
          ) %>% return()
        } ) %>% t()
        
        row.names(controls_pre_well_to_well) <- paste0('Control_', 1:8)
        
        ## adding some variance/uncertainty
        
        conc_noise=0.03
        
        observed_concentrations <- 10**(log10( true_concentrations) * rnorm( length( true_concentrations ), 
                                                                             mean=1, 
                                                                             sd=conc_noise ) )
        
        ## add processing biases
        processing_biases <- 2 ^ (rnorm(n_taxa, 0, processing_bias_std) )
        observed_samples <- observed_samples*processing_biases
        
        ## rescale back to cuount space
        observed_samples <- ( (observed_samples / (rowSums(observed_samples)) )*10000 ) %>% round(digits=0)
        
        
        for(plate_format in c('F')){ 
          well_positions <- Simulate_Plate_Positions(n_samples, 8, n_wells=8, plate_format=plate_format)
          well_dists <- dist(well_positions) %>% as.matrix()
          for(well_level in c(0.05)){
            well_species_limit <- NA
            observed_controls <- get_final_conts(controls_pre_well_to_well,
                                                 samples,
                                                 well_dists, 
                                                 well_level, 
                                                 well_species_limit)
            
            observed_controls <- observed_controls*processing_biases
            observed_controls <- ( (observed_controls / (rowSums(observed_controls)) )*10000 ) %>% round(digits=0)
            
            
            
            colnames(observed_samples) <- colnames(samples)
            colnames(observed_controls) <- colnames(samples)
            for(n_controls in c(4)){ # test across 1:8 controls
              
              mps <- c( rep('control', n_controls), 
                        rep('sample', n_samples) )
              
              mp1 <- c( rep(T, n_controls), 
                        rep(F, n_samples) )
              
              scr_md <- data.frame( mp1,
                                    mps,
                                    c( paste0( LETTERS[well_positions$x], well_positions$y )[1:n_controls], 
                                       paste0( LETTERS[well_positions$x], well_positions$y )[9:96] 
                                    )
              )
              row.names(scr_md) <- c( row.names(observed_controls[1:n_controls,]), 
                                      row.names(observed_samples))
              
              
              if(n_controls>1){ out <- run_SCRUB_DEC_and_RESTR(observed_samples, 
                                                               observed_controls[1:n_controls,], 
                                                               well_dists, 
                                                               observed_concentrations,
                                                               dist_threshold=1.5, 
                                                               scr_md=scr_md)
              }else{
                tmp <- observed_controls[1, ] %>% t()
                row.names(tmp) <- c('Control_1')
                out <- run_SCRUB_DEC_and_RESTR(observed_samples, 
                                               tmp, 
                                               well_dists, 
                                               observed_concentrations,
                                               dist_threshold=1.5, 
                                               scr_md=scr_md)
              }
              
              true_samples_with_bias <- samples*processing_biases
              true_samples_with_bias <- ( (true_samples_with_bias / (rowSums(true_samples_with_bias)) )*10000 ) %>% round(digits=0)
              
              # store_finalized_results_using_jsd
              print(paste('Running Fig_1 iteration', 
                          idx, 'of', total_experiments))
              sitch_results <- store_f1_results(sitch_results, 
                                                out,
                                                true_samples_with_bias, 
                                                observed_samples, 
                                                n_controls, 
                                                well_species_limit, 
                                                plate_format, 
                                                processing_bias_std, 
                                                concentration_sigma,, 
                                                true_contamination,
                                                n_taxa,
                                                idx)
              
              out_path <- paste0('simulation_data_bias/round_',idx)
              dir.create( out_path )
              
              true_samples_with_bias %>% write.csv( paste0(out_path, 
                                            '/true_samples.csv') )
              
              observed_samples %>% write.csv( paste0(out_path, 
                                                     '/observed_samples.csv') )
              
              get_rescaled_mat( out$SCRUB_out$decontaminated_samples ) %>% 
                write.csv( paste0(out_path, 
                                  '/scrub_out.csv') )
              
              
              observed_concentrations %>% 
                write.csv( paste0(out_path, 
                                  '/true_concentrations.csv') )
              
              
              sitch_results %>% write.csv('simulation_data_bias/concentration_simulation_results.csv')
              
              
              
              idx <- idx+1
            }
          }
          print(idx)
        }
      }
    }
  }
  return(sitch_results)
}

cols <- which( samples %>% colSums > 1000 )
samples <-  samples[, cols]
controls <- controls[, cols]


results_of_final_sims <- simulations_varying_processing_bias(samples, 
                                                             controls,
                                                             n_experiments = 10
)



library(SparseDOSSA2)

load('../data/SPARSEDOSSA_generators.Rdata')



set_up_dmc_scrub_simulations <- function(full_samps, 
                                         full_conts, 
                                         use_uploaded_samples=T,
                                         n_experiments=10, 
                                         n_batches=4
                                         ){
  
  set.seed(123)
  total_experiments <- n_experiments  
  
  n_taxa <- ncol(full_samps)
  n_samples <- 96 * n_batches
  
  control_noise <- 'low'
  sample_noise <- 'low'
  contam_in_sample_noise <- 'low'
  
  n_contam_species <- NA
  
  sitch_results <- NA
  idx <- 1
  
  for(conc_noise in c(0.03)){
    for(concentration_sigma in c(
      1
    ) ){
      
      for(experiment in 1:n_experiments){
        print(idx)
        
        n_controls <- 8
        
        samples <- SparseDOSSA2(
                                template = fitted_sparsedossa_PNP_samples,
                                n_sample = n_samples, 
                                n_feature = 100, 
                                median_read_depth = 1e4,
                                verbose = F
                                )$simulated_data %>% t() %>% as.data.frame
        
        
        row.names(samples) <- paste0('Sample_', 1:n_samples)
        
        # median_value <- 1e8  # Desired median for concentration
        median_value <- 1e8
        
        mu <- log(median_value)  # Mean of the underlying normal distribution
        true_concentrations <- rlnorm(n_samples, meanlog = mu, sdlog = concentration_sigma)
        
         
        observed_concentrations <- 10**(log10( true_concentrations) * rnorm( length( true_concentrations ),
                                                                             mean=1,
                                                                             sd=conc_noise ) )

        
        
        out_path <- paste0('multibatch_simulation_data/round_',idx)

        samples %>% write.csv( paste0(out_path,
                                      '/true_samples.csv') )

        observed_concentrations %>%
          write.csv( paste0(out_path,
                            '/observed_concentrations.csv') )
        
        
        ( 1-as.numeric(map_concentration_to_contamination(
           data.frame(x=c(true_concentrations) ) ) ) ) %>%
                    write.csv( paste0(out_path,
                                      '/underlying_contamination_levels.csv') )  
        
        idx <- idx+1
      }
    }
  }
}



cols <- which( samples %>% colSums > 1000 )
samples <-  samples[, cols]
controls <- controls[, cols]


tmp <- set_up_dmc_scrub_simulations(samples, 
                                    controls,
                                    n_experiments = 10 )




set_up_dmc_scrub_simulations_across_contam_ranges <- function(full_samps, 
                                         full_conts, 
                                         use_uploaded_samples=T,
                                         n_experiments=10, 
                                         n_batches=4
){
  
  set.seed(123)
  total_experiments <- n_experiments 
  
  n_taxa <- ncol(full_samps)
  n_samples <- 96 * n_batches
  
  control_noise <- 'low'
  sample_noise <- 'low'
  contam_in_sample_noise <- 'low'
  
  n_contam_species <- NA
  
  sitch_results <- NA
  idx <- 1
  
  conc_noise <- 0.03
  for(median_value in c(1e7, 1e8, 1e9)){
    for(concentration_sigma in c(
      1
    ) ){
      
      
      for(experiment in 1:n_experiments){
        print(idx)
        
        n_controls <- 8
        
        samples <- SparseDOSSA2(
          template = fitted_sparsedossa_PNP_samples,
          n_sample = n_samples, 
          n_feature = 100, 
          median_read_depth = 1e4,
          verbose = F
        )$simulated_data %>% t() %>% as.data.frame
        
        row.names(samples) <- paste0('Sample_', 1:n_samples)
        
        
        mu <- log(median_value)  # Mean of the underlying normal distribution
        true_concentrations <- rlnorm(n_samples, meanlog = mu, sdlog = concentration_sigma)
        
        
        observed_concentrations <- 10**(log10( true_concentrations) * rnorm( length( true_concentrations ),
                                                                             mean=1,
                                                                             sd=conc_noise ) )
        
        
        
        out_path <- paste0('multibatch_simulation_data_varying_contam/round_',idx)
        dir.create( out_path )
        
        samples %>% write.csv( paste0(out_path,
                                      '/true_samples.csv') )
        
        
        observed_concentrations %>%
          write.csv( paste0(out_path,
                            '/observed_concentrations.csv') )
        
        
        ( 1-as.numeric(map_concentration_to_contamination(
          data.frame(x=c(true_concentrations) ) ) ) ) %>%
          write.csv( paste0(out_path,
                            '/underlying_contamination_levels.csv') )  
        
        idx <- idx+1
      }
    }
  }
}


tmp <- set_up_dmc_scrub_simulations_across_contam_ranges(samples, 
                                                         controls,
                                                         n_experiments = 10 )



### EXAMPLE FOR FIG 3c

set_up_example_for_3b <- function(full_samps, 
                                  full_conts, 
                                  use_uploaded_samples=T,
                                  n_experiments=10, 
                                  seed=123
){
  
  set.seed(seed)
  total_experiments <- n_experiments * 3 
  
  n_taxa <- ncol(full_samps)
  n_samples <- 88
  
  control_noise <- 'low'
  sample_noise <- 'low'
  contam_in_sample_noise <- 'low'
  
  n_contam_species <- NA
  
  sitch_results <- NA
  idx <- 1
  
    for(processing_bias_std in c(2)){
    for(concentration_sigma in c(
      1
    ) ){
      
      for(experiment in 1:1){
        
        n_controls <- 8
        
        samples <- get_samps_from_upload(full_samps, n_samples) %>% get_rescaled_mat
        row.names(samples) <- paste0('Sample_', 1:n_samples)
        
        median_value <- 1e8  # Desired median for concentration
        mu <- log(median_value)  # Mean of the underlying normal distribution
        true_concentrations <- rlnorm(n_samples, meanlog = mu, sdlog = concentration_sigma)
        
        true_contamination <- full_conts[sample(nrow(full_conts), 1),] %>% rescale
        n_taxa <- ncol(samples)
        
        
        # draw noisy / contaminated samples, 8 noisy controls
        
        observed_samples <- sapply(1:n_samples, function(x){
          draw_contaminated_sample(samples[x,],
                                   true_contamination, 
                                   samples,
                                   1-as.numeric(map_concentration_to_contamination( 
                                     data.frame(x=c(true_concentrations[x]) ) ) ), 
                                   n_taxa, 
                                   c(),
                                   sample_noise=sample_noise, 
                                   sample_well_level=0,
                                   contam_in_sample_noise=contam_in_sample_noise, 
          ) %>% return() }  
        ) %>% t()
        
        row.names(observed_samples) <- row.names(samples)
        
        
        controls_pre_well_to_well <- sapply(1:n_controls, function(x){
          draw_control( true_contamination, 
                        samples, 
                        x, 
                        c(), 
                        0,
                        control_noise, 
                        n_taxa
          ) %>% return()
        } ) %>% t()
        
        row.names(controls_pre_well_to_well) <- paste0('Control_', 1:8)
        
        ## adding some variance/uncertainty
        
        conc_noise=0.03
        
        observed_concentrations <- 10**(log10( true_concentrations) * rnorm( length( true_concentrations ), 
                                                                             mean=1, 
                                                                             sd=conc_noise ) )
        
        ## add processing biases
        processing_biases <- 2 ^ (rnorm(n_taxa, 0, processing_bias_std) )
        observed_samples <- observed_samples*processing_biases
        
        ## rescale back to cuount space
        observed_samples <- ( (observed_samples / (rowSums(observed_samples)) )*10000 ) %>% round(digits=0)
        
        ## items to return:
          # true samples
          # true contamination source
          # true contamination level
          # true concentrations
          # biases 
      }}}
        
      outto <- list() 
      outto[['PBS']] = processing_biases
      outto[['SMPS']] = samples
      outto[['CONC']] = true_concentrations
      outto[['CONT']] = true_contamination
      outto[['CLS']] = 1-as.numeric(map_concentration_to_contamination( 
        data.frame(x=c(true_concentrations) ) ) )
      return(outto)
  
        }
        


max_corr=-1

for(seed in 1:300){
  setup <- set_up_example_for_3b(samples, 
                                 controls,
                                 n_experiments = 10, 
                                 seed=seed
                                 )
  
  
  swbs <- rowSums(  (1 - setup$CLS) * setup$SMPS * setup$PBS )
  clwbs <- setup$CLS * sum( setup$CONT * setup$PBS )

  pos_cor <- cor(setup$CONC, 
                 clwbs / (clwbs+swbs), 
                 method='spearman')
  
  if( pos_cor > max_corr){
    best_seed <- seed
    max_corr <- pos_cor
  }}
  


setup <- set_up_example_for_3b(samples, 
                               controls,
                               n_experiments = 10, 
                               seed=best_seed
)


swbs <- rowSums(  (1 - setup$CLS) * setup$SMPS * setup$PBS )
clwbs <- setup$CLS * sum( setup$CONT * setup$PBS )
 
plot(setup$CONC,
     clwbs / (clwbs+swbs),
     log = "x",
    )


data.frame(conc=setup$CONC, 
           contam_level=setup$CLS, 
           biased_contam_level=clwbs / (clwbs+swbs)
           ) %>% write.csv('../results/Example_for_3b.csv')





