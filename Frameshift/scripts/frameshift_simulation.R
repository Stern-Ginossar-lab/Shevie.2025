#!/usr/bin/env Rscript

#
# 	For each percentage of framshit  unique( as.integer( 10^seq( 0, 2, .1 ) ) ) , we should do the nof_iters on the set of windows as we
#	found them on the spike / fluc, and then for each window to find the Z scores 
#	

{
    library( tidyr )
    library( dplyr )

    root = paste0( my_home, "/People/Shevie/2025.08/" )
}

{
    if ( Sys.getenv( "RSTUDIO" ) == 1 )
    {
        file_name = "minU"
        modif = "n1m"
    }
    else
    {
        params = commandArgs( TRUE )
        file_name = as.character( params[ 1 ] )
        modif = as.character( params[ 2 ] )
    }
    nof_iters = 10000
}


{
#    peaks.in = "SS" ; peaks.out = "SS.wT"
#    peaks.in = "SS" ; peaks.out = "SS.woT"
#    peaks_are = "SS" ; peaks.out = "SS.all"
#    peaks.in = "DE" ; peaks.out = "DE.wT"
#    peaks.in = "DE" ; peaks.out = "DE.woT"
    peaks_are = "UUUC" ; peaks.out = "UUUC.all"
    
    left_trim = 3 * 0
    right_trim = 3 * 0    
}

{
    data_file = paste0( root, "/inps/", file_name, ".nucs.tsv" )
    sim_file = paste0( root, "outs/frameshift/NO_FS_simulations/csvs/", file_name, "_", modif, ".10000.all.tsv" )
    win_file = paste0( root, "/outs/frameshift/frameshift/", peaks_are, "/csvs/", file_name, "_", modif, ".tsv" )

    frames = c( "F", "Fp1", "Fm1" )
}

#   The raw data of the reads
{
    data.o = my_read_table( data_file, my_row.names = -1 )
    data.cds = data.o[ data.o$region == "cds", ]
    data.trim = data.cds[ ( left_trim + 1 ) : ( nrow ( data.cds ) - right_trim ), ]
    data.trim[ , modif ] = rowSums( data.trim[ , grep( modif, colnames( data.trim ) ) ] )
    vals = data.trim[ , modif ]
    
    df.vals = data.frame( matrix( vals, nrow = nrow( data.trim ) / 3, 3, byrow = TRUE ) )
    colnames( df.vals ) = frames
}

#   The simulation data (for each window size we have mean and std that will be used to calculate the Z-score)
{
    sim.o = my_read_table( sim_file )
    
#    > colnames( sim.o )
#    [1]  window_size F_0.05      F_0.5       F_0.95      F_mean     
#    [6]  F_sd        Fp1_0.05    Fp1_0.5     Fp1_0.95    Fp1_mean   
#    [11] Fp1_sd      Fm1_0.05    Fm1_0.5     Fm1_0.95    Fm1_mean   
#    [16] Fm1_sd      Fp1toF_0.05 Fp1toF_0.5  Fp1toF_0.95 Fp1toF_mean
#    [21] Fp1toF_sd
}

#   The windows data 
{
    win.o = my_read_table( win_file, my_row.names = -1 )

#    > colnames( win.o)
#    [1]  trans.pos PA.seq    nof_T     dist_next dist_stop codon_pos
#    [7]  length    win_start win_end   F         Fp1       Fm1      
#    [13] total     F.p       Fp1.p     Fm1.p     Fp1toF    F.p.z    
#    [19] Fp1.p.z   Fp1toF.z 

	if ( length( grep( "wT", peaks.out ) == 1 ) )
		win = win.o %>% filter( nof_T > 0 )
	if ( length( grep( "woT", peaks.out ) == 1 ) )
		win = win.o %>% filter( nof_T == 0 )
	if ( length( grep( "wT", peaks.out ) == 0 ) & length( grep( "woT", peaks.out ) == 0 ) )
		win = win.o %>% filter( nof_T > 0 )
	if ( length( grep( "all", peaks.out ) == 1 ) )
		win = win.o
    win_lens = win$length
}

{
    sim_data_sampling <- function( frame, nof_vals = window_size )
    {
        sample( df.vals[ , frame ], nof_vals, replace = TRUE )
    }
    
    my_summ.titles = c( "mean", "median", "min", "max" )
    my_summ <- function( arr )
    {
        c( mean( arr ), median( arr ), min( arr ), max( arr ) )
    }
}

{
    res_df = data.frame()
    for ( frameshift_perc in unique( as.integer( 10^seq( 0, log10( 26 ), .1 ) ) ) )
    {
        frac = frameshift_perc / 100.0
        for ( iter in 1 : nof_iters )
        {
            Z.F = Z.Fp1 = Z.Fp1toF = c()
            for ( window_size in win_lens )
            {
                sum.F = sum( ( 1 - frac ) * sim_data_sampling( "F" ) + frac * sim_data_sampling( "Fm1" ) )
                sum.Fp1 = sum( ( 1 - frac ) * sim_data_sampling( "Fp1" ) + frac * sim_data_sampling( "F" ) )
                sum.Fm1 = sum( ( 1 - frac ) * sim_data_sampling( "Fm1" ) + frac * sim_data_sampling( "Fp1" ) )
                
                frame_fracs = 100.0 * c( sum.F, sum.Fp1, sum.Fm1 ) / sum( c( sum.F, sum.Fp1, sum.Fm1 ) )
                names( frame_fracs ) = frames
                
                sim_vals = sim.o[ as.character( window_size ), , drop = FALSE ]
                Z.F = c( Z.F, ( frame_fracs[ "F" ] - sim_vals$F_mean ) / sim_vals$F_sd ) 
                Z.Fp1 = c( Z.Fp1, ( frame_fracs[ "Fp1" ] - sim_vals$Fp1_mean ) / sim_vals$Fp1_sd )
                Z.Fp1toF = c( Z.Fp1toF, ( frame_fracs[ "Fp1" ] / frame_fracs[ "F" ] - sim_vals$Fp1toF_mean ) / sim_vals$Fp1toF_sd )
            }
            stats = c( frameshift_perc, my_summ( Z.F ), my_summ( Z.Fp1 ), my_summ( Z.Fp1toF ) )
            res_df = rbind( res_df, stats )
        }
    }
    colnames( res_df ) = c( "FS_perc", paste.cross( c( "F", "Fp1", "Fp1toF" ), my_summ.titles ) )
}

{
    out_file = paste0( root, "outs/frameshift/FS_simulations/", peaks_are, "/", file_name, "_", modif, ".", nof_iters, ".tsv" )
    write.table( res_df, out_file, sep = "\t", row.names = FALSE, quote = FALSE )
}