#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------#
# SET ENVIRONMENT #
#-----------------#
# You will need to modify the path to `Site_models` if you run this
# script elsewhere!
wd <- paste( "~/my_session/day3/01_site_models/" )
setwd( wd )

#-----------#
# LOAD DATA #
#-----------#
# 1. Load our text file with the lnL values
lnL_vals <- read.table( file = "Site_models/lnL_sites.txt", sep= " ",
                        stringsAsFactors = FALSE, header = TRUE )

# 2. We can now compute the LRT statistic, but this can 
# only be done between models that are nested:
#     M0 < M1a < M2a 
#     M7 < M8 

# First, we will compare M1a (alternative) to the M0 (null). 
# If the null is rejected, then we can compare M1a (null) to 
# M2a (alternative).
## M0 vs M1a ##
diff_M0vsM1a <- 2*( lnL_vals$NSsites_1[1] - lnL_vals$NSsites_0[1] )
diff_M0vsM1a
# diff = 912.0554
# Find degrees of freedom
df_M0vsM1a <- lnL_vals$NSsites_1[2] - lnL_vals$NSsites_0[2]
df_M0vsM1a
# 1 degree of freedom
pchisq( diff_M0vsM1a, df = df_M0vsM1a, lower.tail=F )
# p-val = 2.350442e-200 < 0.05
Chisq.crit.M0vsM1a <- qchisq( p = 0.95, df = df_M0vsM1a )
Chisq.crit.M0vsM1a
# alpha critical value at 5% = 3.841459
Chisq.crit.M0vsM1a_2 <- qchisq( p = 0.99, df = df_M0vsM1a )
Chisq.crit.M0vsM1a_2
# alpha critical value at 1% = 6.634897

# As M1a is a better fit to the data than M0, we can compare M1a 
# (Nearly Neutral) against M2a (Positive Selection). 
## M1 vs M2a ##
diff_M1avsM2a <- 2*( lnL_vals$NSsites_2[1] - lnL_vals$NSsites_1[1] )
diff_M1avsM2a
# diff = 0
# Find degree of freedom
df_M1avsM2a <- lnL_vals$NSsites_2[2] - lnL_vals$NSsites_1[2]
df_M1avsM2a
# 2 degree of freedom
pchisq( diff_M1avsM2a, df = df_M1avsM2a, lower.tail=F )
# p-val = 1 > 0
Chisq.crit.M1vsM2a <- qchisq( p = 0.95, df = df_M1avsM2a )
Chisq.crit.M1vsM2a
# alpha critical value at 5% level = 5.991465
Chisq.crit.M1vsM2a_2 <- qchisq( p = 0.99, df = df_M1avsM2a )
Chisq.crit.M1vsM2a_2
# alpha critical value at 1% level = 9.21034

# In addition, we can run an additional comparison between 
# M7 (beta) and M8 (beta&omega).
## M7 vs M8 ##
diff_M7vsM8 <- 2*( lnL_vals$NSsites_8[1] - lnL_vals$NSsites_7[1] )
diff_M7vsM8
# diff = 39.14911
# Find degree of freedom
df_M7vsM8 <- lnL_vals$NSsites_8[2] - lnL_vals$NSsites_7[2]
df_M7vsM8
# 2 degree of freedom
pchisq( diff_M7vsM8, df = df_M7vsM8, lower.tail=F )
# p-val = 3.154118e-09 < 0.05
Chisq.crit.M7vsM8 <- qchisq( p = 0.95, df = df_M7vsM8 )
# alpha critical value at 5% level = 5.991465
Chisq.crit.M7vsM8_2 <- qchisq( p = 0.99, df = df_M7vsM8 )
# alpha critical value at 1% level = 9.21034

# 3. Plot results 
par( mfrow = c(1, 3) )

# M0 vs M1a
curve( dchisq( x, df = 1 ), from = 0, to =  930 )
abline( v = c( Chisq.crit.M0vsM1a, Chisq.crit.M0vsM1a_2, diff_M0vsM1a ), col = c( "darkgray", "brown", "red" ) )
coords_dev    <- c( 407, 0.0012 )
coords_pval   <- c( 410.4, 0.0011 )
coords_alphac <- c( 390, 0.0010 )
coords_alphac2 <- c( 390, 0.0009 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = 912.0554', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["1,0.05"]^"2", "= 3.84", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_alphac2[1], y = coords_alphac2[2],
      labels = expression( atop( paste( chi["1,0.01"]^"2", "= 6.63", sep = " " ) ) ),
      cex = 1.2, col = "brown" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 2.35e-200e-123' ) ),
      cex = 1.2, col = "black" )
title( "A) M0 vs M1a" )

#M1a vs M2a 
curve( dchisq( x, df = 2 ), from = 0, to =  10 )
abline( v = c( Chisq.crit.M1vsM2a, Chisq.crit.M1vsM2a_2, diff_M1avsM2a ), col = c( "darkgray", "brown", "red" ) )
coords_dev     <- c( 7.85, 0.48 )
coords_pval    <- c( 7.8, 0.45 )
coords_alphac  <- c( 8, 0.41 )
coords_alphac2 <- c( 8, 0.38 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = 0', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["2,0.05"]^"2", "= 5.99", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_alphac2[1], y = coords_alphac2[2],
      labels = expression( atop( paste( chi["2,0.01"]^"2", "= 9.21", sep = " " ) ) ),
      cex = 1.2, col = "brown" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 1' ) ),
      cex = 1.2, col = "black" )
title( "B) M1a vs M2a" )

# M7 vs M8
curve( dchisq( x, df = 2 ), from = 0, to =  42 )
abline( v = c( Chisq.crit.M7vsM8, Chisq.crit.M7vsM8_2, diff_M7vsM8 ), col = c( "darkgray", "brown", "red" ) )
coords_dev     <- c( 20, 0.48 )
coords_pval    <- c( 19, 0.45 )
coords_alphac  <- c( 20, 0.41 )
coords_alphac2 <- c( 20, 0.38 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = 39.14911', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["2,0.05"]^"2", "= 5.99", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_alphac2[1], y = coords_alphac2[2],
      labels = expression( atop( paste( chi["2,0.01"]^"2", "= 9.21", sep = " " ) ) ),
      cex = 1.2, col = "brown" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 3.15e-09' ) ),
      cex = 1.2, col = "black" )
title( "C) M7 vs M8" )
