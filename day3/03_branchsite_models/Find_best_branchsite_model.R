#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------#
# SET ENVIRONMENT #
#-----------------#
# You will need to modify the path to `Site_models` if you run this
# script elsewhere!
wd <- paste( "~/my_session/day3/03_branchsite_models/" )
setwd( wd )

#-----------#
# LOAD DATA #
#-----------#
# 1. Load our text file with the lnL values
#    Order in columns: 1) BS-Bird-Alt | 2) BS-Bird-Null
#                      3) BS-DuckChicken-Alt | 4) BS-Duckchicken-Null
lnL_vals <- read.table( file = "lnL_branchsite_mods.txt", sep= " ",
                        stringsAsFactors = FALSE, header = FALSE )
colnames( lnL_vals ) <- c( "BS-bird-alt", "BS-bird-null",
                           "BS-duckchicken-alt", "BS-duckchicken-null" )

# 2. We can now compute the LRT statistic.

## 2.1. BS-duckchicken-null vs BS-duckchicken-alt ##
diff_BS_duckchick <- 2*( lnL_vals$`BS-duckchicken-alt`[1] - lnL_vals$`BS-duckchicken-null`[1] )
diff_BS_duckchick
# diff = 15.81525
# Find degree of freedom
df_BS_duckchick <- lnL_vals$`BS-duckchicken-alt`[2] - lnL_vals$`BS-duckchicken-null`[2]
df_BS_duckchick
# 1 degree of freedom
pchisq( diff_BS_duckchick, df = df_BS_duckchick, lower.tail=F )
# p-val = 6.98375e-05 < 0.05

## 2.2. BS-bird-null vs BS_bird-alt ##
diff_BS_bird <- 2*( lnL_vals$`BS-bird-alt`[1] - lnL_vals$`BS-bird-null`[1] )
diff_BS_bird
# diff = 11.07439
# Find degree of freedom
df_BS_bird <- lnL_vals$`BS-bird-alt`[2] - lnL_vals$`BS-bird-null`[2]
df_BS_bird
# 1 degree of freedom
pchisq( diff_BS_bird, df = df_BS_bird, lower.tail=F )
# p-val = 0.0008752812 < 0.05

# NOTE: The LRT statistic should be compared with the 50:50 mixture of
# point mass 0 and 1,5%2=2.71 and 1,1%2=5.41(Self and Liang 1987). In that
# way, we would need to divide the value obtained from `pchisq` into 2
# and use this as a p-value.
# Nevertheless, we use critical values \chi_{1,5%}^2=3.84 and
# \chi_{1,1%}^2=5.99 to guide against violations of model assumptions as
# recommended in the PAML documentation
Chisq.crit1 <- 3.84
Chisq.crit2 <- 5.99

# 3. Plot results
par( mfrow = c( 1,2 ) )

# DuckChicken branch
curve( dchisq( x, df = 1 ), from = 0, to =  16 )
abline( v = c( Chisq.crit1, Chisq.crit2, diff_BS_duckchick ), col = c( "darkgray", "brown", "red" ) )
coords_dev     <- c( 12.4, 0.85 )
coords_pval    <- c( 12.15, 0.78 )
coords_alphac  <- c( 12, 0.67 )
coords_alphac2 <- c( 12, 0.60 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = 15.82', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["1,0.05"]^"2", "= 3.84", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_alphac2[1], y = coords_alphac2[2],
      labels = expression( atop( paste( chi["1,0.01"]^"2", "= 5.99", sep = " " ) ) ),
      cex = 1.2, col = "brown" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 6.98e-5' ) ),
      cex = 1.2, col = "black" )
title( expression( 'A) '*italic(DuckChicken)*': branch-site model A VS branch-site model A with '*omega*'=1' ) )

# Bird branch
curve( dchisq( x, df = 1 ), from = 0, to =  16 )
abline( v = c( Chisq.crit1, Chisq.crit2, diff_BS_bird ), col = c( "darkgray", "brown", "red" ) )
coords_dev     <- c( 12.4, 0.85 )
coords_pval    <- c( 12.15, 0.78 )
coords_alphac  <- c( 12, 0.67 )
coords_alphac2 <- c( 12, 0.60 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = 11.07', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["1,0.05"]^"2", "= 3.84", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_alphac2[1], y = coords_alphac2[2],
      labels = expression( atop( paste( chi["1,0.01"]^"2", "= 5.99", sep = " " ) ) ),
      cex = 1.2, col = "brown" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 8.8e-4' ) ),
      cex = 1.2, col = "black" )
title( expression( 'B) '*italic(Bird)*': branch-site model A VS branch-site model A with '*omega*'=1' ) )