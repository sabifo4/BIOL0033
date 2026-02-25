#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------#
# SET ENVIRONMENT #
#-----------------#
# You will need to modify the path to `Site_models` if you run this
# script elsewhere!
wd <- paste( "~/my_session/day3/02_branch_models/" )
setwd( wd )

#-----------#
# LOAD DATA #
#-----------#
# 1. Load our text file with the lnL values
#    Order in columns: Branch-ChickenDuck | Branch-Bird | M0
lnL_vals <- read.table( file = "lnL_branch_mods.txt", sep= " ",
                        stringsAsFactors = FALSE, header = FALSE )
colnames( lnL_vals ) <- c( "Bird-branch", "DuckChicken-branch", "M0-nobranch" )

# 2. We can now compute the LRT statistic

## 2.1. DuckChicken-branch vs M0-nobranch ##
diff_duckchickVSM0nob <- 2*( lnL_vals$`DuckChicken-branch`[1] - lnL_vals$`M0-nobranch`[1] )
diff_duckchickVSM0nob
# diff = 62.10233
# Find degree of freedom
df_duckchickVSM0nob <- lnL_vals$`DuckChicken-branch`[2] - lnL_vals$`M0-nobranch`[2]
df_duckchickVSM0nob
# 1 degree of freedom
pchisq( diff_duckchickVSM0nob, df = df_duckchickVSM0nob, lower.tail=F )
# p-val = 3.260654e-15 < 0.05
Chisq.crit.duckchick <- qchisq( p = 0.95, df = df_duckchickVSM0nob )
Chisq.crit.duckchick
# alpha critical value = 3.841459

## 2.2. Bird-branch vs M0-nobranch ##
diff_birdVSM0nob <- 2*( lnL_vals$`Bird-branch`[1] - lnL_vals$`M0-nobranch`[1] )
diff_birdVSM0nob
# diff = 62.1022
# Find degree of freedom
df_birdVSM0no <- lnL_vals$`Bird-branch`[2] - lnL_vals$`M0-nobranch`[2]
df_birdVSM0no
# 2 degree of freedom
pchisq( diff_birdVSM0nob, df = df_birdVSM0no, lower.tail=F )
# p-val = 3.270979e-14 < 0.05
Chisq.crit.bird <- qchisq( p = 0.95, df = df_birdVSM0no )
Chisq.crit.bird
# alpha critical value = 5.991465

# 3. Plot results
par( mfrow = c( 1, 2 ) )

# DuckChicken branch
curve( dchisq( x, df = 1 ), from = 0, to =  65 )
abline( v = c( Chisq.crit.duckchick, diff_duckchickVSM0nob ), col = c( "darkgray", "red" ) )
coords_dev    <- c( 28.6, 0.3 )
coords_pval   <- c( 29.4, 0.25 )
coords_alphac <- c( 28, 0.20 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = 62.10', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["1,0.05"]^"2", "= 3.84", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 3.26-15' ) ),
      cex = 1.2, col = "black" )
title( expression( 'A) '*italic(Duck-Chicken)*': Branch model VS M0 model ' ) )

# Bird branch
curve( dchisq( x, df = 1 ), from = 0, to =  65 )
abline( v = c( Chisq.crit.bird, diff_birdVSM0nob ), col = c( "darkgray", "red" ) )
coords_dev    <- c( 28.6, 0.3 )
coords_pval   <- c( 29.4, 0.25 )
coords_alphac <- c( 28, 0.20 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = 62.10', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["2,0.05"]^"2", "= 5.99", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 3.27-14' ) ),
      cex = 1.2, col = "black" )
title( expression( 'B) '*italic(Bird)*': Branch model VS M0 model ' ) )
