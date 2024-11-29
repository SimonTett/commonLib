Code from Simon Brown (Met Office). This code filters out spurious extremes from the data
Sub_filter_hourly.bash calls
                Batch.filter_hourly.R calls
                                filt_point_vbar_hbar_wnorm() in filter_fn.R

Key numbers in Batch.filter_hourly.R are:

min_pr_th <- 10    # dont bother with fields who’s max rain rate is below this value, mm/h
rim       <- 54   # number of grid points to ignore around the rim
badb      <- 2     # no grid points to remove near bad data – long dimension
bads      <- 1     # no grid points to remove near bad data – short dimension

# thresholds for hbar, vbar and point events
if(stseas=='djf') {
    badth   <- 180 # hbar or vbar over this is bad dta     
    badthp  <- 12  # point over this is bad dta            
    }
if(stseas=='jja') {
    badth   <- 160 #  hbar or vbar over this is bad dta  
    badthp  <- 31  #  point over this is bad dta         
}

The normalisation I do is to take the square root of the returned convolution value – this seems to work quite well at getting some of the weaker events.

NB if you do try and run the R code note that the filtered data is written into the original files OVERWRITING the original.