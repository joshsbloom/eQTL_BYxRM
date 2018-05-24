library(reshape)
library(gdata)

# total.DNA = target OD in final plate
# rxn.vol = target final vol in final plate
# pre.dilun = whole plate dilution factor (i.e. 5 = 1:5 dilution before running on plate reader)
parameters <- list(total.DNA = .05, rxn.vol = 1000, pre.dilun=3)

#DNA = OD of overnight culture
#total.DNA = target OD of new culture
#DNAneeded = volumne

# put xlsx files (8 rows by 12 columns each) of 96-well plates in a directory
dir.in = '/data/eQTL/OD_dilution/12_13/'
fs = list.files(dir.in, pattern='xlsx')
blank=0.07875


attach(parameters)
for (filein in fs) { 

    DNAfile = paste(dir.in, filein, sep='')
    print(DNAfile)
    outfile = gsub('xls', 'txt', DNAfile)


    intensity = read.xls(DNAfile, header=F,sheet=1)
    DNA  = intensity[,1:12] - blank
    DNA[DNA==0]=0.000001

    DNAneeded <- (total.DNA*rxn.vol)/(DNA*pre.dilun)
    DNAneeded[DNAneeded>rxn.vol]=rxn.vol
    
    # Reconfigure info into lists compatible with the robot
    positions <- as.character(melt(
                                   sapply(1:12, function(num) {
                                          sapply( toupper(letters[1:8]),
                                            function(letter) {paste(letter, num, sep = "")}, USE.NAMES = FALSE)}, 
                                                                USE.NAMES = FALSE))$value)
    DNAneeded <- melt(DNAneeded)$value
    H2Oneeded <- (rxn.vol - DNAneeded)
    robot.df <- data.frame(cultureplate = "Culturesource", mediasource = "YNBsource", newplate = "NewCulture", 
            startpos = positions, Cultvol = DNAneeded, YNBvol = H2Oneeded, endpos = positions, stringsAsFactors = FALSE)
    write.table(robot.df, outfile, col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
    
    #watch out, not installed by default on OS X
    system(paste('unix2dos', outfile))
}

