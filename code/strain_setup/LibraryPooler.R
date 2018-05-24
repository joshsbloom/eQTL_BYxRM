library(reshape)
library(gdata)

#target amount of DNA in ng in final pool for each sample
total.DNA=100
target.pos = 'A1'
# put xlsx files (8 rows by 12 columns each) of 96-well plates in this directory
dir.in = '/data/eQTL/dilution/02/'

# standards in ng
standardNGin = c(0,5,10,20,40,60,80,100)

# note excel extension, change to xls if necessary
fs = list.files(dir.in, pattern='xlsx')
for (filein in fs) { 

    DNAfile = paste(dir.in, filein, sep='')
    print(DNAfile)

    #output file will have same name as input but with txt extension
    outfile = gsub('xlsx', 'txt', DNAfile)

    #8 x 12 matrix of readings of sample fluorescence from plate reader
    intensity = read.xls(DNAfile, header=F,sheet=1)
    # 8x12 matrix of readings of standards fluorescence from plate reader
    standard  = read.xls(DNAfile, header=F,sheet=2)
    rownames(standard)=standardNGin
    
    #average the standard intensities
    # assuming standards are in first three columns, adjust as necessary
    g1=apply(standard[,c(1,3)], 1, mean)
    #fit line, toss highest standard
    g1model = lm(standardNGin[-c(8)]~g1[-c(8)]-1)

    #diagnose standard curve ---------------------------------------------------
    #plot(g1,standardNGin)
    #abline(g1model)
    #---------------------------------------------------------------------------
    

    DNA  = cbind( intensity[,1:12]*as.numeric(g1model$coeff[1]))
    # if less than 1 ng/ul then NA

    #diagnostics  ---------------------------------------------------------------
    #hist(unlist(DNA), breaks=40)
    #----------------------------------------------------------------------------

    DNA[DNA<1]=NA    
    DNAneeded = total.DNA/DNA 
    # Reconfigure info into lists compatible with the robot
    positions <- as.character(melt(
                                   sapply(1:12, function(num) {
                                          sapply( toupper(letters[1:8]),
                                            function(letter) {paste(letter, num, sep = "")}, USE.NAMES = FALSE)}, 
                                                                USE.NAMES = FALSE))$value)
    DNAneeded = melt(DNAneeded)$value
    DNAneeded[is.na(DNAneeded)]=0
    
    #this is crude ... check sanity of values for each pool
    
    # if volume of any one well is greater than 18 then break for loop and error
    if(sum(DNAneeded>18)) { print('ERROR, reduce total DNA requirement'); break;   }

    # if total volume greater than 1.3 mL then break for loop and error    
    if(sum(DNAneeded)>1300) { print('ERROR, reduce total DNA requirement'); break; }

    robot.df = data.frame(sourceplate = "DNAsource", receivingplate = "DNAreceiving", 
            startpos = positions, DNAvol = DNAneeded, endpos =target.pos, stringsAsFactors = FALSE)
    robot.df$DNAvol=round(robot.df$DNAvol)

    rand.ind=sort(sample(1:length(robot.df$DNAvol) ,45))
    robot.df$DNAvol[rand.ind]=robot.df$DNAvol[rand.ind]+.5
    outfile = "/data/eQTL/dilution/02/eqtl_02_half_ul_precision.txt"
    outfile = "/data/eQTL/dilution/02/eqtl_02_rounded.txt"
    # set eol to output DOS new line characters
    write.table(robot.df, outfile, col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE, eol='\r\n')
}

