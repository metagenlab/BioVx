# Documentation for script: _tmsplit_

## Summary
This program cuts protein sequences under 4 different criteria: 1) cut a residue range 


## Dependencies
The following Python module needs to be available to the script: 

1. **tmsFunction**  
This module is part of the distribution.  


## Command line options
The following options are available. You can also run the 
script without arguments to display the options:


    -h, --help            show this help message and exit.  
    -if INPUTFILE         The path of the file containing the FASTA sequences 
                          that are to be processed.  
    -od OUTPUTDIR         The path to the directory for the output. If the 
                          directory is not present,it will be created.  
    -of OUTFILE           Name of the file in which the processed sequence will 
                          be placed. NOTE: Do not add any extensions to the  
                          name. ".faa" will automaticaly be assigned.  
    -rangeCut             Extract a segment of variable size from a protein  
                          sequence. Start(-s or --start) and End(-e or --end)  
                          required for range cut.  
    -equal                Cut protein(s) into segments of equal length. Parts 
                          (-p or --parts) required for equal cut.  
    -tmsCut               Cut multiple sequences between any two TMS's. Start(-s 
                          or --start) and End(-e or --end) required for TMS 
                          cut.The tail argument is optional, but the default is 3.  
    -split                Split protein(s) into groups of TMS. Choose  
                          overlapping(-over) or non-overlapping(-non). DEFAULT  
                          is non-overlapping. TMS(-tms) argument required for  
                          TMS split.The tail argument is optional, but the  
                          default is 3.  
    -t TAIL, --tail TAIL  The number of residues left before and after the  
                          sequence segments. DEFAULT is 3.  
    -p PARTS, --parts PARTS  
                          Argument for equal cut. The number of segments of  
                          equal length the sequence(s) should be cut into.  
    -s START, --start START 
                          Argument for TMS cut. The number of the TMS or Residue  
                          that will be the first in the segment.  
    -e END, --end END     Argeument for TMS cut. The number of the TMS or  
                          Residue that will be the last in the segment.  
    -tms TMS              Argument for TMS Split. The number of TMS per group.  
    -over                 Argument for TMS Split. TMS are grouped by overlapping  
                          sections. Ex [1,2,3],[2,3,4] etc...  
    -non                  DEFAULT Argument for TMS Split. TMSs are grouped by  
                          non-overlapping sections. Ex [1,2,3],[4,5,6] etc...  
                          Also requires a condition -ignore, -append, or -new  
    -ignore               Argument for TMS Split Non-Overlapping. Any TMS that  
                          are not able to be formed into a full group are  
                          ignored.  
    -append               DEFAULT Argument for TMS Split Non-Overlapping.Any TMS  
                          that are not able to be formed into a full group are  
                          added to the last possible segment.  
    -new                  Argument for TMS Split Non-Overlapping. Any TMS that  
                          are not able to be formed into a full group are added  
                          to another segment


