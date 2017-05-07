# Welcome to the Extended repository of the BioV suite of programs

Most of the original programs in the suite were originally written by Vamsee Reddy, many have been updated and new programs have been added by different members of the lab and collaborators.

Find below the list of programs available. Scripts developed by our group as dependencies are not listed here but are included in the distribution:  

<ul>
  <li>protocol1.py</li>  
  <li>protocol2.py</li>  
  <li>gsat.py</li>  
  <li>gblast3.py</li>  
  <li>ancient.py</li>  
  <li>quod.py</li>  
  <li>hvordan.py</li>  
  <li>tmsplit.py</li>  
</ul>

---  

**Note 1:**  
> The manual is not updated so far.  

**Note 2:**  
> Some modifications to protocol 2 and other programs by Gabriel Moreno-Hagelsieb.   

**Note 3:**  
&nbsp;&nbsp;&nbsp;Pranav Iddamsetty and Arturo Medrano  
&nbsp;&nbsp;&nbsp;added the scripts tmsplit and tmsFunction.py.  
&nbsp;&nbsp;&nbsp;These programs provide a command-line interface  
&nbsp;&nbsp;&nbsp;to cut transporter sequences under the same   
&nbsp;&nbsp;&nbsp;criteria as the website:  

&nbsp;&nbsp;&nbsp;http://biotools.tcdb.org/bartms_split.html  

**Note 4:**  
&nbsp;&nbsp;&nbsp;The program gblast3.py was modified by  
&nbsp;&nbsp;&nbsp;Pranav Iddamsetty and Arturo Medrano to  
&nbsp;&nbsp;&nbsp;generate a tabulated file (results.tsv) with  
&nbsp;&nbsp;&nbsp;the same output as results.html. This will  
&nbsp;&nbsp;&nbsp;allow further automatic processing of gblast  
&nbsp;&nbsp;&nbsp;output (e.g. compare different outputs of gblast  
&nbsp;&nbsp;&nbsp;to identify commonalities/differences between  
&nbsp;&nbsp;&nbsp;trasporter systems in different genomes).   
  This program was also modified to add annotated  
  substrates in tcdb to the output. Currently the  
  substrates are read from flat files but as soon  
  as the substrates are uploaded to TCDB we will  
  modify the script to read the substrates directly  
  from the database.  

**Note 5:**  
  The programs quod.py, tcblast.py and hvordan.py were  
  added by Kevin Hendargo and Arturo Medrano. The script  
  quod.py runs WHAT from the command line on one or more  
  sequences and is able to generate the plots in different  
  formats and qualities. The script tcblast.py provides  
  functions to run blast against tcdb from the command  
  line and provides graphical display that will be used  
  by the last script hvordan.py, which runs blast and WHAT  
  in order to generated an html ouput file that will help  
  the user to make the biological interpretation of  
  protocol2 top hits.  

**For developers:**  
  Please note: There is no need to add a link for your python   
  binary pointing to '/usr/local/bin/python'  
