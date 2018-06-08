
# Genomer (Genome Processor)

Used to make simple modifications to FASTA formatted genome files.

## Modes

### Collapse

For collapsing multiline files into one. A FASTA comment will cause a new block to start.

### Extract

For extracting a subsequence within a file at specified start and end points. Ignores FASTA comments.

### Analyse

Provides some basic statistics about a FASTA file.

### Refadjust

Takes the UCSC refGene MySQL table (as a tab separated file) and modifies the start and end points according to an offset (useful for after running 'extract' mode).

This realigns the annotation with modified genome files.

RefAdjust calculates the new position by: `original position - offset`

`offset`: generally, the start position used in extract mode

## Verification method

Note: this to verify that collapsing, modifying an annotation and extracting DNA sequences is done correctly.

 1) Visit this URL to view GRCm38/mm10 chr19:16,767,421-20,818,303 in the UCSC
       Genome Browser. It will show Foxb2 gene. Follow link below.
       
[UCSC Browser - GRCm38/mm10 Chr19](https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A16767421-20818303&hgsid=674918601_YLI5ZKYBVVEjC29faHGp1IVNABk8)
 
 2) Right click on Foxb2 and:
 
       - Get DNA for Foxb2, AND
       
       - Open details page in new window; OR
       
       - Follow the provided links below
       
[Get DNA for Foxb2](https://genome.ucsc.edu/cgi-bin/hgc?hgsid=674918601_YLI5ZKYBVVEjC29faHGp1IVNABk8&g=htcGetDna2&table=&i=mixed&l=16872315&r=16873830&getDnaPos=chr19%3A16%2C872%2C316-16%2C873%2C830&db=mm10&hgSeq.cdsExon=1&hgSeq.padding5=0&hgSeq.padding3=0&hgSeq.casing=upper&boolshad.hgSeq.maskRepeats=0&hgSeq.repMasking=lower&boolshad.hgSeq.revComp=0&submit=get+DNA)

[Open details page in new window](https://genome.ucsc.edu/cgi-bin/hgGene?hgg_gene=uc008gxc.1&hgg_prot=uc008gxc.1&hgg_chrom=chr19&hgg_start=16872315&hgg_end=16873830&hgg_type=knownGene&db=mm10&c=chr19&l=16767420&r=20818303)


 3) Run genomer.py in collapse mode to collapse chr19 of mm10.
       
       Obtain chr19 data:
           
           UCSC data: `ftp://hgdownload.cse.ucsc.edu`
           
           Navigate to: `/apache/htdocs/goldenPath/mm10/chromosomes`
           
           Download: `chr19.fa.gz`
       
       ```python genomer.py -m collapse -f "C:\genomes\mm10-ucsc-mod\chr19\chr19.fa"```

       You can verify the file was collapsed using the analyse mode

       ```python genomer.py -m analyse -f  "C:\genomes\mm10-ucsc-mod\chr19\chr19.fa"```

 4) Run genomer.py in extract mode and extract 10m to 20m BP from the chr19 UCSC genome data

       ```python genomer.py -m extract -s 10000000 -e 20000000 -f "C:\genomes\mm10-ucsc-mod\chr19\chr19.fa.collapse"```

 5) Run genome.py in refadjust mode to adjust the UCSC genome annotation file

       Obtain annotation:
           
           MySQL server: genome-mysql.soe.ucsc.edu (username: genome, password: <your-email>)
           
           Find `mm10` database. Find `refGene` table. Export as tab separated file. 
           
           This is a simple task using MySQL Workbench (use Table Data Export Wizard)
           
           In our example, we have named this TSV file `refGene.txt`

       ```python genomer.py -m refadjust -f "C:\genomes\mm10-ucsc-mod\chr19\refGene.txt" -o 10000000 -l 50```
           
           -o flag (offset flag) should be equal to the -s flag in step 3
           
           -l flag indicates length of each line in original UCSC genome data file
       
 6) Find Foxb2 in the adjusted annotation file.
      
       Note: each numerical value in "adjusted" is 10m less than in "original".

       Original:
       
       `713	NM_008023	chr19	-	16872315	16873830	16872353	16873640	1	16872315,	16873830,	0	Foxb2	cmpl	cmpl	0,`
       
       Adjusted:
       
       `713	NM_008023	chr19	-	6872315	6873830	6872353	6873640	1	6872315,	6873830,	0	Foxb2	cmpl	cmpl	0,`

 7) Run genomer.py in extract mode and extract Foxb2 from the output generated in step 4
       
       Note: -s and -e flags match the adjusted gene start and end values from step 6
           
           Refer to the MySQL table structure for what each column represents

       ```python genomer.py -m extract -f "C:\genomes\mm10-ucsc-mod\chr19\chr19.fa.collapse.extract" -s 6872315 -e 6873830```    

 8) Verify the sequence provided in "Get DNA for Foxb2" (step 2) from UCSC Genome Browser matches the extract produced in step 7.

       UCSC "Get DNA for Foxb2":
           `TCTCTCGACA`...`GGTCCCCGCA`

       Genomer.py process:
           
           Note: if you followed this process precisely, the final file you should look at is called:
           
               `<dirs>/chr19.fa.collapse.extract.extract`
               
           `TCTCTCGACA`...`GTCCCCGCAA`

 9) Verification complete
