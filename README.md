
# FASTAr (FASTA Processor)

Used to make simple modifications to FASTA formatted FASTA files (generally provided by [UCSC](https://FASTA.ucsc.edu/))

Written by Jake Bradford. Checkout our lab website! [biomedicaldatascience.com](http://biomedicaldatascience.com)

## Modes

#### Collapse

For collapsing multi-line FASTA blocks files into single-line blocks. A FASTA comment will cause a new block to start.

#### Extract

For extracting a sequence within a file at specified start and end points. Ignores FASTA comments.

#### Analyse

Provides some basic statistics about a FASTA file.

#### Refadjust

Realigns UCSC refGene annotations according to an offset. Useful for after extracting a sequence from a FASTA file and maintaining the correct annotation for that region.

#### Refgeneextract

Extract gene names between start and end position using the given reference.

#### Singularise

Converts a single multi-FASTA formatted file into multiple single FASTA files.

## Files Required

#### fastar.py

By cloning this repository: `git clone https://github.com/jakeb1996/FASTAr.git`

#### FASTA sequence

Obtained from UCSC FTP server ([help](https://FASTA.ucsc.edu/goldenpath/help/ftp.html)).

Server: `ftp://hgdownload.cse.ucsc.edu` (user: `anonymous` pass: `<your-email>`)

Chromosome example: `/apache/htdocs/goldenPath/mm10/chromosomes/chr19.fa.gz`

#### Annotation file

Obtained from UCSC MySQL server ([help](https://FASTA.ucsc.edu/goldenpath/help/mysql.html))

MySQL server: FASTA-mysql.soe.ucsc.edu (username: FASTA, password: <your-email>)
           
Find `mm10` database. Find `refGene` table. Export as tab separated file. 
           
This is a simple task using the Table Data Export Wizard in MySQL Workbench ([help](https://dev.mysql.com/doc/workbench/en/wb-admin-export-import-table.html)).

## Verification

**1) Understand the context of this method**

View GRCm38/mm10 chr19:16,767,421-20,818,303 in the UCSC FASTA Browser. It will show Foxb2 gene. [UCSC Browser - GRCm38/mm10 Chr19](https://FASTA.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A16767421-20818303&hgsid=674918601_YLI5ZKYBVVEjC29faHGp1IVNABk8)
 
This verification method provides an example of:

- Collapsing a FASTA file into a single line (this way line breaks are ignored; as *some* CRISPR related tools do anyway)

- Extracting a sequence from the collapsed file (eg: you may only want a sample file that is 10 million base-pairs in length)

- Updating an annotation file so that the annotations still align with your extract

*Note: this method uses the Foxb2 gene as an example (the gene is very short)*
 
**2) Obtain Foxb2 gene info**

In the UCSC FASTA Browser, right click on Foxb2 and:
 
- [Get DNA for Foxb2](https://FASTA.ucsc.edu/cgi-bin/hgc?hgsid=674918601_YLI5ZKYBVVEjC29faHGp1IVNABk8&g=htcGetDna2&table=&i=mixed&l=16872315&r=16873830&getDnaPos=chr19%3A16%2C872%2C316-16%2C873%2C830&db=mm10&hgSeq.cdsExon=1&hgSeq.padding5=0&hgSeq.padding3=0&hgSeq.casing=upper&boolshad.hgSeq.maskRepeats=0&hgSeq.repMasking=lower&boolshad.hgSeq.revComp=0&submit=get+DNA); AND
       
- [Open details page in new window](https://FASTA.ucsc.edu/cgi-bin/hgGene?hgg_gene=uc008gxc.1&hgg_prot=uc008gxc.1&hgg_chrom=chr19&hgg_start=16872315&hgg_end=16873830&hgg_type=knownGene&db=mm10&c=chr19&l=16767420&r=20818303)

**3) Collapse chr19 data**

Run `fastar.py` in collapse mode to collapse chr19 of mm10 (collapsed multi-line DNA sequence into single-line DNA sequence).
       
If you need CHR19:
           
	1) UCSC data: `ftp://hgdownload.cse.ucsc.edu`

	2) Navigate to: `/apache/htdocs/goldenPath/mm10/chromosomes`

	3) Download: `chr19.fa.gz`
       
```$ python fastar.py -m collapse -f "C:\FASTAs\mm10-ucsc-mod\chr19\chr19.fa"```

You can verify the file was collapsed using the analyse mode

```$ python fastar.py -m analyse -f  "C:\FASTAs\mm10-ucsc-mod\chr19\chr19.fa"```

**4) Extract large portion of Chr19**

Run fastar.py in extract mode and extract 10m to 20m BP from the chr19 UCSC FASTA data

```$ python fastar.py -m extract -s 10000000 -e 20000000 -f "C:\FASTAs\mm10-ucsc-mod\chr19\chr19.fa.collapse"```

**5) Adjust the annotation file**

Seeing that we just took an extract from the chromosome file, the annotations will now be out of alignment.

Run `fastar.py` in refadjust mode to adjust the UCSC FASTA annotation file

Obtain annotation:
           
	1) MySQL server: `FASTA-mysql.soe.ucsc.edu` (username: `FASTA`, password: `<your-email>`)
           
	2) Find `mm10` database. Find `refGene` table. Export as tab separated file. 
           
	3) This is a relatively simple task using MySQL Workbench (use Table Data Export Wizard)
           
In our example, we have named this TSV file `refGene.txt`

```$ python fastar.py -m refadjust -f "C:\FASTAs\mm10-ucsc-mod\chr19\refGene.txt" -o 10000000 -l 50```
           
`-o` flag (offset flag) should be equal to the `-s` flag in step 4
           
`-l` flag indicates length of each line in original UCSC FASTA data file
       
**6) Find Foxb2 in the adjusted annotation file**
      
Note: each numerical value in "adjusted" is 10m less than in "original".

Original: `713	NM_008023	chr19	-	16872315	16873830	16872353	16873640	1	16872315,	16873830,	0	Foxb2	cmpl	cmpl	0,`
       
Adjusted: `713	NM_008023	chr19	-	6872315	6873830	6872353	6873640	1	6872315,	6873830,	0	Foxb2	cmpl	cmpl	0,`

**7) Extract gene from custom chr19 file**

Run `fastar.py` in extract mode and extract Foxb2 from the output generated in step 4
       
Note: `-s` and `-e` flags match the adjusted gene start and end values from step 6
           
Refer to the MySQL table structure for what each column represents

```$ python fastar.py -m extract -f "C:\FASTAs\mm10-ucsc-mod\chr19\chr19.fa.collapse.extract" -s 6872315 -e 6873830```    

**8) Verify sequences match**

Verify the sequence provided in "Get DNA for Foxb2" (step 2) from UCSC FASTA Browser matches the extract produced in step 7.

UCSC "Get DNA for Foxb2": `TCTCTCGACA`...`GGTCCCCGCA`

fastar.py process: `TCTCTCGACA`...`GTCCCCGCAA`
           
Note: if you followed this process precisely, the final file you should look at is called: `<dirs>/chr19.fa.collapse.extract.extract`

**9) Verification complete**

Congratulations!

# License

* See LICENSE *