This is Hapler, version 1.61. Hapler is a tool for producing 
robust haplotype consensus regions given multiple alignments (assemblies, 
really) of genetically diverse sequence data. Hapler compares each sequence 
to every other, and groups sequences together into sets that don't have any 
conflicts (minimum coloring of the sequence 'conflict graph'). This can be 
done in O(n^3) time, because Hapler assumes that sequences contain no gaps 
(e.g., it ignores mate-pair information. Future versions may allow gaps at 
the risk of producing incorrect haplotypes.) Because such a minimum coloring 
is usually not unique, Hapler allows the user to select many pseudo-random 
colorings (using the --random-repetitions parameter, default 10) and only 
keep haplotype groupings which are common to all. In practice this drastically 
increases the correctness of results.

 Example Usage: cat TIGR_alignment.tigr | java -jar Hapler.jar
 Example Usage 2: java -jar Hapler.jar --input SAM_alignment.sam --alignment-type sam --snp-caller 454 --random-repetitions 100

Hapler works on any alphabet (RNA, DNA, Protein, ...). '-' characters are 
treated as conflict-causing alleles, '~' characters are treated as unknown 
gaps (that are not, by default, allowed.)


Option                                  Description                            
------                                  -----------                            
--alignment-type <Multiple Alignment    (default: tigr)                        
  type: tigr or ...>                                                           
--allow-gaps [One of false, split,      (default: false)                       
  true. The default option, false,                                             
  ensures that no gaps (unknown bases,                                         
  such as those between mate pairs)                                            
  are allowed. If sequences do have                                            
  haps, you may either split sequences                                         
  at gaps into seperate sequences (e.                                          
  g., ignoring mate-pair information),                                         
  or allow gaps in sequences. If gaps                                          
  area allowed, the algorithm _may_                                            
  produce incorrect results, in the                                            
  form of internally inconsistent                                              
  haplotypes. These will be marked, in                                         
  this case.]                                                                  
--binomial-alpha [The alpha value to    (default: 0.05)                        
  use when calling SNPs with the                                               
  binomial SNP caller. Significance is                                         
  determined on a per-SNP basis, and                                           
  is multiply corrected by the number                                          
  of tests run (which equals the                                               
  length of the alignment).]                                                   
--binomial-error-rate [The sequencing   (default: 0.005)                       
  error rate to assume when calling                                            
  SNPs with the binomial SNP caller.]                                          
--compute-reconstructions [If this is   (default: true)                        
  set to 'false,' no consensus                                                 
  reconstructions will be computed.                                            
  This can save computational time.]                                           
--evaluate-contigs [Given a fasta file                                         
  of contigs with ids matching                                                 
  multiple alignment names in the                                              
  input, reports how 'good' that                                               
  contig is, similar to how the                                                
  consensus is reconstructed and                                               
  evaluated. It is assumed the start                                           
  of the contigs corresponds to the                                            
  start of the alignment (position 0).]                                        
--ground-truth [Given a fasta file of                                          
  the actual, full haplotypes the data                                         
  was drawn from, hapler will output                                           
  extra data describing how each                                               
  haplotype reconstruction matches                                             
  them.]                                                                       
--help [Show help text.]                                                       
--human-readable [If this is true, the  (default: true)                        
  output is more or less human                                                 
  readable (haplotype regions are                                              
  aligned against each other), if this                                         
  is false, this is not the case,                                              
  which can drastically reduce the                                             
  amount of characters that need to be                                         
  output]                                                                      
--input <Input file to read. If not     (default: -)                           
  specified, or specified as -, read                                           
  from standard input.>                                                        
--max-repetitions [Maximum number of    (default: 100)                         
  repetitions that will be run,                                                
  regardless of --random-repetitions.]                                         
--maximize-one-read-haps [Maximize the  (default: true)                        
  number of haplotypes containing only                                         
  one non-redundant sequence: true or                                          
  false. Defaults to true.]                                                    
--pos-info-tr-reduction [With this      (default: false)                       
  option, the overlap graph will only                                          
  allow overlaps if they share and                                             
  agree at some SNP; further, we'll                                            
  compute the transitive reduction of                                          
  it before path covering (this no                                             
  longer corresponds to coloring, but                                          
  should produce faster and more                                               
  accurate results!)]                                                          
--random-repetitions [Number of times   (default: auto)                        
  to randomize sequence order and re-                                          
  run the bipartite matching/coloring;                                         
  the comonalities amongst colorings                                           
  are kept, thus this option controls                                          
  the conservativeness of the                                                  
  haplotype predictions. If 'auto'                                             
  (the default) is used, the coloring                                          
  will be repeated until no                                                    
  functionally new colorings are                                               
  discovered for 10 repetitions or 10%                                         
  of repetitions, whichever is                                                 
  smaller.]                                                                    
--show-alignments [Outputs the                                                 
  alignments in human-friendly format.]                                        
--snp-caller <SNP caller to use:        (default: simple)                      
  simple, simplestrict, 454 or                                                 
  binomial. Defaults to binomial.>                                             
--snp-list <File containing a column                                           
  of integers, where each integer                                              
  specifies a position which should be                                         
  considered a SNP in the multiple                                             
  alignment (0 indexed). Overrides --                                          
  snp_caller.>                                                                 
--version [Show version information.]                                          
