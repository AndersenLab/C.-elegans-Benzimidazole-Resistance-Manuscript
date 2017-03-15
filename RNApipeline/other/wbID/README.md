##Wormbase ID converter version 0.02  
April 2015, Damien O'Halloran, The George Washington University  

##Getting Started
1. Place a file containing a list of _C. elegans_ gene identifiers to be converted in the appropriate dir  
2. File should have one identifier per line (for example see: _in_names.txt_)  
3. _names.csv_ contains all genes identified as: WBID,gene_name,transcript_name from Wormbase Ver.WS240  
4. Converted file will have only unique entries

##Usage  
  `perl wormbase_ID_converter.pl -a gene -b WBID`  

##GetOpts
```per
  -a input_type # 'gene' OR 'transcript' OR 'WBID'  
  -b output_type # 'gene' OR 'transcript' OR 'WBID' 
```
##Contemporize gene names  
To get the most up to date names from WormBase run "get_updated_names.pl" script  
You will be asked to enter the version of WormBase names you want e.g. WS237 or WS240  
Visit here to get available versions: ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/annotation/geneIDs/ 

##Usage to contemporize gene names  
  `perl get_updated_names.pl`  
  
## Contributing
All contributions are welcome.

## Support
If you have any problem or suggestion please open an issue [here](https://github.com/dohalloran/wormbase_ID_converter/issues).

## License 
The MIT License
  
