# Column information for the mapping files "16s_map_rainx.txt" and "ITS_map_rainx.txt"

```
SampleID(2)- Treatment and rep, abbreviated by letter codes. C/S=Clay or Sand, R/S=Rain or Sterile treatment, R/A=Reduced or Ambient rainfall
Region- Region of the DNA sequenced
BarcodeSequence- Short unique barcode sequence used for multiplexing
LinkerPrimerSequence- Linker and primer sequence used to target and amplify the region of interest	
InputFastaFileName- Name of fasta file for each individual sample
ProjectName- Name of project that the sequencing was conducted for	
Treatments- Treatments combined into one string with "space" a delimiter
SoilType- Two soil types (Clay- higher percentage of clay, Sand- higher percentage of sand, Rain- rain sample)
SampleType- Identifier of soil versus rain sample
RainLevel- How much water soils received (Ambient= ambient amount of rainfall. Reduced= drought treatment, 25% of ambient)
RainCom- Whether the soil received rain water with microbes (NonSterile- unfiltered rain, Sterile- rain filtered with 0.22 µm filter, Rain- rain sample)
SampleTime- Time of sample collection (End- Soils collected at the end of the experiment, STL- Autoclave sterilized soils collected at the end of the experiment, Start- Soil samples collected at beginning of experiment and archived to represent initial community, August through February- month of rain collection)
SampleTimeTrunc- Time of sample collection with months truncated to three letters
SampleDate-	Date of rain collection in serial number format
SampleDate2- Date of rain collection in mm/dd/yyyy format
precip-	Amount of rain that fell on the day of collection from MET data (https://lter.kbs.msu.edu/datatables/392)
amt_filter- Estimate amount of rain filtered for DNA extraction
com.rain.time- String combination of RainCom, RainLevel, and SampleTime treatments
PW- Identifier of Blank Power Water extraction kit samples used as controls
```