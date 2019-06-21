#Column information for the mapping files "16s_map_rainx.txt" and "ITS_map_rainx.txt"

```
SampleID(2)- Unique idea for every sample
Region- Region of the DNA sequenced
BarcodeSequence- Short unique barcode sequence used for multiplexing
LinkerPrimerSequence- Linker and primer sequence used to target and amplify the region of interest	
InputFastaFileName- Name of fasta file for each individual sample
ProjectName- Name of project that the sequencing was conducted for	
Treatments- Treatments combined into one string with "space" a deliminator
SoilType- Soil characteristic of the soil core under analyses (Clay- higher percentage of clay, Sand- higer percentage of sand, Rain- rain sample)
SampleType- Identifier of soil versus rain sample
RainLevel- Watering level soil cores received (Ambient- ambient level of collected rain, Drought- 25% of ambient, Rain- rain sample)
RainCom- Disperal treatment each core received (NonSterile- unfiltered rain, Sterile- rain filtered with 0.22 µm filter, Rain- rain sample)
SampleTime- Time of sample collection (End- Soils collected at the end of the experiment, STL- Autoclave sterilized soils collected at the end of the experiment, Start- soil samples collected at beginning of experiment and archived to represent intial community, August through Febuary- month of rain collection)
SampleTimeTrunc- Time of sample collection with months truncated to three letters
SampleDate-	Date of rain collection in serial number format
SampleDate2- Date of rain collection in mm/dd/yyyy format
precip-	amount of rain that fell on the day of collection from MET data (https://lter.kbs.msu.edu/datatables/392)
amt_filter- estimate amount of rain filtered for DNA extraction
com.rain.time- string combination of RainCom, RainLevel, and SampleTime treatments
PW- indentifier of Blank power water samples used as controls
```