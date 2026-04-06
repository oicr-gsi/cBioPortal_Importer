# cBioPortal Importer

Script pycbio.py is used to generate an import folder with all the data and metadata files required for uploading data to cBioPortal.

The script is available as a module:

```module load cbioportal-importer```

The module will also load accessory tools in the environment required for processing and annotating mutations.


## Mapping file

Currently, data accepted for CbioPortal uploads are maf files from the VEP workflow, .seg files from sequenza or purple, .genes.results files from the rsem workflow and .tab from the mavis workflow.

The data should be organized in a comma-separated map.csv file with the following information:

`patient_id,sample_id,maf_file.maf.gz,seg_file.seg,rsem.genes.results,mavis.tab`

The mapping file can be generated using the cbioportal importer with an input json that includes donor and sample Ids and workflow output files. This json can also be downloaded from Waterzooi.


Parameters of the ```map``` function

| argument | purpose | required/optional                                    |
| ------- | ------- | ------------------------------------------ |
| -o | Path to the output directory  | required              |
| -i | Path to the json file with input data   | required              |
| -g | Gamma value of the sequenza output. Default is 500 | Required with default value              |


```cbio_importer map -o path/to/the/output/directory -i path/to/the/json/file/with/input/data -g gamma/value/of/the/sequenza/output(=500 by default)```


## Generating the import folder

Options, including path the output directory, path to the mapping file and filters can be specified in the config file.

Parameters of the ```generate``` function

| argument | purpose | required/optional                                    |
| ------- | ------- | ------------------------------------------ |
| -cf | Path to configuration file  | required              |
| -cl | Path to sample clinical information file   | optional              |



Generate the import folder with the following command: 

```cbio_importer generate -cf /path/to/the/config/file```


By default the only required clinical sample information are the patient and sample identifiers, which originate from the mapping file map.csv. It is possible however to add user-defined clinical sample information. This information is displayed in the Clinical Data tab of the project dashboard in cBioPortal. Clinical information should be organized in a tab-delimited table. The first two columns should be labeled Patient and Sample and must contain the same patient and sample identifiers as in the map.csv file. The Patient and Sample IDs are used to map clinical information. Any other column names are valid but each column must contain a single data type (eg, boolean, string or number). The column names provided in the clinical information file will be displayed in cBioPortal. 


```cbio_importer generate -cf /path/to/the/configfile -cl /path/to/the/clinicalfile```


## Merging import folders


Data uploaded to cBioPortal will replace any data already on the server. Import folders can be generated in batches and subsequently merged when data becomes available to add new data to cBioPortal without replacing the uploaded data. One important assumption is that all merging import folders were processed using the same parameters described in the config file. Nevertheless, it is possible to force merging when some of these parameters have differences. The only parameter required to be consistent is the genome reference.  

The tool takes a list of config files and then compares the parameters used to generate the import folders. If differences are found, it will not merge the import folders and will emit a warning with a list of parameters that show different values. The user can then review the parameters and force merging using the flag ```--merge_with_differences``` if these differences are estimated to be negligeable and/or the original data does not exist. However, this flag does not force merging if the genome reference differs. 

The paths to the config files and import folders can be specified using parameters ```-configs``` and ```-importFolders``` respectively. For convenience, paths to the config files and import folders can be specified in a file when several import folders need to be merged using options ```-configPaths``` and ```-importPaths```. 


Parameters of the ```merge``` function

| argument | purpose | required/optional                                    |
| ------- | ------- | ------------------------------------------ |
| -configs | White space separated list of config file paths | required              |
| -importFolders | White space separated list of of paths to import folders | required              |
| -configPaths | File with paths to config files   | required              |
| -importPaths | File with paths to import folders   | required              |
| -cancerCode | Cancer code. see  http://oncotree.mskcc.org    | required              |
| -genome | Genome reference hg19 or hg38   | required              |
| -outdir | Path to the output directory   | required              |
| -center | Genomic center   | required              |
| -project | Project name   | required              |
| -description | Project description   | required              |
| -study | Study name in the format: Project: Top-Level-OncoTree, Concept (PI, Centre)   | required              |
| -excludeSamples | White space separated list of excluded samples   | required              |
| -excludeFile | Path to the file with samples to exclude   | required              |
| --merge_with_differences | Merge import folders with parameter differences if activated   | optional              |


Example commands to merge import folders:

- merge import folders 

```
module load cbioportal-importer;
cbio_importer merge -configs path/to/config1 path/to/config2 path/to/config3 -importFolders path/to/import_folder1 path/to/import_folder2 path/to/import_folder3 -cancerCode MIXED -genome hg38 -outdir path/to/the/output/directory -center -project ProjectName -description "short description of the project" -study "ProjectName: Top-Level-OncoTree, Concept (PI, OICR)" 
```

- merge import folders specified in a file and enabling differences found in the config files

```
module load cbioportal-importer;
cbio_importer merge -configPaths path/to/configfiles -importPaths path/to/importfolders -cancerCode MIXED -genome hg38 -outdir path/to/the/output/directory -center -project ProjectName -description "short description of the project" -study "ProjectName: Top-Level-OncoTree, Concept (PI, OICR)"  --merge_with_differences
```

## Excluding samples when merging import folders


Samples previously included in an import folder can also be removed during merging if needed, for instance because of a sample swap or other considerations. The list of samples to be excluded can either be specified in the command using ```-excludeSamples``` or in a file using ```-excludeSamples```. The samples to exclude must be in some of the data and metadafiles of the import folders or will be ignored. 

merge import folders while excluding samples and forcing merging because of parameter differences

```
module load cbioportal-importer;
cbio_importer merge -configPaths path/to/configfiles -importPaths path/to/importfolders -cancerCode MIXED -genome hg38 -outdir path/to/the/output/directory -center -project ProjectName -description "short description of the project" -study "ProjectName: Top-Level-OncoTree, Concept (PI, OICR)"  --merge_with_differences -excludeSamples sample1 sample2 sample3
```


