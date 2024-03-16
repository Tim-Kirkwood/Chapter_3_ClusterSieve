# Chapter_3_ClusterSieve
Set up and start pip_env environment as described for Chapter_3_ModuleMapper.  Copy this directory to your local computer.  Choose a BGC query, run using CBLASTER (locally or with CAGECAT).  Do not find intermediate genes (makes no difference and will increase run time).  Under "Advanced", "Binary table", set "Delimeter" to "," (comma).  Can alter other settings acording to your use case.  Run ClusterSieve_thesis.py - before running, set neighborhood_size (line 108), similarity_filter (line 123) and folders (line 126) variables.  Other variables can be left as default (or changed if you want to get all cblaster files written in gbk format for example - see code comments in lines 102 to 126). 
