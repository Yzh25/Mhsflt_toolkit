# Mhsflt(MGIEasy high sensitive full-length transcriptome) analysis toolkit  
## Usage:  
###    Step1: Prerequisites  
####        Make sure you are using Python3 and have Pyyaml installed. If not, use the following command to install it:  
            pip install PyYAML  
###    Step2: Clone repository   
####        Use the following command to clone the repository:  
            git clone https://github.com/Yzh25/Mhsflt_toolkit.git  
###    Step3: Release environment   
####        Once you have cloned the repository, navigate to the root directory of the project and run following command to release environment:  
            sh release_env.sh  
###    Step4: Build reference index  
####        Download genome sequence fasta file (Primary assembly version recommend) and gene annotation gtf file （Basic gene annotation version recommend）from GENCODE (www.gencodegenes.org) and build genome index use STAR in following path:   
            /path/to/local/repository/zUMIs-env/bin/STAR  
###    Step5: Build config file  
####       Create a YAML format configuration file for your project, and the content of your configuration file can refer to 'example_config.yaml'.  
###    Step6: Run analysis  
####        Use following command to run analysis pipeline:  
            python3 /path/to/local/repository/run_analysis_pipeline.py -y your_config.yaml  
###    Step6: Check results  
####        Check analysis results in 'out_dir(defind in your_config.yaml)/results' directory.