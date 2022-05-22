# sxtARFModule-v03
### ARF file for custom source region in SXT data analysis 
>Please uncompress the downloaded module. You may use the following depending on archive type.  
### tar xzvf SXTEEFModule_v03.tar.gz 
### OR 
### gunzip SXTEEFModule_v03.zip 
>In order to install this tool you need to run either of the following:
### sh install_sxtEEFModule.sh
### OR
### csh install_sxtEEFModule.csh
>This step will created a directory 'auxpyscrpt/' in $HOME directory and copy files with appropriate links... 
This way a keyword 'sxtARFModule' in the .bashrc or .cshrc or .profile file of your computer
You may delete the downloaded directory now
This code can be used from any terminal using:
### sxtARGModule -h [--help] to access the options
### sxtARFModule --evtfile=EVENTS_FILE_NAME --sxtpha=PHA_FILE_NAME --outstem=NAME_OUT_FILE_STEM --coordtyp=sky --vigcorrflag=yes --pltflag=yes
>if you do want to use some other verions of ARF file i.e., other than 'sxt_pc_excl00_v04_20190608.arf' please provide it in command line options 
### sxtARFModule --evtfile=EVENTS_FILE_NAME --sxtpha=PHA_FILE_NAME --outstem=NAME_OUT_FILE_STEM --coordtyp=sky --vigcorrflag=yes --pltflag=yes --sxtarf=FULL_PATH_OF_ARF_FILE
