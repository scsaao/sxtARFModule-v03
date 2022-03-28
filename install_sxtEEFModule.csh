#! /bin/csh -f


#Write and decide the default python version in your computer
set verdist="v03"
set scrptname=`ls sxt_*ARF*.py`

set scrptname2=`ls make_*EEF*.py`

mkdir -p $HOME/auxpyscrpt/

set readmefile="readme.pdf"
set readmefilehtml="readme.html"
set docdir="$HOME/auxpyscrpt/docs/"

#Adding pdf readme file...
if (-f "$readmefile") then

    echo $readmefile;
    mkdir -p $docdir
    cp $readmefile $docdir/sxtEEFModule_help_${verdist}.pdf;
    echo "The help file sxtEEFModule_help_${verdist}.pdf is added to the directory $docdir... \n Please read it carefully before using this tool";

endif

#Adding html readme file...
if (-f "$readmefilehtml") then

    echo $readmefilehtml;
    mkdir -p ${docdir}
    cp $readmefilehtml ${docdir}/sxtEEFModule_help_${verdist}.html;
    echo "The help file sxtEEFModule_help_${verdist}.html is added to the directory ${docdir}... \n Please read it carefully before using this tool";

endif

#adding/aliasing shortcut to paths....
 
cp ./$scrptname $HOME/auxpyscrpt/
cp ./$scrptname2 $HOME/auxpyscrpt/
    
 
if  (-f $HOME/.cshrc ) then 

    sed -i '/EEF/d' $HOME/.cshrc
    sed -i '/sxteefmaker/d' $HOME/.cshrc
    sed -i '/sxtARFModule/d' $HOME/.cshrc

    echo "#For using SXT EEF/ARF Modules" >>$HOME/.cshrc
    echo 'alias sxtARFModule "python '$HOME/auxpyscrpt/$scrptname' "' >>$HOME/.cshrc
    echo 'alias sxtEEFmaker "python '$HOME/auxpyscrpt/$scrptname2' "' >>$HOME/.cshrc
    echo " The alias commands \'sxtARFModule\' and \'sxtEEFmaker\' are added to file : $HOME/.cshrc "

else if (-f $HOME/.cshrc ) then

    sed -i '/EEF/d' $HOME/.profile
    sed -i '/sxteefmaker/d' $HOME/.profile
    sed -i '/sxtARFModule/d' $HOME/.profile

    echo "#For using SXT EEF/ARF Modules" >>$HOME/.profile
    echo 'alias sxtARFModule "python '$HOME/auxpyscrpt/$scrptname' "' >>$HOME/.profile
    echo 'alias sxtEEFmaker "python '$HOME/auxpyscrpt/$scrptname2' "' >>$HOME/.profile
    echo " The alias commands \'sxtARFModule\' and \'sxtEEFmaker\' are added to file : $HOME/.profile "

else 

    echo "#For using SXT EEF/ARF Modules" >>$HOME/.profile
    echo 'alias sxtARFModule "python '$HOME/auxpyscrpt/$scrptname' "' >>$HOME/.profile
    echo 'alias sxtEEFmaker "python '$HOME/auxpyscrpt/$scrptname2' "' >>$HOME/.profile
    echo " The alias commands \'sxtARFModule\' and \'sxtEEFmaker\' are added to file : $HOME/.profile "

endif
       

