#!/bin/bash

#Decision about your working terminal...

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=Mac;;
    CYGWIN*)    machine=Cygwin;;
    MINGW*)     machine=MinGw;;
    *)          machine="UNKNOWN:${unameOut}"
esac
echo ${machine}


#Write and decide the default python version in your computer

verdist="v03"

scrptname=`ls sxt_*ARF*.py`

scrptname2=`ls make_*EEF*.py`

mkdir -p $HOME/auxpyscrpt/

readmefile="readme.pdf"
readmefilehtml="readme.html"
docdir=$HOME/auxpyscrpt/docs/

#Write PDF help file...
if [ -f "$readmefile" ]; 
then
    echo $readmefile
    mkdir -p $docdir
    cp $readmefile $docdir/sxtEEFModule_help_${verdist}.pdf;
    echo "The help file sxtEEFModule_help_${verdist}.pdf(.html) is added to the directory $docdir... \n Please read it carefully before using this tool";
fi

#Write html help file...
if [ -f "$readmefilehtml" ]; 
then
    echo $readmefilehtml
    mkdir -p $docdir
    cp $readmefilehtml $docdir/sxtEEFModule_help_${verdist}.html;
    echo "The help file sxtEEFModule_help_${verdist}.html is added to the directory $docdir... \n Please read it carefully before using this tool";
fi


#adding/aliasing shortcut to paths....
 
cp ./$scrptname $HOME/auxpyscrpt/
cp ./$scrptname2 $HOME/auxpyscrpt/
    
if [ -f $HOME/.bashrc ]; then 

    sed -i '/EEF/d' $HOME/.bashrc
    sed -i '/sxtEEFmaker/d' $HOME/.bashrc
    sed -i '/sxtARFModule/d' $HOME/.bashrc

    echo "#For using SXT EEF/ARF Modules" >>$HOME/.bashrc
    echo 'alias sxtARFModule="python '$HOME/auxpyscrpt/$scrptname' "' >>$HOME/.bashrc
    echo 'alias sxtEEFmaker="python '$HOME/auxpyscrpt/$scrptname2' "' >>$HOME/.bashrc
    echo " The alias commands \'sxtARFModule\' and \'sxtEEFmaker\' are added to file : $HOME/.bashrc "

elif [ -f $HOME/.profile ]; then

    sed -i '/EEF/d' $HOME/.profile
    sed -i '/sxtEEFmaker/d' $HOME/.profile
    sed -i '/sxtARFModule/d' $HOME/.profile

    echo "#For using SXT EEF/ARF Modules" >>$HOME/.profile
    echo 'alias sxtARFModule="python '$HOME/auxpyscrpt/$scrptname' "' >>$HOME/.profile
    echo 'alias sxtEEFmaker="python '$HOME/auxpyscrpt/$scrptname2' "' >>$HOME/.profile
    echo " The alias commands \'sxtARFModule\' and \'sxtEEFmaker\' are added to file : $HOME/.profile "

else

    echo "#For using SXT EEF/ARF Modules" >>$HOME/.profile
    echo 'alias sxtARFModule="python '$HOME/auxpyscrpt/$scrptname' "' >>$HOME/.profile
    echo 'alias sxtEEFmaker="python '$HOME/auxpyscrpt/$scrptname2' "' >>$HOME/.profile
    echo " The alias commands \'sxtARFModule\' and \'sxtEEFmaker\' are added to file : $HOME/.profile "

fi

 



