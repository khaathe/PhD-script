#!/usr/bin/bash

#Package to install
PACKAGE_MANAGER="apt"
PACKAGE="libcurl4-gnutls-dev libssl-dev default-jre libxml2-dev curl"
TOOLS="rclone rstudio r-base"
SNAP_PACKAGE="code"

#Install Command
INSTALL_COMMAND="$PACKAGE_MANAGER install $PACKAGE $TOOLS"
SNAP_COMMAND="snap install --classic $SNAP_PACKAGE"

#Update system Command
UDPATE_COMMAND="$PACKAGE_MANAGER upgrade"

#Python version
PYTHON3="3.8"

#Update System
$UDPATE_COMMAND

#Install packages, tools and snap
$INSTALL_COMMAND
$SNAP_COMMAND

alias python="/usr/bin/python$PYTHON3"

#echo '#Adding directories into PATH' >> ~/.profile
#echo 'PATH=$PATH:/home/spinicck/Bin/:/home/spinicck/Bin/sratoolkit.2.10.8-ubuntu64/bin' >> ~/.profile
