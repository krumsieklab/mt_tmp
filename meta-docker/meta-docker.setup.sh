#FILE:    meta-docker.setup.sh
#AUTHOR:  Karsten Suhre
#DATE:    Mon Jan 27 10:48:47 +03 2020
#PURPOSE: set-up a docker image for meta tools
#BUGS:    
#MODIF:   Mon Jan 27 16:45:08 +03 2020
#         make it run on Linux subsystem for Windows
#         Thu Jan 30 11:42:27 +03 2020
#         add R packages for next version
#         Jan's to-do list:
#         firefox https://docs.google.com/spreadsheets/d/1P3hEAb_vrABsucvwD7RT0jTztRHw7WYE_TnhDK7F2ng/edit?usp=sharing
#
#         Version 1.2
#         Thu Jan 30 16:45:56 +03 2020
#         make keras work
#
#         Version 1.3::
#         Sun Feb  2 12:47:49 +03 2020
#         add cloning autonomics into /home/rstudio
#
#         Version 1.3.1:
#         Sun Feb 23 14:10:49 +03 2020
#         add missing tensorflow install command
#
#         Version 1.3.2:
#         add more R libs
#         debug install of keras/tensorflow
#
#         Version 1.3.3:
#         installed java developper kit
#         used a dirty trick (ln -s) to arrange some library problems when installing rJava (needed for glmnet)
#         add additional R packages: ggforce car rJava glmnet stabs mboost
#         
#         Version 1.3.4:
#         fix change group for rstudio to users
#         add cloning of snippets.git from Jan's lab
#         add "ln -s mt Mt" fix
#
#         Version 1.3.5
#         install.packages("glmulti")
#         remove git clone of autonomics
#         add Olink R package and dependencies
#         install.packages("ggfortify")
#         install.packages("emmeans")
#         devtools::install_github(repo ='Olink-Proteomics/OlinkRPackage/OlinkAnalyze')
#         add more libs
#         install.packages("kableExtra")
#         BiocManager::install("ppcor", update = FALSE)
#         install.packages("ggpubr")
#
#         Version 1.3.6
#         move to rocker/tidyverse:3.6.3
#         commented out: quick fix for Mt upper-lower case problem (let's see whether the problem disappeared)
#         commented out: pull of qmdiab.git and snippet.git (let's look at that later)
#         add apt-get install python3-venv 
#	  a few minor fixes
#         
#         error installing tensorflow remains:
#          
#          > install_tensorflow()
#          Creating virtual environment '~/.virtualenvs/r-reticulate' ...
#          Using python: /usr/bin/python3.7
#          The virtual environment was not created successfully because ensurepip is not
#          available.  On Debian/Ubuntu systems, you need to install the python3-venv
#          package using the following command.
#          
#              apt-get install python3-venv
#          
#          You may need to use sudo with that command.  After installing the python3-venv
#          package, recreate your virtual environment.
#          
#          Failing command: ['/root/.virtualenvs/r-reticulate/bin/python3.7', '-Im', 'ensurepip', '--upgrade', '--default-pip']
#          
#          ^[[91mError: Error creating virtual environment '~/.virtualenvs/r-reticulate' [error code 1]
#          Execution halted
#          
#         https://stackoverflow.com/questions/39539110/pyvenv-not-working-because-ensurepip-is-not-available
#         ==> that fix does not work
#         
#         Version 1.4
# 
#         10 July - total structural overhaul
#                   metatools as package
#                   logical packing of the different install scripts
#                   inactivated keras/tensorflow for now
#                   inactivated autonomics for now
#
#          13 July - further testing 
#
# personal notes KS regarding Windows:
#            We recommend to convert this distro to WSL 2 and activate
#            the WSL integration in Docker Desktop settings.
#
# install Olink code 
# devtools::install_github(repo ='Olink-Proteomics/OlinkRPackage/OlinkAnalyze', upgrade = FALSE)


set -

# enter the version tag for the docker image here
export VERSION="1.4.0"

# enter your GIT credentials here (needed to access to Jan's repo)
# in the future we should change this to using git tokens
export GIT_CREDENTIALS='ksuhre:4meta-TOOLs'

# define the initial target image from rocker
# see https://hub.docker.com/r/rocker/tidyverse/tags
BASE="rocker/tidyverse:3.6.3" 

# deal with weird calling of docker.exe from Linux Subsystem for Windows
# for some reason it is not possible to create an alias (because the docker command already exists under Windows, I guess)
uname -a | grep Microsoft 
if [ $? -eq 0 ] ; then
  export DOCKER=docker.exe
else 
  export DOCKER=docker
fi

# create files that will be moved to the docker image for install in meta-docker.dir
if [ -d meta-docker.dir ] ; then
  echo "KS: warning, directory meta-docker.dir alreay exists, some files may be overwritten"
else
  mkdir meta-docker.dir
fi


#################################################################################################
# create meta-tools script KSinstall_libs.sh
#################################################################################################
echo "KS: creating meta-tools install script KSinstall_libs.sh"
cat > meta-docker.dir/KSinstall_libs.sh <<EOF2KS1
echo "running \$0"
echo "This is \$0 version $VERSION" >> /meta-docker.version
date >> /meta-docker.version
date
uname -a

echo "installing Ubuntu libs"
apt-get install -y libjpeg-dev  # needed for R package "remote"
apt-get install vim

# install java, used by glmnet which requires rJava 
apt-get install -y default-jdk

# a dirty trick to get rJava installed (circumvent missing library error)
ln -sf  /lib/x86_64-linux-gnu/libbz2.so.1 /usr/local/lib/R/lib/libbz2.so
ln -sf  /lib/x86_64-linux-gnu/liblzma.so.5 /usr/local/lib/R/lib/liblzma.so

EOF2KS1

#################################################################################################
# create meta-tools install script for  keras/tensorflow
#################################################################################################
echo "KS: creating meta-tools install script KSinstall_keras.sh"
cat > meta-docker.dir/KSinstall_keras.sh <<EOF2KS2
echo "running \$0"
echo "This is \$0 version $VERSION" >> /meta-docker.version
date >> /meta-docker.version
date
uname -a

# libs needed for keras
apt-get install -y python-pip  
apt-get install -y python3-pip
apt-get install python3-venv
pip install virtualenv

echo "installing required R-packages"
#------------------------------
R <<EEOOFF
cat("installing keras/tensorflow\n")

# install keras/tensorflow
install.packages("tensorflow")
library(tensorflow)
install_tensorflow()
install.packages("keras")
library(keras)
install_keras()
EEOOFF

# move the virtual environment (created by install_tensorflow) to /home/rstudio (not very elegant, better to run the install under the rstudio user, but I couldn't get this to run, su -l rstudio had no effect)
mv /root/.virtualenvs/ /home/rstudio/

chown rstudio /home/rstudio/.virtualenvs
chgrp users /home/rstudio/.virtualenvs

EOF2KS2


#################################################################################################
# create meta-tools install script for diverse R packages
#################################################################################################
echo "KS: creating meta-tools install script KSinstall_packages.sh"
cat > meta-docker.dir/KSinstall_packages.sh <<EOF2KS3
echo "running \$0"
echo "This is \$0 version $VERSION" >> /meta-docker.version
date >> /meta-docker.version
date
uname -a

R <<EEOOFF
cat("installing R packages\n")

# R packages
install.packages("hash")
install.packages("igraph")
install.packages("d3heatmap")
install.packages("DT")
install.packages("plotly")
install.packages("pROC")
install.packages("caret")
install.packages("car")
install.packages("voronoiTreemap")
install.packages("broom.mixed")
install.packages("lmerTest")
install.packages("ggforce")
install.packages("car")
install.packages("rJava")
install.packages("glmnet")
install.packages("stabs")
install.packages("mboost")
install.packages("glmulti")
install.packages("kableExtra")
install.packages("ggpubr")

# Bioconductor packages, most required by metatools
install.packages('BiocManager')
BiocManager::install("graphite", update = FALSE)
BiocManager::install('SummarizedExperiment', update = FALSE)
BiocManager::install('gage', update = FALSE)
BiocManager::install('tictoc', update = FALSE)
BiocManager::install('logging', update = FALSE)
BiocManager::install('gdata', update = FALSE)
BiocManager::install('ggrepel', update = FALSE)
BiocManager::install('GGally', update = FALSE)
BiocManager::install('sna', update = FALSE)
BiocManager::install('ggnetwork', update = FALSE)
BiocManager::install('visNetwork', update = FALSE)
BiocManager::install('dils', update = FALSE)
BiocManager::install('ggbeeswarm', update = FALSE)
BiocManager::install('sva', update = FALSE)
BiocManager::install('limma', update = FALSE)
BiocManager::install('GeneNet', update = FALSE)
BiocManager::install('formula.tools', update = FALSE)
BiocManager::install("rpubchem", update = FALSE)
BiocManager::install('pathview', update = FALSE)
BiocManager::install("ppcor", update = FALSE)

# autonomics
# https://github.com/bhagwataditya/autonomics
BiocManager::install('mixOmics', update = FALSE)   # CRAN -> BioC, requires explicit installation

# Install autonomics (drop ref = 'dev' to install older autonomics stable)
install.packages('remotes')
remotes::install_github('bhagwataditya/autonomics/autonomics.data',       ref = 'dev', upgrade = FALSE)
remotes::install_github('bhagwataditya/autonomics/autonomics.support',    ref = 'dev', upgrade = FALSE)
remotes::install_github('bhagwataditya/autonomics/autonomics.annotate',   ref = 'dev', upgrade = FALSE)
remotes::install_github('bhagwataditya/autonomics/autonomics.import',     ref = 'dev', upgrade = FALSE)
remotes::install_github('bhagwataditya/autonomics/autonomics.preprocess', ref = 'dev', upgrade = FALSE)
remotes::install_github('bhagwataditya/autonomics/autonomics.plot',       ref = 'dev', upgrade = FALSE)
remotes::install_github('bhagwataditya/autonomics/autonomics.find',       ref = 'dev', upgrade = FALSE)
remotes::install_github('bhagwataditya/autonomics/autonomics.ora',        ref = 'dev', upgrade = FALSE)
remotes::install_github('bhagwataditya/autonomics/autonomics',            ref = 'dev', upgrade = FALSE)

EEOOFF

# pass ownership of /home/rstudio to the rstudio user so that a user can modify the libraries if required
chown -R rstudio /home/rstudio/mt
chgrp -R users /home/rstudio/mt
EOF2KS3

#################################################################################################

#################################################################################################
# create meta-tools install script KSinstall_mt.sh
#################################################################################################
echo "KS: creating meta-tools install script KSinstall_mt.sh"
cat > meta-docker.dir/KSinstall_mt.sh <<EOF2KS4
echo "running \$0"
echo "This is \$0 version $VERSION" >> /meta-docker.version
date >> /meta-docker.version
date
uname -a

# change into the home dir of the rstudio user (although the script will run as root)
cd /home/rstudio

# we should replace the below using tokens 
echo "cloning meta-tools"
git clone  https://'$GIT_CREDENTIALS'@gitlab.com/krumsieklab/mt.git
git clone  https://'$GIT_CREDENTIALS'@gitlab.com/krumsieklab/snippets.git

# a quick fix to work around different versions of mt/Mt/MT (may only concern linux systems)
ln -s /home/rstudio/mt /home/rstudio/MT
ln -s /home/rstudio/mt /home/rstudio/Mt

R <<EEOOFF1
cat("installing meta-tools\n")

# the below command would install from the GitLab - but we want an install it from a local directory for now
# devtools::install_gitlab(repo="krumsieklab/mt", subdir = "MetaboTools", auth_token = "Eu8aCTCZn7grLyrT3xUm"

# the below should not be needed, but I keep it as a reminder for now
# codes.makepath <- function(p){file.path("/home/rstudio",p)}
# source(codes.makepath("snippets/mt.checkout.R"))
# mt.quickload <- function()source(codes.makepath("mt/quickload.R"))
# mt.quickload()

devtools::install("mt/MetaboTools", keep_source = T)

cat("installing meta-tools done\n")

EEOOFF1

# pass ownership of /home/rstudio to the rstudio user so that a user can modify the libraries if required
chown -R rstudio /home/rstudio/mt
chgrp -R users /home/rstudio/mt
chown -R rstudio /home/rstudio/snippets
chgrp -R users /home/rstudio/snippets
chown rstudio /home/rstudio/Mt
chgrp users /home/rstudio/Mt
chown rstudio /home/rstudio/MT
chgrp users /home/rstudio/MT

EOF2KS4

#################################################################################################
# create Dockerfile
#################################################################################################
echo "KS: creating Dockerfile"
cat > meta-docker.dir/Dockerfile << EOF
FROM $BASE
COPY . /app
#RUN sh /app/KSinstall_libs.sh
#RUN sh /app/KSinstall_keras.sh
#RUN sh /app/KSinstall_packages.sh
RUN sh /app/KSinstall_mt.sh
EOF

#################################################################################################
# build image (use --no-cache)
#################################################################################################
echo "KS: to build the image"
echo ${DOCKER} build -t registry.gitlab.com/krumsieklab/mt/meta-docker:${VERSION}  --no-cache meta-docker.dir 

# list existing images
echo "KS: to list existing images:"
echo ${DOCKER} image ls

echo "KS: to push the image to the GitLab registry:"
echo ${DOCKER} push registry.gitlab.com/krumsieklab/mt/meta-docker:${VERSION} 

echo "to start the image (adapt the work directory you wish to mount using -v):"
echo "${DOCKER} run -v/tmp:/home/rstudio/home -e PASSWORD=pwd -p 8787:8787 --detach --name meta registry.gitlab.com/krumsieklab/mt/meta-docker:$VERSION"

echo "KS: to open a shell inside the container:"
echo "${DOCKER} exec -it meta /bin/bash"

echo "KS: to stop the container:"
echo "${DOCKER} stop meta"

echo "KS: to start the container again:"
echo "${DOCKER} start meta"

echo "KS: to delete the container:"
echo "${DOCKER} rm meta"

echo "KS: to list container and images"
echo "${DOCKER} container ls"
echo "${DOCKER} image ls"

