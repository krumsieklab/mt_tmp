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

set -x

# enter the version tag for the docker image here
export VERSION="1.3.6"

# enter your GIT credentials here (needed to access to Jan's repo)
export GIT_CREDENTIALS='ksuhre:4meta-TOOLs'

# define the meta-dockerinitial target image
# https://hub.docker.com/r/rocker/tidyverse/tags
BASE="rocker/tidyverse:3.6.3" 

# deal with weird calling of docker.exe from Linux Subsystem for Windows
uname -a | grep Microsoft 
if [ $? -eq 0 ] ; then
  export DOCKER=docker.exe
else 
  export DOCKER=docker
fi

echo "KS: building meta-docker"

# check whether the directory already exists
# files may be overwritten
if [ -d meta-docker.dir ] ; then
  echo "KS: warning, directory meta-docker.dir alreay exists, some files may be overwritten"
else
  mkdir meta-docker.dir
fi

#################################################################################################
# create meta-tools install script 1
#################################################################################################
echo "KS: creating meta-tools install script 1"
cat > meta-docker.dir/KSinstall1.sh <<EOF2KS1
echo "running $0"
echo "This is meta-docker version $VERSION" > /meta-docker.version
date >> /meta-docker.version
date
uname -a
echo "installing Ubuntu libs"
apt-get -y update
apt-get install -y libjpeg-dev  # needed for R package "remote"

# install java, used by glmnet which requires rJava 
apt-get install -y default-jdk

# a dirty trick to get rJava installed (circumvent missing library error)
ln -sf  /lib/x86_64-linux-gnu/libbz2.so.1 /usr/local/lib/R/lib/libbz2.so
ln -sf  /lib/x86_64-linux-gnu/liblzma.so.5 /usr/local/lib/R/lib/liblzma.so

# libs needed for keras
apt-get install -y python-pip  
apt-get install -y python3-pip
apt-get install python3-venv
pip install virtualenv

echo "cloning meta-tools"
# git clone  https://'$GIT_CREDENTIALS'@gitlab.com/krumsieklab/qmdiab.git
# mv qmdiab /home/rstudio
git clone  https://'$GIT_CREDENTIALS'@gitlab.com/krumsieklab/mt.git
mv mt /home/rstudio
# git clone  https://'$GIT_CREDENTIALS'@gitlab.com/krumsieklab/snippets.git
# mv snippets /home/rstudio

# echo "cloning autonomics files"
# git clone https://github.com/bhagwataditya/autonomics
# mv autonomics /home/rstudio

# a quick-fix for upper-lower case mixup in metatools
# ln -s /home/rstudio/mt /home/rstudio/MT
# ln -s /home/rstudio/mt /home/rstudio/Mt

# pass ownership of /home/rstudio to the rstudio user so that a user can modify the libraries if required
chown -R rstudio /home/rstudio
chgrp -R users /home/rstudio
EOF2KS1

#################################################################################################
# create meta-tools install script 2
#################################################################################################
echo "KS: creating meta-tools install script 2"
cat > meta-docker.dir/KSinstall2.sh <<EOF2KS2
echo "running $0"
date
uname -a

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

# move the virtual environment (create by install_tensorflow) to /home/rstudio (not very elegant, better to run the install under the rstudio user, but I couldn't get this to run, su -l rstudio had no effect)
mv /root/.virtualenvs/ /home/rstudio/

# pass ownership of /home/rstudio to the rstudio user so that a user can modify the libraries if required
chown -R rstudio /home/rstudio
chgrp -R users /home/rstudio
EOF2KS2


#################################################################################################
# create meta-tools install script 3
#################################################################################################
echo "KS: creating meta-tools install script 3"
cat > meta-docker.dir/KSinstall3.sh <<EOF2KS3
echo "running $0"
date
uname -a

echo "installing required R-packages"
R <<EEOOFF
cat("installing all required packages\n")

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

# install Olink code
devtools::install_github(repo ='Olink-Proteomics/OlinkRPackage/OlinkAnalyze', upgrade = FALSE)

EEOOFF

# pass ownership of /home/rstudio to the rstudio user so that a user can modify the libraries if required
chown -R rstudio /home/rstudio
chgrp -R users /home/rstudio
EOF2KS3

#################################################################################################
# create Dockerfile
#################################################################################################
echo "KS: creating Dockerfile"
cat > meta-docker.dir/Dockerfile << EOF
FROM $BASE
COPY . /app
RUN sh /app/KSinstall1.sh
RUN sh /app/KSinstall2.sh
RUN sh /app/KSinstall3.sh
EOF

#################################################################################################
# build image (use --no-cache)
#################################################################################################
echo "KS: building image"
${DOCKER} build -t registry.gitlab.com/krumsieklab/mt/meta-docker:${VERSION}  --no-cache meta-docker.dir 

# list existing images
echo "KS: listing image"
${DOCKER} image ls

# list running containers
echo "KS: listing running container (including stopped ones)"
${DOCKER} container ls -a

echo "KS: run the following to push the image"
echo ${DOCKER} push registry.gitlab.com/krumsieklab/mt/meta-docker:${VERSION} 

echo "now run the following:"
echo "docker run -v\`pwd\`:/home/rstudio/home -e PASSWORD=pwd -p 8787:8787 --detach --name meta registry.gitlab.com/krumsieklab/mt/meta-docker:$VERSION"
