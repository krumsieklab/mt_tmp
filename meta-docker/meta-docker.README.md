# meta-docker

**Integration of meta-tools and autonomics in an rstudio docker image from the rocker project**

FILE:    README.md<br>
DATE:    29 January 2020<br>
AUTHOR:  Karsten Suhre<br>
PURPOSE: Readme file for the meta-rocker project<br>
<p>

This project will install autonomics ([maintained by Aditya Bhagwat](https://github.com/bhagwataditya/autonomics)) and metatools ([Jan Krumsiek lab](https://gitlab.com/krumsieklab/mt)) into a rocker/tidyverse environment and make the docker image available via GitLab.

The user can mount their own working directory using the docker -v option and then work with metatools and autonomics in an rstudio window via any web browser.

Planned extensions are routines to ease the cross-usage of metatools and autonomics.

**HOWTO use the meta-docker image**<p>

**establish your credentials (needed only once):**<br>
docker login code.qatar-med.cornell.edu

**Under Windows, open a CMD terminal and run the following command (adapt the directory to mount using the -v option)**<br>
docker run -v%USERPROFILE%:/home/rstudio/home -e PASSWORD=pwd -p 8787:8787 --detach --name meta code.qatar-med.cornell.edu/kas2049/meta-docker/meta-docker:1.3.4

**then open a browser and navigate to localhost:8787**

    username: rstudio
    pwd: pwd

**Notes for developpers**<br>
The following is needed to create a docker image from scratch

**test connection to git:**<br>
ssh -T -p 10022  git@code.qatar-med.cornell.edu

**clone setup script from WCM-Q Gitlab:**<br>
git clone https://code.qatar-med.cornell.edu/kas2049/meta-docker

**run setup script to create a new image:**<br>
./meta-docker.setup.sh<br>
This will generate all auxillary files needed and build the final docker image. In particular, it will clone metatools and autonomics into the docker image, and install all R packages that are needed by these.


**REFERENCES**

**Rocker project:**<br>
https://www.rocker-project.org/images/<br>
https://hub.docker.com/u/rocker<br>
https://hub.docker.com/r/rocker/tidyverse/tags

**Docker documentations:**<br>
https://docs.docker.com/

**autonomics**<br>
https://github.com/bhagwataditya/autonomics

**metatools**<br>
https://gitlab.com/krumsieklab/mt<br>
(requires access to Jan's gitlab, note that my personal credentials are presently in the setup script, so please do not distribute)


**Notes**<br>
docker image ls<br>
docker container ls -a<br>
docker container stop meta<br>
docker container start meta<br>
docker container rm meta<br>
docker exec meta df<br>
docker exec -it meta /bin/bash<br>

