# meta-docker

**Integration of meta-tools and autonomics in an rstudio docker image from the rocker project**

FILE:    meta-docker.README.md<br>
DATE:    29 January 2020<br>
AUTHOR:  Karsten Suhre<br>
PURPOSE: Readme file for the meta-rocker project<br>
MODIF:   02 July 2020 - port from WCM-Q code.qatar-med.cornell.edu/kas2049/meta-docker to krumsieklab GitLab<br>
<p>


**HOWTO use the meta-docker image**<p>

**establish connection with the GitLab docker registry using your GitLab credentials (needs to be run only once):**<br>
docker login registry.gitlab.com

**Under Windows, open a CMD terminal and run the following command:**<br>
docker run -v%USERPROFILE%:/home/rstudio/home -e PASSWORD=pwd -p 8787:8787 --detach --name meta code.qatar-med.cornell.edu/kas2049/meta-docker/meta-docker:1.3.4

**Under Linux, run the following command in a shell (adapt the directory to mount using the -v option):**<br>
docker run -v$HOME:/home/rstudio/home -e PASSWORD=pwd -p 8787:8787 --detach --name meta code.qatar-med.cornell.edu/kas2049/meta-docker/meta-docker:1.3.5

**Then open a browser and navigate to localhost:8787**

    username: rstudio
    pwd: pwd

<br><br><br>

The docker -v option can be used to mount any directory that is needed. In the example above, the user home directory is mounted.


**HOWTO create a new meta-docker image:**<br>
./meta-docker.setup.sh<br>

Please note: this step needs only be performed when a new version of the docker image is created.
The "normal" user shall not have to do this.
The shell script meta-docker.setup.sh will generate all auxillary files needed and build the final docker image.
In particular, it will install autonomics ([maintained by Aditya Bhagwat](https://github.com/bhagwataditya/autonomics)) 
and metatools ([Jan Krumsiek lab](https://gitlab.com/krumsieklab/mt)) into a rocker/tidyverse environment.
Keras, Tensorflow and rJava is also included.<br>

<br><br><br>

**REFERENCES & NOTES**

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

**meta-docker at WCM-Q**<br>
https://code.qatar-med.cornell.edu/kas2049/meta-docker


**Cheat-sheet for docker**<br>
docker image ls<br>
docker container ls -a<br>
docker container stop meta<br>
docker container start meta<br>
docker container rm meta<br>
docker exec meta df<br>
docker exec -it meta /bin/bash<br>

