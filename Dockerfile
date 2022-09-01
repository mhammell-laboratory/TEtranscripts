FROM python

#PREAMBLE

WORKDIR /home/genomics
COPY . /home/genomics
RUN cd /home/genomics

RUN apt-get --assume-yes update \
	&& apt-get --assume-yes upgrade

#MAIN

RUN apt-get --assume-yes install r-base

RUN R -e "install.packages('locfit', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
  && R -e "install.packages('BiocManager', dependencies=TRUE, repos='http://cran.rstudio.com/'); BiocManager::install()" \
	&& R -e "BiocManager::install(\"DESeq2\")" \
	&& pip install pysam \
	&& pip install TEtranscripts \
	&& rm -rf *.tgz *.tar *.zip \
	&& rm -rf /var/cache/apk/* \
	&& rm -rf /tmp/*
