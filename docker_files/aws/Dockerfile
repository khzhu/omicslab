# Using OpenJDK 8
# This Dockerfile does not require any files that are in the GATK4 repo.
FROM ubuntu:18.04

#### Basic image utilities
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
    python \
    wget \
    curl \
    bc \
    unzip \
    bzip2 \
    less \
    bedtools \
    samtools \
    openjdk-8-jdk \
    tabix \
    gpg-agent \
    build-essential \
    software-properties-common && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove


ENV JAVA_LIBRARY_PATH /usr/lib/jni
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/

RUN java -version

RUN ls .
ADD . /gatk
WORKDIR /gatk

COPY ./nio-spi-for-s3-1.1.0-all.jar $JAVA_HOME/jre/lib/ext/

RUN wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip -O /tmp/gatk-4.2.6.1.zip\
	&& unzip /tmp/gatk-4.2.6.1.zip -d / \
	&& rm /tmp/gatk-4.2.6.1.zip -f \
	&& cd / \
	&& mv gatk-4.2.6.1 gatk \
	&& /gatk/gatk --list
	
# Deleting unneeded caches
RUN rm -rf /var/lib/apt/lists/*

#Setup linked jars that may be needed for running gatk
RUN ln -s $( find /gatk -name "gatk*local.jar" ) gatk.jar
RUN ln -s $( find /gatk -name "gatk*local.jar" ) /root/gatk.jar
RUN ln -s $( find /gatk -name "gatk*spark.jar" ) gatk-spark.jar

WORKDIR /root

# Make sure we can see a help message
RUN java -jar gatk.jar -h

# Set environment variables.
ENV HOME /gatk
ENV CLASSPATH /gatk/gatk.jar:$CLASSPATH

# Define working directory.
WORKDIR /gatk
RUN chmod -R a+rw /gatk

# Define default command.
CMD ["bash"]
