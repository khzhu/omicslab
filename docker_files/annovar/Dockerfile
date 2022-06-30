FROM ubuntu:16.04

##############################################################
# Dockerfile Version:   1.0
# Software:             ANNOVAR
# Software Version:     current
# Software Website:     http://www.openbioinformatics.org/annovar/
# Description:          ANNOVAR
##############################################################

#Updating ubuntu and installing other necessary software
RUN apt-get clean all \
	&& apt-get -y update  \
	&& apt-get upgrade -y \
	&& apt-get install -y \
		unzip \
		wget \
		build-essential \
		zlib1g-dev \
		libncurses5-dev \
	&& apt-get clean


###########
# ANNOVAR #
###########

ENV TOOLS=/home/TOOLS/tools
ENV TOOL_NAME=annovar
ENV TOOL_VERSION=current
ENV TARBALL_LOCATION=http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/
ENV TARBALL=annovar.latest.tar.gz
ENV TARBALL_FOLDER=$TOOL_NAME
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz

# INSTALL
RUN wget $TARBALL_LOCATION/$TARBALL ; \
    tar xf $TARBALL --wildcards *pl ; \
    rm -rf $TARBALL ; \
    cd $TARBALL_FOLDER ; \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
    cp *pl $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/ -R ; \
    mkdir /databases ; \
    ln -s /databases/ $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/humandb ; \
	perl $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/humandb/ ;\
	perl $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar knownGene $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/humandb/ ;\
	perl $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp35c $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/humandb/ ;\
	perl $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar cosmic70 $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/humandb/ ;\
	perl $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar esp6500siv2_all $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/humandb/ ;\
	perl $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/humandb/ ;\
	perl $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03nontcga $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/humandb/ ;\
	perl $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp150 $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/humandb/ ;\
	perl $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20180603 $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/humandb/ ;\
    cd ../ ; \
    rm -rf $TARBALL_FOLDER ;

WORKDIR "$TOOLS/$TOOL_NAME/current/bin"

CMD ["/bin/bash"]
