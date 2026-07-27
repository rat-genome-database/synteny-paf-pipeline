#!/usr/bin/env bash
#
# as an cmdline parameter, pass the species type key (f.e. 3), or common species name (f.e. rat)
#  to run the pipeline for given species

. /etc/profile

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

APPNAME=synteny-paf-pipeline
APPDIR=/home/rgddata/pipelines/$APPNAME

cd $APPDIR

rm -rf "${APPDIR}/out/*"

java -Dspring.config=$APPDIR/../properties/default_db2.xml \
    -Dlog4j.configurationFile=file://$APPDIR/properties/log4j2.xml \
    -jar lib/$APPNAME.jar "${APPDIR}/out" 2>&1

# publish PAF files to the JBrowse2 server only on the production server (REED)
if [ "$SERVER" == "REED" ]; then
  scp "${APPDIR}/out/*" rgdpub@pipelines.rgd.mcw.edu:/data/data/jbrowse2/paf/
fi
