#!/usr/bin/env bash
#
# generates PAF files for assembly pairs having synteny data into $APPDIR/out
#  and, on the production server, publishes them to the JBrowse2 server

. /etc/profile

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

APPNAME=synteny-paf-pipeline
APPDIR=/home/rgddata/pipelines/$APPNAME

cd $APPDIR || exit 1

# note: the glob must stay outside the quotes, otherwise it is not expanded
rm -rf "${APPDIR:?}/out/"*

java -Dspring.config=$APPDIR/../properties/default_db2.xml \
    -Dlog4j.configurationFile=file://$APPDIR/properties/log4j2.xml \
    -jar lib/$APPNAME.jar "${APPDIR}/out" 2>&1
STATUS=$?
if [ $STATUS -ne 0 ]; then
  echo "$APPNAME FAILED (exit code $STATUS); PAF files not published"
  exit $STATUS
fi

# publish PAF files to the JBrowse2 server only on the production server (REED)
if [ "$SERVER" == "REED" ]; then
  scp "${APPDIR}/out/"* rgdpub@pipelines.rgd.mcw.edu:/data/data/jbrowse2/paf/
fi
