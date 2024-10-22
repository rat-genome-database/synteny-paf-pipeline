#!/usr/bin/env bash
#
# as an cmdline parameter, pass the species type key (f.e. 3), or common species name (f.e. rat)
#  to run the pipeline for given species

. /etc/profile

APPNAME=synteny-paf-pipeline
APPDIR=/home/rgddata/pipelines/$APPNAME

cd $APPDIR

rm -rf "${APPDIR}/out/*

java -Dspring.config=$APPDIR/../properties/default_db2.xml \
    -Dlog4j.configurationFile=file://$APPDIR/properties/log4j2.xml \
    -jar lib/$APPNAME.jar "${APPDIR}/out" 2>&1

