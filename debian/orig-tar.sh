#!/bin/sh -e

# called by uscan with '--upstream-version' <version> <file>
SOURCE_NAME=giac
VERSION=$2
DEBIAN_VERSION=$VERSION-dfsg1
UPSTREAM_SOURCE_DIR=${SOURCE_NAME}-$VERSION
DEBIAN_SOURCE_DIR=${SOURCE_NAME}-$DEBIAN_VERSION
TAR=../${SOURCE_NAME}_$DEBIAN_VERSION.orig.tar.gz

# clean up the upstream tarball
tar xzf $3 
# rename upstream source dir
mv ${UPSTREAM_SOURCE_DIR} ${DEBIAN_SOURCE_DIR}
# Remove doc/fr documentation (C) Renée de Graeve
tar -c -z -X debian/orig-tar.exclude -f $TAR ${DEBIAN_SOURCE_DIR}/
rm -rf ${DEBIAN_SOURCE_DIR} $3
echo "$SOURCE_NAME: Applied DFSG removals and renamed tarball to `basename ${TAR}`"
