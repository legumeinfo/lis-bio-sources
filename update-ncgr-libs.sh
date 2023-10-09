#!/bin/bash
for dir in `ls -1d lis-*`
do
    cp /home/shokin/java/ncgr/datastore/build/libs/ncgr-datastore.jar $dir/libs
done
