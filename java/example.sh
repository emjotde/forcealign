#!/bin/bash
SD=`dirname $0`
java -Djava.library.path=$SD/build -cp $SD/build example $SD/../data fr en
