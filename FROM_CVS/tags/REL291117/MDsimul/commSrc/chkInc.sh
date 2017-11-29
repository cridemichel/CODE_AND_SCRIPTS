#! /bin/bash
if [ ! -e measures.h ] 
then
    echo "ERROR: I need measures.h file!"
    echo "Create it using the template in commSrc/templates dir"
fi
if [ ! -e decls.h ] 
then
    echo "ERROR: I need the decls.h file"
    echo "Create it using the template in commSrc/template.h dir"
fi
if [ ! -e defs.h ]
then
    echo "ERROR: I need the defs.h file"
    echo "Create it using the template in commSrc/template.h"
fi
