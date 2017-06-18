#!/bin/sh
ls -1 Cnf* > cnf.list
ls -1 Qnc* > qnc.list
mancanti cnf.list qnc.list > todo.list
echo `pwd`/ > lista ; cat todo.list >> lista
