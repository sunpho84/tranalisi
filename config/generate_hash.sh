#!/bin/bash

rem_newl ()
{
    tr -d "\n"
}

requote ()
{
    for i in $*
    do
	printf "%q " $i
    done
}

echo '#ifndef GIT_HASH'

echo -n ' #define GIT_HASH "'
requote $(git rev-parse HEAD|rem_newl)
echo '"'

echo -n ' #define GIT_TIME "'
requote $(git log -1 --pretty=%ad|rem_newl)
echo '"'

echo -n ' #define GIT_LOG "'
requote $(git log -1 --pretty=%B|rem_newl)
echo '"'

echo '#endif'
