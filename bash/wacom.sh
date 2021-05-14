#!/bin/bash

STYLUS_ID=$(xsetwacom --list dev | grep stylus | cut -f 2 | sed -e 's/\(id:\)//g' -e 's/ //g')
PAD_ID=$(xsetwacom --list dev | grep pad | cut -f 2| sed -e 's/\(id:\)//g' -e 's/ //g')


xsetwacom set "$STYLUS_ID" Mode relative
xsetwacom set "$STYLUS_ID" Suppress 10
#xsetwacom set "$PAD_ID" Button 1 "key Ctrl Z"
#xsetwacom set "$PAD_ID" Button 2 "key Ctrl Y"
#xsetwacom set "$PAD_ID" Button 3 "key Page Up"
#xsetwacom set "$PAD_ID" Button 8 "key Ctrl Y"

