#!/bin/sh


analyze -f /usr/local/share/freeling/config/es-old.cfg --tlevel 10 --nortk --nortkcon --input freeling --inplv splitted < $1 > $2

