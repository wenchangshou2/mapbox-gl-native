#!/bin/bash

tests=`fgrep mean $1 | cut -c1-85 | grep 91m`

if [ -z "$tests" ]
then
      echo "" >> $1
      echo "PASSED, all benchmarks within 5% threshold" >> $1

      ansi2html < $1 > $2
      exit 0
else
      test_names=`fgrep mean $1 | cut -c1-85 | grep 91m | cut -f1 -d' '`

      echo "" >> $1
      echo "FAILED, following benchmarks are not within 5% threshold:" >> $1
      echo "" >> $1
      for t in $test_names;
      do
        echo $t >> $1
      done

      ansi2html < $1 > $2
      exit 1
fi
