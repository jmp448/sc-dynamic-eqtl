#!/bin/bash

all_tests="file1.txt"
dqtl_tests="file2.txt"

comm -12 $all_tests $dqtl_tests
