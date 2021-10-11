#!/bin/bash
yade-2020.01a-batch $3 --log ./logs/$.%.log $1 $2 2>&1 | tee batch-log.log