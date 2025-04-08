#!/bin/bash

nextflow run -c configs/resources.config start.nf -resume -with-trace -with-report nextflow_report.html
