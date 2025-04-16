#!/usr/bin/env nextflow

// Script parameters
nextflow.enable.dsl=2

process bamCopy {

    input:
    val BAM_File

    output:
    file "*.bam"

    """
	ln -s $BAM_File .
    """
}

process createMatrix {

    input:
    file countSortFile
    val pref
    val TMPDIR
    val PUBDIR
    val TOMTX

    output:
    file '*mtx'
    file '*gns.txt'
    file '*barc.txt'

    """
	module load R/4.0.2-gccmkl

	~/workdir/PROGRAMS/Drop-seq_tools-2.5.1/./DigitalExpression I=$countSortFile O=$pref"_count.mat" NUM_CORE_BARCODES=80000 CELL_BARCODE_TAG='XC' MOLECULAR_BARCODE_TAG='XM'

	Rscript $TOMTX $pref"_count.mat" $pref"_count"

	cp $pref*".mtx" $PUBDIR
	cp $pref*".txt" $PUBDIR

    """
}

workflow {
	bamCopy(params.bamfile)
	createMatrix(bamCopy.out,params.pref,params.tmpdir,params.pubdir,params.tomtx)
}





