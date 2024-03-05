$HOSTNAME = ""
params.outdir = 'results'  

//* autofill
if ($HOSTNAME == "default"){
    $DOCKER_IMAGE = "immcantation/suite:4.3.0"
    $DOCKER_OPTIONS = "-v /work:/work"

}

//* platform
//* platform


//* autofill

// Process Parameters:

// Process Parameters for Assemble_pairs_assemble_pairs:
//* params.Assemble_pairs_assemble_pairs.method =  "align"  //* @dropdown @options:"align, sequential, reference, join" @description:"Assembly method. Default align (See https://presto.readthedocs.io/ for more details.)"
//* params.Assemble_pairs_assemble_pairs.coord =  "sra"  //* @dropdown @options:"illumina, solexa, sra, 454, presto" @description:"The format of the sequence identifier which defines shared coordinate information across mate pairs. Default presto" @title:"General params"
//* params.Assemble_pairs_assemble_pairs.rc =  "tail"  //* @dropdown @options:"tail, head, both, none" @description: "Specify which read to reverse complement before stitching. Default tail"
//* params.Assemble_pairs_assemble_pairs.head_fields_R1 =  ""  //* @input @description:"Annotation fields to copy from R1 file head records into assembled record. Input should be seperated by space. Default CONSCOUNT"
//* params.Assemble_pairs_assemble_pairs.head_fields_R2 =  ""  //* @input @description:"Annotation fields to copy from R2 file head records into assembled record. Input should be seperated by space. Default CONSCOUNT PRCONS"
//* params.Assemble_pairs_assemble_pairs.failed =  "false"  //* @checkbox @description:"Check the box to activate failed option. Default false" @tooltip:"Specify to output the failed sequences as well."
//* params.Assemble_pairs_assemble_pairs.fasta =  "false"  //* @checkbox @description:"Check the box to get fasta file as output. Default false"
//* params.Assemble_pairs_assemble_pairs.nproc =  "40"  //* @input @description:"Number of nproc to use for running FilterSeq. Default value 1."
//* params.Assemble_pairs_assemble_pairs.alpha =  0.00001  //* @input @description:"Significance threshold for de novo paired-end assembly. Default 1e-05" @title:"De novo assembly arguments"
//* params.Assemble_pairs_assemble_pairs.maxerror =  0.2  //* @input @description:"Maximum allowable error rate for de novo assembly. Default 0.3"
//* params.Assemble_pairs_assemble_pairs.minlen =  20  //* @input @description:"Minimum sequence length to scan for overlap in de novo assembly. Default 8"
//* params.Assemble_pairs_assemble_pairs.maxlen =  1000  //* @input @description:"Maximum sequence length to scan for overlap in de novo assembly. Default 1000"
//* params.Assemble_pairs_assemble_pairs.scanrev =  "false"  //* @checkbox @description:"If specified, scan past the end of the tail sequence in de novo assembly to allow the head sequence to overhang the end of the tail sequence. Default false"
//* params.Assemble_pairs_assemble_pairs.minident =  0.5  //* @input @description:"Minimum identity of the assembled sequence required to call a valid reference guided assembly (between 0 and 1). Default 0.5"
//* params.Assemble_pairs_assemble_pairs.evalue =  0.00001  //* @input @description:"Minimum E-value for reference alignment for both the head and tail sequence. Default 1e-05"
//* params.Assemble_pairs_assemble_pairs.maxhits =  100  //* @input @description:"Maximum number of hits from the reference alignment to check for matching head and tail sequence assignments. Default 100"
//* params.Assemble_pairs_assemble_pairs.fill =  "false"  //* @checkbox @description:"Check the box to change the behavior of inserted characters when the head and tail sequences do not overlap during reference guided assembly. Default: False" @tolltip:"If specified, this will result in inserted of the V region reference sequence instead of a sequence of Ns in the non-overlapping region. Warning: you could end up making chimeric sequences by using this option."
//* params.Assemble_pairs_assemble_pairs.aligner =  "blastn"  //* @dropdown @options:"blastn, usearch" @description:"The local alignment tool to use. Must be one blastn (blast+ nucleotide) or usearch (ublast algorithm). Default blastn"
//* params.Assemble_pairs_assemble_pairs.// align_exec =   ""    //* @input @description:"The name or location of the aligner executable file (blastn or usearch). Defaults to the name specified by the aligner argument. Default: None"
//* params.Assemble_pairs_assemble_pairs.// dbexec =   ""    //* @input @description:"The name or location of the executable file that builds the reference database. This defaults to makeblastdb when blastn is specified to the aligner argument, and usearch when usearch is specified. Default None"
//* params.Assemble_pairs_assemble_pairs.gap =  0  //* @input @description:"Number of N characters to place between ends. Default 0" @title:"join assembly arguments"
//* params.Assemble_pairs_assemble_pairs.usearch_version =  "11.0.667"  //* @input @description:"The usearch version to download and run. Default 11.0.667"

// Process Parameters for Assemble_pairs_parse_log_AP:
//* params.Assemble_pairs_parse_log_AP.field_to_parse =  "ID REFID LENGTH OVERLAP GAP ERROR IDENTITY PVALUE EVALUE1 EVALUE2"  //* @input @description:"List of fields to collect, the fields should be seperated by space. Default ID REFID LENGTH OVERLAP GAP ERROR IDENTITY" @tooltip:"The sequence identifier may be specified using the hidden field name <ID>."

// Process Parameters for Filter_Sequence_Quality_filter_seq_quality:
//* params.Filter_Sequence_Quality_filter_seq_quality.method =  "quality"  //* @dropdown @options:"quality, length, missing" @description: "Filtering operation" @tooltip:"See https://presto.readthedocs.io/ for more details." @title:"General parametrs"
//* params.Filter_Sequence_Quality_filter_seq_quality.nproc =  "40"  //* @input @description:"Number of nproc to use for running FilterSeq. Default value 1."
//* params.Filter_Sequence_Quality_filter_seq_quality.q =  "20"  //* @input @description:"Discards reads where the Phred score below this threshold. Default value 20." @title:"Quality params"
//* params.Filter_Sequence_Quality_filter_seq_quality.n_length =  "35"  //* @input @description:"Minimum sequence length to retain. Default value 35." @title:"Length params"
//* params.Filter_Sequence_Quality_filter_seq_quality.n_missing =  "10"  //* @input @description:"Threshold for fraction of gap or N nucleotides. Default value 10." @title:"Missing params"

// Process Parameters for Mask_Primer_global_primer_MaskPrimers:
//* params.Mask_Primer_global_primer_MaskPrimers.method =  ["align"]  //* @dropdown @options:"score, align, extract" @description: "MaskPrimer primer identification options. Default score " @tooltip:"See https://presto.readthedocs.io/ for more details."
//* params.Mask_Primer_global_primer_MaskPrimers.barcode_field =  ["BARCODE"]  //* @input @description:"Name of the annotation field containing the barcode name. Default BARCODE"
//* params.Mask_Primer_global_primer_MaskPrimers.primer_field =  ["SPRIMER"]  //* @input @description:"Name of the annotation field containing the primer name. Default PRIMER"
//* params.Mask_Primer_global_primer_MaskPrimers.barcode =  ["false"]  //* @checkbox @description:"Check the box to remove the sequence preceding the extracted region and annotate the read with that sequence. Default false"
//* params.Mask_Primer_global_primer_MaskPrimers.revpr =  ["true"]  //* @checkbox @description:"Check the box to activate revpr option. Default false." @tooltip:"Specify to match the tail-end of the sequence against the reverse complement of the primers. This also reverses the behavior of the <maxlen> argument, such that the search window begins at the tail-end of the sequence."
//* params.Mask_Primer_global_primer_MaskPrimers.mode =  ["cut"]  //* @dropdown @options:"cut, mask, trim, tag" @description: "Which action to take with the primer sequence. Default cut." @tooltip:"The *cut* mode will remove both the primer region and the preceding sequence. The *mask* mode will replace the primer region with Ns and remove the preceding sequence. The *trim* mode will remove the region preceding the primer, but leave the primer region intact. The *tag* mode will leave the input sequence unmodified."
//* params.Mask_Primer_global_primer_MaskPrimers.failed =  "true"  //* @checkbox @description:"Check the box to output the failed sequences. Default false"
//* params.Mask_Primer_global_primer_MaskPrimers.nproc =  "20"  //* @input @description:"Number of nproc to use for running MaskPrimers. Default value 1."
//* params.Mask_Primer_global_primer_MaskPrimers.maxerror =  ["0.5"]  //* @input @description:"Maximum allowable error rate. Default value 0.2."
//* params.Mask_Primer_global_primer_MaskPrimers.umi_length =  ["0"]  //* @input @description:"The UMI length. Default value 0." @tooltip:"In the score and extract methods, setting a  <umi_length> will be added to the set <start> primer position. Such that the primer will start at the end of the UMI."
//* params.Mask_Primer_global_primer_MaskPrimers.start =  ["0"]  //* @input @description:"The starting position of the primer. Default 0"
//* params.Mask_Primer_global_primer_MaskPrimers.extract_length =  [""]  //* @input @description:"The sequence length to extract, only applicable for method extract. Default value 0."
//* params.Mask_Primer_global_primer_MaskPrimers.maxlen =  ["20"]  //* @input @description:"Length of the sequence window to scan for primers. Default value 50."
//* params.Mask_Primer_global_primer_MaskPrimers.skiprc =  ["false"]  //* @checkbox @description:"Check the box to prevent checking of sample reverse complement sequences. Default false"
params.Mask_Primer_global_primer_MaskPrimers.R1_primers = "${projectDir}/primers/IGHC_primers_universal.fasta"

// Process Parameters for Mask_Primer_constant_region_MaskPrimers:
//* params.Mask_Primer_constant_region_MaskPrimers.method =  ["score"]  //* @dropdown @options:"score, align, extract" @description: "MaskPrimer primer identification options. Default score " @tooltip:"See https://presto.readthedocs.io/ for more details."
//* params.Mask_Primer_constant_region_MaskPrimers.barcode_field =  ["BARCODE"]  //* @input @description:"Name of the annotation field containing the barcode name. Default BARCODE"
//* params.Mask_Primer_constant_region_MaskPrimers.primer_field =  ["CPRIMER"]  //* @input @description:"Name of the annotation field containing the primer name. Default PRIMER"
//* params.Mask_Primer_constant_region_MaskPrimers.barcode =  ["true"]  //* @checkbox @description:"Check the box to remove the sequence preceding the extracted region and annotate the read with that sequence. Default false"
//* params.Mask_Primer_constant_region_MaskPrimers.revpr =  ["true"]  //* @checkbox @description:"Check the box to activate revpr option. Default false." @tooltip:"Specify to match the tail-end of the sequence against the reverse complement of the primers. This also reverses the behavior of the <maxlen> argument, such that the search window begins at the tail-end of the sequence."
//* params.Mask_Primer_constant_region_MaskPrimers.mode =  ["cut"]  //* @dropdown @options:"cut, mask, trim, tag" @description: "Which action to take with the primer sequence. Default cut." @tooltip:"The *cut* mode will remove both the primer region and the preceding sequence. The *mask* mode will replace the primer region with Ns and remove the preceding sequence. The *trim* mode will remove the region preceding the primer, but leave the primer region intact. The *tag* mode will leave the input sequence unmodified."
//* params.Mask_Primer_constant_region_MaskPrimers.failed =  "true"  //* @checkbox @description:"Check the box to output the failed sequences. Default false"
//* params.Mask_Primer_constant_region_MaskPrimers.nproc =  "20"  //* @input @description:"Number of nproc to use for running MaskPrimers. Default value 1."
//* params.Mask_Primer_constant_region_MaskPrimers.maxerror =  ["0.2"]  //* @input @description:"Maximum allowable error rate. Default value 0.2."
//* params.Mask_Primer_constant_region_MaskPrimers.umi_length =  ["15"]  //* @input @description:"The UMI length. Default value 0." @tooltip:"In the score and extract methods, setting a  <umi_length> will be added to the set <start> primer position. Such that the primer will start at the end of the UMI."
//* params.Mask_Primer_constant_region_MaskPrimers.start =  ["0"]  //* @input @description:"The starting position of the primer. Default 0"
//* params.Mask_Primer_constant_region_MaskPrimers.extract_length =  [""]  //* @input @description:"The sequence length to extract, only applicable for method extract. Default value 0."
//* params.Mask_Primer_constant_region_MaskPrimers.maxlen =  ["50"]  //* @input @description:"Length of the sequence window to scan for primers. Default value 50."
//* params.Mask_Primer_constant_region_MaskPrimers.skiprc =  ["false"]  //* @checkbox @description:"Check the box to prevent checking of sample reverse complement sequences. Default false"
params.Mask_Primer_constant_region_MaskPrimers.R1_primers = "${projectDir}/primers/IGHC_primers_v1.fasta"

// Process Parameters for Mask_Primer_v_region_MaskPrimers:
//* params.Mask_Primer_v_region_MaskPrimers.method =  ["align"]  //* @dropdown @options:"score, align, extract" @description: "MaskPrimer primer identification options. Default score " @tooltip:"See https://presto.readthedocs.io/ for more details."
//* params.Mask_Primer_v_region_MaskPrimers.barcode_field =  ["BARCODE"]  //* @input @description:"Name of the annotation field containing the barcode name. Default BARCODE"
//* params.Mask_Primer_v_region_MaskPrimers.primer_field =  ["VPRIMER"]  //* @input @description:"Name of the annotation field containing the primer name. Default PRIMER"
//* params.Mask_Primer_v_region_MaskPrimers.barcode =  ["false"]  //* @checkbox @description:"Check the box to remove the sequence preceding the extracted region and annotate the read with that sequence. Default false"
//* params.Mask_Primer_v_region_MaskPrimers.revpr =  ["true"]  //* @checkbox @description:"Check the box to activate revpr option. Default false." @tooltip:"Specify to match the tail-end of the sequence against the reverse complement of the primers. This also reverses the behavior of the <maxlen> argument, such that the search window begins at the tail-end of the sequence."
//* params.Mask_Primer_v_region_MaskPrimers.mode =  ["cut"]  //* @dropdown @options:"cut, mask, trim, tag" @description: "Which action to take with the primer sequence. Default cut." @tooltip:"The *cut* mode will remove both the primer region and the preceding sequence. The *mask* mode will replace the primer region with Ns and remove the preceding sequence. The *trim* mode will remove the region preceding the primer, but leave the primer region intact. The *tag* mode will leave the input sequence unmodified."
//* params.Mask_Primer_v_region_MaskPrimers.failed =  "true"  //* @checkbox @description:"Check the box to output the failed sequences. Default false"
//* params.Mask_Primer_v_region_MaskPrimers.nproc =  "20"  //* @input @description:"Number of nproc to use for running MaskPrimers. Default value 1."
//* params.Mask_Primer_v_region_MaskPrimers.maxerror =  ["0.2"]  //* @input @description:"Maximum allowable error rate. Default value 0.2."
//* params.Mask_Primer_v_region_MaskPrimers.umi_length =  ["0"]  //* @input @description:"The UMI length. Default value 0." @tooltip:"In the score and extract methods, setting a  <umi_length> will be added to the set <start> primer position. Such that the primer will start at the end of the UMI."
//* params.Mask_Primer_v_region_MaskPrimers.start =  ["0"]  //* @input @description:"The starting position of the primer. Default 0"
//* params.Mask_Primer_v_region_MaskPrimers.extract_length =  [""]  //* @input @description:"The sequence length to extract, only applicable for method extract. Default value 0."
//* params.Mask_Primer_v_region_MaskPrimers.maxlen =  ["50"]  //* @input @description:"Length of the sequence window to scan for primers. Default value 50."
//* params.Mask_Primer_v_region_MaskPrimers.skiprc =  ["false"]  //* @checkbox @description:"Check the box to prevent checking of sample reverse complement sequences. Default false"
params.Mask_Primer_constant_region_MaskPrimers.R1_primers = "${projectDir}/primers/BIOMED2_VHFR1.fasta"

// Process Parameters for Cluster_Sets_BARCODE_cluster_sets:
//* params.Cluster_Sets_BARCODE_cluster_sets.method =  ["barcode"]  //* @dropdown @options:"set,all,barcode" @description:"Clustering method." @tooltip:"Set - Cluster sequences within annotation sets.\nAll - Cluster all sequences regardless of annotation.\nBarcode - Cluster reads by clustering barcode sequences"
//* params.Cluster_Sets_BARCODE_cluster_sets.failed =  ["false"]  //* @checkbox @description:"Check the box to activate failed option. " @tooltip:"Specify to output the failed sequences as well."
//* params.Cluster_Sets_BARCODE_cluster_sets.nproc =  "20"  //* @input @description:"Number of nproc to use for running MaskPrimers. Default value 1."
//* params.Cluster_Sets_BARCODE_cluster_sets.cluster_field =  ["INDEX_UMI"]  //* @input @description:"The name of the output annotation field to add with the cluster information for each sequence. Default CLUSTER."
//* params.Cluster_Sets_BARCODE_cluster_sets.ident =  ["0.9"]  //* @input @description:"The sequence identity threshold to use for clustering. Default 0.9" @tooltip:" Note, how identity is calculated is specific to the clustering application used."
//* params.Cluster_Sets_BARCODE_cluster_sets.length =  ["0"]  //* @input @description:"The minimum allowed shorter/longer sequence length ratio allowed within a cluster. Default 0" @tooltip:"Setting this value to 1.0 will require identical length matches within clusters. A value of 0.0 will allow clusters containing any length of substring."
//* params.Cluster_Sets_BARCODE_cluster_sets.prefix =  [""]  //* @input @description:"A string to use as the prefix for each cluster identifier. By default, cluster identifiers will be numeric values only. Default none"
//* params.Cluster_Sets_BARCODE_cluster_sets.cluster_tool =  ["usearch"]  //* @dropdown @options:"usearch,vsearch,cd-hit-est" @description:"The clustering tool to use for assigning clusters. Default usearch" @tooltip:"Must be one of usearch, vsearch or cd-hit-est. Note, for cd-hit-est the maximum memory limit is set to 3GB."
//* params.Cluster_Sets_BARCODE_cluster_sets.cluster_exec =  [""]  //* @input @description:"The name or path of the usearch, vsearch or cd-hit-est executable."
//* params.Cluster_Sets_BARCODE_cluster_sets.set_field =  ["BARCODE"]  //* @input @description:"The annotation field containing annotations, such as UMI barcode, for sequence grouping. Default BARCODE"
//* params.Cluster_Sets_BARCODE_cluster_sets.start =  ["0"]  //* @input @desciption:"The start of the region to be used for clustering. Together with end, this parameter can be used to specify a subsequence of each read to use in the clustering algorithm. Default 0"
//* params.Cluster_Sets_BARCODE_cluster_sets.end =  [""]  //* @input @description:"The end of the region to be used for clustering. Default none"
//* params.Cluster_Sets_BARCODE_cluster_sets.barcode_field =  ["BARCODE"]  //* @input @description:"The annotation field containing barcode sequences. Default BARCODE"

// Process Parameters for Cluster_Sets_VDJ_cluster_sets:
//* params.Cluster_Sets_VDJ_cluster_sets.method =  ["set"]  //* @dropdown @options:"set,all,barcode" @description:"Clustering method." @tooltip:"Set - Cluster sequences within annotation sets.\nAll - Cluster all sequences regardless of annotation.\nBarcode - Cluster reads by clustering barcode sequences"
//* params.Cluster_Sets_VDJ_cluster_sets.failed =  ["false"]  //* @checkbox @description:"Check the box to activate failed option. " @tooltip:"Specify to output the failed sequences as well."
//* params.Cluster_Sets_VDJ_cluster_sets.nproc =  "20"  //* @input @description:"Number of nproc to use for running MaskPrimers. Default value 1."
//* params.Cluster_Sets_VDJ_cluster_sets.cluster_field =  ["INDEX_SEQ"]  //* @input @description:"The name of the output annotation field to add with the cluster information for each sequence. Default CLUSTER."
//* params.Cluster_Sets_VDJ_cluster_sets.ident =  ["0.8"]  //* @input @description:"The sequence identity threshold to use for clustering. Default 0.9" @tooltip:" Note, how identity is calculated is specific to the clustering application used."
//* params.Cluster_Sets_VDJ_cluster_sets.length =  ["0"]  //* @input @description:"The minimum allowed shorter/longer sequence length ratio allowed within a cluster. Default 0" @tooltip:"Setting this value to 1.0 will require identical length matches within clusters. A value of 0.0 will allow clusters containing any length of substring."
//* params.Cluster_Sets_VDJ_cluster_sets.prefix =  [""]  //* @input @description:"A string to use as the prefix for each cluster identifier. By default, cluster identifiers will be numeric values only. Default none"
//* params.Cluster_Sets_VDJ_cluster_sets.cluster_tool =  ["usearch"]  //* @dropdown @options:"usearch,vsearch,cd-hit-est" @description:"The clustering tool to use for assigning clusters. Default usearch" @tooltip:"Must be one of usearch, vsearch or cd-hit-est. Note, for cd-hit-est the maximum memory limit is set to 3GB."
//* params.Cluster_Sets_VDJ_cluster_sets.cluster_exec =  [""]  //* @input @description:"The name or path of the usearch, vsearch or cd-hit-est executable."
//* params.Cluster_Sets_VDJ_cluster_sets.set_field =  ["INDEX_UMI"]  //* @input @description:"The annotation field containing annotations, such as UMI barcode, for sequence grouping. Default BARCODE"
//* params.Cluster_Sets_VDJ_cluster_sets.start =  ["0"]  //* @input @desciption:"The start of the region to be used for clustering. Together with end, this parameter can be used to specify a subsequence of each read to use in the clustering algorithm. Default 0"
//* params.Cluster_Sets_VDJ_cluster_sets.end =  [""]  //* @input @description:"The end of the region to be used for clustering. Default none"
//* params.Cluster_Sets_VDJ_cluster_sets.barcode_field =  [""]  //* @input @description:"The annotation field containing barcode sequences. Default BARCODE"

// Process Parameters for Parse_header_parse_headers:
//* params.Parse_header_parse_headers.method =  "merge"  //* @dropdown @options:"collapse, add, copy, delete, expand, merge, rename, table" @description: "Parse method. Default collapse (See https://presto.readthedocs.io/ for more details.)"
//* params.Parse_header_parse_headers.act =  "cat"  //* @dropdown @options:"min, max, sum, first, last, set, cat" @description: "Only applicable for methods collapse and add. List of actions to take for each field defining how each annotation will be combined into a single value. Default min (See https://presto.readthedocs.io/ for more details.)"
//* params.Parse_header_parse_headers.args =  "-f BARCODE INDEX_UMI CPRIMER INDEX_SEQ -k INDEX_MERGE"  //* @input @description: "Additional arrguments for ParseHeader function. Defualt is '-f CONSCOUNT' for method collapse."

// Process Parameters for Build_Consensus_build_consensus:
//* params.Build_Consensus_build_consensus.failed =  "false"  //* @checkbox @description:"Check the box to activate failed option. " @tooltip:"Specify to output the failed sequences as well." @title:"General params"
//* params.Build_Consensus_build_consensus.nproc =  "20"  //* @input @description:"Number of nproc to use for running MaskPrimers. Default value 1."
//* params.Build_Consensus_build_consensus.barcode_field =  ["INDEX_MERGE"]  //* @input @description:"Position of description barcode field to group sequences by. Default BARCODE." @title:"Consensus generation copy fields and actions"
//* params.Build_Consensus_build_consensus.primer_field =  ["CPRIMER"]  //* @input @description:"Specifies the field name of the primer annotations. Default is none." @tooltip:"In most processing pipeline this parameter is set to PRIMER"
//* params.Build_Consensus_build_consensus.act =  ["set"]  //* @dropdown @options:"none,min,max,sum,set,majority" @description:"List of actions to take for each copy field which defines how each annotation will be combined into a single value. Default none." @tooltip:"The actions “min”, “max”, “sum” perform the corresponding mathematical operation on numeric annotations. The action “set” combines annotations into a comma delimited list of unique values and adds an annotation named <FIELD>_COUNT specifying the count of each item in the set. The action “majority” assigns the most frequent annotation to the consensus annotation and adds an annotation named <FIELD>_FREQ specifying the frequency of the majority value."
//* params.Build_Consensus_build_consensus.copy_field =  ["VPRIMER"]  //* @input @description:"Specifies a set of additional annotation fields to copy into the consensus sequence annotations. Default None" @tooltip:"If an action is specified under the <act> paramter, a copy field is needed as well."
//* params.Build_Consensus_build_consensus.mincount =  [1]  //* @input @description:"The minimum number of sequences needed to define a valid consensus. Default is 1" @title:"Consensus generation groups params"
//* params.Build_Consensus_build_consensus.minqual =  [0]  //* @input @description:"Consensus quality score cut-off under which an ambiguous character is assigned. Default value 0." @tooltip:"Does not apply when quality scores are unavailable."
//* params.Build_Consensus_build_consensus.minfreq =  [0.6]  //* @input @description:"Fraction of character occurrences under which an ambiguous character is assigned. Default value 0.6."
//* params.Build_Consensus_build_consensus.maxerror =  ["none"]  //* @input @description:"Maximum allowable error rate. Default is none (A numeric field from 0 to 1)." @tooltip:"Specify to calculate the error rate of each read group (rate of mismatches from consensus) and remove groups exceeding the given error threshold. The error rate is calculated against the final consensus sequence, which may include masked positions due to the <minqual> and <minfreq> arguments and may have deleted positions due to the <maxgap> argument. Mutually exclusive with <maxdiv>."
//* params.Build_Consensus_build_consensus.prcons =  ["none"]  //* @input @description:"Minimum primer frequency required. Default is none (A numeric field from 0 to 1)." @tooltip:"Specify to define a minimum primer frequency required to assign a consensus primer, and filter out sequences with minority primers from the consensus building step."
//* params.Build_Consensus_build_consensus.maxgap =  ["none"]  //* @input @description:"A cut-off for the frequency allowed gao values for each position. Default is none (A numeric field from 0 to 1)." @tooltip:"If specified, this defines a cut-off for the frequency of allowed gap values for each position. Positions exceeding the threshold are deleted from the consensus. If not defined, positions are always retained. "
//* params.Build_Consensus_build_consensus.maxdiv =  ["none"]  //* @input @description:"Maximum allowable diversity rate. Default is none (A numeric field from 0 to 1)" @tooltip:"Specify to calculate the nucleotide diversity of each read group (average pairwise error rate) and remove groups exceeding the given diversity threshold. Diversity is calculate for all positions within the read group, ignoring any character filtering imposed by the <minqual>, <minfreq> and <maxgap> arguments. Mutually exclusive with <maxerror>."
//* params.Build_Consensus_build_consensus.dep =  ["false"]  //* @checkbox @description:"Check the box to calculate consensus quality with a non-independence assumption. Default false"

// Process Parameters for collapse_sequences_collapse_seq:
//* params.collapse_sequences_collapse_seq.max_missing =  20  //* @input @description:"Maximum number of missing nucleotides to consider for collapsing sequences. A sequence will be considered undetermined if it contains too many missing nucleotides. Default is 0"
//* params.collapse_sequences_collapse_seq.inner =  "true"  //* @checkbox @description:"Exclude consecutive missing characters at either end of the sequence. Default is false."
//* params.collapse_sequences_collapse_seq.fasta =  "true"  //* @checkbox @description:"Specify to force output as FASTA rather than FASTQ. Default is false."
//* params.collapse_sequences_collapse_seq.act =  "set"  //* @dropdown @options:"none, min, max, sum, set" @description:"Only applicable for methods collapse and add. List of actions to take for each field defining how each annotation will be combined into a single value. Default none"
//* params.collapse_sequences_collapse_seq.uf =  "CLUSTER"  //* @input @description:"Specifies a set of annotation fields that must match for sequences to be considered duplicates. Default none"
//* params.collapse_sequences_collapse_seq.cf =  "VPRIMER"  //* @input @description:"Specifies a set of annotation fields to copy into the unique sequence output. Default none"
//* params.collapse_sequences_collapse_seq.nproc =  "40"  //* @input @description:"Number of nproc to use for running FilterSeq. Default value 1."
//* params.collapse_sequences_collapse_seq.failed =  "true"  //* @checkbox @description:"Check the box to activate failed option. " @tooltip:"Specify to output the failed sequences as well." @title:"General params"

if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 
if (!params.mate_single){params.mate_single = ""} 

if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_0_0_g2_12}
 } else {  
	g_0_0_g2_12 = Channel.empty()
 }

Channel.value(params.mate).into{g_1_1_g2_15;g_1_1_g2_19;g_1_1_g2_12}
Channel.value(params.mate_single).into{g_4_1_g32_14;g_4_1_g35_14;g_4_1_g37_0;g_4_1_g37_5;g_4_0_g37_7;g_4_0_g6_9;g_4_1_g6_12;g_4_0_g6_11;g_4_0_g26_9;g_4_1_g26_12;g_4_0_g26_11;g_4_0_g5_9;g_4_1_g5_12;g_4_0_g5_11;g_4_1_g15_14;g_4_1_g15_12;g_4_1_g15_10;g_4_1_g21_16}


process Assemble_pairs_assemble_pairs {

input:
 set val(name),file(reads) from g_0_0_g2_12
 val mate from g_1_1_g2_12

output:
 set val(name),file("*_assemble-pass.f*")  into g2_12_reads00_g37_0
 set val(name),file("AP_*")  into g2_12_logFile10_g2_15
 set val(name),file("*_assemble-fail.f*") optional true  into g2_12_reads_failed22
 set val(name),file("out*")  into g2_12_logFile33

script:
method = params.Assemble_pairs_assemble_pairs.method
coord = params.Assemble_pairs_assemble_pairs.coord
rc = params.Assemble_pairs_assemble_pairs.rc
head_fields_R1 = params.Assemble_pairs_assemble_pairs.head_fields_R1
head_fields_R2 = params.Assemble_pairs_assemble_pairs.head_fields_R2
failed = params.Assemble_pairs_assemble_pairs.failed
fasta = params.Assemble_pairs_assemble_pairs.fasta
nproc = params.Assemble_pairs_assemble_pairs.nproc
alpha = params.Assemble_pairs_assemble_pairs.alpha
maxerror = params.Assemble_pairs_assemble_pairs.maxerror
minlen = params.Assemble_pairs_assemble_pairs.minlen
maxlen = params.Assemble_pairs_assemble_pairs.maxlen
scanrev = params.Assemble_pairs_assemble_pairs.scanrev
minident = params.Assemble_pairs_assemble_pairs.minident
evalue = params.Assemble_pairs_assemble_pairs.evalue
maxhits = params.Assemble_pairs_assemble_pairs.maxhits
fill = params.Assemble_pairs_assemble_pairs.fill
aligner = params.Assemble_pairs_assemble_pairs.aligner
// align_exec = params.Assemble_pairs_assemble_pairs.// align_exec
// dbexec = params.Assemble_pairs_assemble_pairs.// dbexec
gap = params.Assemble_pairs_assemble_pairs.gap
usearch_version = params.Assemble_pairs_assemble_pairs.usearch_version
assemble_reference = params.Assemble_pairs_assemble_pairs.assemble_reference
head_seqeunce_file = params.Assemble_pairs_assemble_pairs.head_seqeunce_file
//* @style @condition:{method="align",alpha,maxerror,minlen,maxlen,scanrev}, {method="sequential",alpha,maxerror,minlen,maxlen,scanrev,ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="reference",ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="join",gap} @multicolumn:{method,coord,rc,head_fields_R1,head_fields_R2,failed,nrpoc,usearch_version},{alpha,maxerror,minlen,maxlen,scanrev}, {ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec}, {gap} 

// args
coord = "--coord ${coord}"
rc = "--rc ${rc}"
head_fields_R1 = (head_fields_R1!="") ? "--1f ${head_fields_R1}" : ""
head_fields_R2 = (head_fields_R2!="") ? "--2f ${head_fields_R2}" : ""
failed = (failed=="false") ? "" : "--failed"
fasta = (fasta=="false") ? "" : "--fasta"
nproc = "--nproc ${nproc}"

scanrev = (scanrev=="false") ? "" : "--scanrev"
fill = (fill=="false") ? "" : "--fill"

// align_exec = (align_exec!="") ? "--exec ${align_exec}" : ""
// dbexec = (dbexec!="") ? "--dbexec ${dbexec}" : ""


ref_file = (assemble_reference!='') ? "-r ${assemble_reference}" : ""



args = ""

if(method=="align"){
	args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev}"
}else{
	if(method=="sequential"){
		args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev} ${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
	}else{
		if(method=="reference"){
			args = "${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
		}else{
			args = "--gap ${gap}"
		}
	}
}


readArray = reads.toString().split(' ')	


if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	
	if(R1.contains("."+head_seqeunce_file)){
		R1 = readArray[0]
		R2 = readArray[1]
	}else{
		R2 = readArray[0]
		R1 = readArray[1]
	}
	
	"""
	if [ "${method}" != "align" ]; then
		if  [ "${aligner}" == "usearch" ]; then
			wget -q --show-progress --no-check-certificate https://drive5.com/downloads/usearch${usearch_version}_i86linux32.gz
			gunzip usearch${usearch_version}_i86linux32.gz
			chmod +x usearch${usearch_version}_i86linux32
			mv usearch${usearch_version}_i86linux32 /usr/local/bin/usearch2
			align_exec="--exec /usr/local/bin/usearch2"
			dbexec="--dbexec /usr/local/bin/usearch2"
		else
			align_exec="--exec /usr/local/bin/blastn"
			dbexec="--dbexec /usr/local/bin/makeblastdb"
		fi
	else
		align_exec=""
		dbexec=""
	fi

	AssemblePairs.py ${method} -1 ${R1} -2 ${R2} ${coord} ${rc} ${head_fields_R1} ${head_fields_R2} ${args} \$align_exec \$dbexec ${fasta} ${failed} --log AP_${name}.log ${nproc}  2>&1 | tee out_${R1}_AP.log
	"""

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}

}


process Assemble_pairs_parse_log_AP {

input:
 set val(name),file(log_file) from g2_12_logFile10_g2_15
 val mate from g_1_1_g2_15

output:
 set val(name),file("*.tab")  into g2_15_logFile00_g2_19

script:
field_to_parse = params.Assemble_pairs_parse_log_AP.field_to_parse
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ${field_to_parse}
"""


}


process Assemble_pairs_report_assemble_pairs {

input:
 set val(name),file(log_files) from g2_15_logFile00_g2_19
 val matee from g_1_1_g2_19

output:
 file "*.rmd"  into g2_19_rMarkdown00_g2_22



shell:

if(matee=="pair"){
	readArray = log_files.toString().split(' ')
	assemble = readArray[0]
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("assemble_length", "Histogram showing the distribution assembled sequence lengths in 
	                            nucleotides for the Align step (top) and Reference step (bottom).")
	figures("assemble_overlap", "Histogram showing the distribution of overlapping nucleotides between 
	                             mate-pairs for the Align step (top) and Reference step (bottom).
	                             Negative values for overlap indicate non-overlapping mate-pairs
	                             with the negative value being the number of gap characters between
	                             the ends of the two mate-pairs.")
	figures("assemble_error", "Histograms showing the distribution of paired-end assembly error 
	                           rates for the Align step (top) and identity to the reference germline 
	                           for the Reference step (bottom).")
	figures("assemble_pvalue", "Histograms showing the distribution of significance scores for 
	                            paired-end assemblies. P-values for the Align mode are shown in the top
	                            panel. E-values from the Reference step's alignment against the 
	                            germline sequences are shown in the bottom panel for both input files
	                            separately.")
	```
	
	```{r, echo=FALSE, warning=FALSE}
	assemble_log <- loadLogTable(file.path(".", "!{assemble}"))
	
	# Subset to align and reference logs
	align_fields <- c("ERROR", "PVALUE")
	ref_fields <- c("REFID", "GAP", "EVALUE1", "EVALUE2", "IDENTITY")
	align_log <- assemble_log[!is.na(assemble_log$ERROR), !(names(assemble_log) %in% ref_fields)]
	ref_log <- assemble_log[!is.na(assemble_log$REFID), !(names(assemble_log) %in% align_fields)]
	
	# Build log set
	assemble_list <- list()
	if (nrow(align_log) > 0) { assemble_list[["Align"]] <- align_log }
	if (nrow(ref_log) > 0) { assemble_list[["Reference"]] <- ref_log }
	plot_titles <- names(assemble_list)
	```
	
	# Paired-End Assembly
	
	Assembly of paired-end reads is performed using the AssemblePairs tool which 
	determines the read overlap in two steps. First, de novo assembly is attempted 
	using an exhaustive approach to identify all possible overlaps between the 
	two reads with alignment error rates and p-values below user-defined thresholds. 
	This method is denoted as the `Align` method in the following figures. 
	Second, those reads failing the first stage of de novo assembly are then 
	mapped to the V-region reference sequences to create a full length sequence, 
	padding with Ns, for any amplicons that have insufficient overlap for 
	de novo assembly. This second stage is referred to as the `Reference` step in the
	figures below.
	
	## Assembled sequence lengths
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="length", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_length")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="overlap", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_overlap")`
	
	## Alignment error rates and significance
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="error", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_error")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="pvalue", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```

	`r figures("assemble_pvalue")`

	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}
}


process Assemble_pairs_render_rmarkdown {

input:
 file rmk from g2_19_rMarkdown00_g2_22

output:
 file "*.html"  into g2_22_outputFileHTML00

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Filter_Sequence_Quality_filter_seq_quality {

input:
 set val(name),file(reads) from g2_12_reads00_g37_0
 val mate from g_4_1_g37_0

output:
 set val(name), file("*_${method}-pass.fast*")  into g37_0_reads01_g6_11
 set val(name), file("FS_*")  into g37_0_logFile10_g37_5
 set val(name), file("*_${method}-fail.fast*") optional true  into g37_0_reads22
 set val(name),file("out*") optional true  into g37_0_logFile33

script:
method = params.Filter_Sequence_Quality_filter_seq_quality.method
nproc = params.Filter_Sequence_Quality_filter_seq_quality.nproc
q = params.Filter_Sequence_Quality_filter_seq_quality.q
n_length = params.Filter_Sequence_Quality_filter_seq_quality.n_length
n_missing = params.Filter_Sequence_Quality_filter_seq_quality.n_missing
fasta = params.Filter_Sequence_Quality_filter_seq_quality.fasta
//* @style @condition:{method="quality",q}, {method="length",n_length}, {method="missing",n_missing} @multicolumn:{method,nproc}

if(method=="missing"){
	q = ""
	n_length = ""
	n_missing = "-n ${n_missing}"
}else{
	if(method=="length"){
		q = ""
		n_length = "-n ${n_length}"
		n_missing = ""
	}else{
		q = "-q ${q}"
		n_length = ""
		n_missing = ""
	}
}

readArray = reads.toString().split(' ')	

fasta = (fasta=="true") ? "--fasta" : ""

if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R1_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R2_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	"""
}else{
	R1 = readArray[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	"""
}


}


process Mask_Primer_global_primer_MaskPrimers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_primers-fail.fastq$/) "MP_fail_reads_s_primer/$filename"}
input:
 val mate from g_4_0_g6_11
 set val(name),file(reads) from g37_0_reads01_g6_11

output:
 set val(name), file("*_primers-pass.fastq")  into g6_11_reads01_g26_11
 set val(name), file("*_primers-fail.fastq") optional true  into g6_11_reads_failed11
 set val(name), file("MP_*")  into g6_11_logFile21_g6_9
 set val(name),file("out*")  into g6_11_logFile33

script:
method = params.Mask_Primer_global_primer_MaskPrimers.method
barcode_field = params.Mask_Primer_global_primer_MaskPrimers.barcode_field
primer_field = params.Mask_Primer_global_primer_MaskPrimers.primer_field
barcode = params.Mask_Primer_global_primer_MaskPrimers.barcode
revpr = params.Mask_Primer_global_primer_MaskPrimers.revpr
mode = params.Mask_Primer_global_primer_MaskPrimers.mode
failed = params.Mask_Primer_global_primer_MaskPrimers.failed
nproc = params.Mask_Primer_global_primer_MaskPrimers.nproc
maxerror = params.Mask_Primer_global_primer_MaskPrimers.maxerror
umi_length = params.Mask_Primer_global_primer_MaskPrimers.umi_length
start = params.Mask_Primer_global_primer_MaskPrimers.start
extract_length = params.Mask_Primer_global_primer_MaskPrimers.extract_length
maxlen = params.Mask_Primer_global_primer_MaskPrimers.maxlen
skiprc = params.Mask_Primer_global_primer_MaskPrimers.skiprc
R1_primers = params.Mask_Primer_global_primer_MaskPrimers.R1_primers
R2_primers = params.Mask_Primer_global_primer_MaskPrimers.R2_primers
//* @style @condition:{method="score",umi_length,start,maxerror}{method="extract",umi_length,start},{method="align",maxerror,maxlen,skiprc}, {method="extract",start,extract_length} @array:{method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc} @multicolumn:{method,barcode_field,primer_field,barcode,revpr,mode,failed,nproc,maxerror,umi_length,start,extract_length,maxlen,skiprc}

method = (method.collect().size==2) ? method : [method[0],method[0]]
barcode_field = (barcode_field.collect().size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]
primer_field = (primer_field.collect().size==2) ? primer_field : [primer_field[0],primer_field[0]]
barcode = (barcode.collect().size==2) ? barcode : [barcode[0],barcode[0]]
revpr = (revpr.collect().size==2) ? revpr : [revpr[0],revpr[0]]
mode = (mode.collect().size==2) ? mode : [mode[0],mode[0]]
maxerror = (maxerror.collect().size==2) ? maxerror : [maxerror[0],maxerror[0]]
umi_length = (umi_length.collect().size==2) ? umi_length : [umi_length[0],umi_length[0]]
start = (start.collect().size==2) ? start : [start[0],start[0]]
extract_length = (extract_length.collect().size==2) ? extract_length : [extract_length[0],extract_length[0]]
maxlen = (maxlen.collect().size==2) ? maxlen : [maxlen[0],maxlen[0]]
skiprc = (skiprc.collect().size==2) ? skiprc : [skiprc[0],skiprc[0]]
failed = (failed=="true") ? "--failed" : ""

def args_values = [];
[method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc].transpose().each { m,bf,pf,bc,rp,md,mr,ul,s,el,ml,sk -> {
    
    if(m=="align"){
        s = ""
    }else{
        if(bc=="false"){
            s = "--start ${s}"
        }else{
            s = s + ul
            s = "--start ${s}"
        }
    }
    
    el = (m=="extract") ? "--len ${el}" : ""
    mr = (m=="extract") ? "" : "--maxerror ${mr}" 
    ml = (m=="align") ? "--maxlen ${ml}" : "" 
    sk = (m=="align" && sk=="true") ? "--skiprc" : "" 
    
    PRIMER_FIELD = "${pf}"
    
    // all
    bf = (bf=="") ? "" : "--bf ${bf}"
    pf = (pf=="") ? "" : "--pf ${pf}"
    bc = (bc=="false") ? "" : "--barcode"
    rp = (rp=="false") ? "" : "--revpr"
    args_values.add("${m} --mode ${md} ${bf} ${pf} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk}")
    
    
}}

readArray = reads.toString().split(' ')
if(mate=="pair"){
	args_1 = args_values[0]
	args_2 = args_values[1]
	
  


	R1 = readArray[0]
	R2 = readArray[1]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
	
	
	"""
	
	MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	"""
}else{
	args_1 = args_values[0]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	
	R1 = readArray[0]

	"""
	echo -e "Assuming inputs for R1\n"
	
	MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	"""
}

}


process Mask_Primer_global_primer_parse_log_MP {

input:
 val mate from g_4_0_g6_9
 set val(name), file(log_file) from g6_11_logFile21_g6_9

output:
 set val(name), file("*.tab")  into g6_9_logFile00_g6_12

script:
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ID PRIMER BARCODE ERROR
"""

}


process Mask_Primer_global_primer_try_report_maskprimer {

input:
 set val(name), file(primers) from g6_9_logFile00_g6_12
 val matee from g_4_1_g6_12

output:
 file "*.rmd"  into g6_12_rMarkdown00_g6_16


shell:

if(matee=="pair"){
	readArray = primers.toString().split(' ')	
	primers_1 = readArray.grep(~/.*R1.*/)[0]
	primers_2 = readArray.grep(~/.*R2.*/)[0]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read 1", "Read 2")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom),",
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers_1}"))
	primer_log_2 <- loadLogTable(file.path(".", "!{primers_2}"))
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{

	readArray = primers.toString().split(' ')
	primers = readArray.grep(~/.*.tab*/)[0]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1],
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers}"))
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}
}


process Mask_Primer_global_primer_render_rmarkdown {

input:
 file rmk from g6_12_rMarkdown00_g6_16

output:
 file "*.html"  into g6_16_outputFileHTML00

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Mask_Primer_constant_region_MaskPrimers {

input:
 val mate from g_4_0_g26_11
 set val(name),file(reads) from g6_11_reads01_g26_11

output:
 set val(name), file("*_primers-pass.fastq")  into g26_11_reads01_g5_11
 set val(name), file("*_primers-fail.fastq") optional true  into g26_11_reads_failed11
 set val(name), file("MP_*")  into g26_11_logFile21_g26_9
 set val(name),file("out*")  into g26_11_logFile33

script:
method = params.Mask_Primer_constant_region_MaskPrimers.method
barcode_field = params.Mask_Primer_constant_region_MaskPrimers.barcode_field
primer_field = params.Mask_Primer_constant_region_MaskPrimers.primer_field
barcode = params.Mask_Primer_constant_region_MaskPrimers.barcode
revpr = params.Mask_Primer_constant_region_MaskPrimers.revpr
mode = params.Mask_Primer_constant_region_MaskPrimers.mode
failed = params.Mask_Primer_constant_region_MaskPrimers.failed
nproc = params.Mask_Primer_constant_region_MaskPrimers.nproc
maxerror = params.Mask_Primer_constant_region_MaskPrimers.maxerror
umi_length = params.Mask_Primer_constant_region_MaskPrimers.umi_length
start = params.Mask_Primer_constant_region_MaskPrimers.start
extract_length = params.Mask_Primer_constant_region_MaskPrimers.extract_length
maxlen = params.Mask_Primer_constant_region_MaskPrimers.maxlen
skiprc = params.Mask_Primer_constant_region_MaskPrimers.skiprc
R1_primers = params.Mask_Primer_constant_region_MaskPrimers.R1_primers
R2_primers = params.Mask_Primer_constant_region_MaskPrimers.R2_primers
//* @style @condition:{method="score",umi_length,start,maxerror}{method="extract",umi_length,start},{method="align",maxerror,maxlen,skiprc}, {method="extract",start,extract_length} @array:{method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc} @multicolumn:{method,barcode_field,primer_field,barcode,revpr,mode,failed,nproc,maxerror,umi_length,start,extract_length,maxlen,skiprc}

method = (method.collect().size==2) ? method : [method[0],method[0]]
barcode_field = (barcode_field.collect().size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]
primer_field = (primer_field.collect().size==2) ? primer_field : [primer_field[0],primer_field[0]]
barcode = (barcode.collect().size==2) ? barcode : [barcode[0],barcode[0]]
revpr = (revpr.collect().size==2) ? revpr : [revpr[0],revpr[0]]
mode = (mode.collect().size==2) ? mode : [mode[0],mode[0]]
maxerror = (maxerror.collect().size==2) ? maxerror : [maxerror[0],maxerror[0]]
umi_length = (umi_length.collect().size==2) ? umi_length : [umi_length[0],umi_length[0]]
start = (start.collect().size==2) ? start : [start[0],start[0]]
extract_length = (extract_length.collect().size==2) ? extract_length : [extract_length[0],extract_length[0]]
maxlen = (maxlen.collect().size==2) ? maxlen : [maxlen[0],maxlen[0]]
skiprc = (skiprc.collect().size==2) ? skiprc : [skiprc[0],skiprc[0]]
failed = (failed=="true") ? "--failed" : ""

def args_values = [];
[method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc].transpose().each { m,bf,pf,bc,rp,md,mr,ul,s,el,ml,sk -> {
    
    if(m=="align"){
        s = ""
    }else{
        if(bc=="false"){
            s = "--start ${s}"
        }else{
            s = s + ul
            s = "--start ${s}"
        }
    }
    
    el = (m=="extract") ? "--len ${el}" : ""
    mr = (m=="extract") ? "" : "--maxerror ${mr}" 
    ml = (m=="align") ? "--maxlen ${ml}" : "" 
    sk = (m=="align" && sk=="true") ? "--skiprc" : "" 
    
    PRIMER_FIELD = "${pf}"
    
    // all
    bf = (bf=="") ? "" : "--bf ${bf}"
    pf = (pf=="") ? "" : "--pf ${pf}"
    bc = (bc=="false") ? "" : "--barcode"
    rp = (rp=="false") ? "" : "--revpr"
    args_values.add("${m} --mode ${md} ${bf} ${pf} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk}")
    
    
}}

readArray = reads.toString().split(' ')
if(mate=="pair"){
	args_1 = args_values[0]
	args_2 = args_values[1]
	
  


	R1 = readArray[0]
	R2 = readArray[1]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
	
	
	"""
	
	MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	"""
}else{
	args_1 = args_values[0]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	
	R1 = readArray[0]

	"""
	echo -e "Assuming inputs for R1\n"
	
	MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	"""
}

}


process Mask_Primer_constant_region_parse_log_MP {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "MP_log_C_table/$filename"}
input:
 val mate from g_4_0_g26_9
 set val(name), file(log_file) from g26_11_logFile21_g26_9

output:
 set val(name), file("*.tab")  into g26_9_logFile00_g26_12

script:
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ID PRIMER BARCODE ERROR
"""

}


process Mask_Primer_constant_region_try_report_maskprimer {

input:
 set val(name), file(primers) from g26_9_logFile00_g26_12
 val matee from g_4_1_g26_12

output:
 file "*.rmd"  into g26_12_rMarkdown00_g26_16


shell:

if(matee=="pair"){
	readArray = primers.toString().split(' ')	
	primers_1 = readArray.grep(~/.*R1.*/)[0]
	primers_2 = readArray.grep(~/.*R2.*/)[0]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read 1", "Read 2")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom),",
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers_1}"))
	primer_log_2 <- loadLogTable(file.path(".", "!{primers_2}"))
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{

	readArray = primers.toString().split(' ')
	primers = readArray.grep(~/.*.tab*/)[0]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1],
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers}"))
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}
}


process Mask_Primer_constant_region_render_rmarkdown {

input:
 file rmk from g26_12_rMarkdown00_g26_16

output:
 file "*.html"  into g26_16_outputFileHTML00

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Mask_Primer_v_region_MaskPrimers {

input:
 val mate from g_4_0_g5_11
 set val(name),file(reads) from g26_11_reads01_g5_11

output:
 set val(name), file("*_primers-pass.fastq")  into g5_11_reads00_g35_14
 set val(name), file("*_primers-fail.fastq") optional true  into g5_11_reads_failed11
 set val(name), file("MP_*")  into g5_11_logFile21_g5_9
 set val(name),file("out*")  into g5_11_logFile33

script:
method = params.Mask_Primer_v_region_MaskPrimers.method
barcode_field = params.Mask_Primer_v_region_MaskPrimers.barcode_field
primer_field = params.Mask_Primer_v_region_MaskPrimers.primer_field
barcode = params.Mask_Primer_v_region_MaskPrimers.barcode
revpr = params.Mask_Primer_v_region_MaskPrimers.revpr
mode = params.Mask_Primer_v_region_MaskPrimers.mode
failed = params.Mask_Primer_v_region_MaskPrimers.failed
nproc = params.Mask_Primer_v_region_MaskPrimers.nproc
maxerror = params.Mask_Primer_v_region_MaskPrimers.maxerror
umi_length = params.Mask_Primer_v_region_MaskPrimers.umi_length
start = params.Mask_Primer_v_region_MaskPrimers.start
extract_length = params.Mask_Primer_v_region_MaskPrimers.extract_length
maxlen = params.Mask_Primer_v_region_MaskPrimers.maxlen
skiprc = params.Mask_Primer_v_region_MaskPrimers.skiprc
R1_primers = params.Mask_Primer_v_region_MaskPrimers.R1_primers
R2_primers = params.Mask_Primer_v_region_MaskPrimers.R2_primers
//* @style @condition:{method="score",umi_length,start,maxerror}{method="extract",umi_length,start},{method="align",maxerror,maxlen,skiprc}, {method="extract",start,extract_length} @array:{method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc} @multicolumn:{method,barcode_field,primer_field,barcode,revpr,mode,failed,nproc,maxerror,umi_length,start,extract_length,maxlen,skiprc}

method = (method.collect().size==2) ? method : [method[0],method[0]]
barcode_field = (barcode_field.collect().size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]
primer_field = (primer_field.collect().size==2) ? primer_field : [primer_field[0],primer_field[0]]
barcode = (barcode.collect().size==2) ? barcode : [barcode[0],barcode[0]]
revpr = (revpr.collect().size==2) ? revpr : [revpr[0],revpr[0]]
mode = (mode.collect().size==2) ? mode : [mode[0],mode[0]]
maxerror = (maxerror.collect().size==2) ? maxerror : [maxerror[0],maxerror[0]]
umi_length = (umi_length.collect().size==2) ? umi_length : [umi_length[0],umi_length[0]]
start = (start.collect().size==2) ? start : [start[0],start[0]]
extract_length = (extract_length.collect().size==2) ? extract_length : [extract_length[0],extract_length[0]]
maxlen = (maxlen.collect().size==2) ? maxlen : [maxlen[0],maxlen[0]]
skiprc = (skiprc.collect().size==2) ? skiprc : [skiprc[0],skiprc[0]]
failed = (failed=="true") ? "--failed" : ""

def args_values = [];
[method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc].transpose().each { m,bf,pf,bc,rp,md,mr,ul,s,el,ml,sk -> {
    
    if(m=="align"){
        s = ""
    }else{
        if(bc=="false"){
            s = "--start ${s}"
        }else{
            s = s + ul
            s = "--start ${s}"
        }
    }
    
    el = (m=="extract") ? "--len ${el}" : ""
    mr = (m=="extract") ? "" : "--maxerror ${mr}" 
    ml = (m=="align") ? "--maxlen ${ml}" : "" 
    sk = (m=="align" && sk=="true") ? "--skiprc" : "" 
    
    PRIMER_FIELD = "${pf}"
    
    // all
    bf = (bf=="") ? "" : "--bf ${bf}"
    pf = (pf=="") ? "" : "--pf ${pf}"
    bc = (bc=="false") ? "" : "--barcode"
    rp = (rp=="false") ? "" : "--revpr"
    args_values.add("${m} --mode ${md} ${bf} ${pf} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk}")
    
    
}}

readArray = reads.toString().split(' ')
if(mate=="pair"){
	args_1 = args_values[0]
	args_2 = args_values[1]
	
  


	R1 = readArray[0]
	R2 = readArray[1]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
	
	
	"""
	
	MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	"""
}else{
	args_1 = args_values[0]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	
	R1 = readArray[0]

	"""
	echo -e "Assuming inputs for R1\n"
	
	MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	"""
}

}


process Mask_Primer_v_region_parse_log_MP {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "MP_log_table/$filename"}
input:
 val mate from g_4_0_g5_9
 set val(name), file(log_file) from g5_11_logFile21_g5_9

output:
 set val(name), file("*.tab")  into g5_9_logFile00_g5_12

script:
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ID PRIMER BARCODE ERROR
"""

}


process Mask_Primer_v_region_try_report_maskprimer {

input:
 set val(name), file(primers) from g5_9_logFile00_g5_12
 val matee from g_4_1_g5_12

output:
 file "*.rmd"  into g5_12_rMarkdown00_g5_16


shell:

if(matee=="pair"){
	readArray = primers.toString().split(' ')	
	primers_1 = readArray.grep(~/.*R1.*/)[0]
	primers_2 = readArray.grep(~/.*R2.*/)[0]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read 1", "Read 2")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom),",
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers_1}"))
	primer_log_2 <- loadLogTable(file.path(".", "!{primers_2}"))
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{

	readArray = primers.toString().split(' ')
	primers = readArray.grep(~/.*.tab*/)[0]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1],
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers}"))
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}
}


process Mask_Primer_v_region_render_rmarkdown {

input:
 file rmk from g5_12_rMarkdown00_g5_16

output:
 file "*.html"  into g5_16_outputFileHTML00

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Filter_Sequence_Quality_parse_log_FS {

input:
 set val(name), file(log_file) from g37_0_logFile10_g37_5
 val mate from g_4_1_g37_5

output:
 file "*table.tab"  into g37_5_logFile01_g37_7, g37_5_logFile01_g37_16

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray}  -f ID QUALITY
"""

}


process Filter_Sequence_Quality_report_filter_Seq_Quality {

input:
 val mate from g_4_0_g37_7
 file log_files from g37_5_logFile01_g37_7

output:
 file "*.rmd"  into g37_7_rMarkdown00_g37_16


shell:

if(mate=="pair"){
	readArray = log_files.toString().split(' ')	
	R1 = readArray[0]
	R2 = readArray[1]

	name = R1 - "_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- c("Read 1", "Read 2")
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("quality", 
	        paste("Mean Phred quality scores for",  plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The dotted line indicates the average quality score under which reads were removed."))
	```
	
	```{r, echo=FALSE}
	quality_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	quality_log_2 <- loadLogTable(file.path(".", "!{R2}"))
	```
	
	# Quality Scores
	
	Quality filtering is an essential step in most sequencing workflows. pRESTO’s
	FilterSeq tool remove reads with low mean Phred quality scores. 
	Phred quality scores are assigned to each nucleotide base call in automated 
	sequencer traces. The quality score (`Q`) of a base call is logarithmically 
	related to the probability that a base call is incorrect (`P`): 
	$Q = -10 log_{10} P$. For example, a base call with `Q=30` is incorrectly 
	assigned 1 in 1000 times. The most commonly used approach is to remove read 
	with average `Q` below 20.
	
	```{r, echo=FALSE}
	plotFilterSeq(quality_log_1, quality_log_2, titles=plot_titles, sizing="figure")
	```
	
	`r figures("quality")`
		
	EOF
	
	open OUT, ">FSQ_!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{

	readArray = log_files.toString().split(' ')
	R1 = readArray[0]
	name = R1 - "_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- c("Read")#params$quality_titles
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("quality", 
	        paste("Mean Phred quality scores for",  plot_titles[1],
	              "The dotted line indicates the average quality score under which reads were removed."))
	```
	
	```{r, echo=FALSE}
	quality_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	```
	
	# Quality Scores
	
	Quality filtering is an essential step in most sequencing workflows. pRESTO’s
	FilterSeq tool remove reads with low mean Phred quality scores. 
	Phred quality scores are assigned to each nucleotide base call in automated 
	sequencer traces. The quality score (`Q`) of a base call is logarithmically 
	related to the probability that a base call is incorrect (`P`): 
	$Q = -10 log_{10} P$. For example, a base call with `Q=30` is incorrectly 
	assigned 1 in 1000 times. The most commonly used approach is to remove read 
	with average `Q` below 20.
	
	```{r, echo=FALSE}
	plotFilterSeq(quality_log_1, titles=plot_titles[1], sizing="figure")
	```
	
	`r figures("quality")`
	
	EOF
	
	open OUT, ">FSQ_!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}
}


process Filter_Sequence_Quality_presto_render_rmarkdown {

input:
 file rmk from g37_7_rMarkdown00_g37_16
 file log_file from g37_5_logFile01_g37_16

output:
 file "*.html"  into g37_16_outputFileHTML00
 file "*csv" optional true  into g37_16_csvFile11

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Cluster_Sets_BARCODE_cluster_sets {

input:
 set val(name),file(reads) from g5_11_reads00_g35_14
 val mate from g_4_1_g35_14

output:
 set val(name),file("*_cluster-pass.fastq")  into g35_14_reads00_g32_14
 set val(name),file("*_cluster-fail.fastq") optional true  into g35_14_reads_failed11

script:
method = params.Cluster_Sets_BARCODE_cluster_sets.method
failed = params.Cluster_Sets_BARCODE_cluster_sets.failed
nproc = params.Cluster_Sets_BARCODE_cluster_sets.nproc
cluster_field = params.Cluster_Sets_BARCODE_cluster_sets.cluster_field
ident = params.Cluster_Sets_BARCODE_cluster_sets.ident
length = params.Cluster_Sets_BARCODE_cluster_sets.length
prefix = params.Cluster_Sets_BARCODE_cluster_sets.prefix
cluster_tool = params.Cluster_Sets_BARCODE_cluster_sets.cluster_tool
cluster_exec = params.Cluster_Sets_BARCODE_cluster_sets.cluster_exec
set_field = params.Cluster_Sets_BARCODE_cluster_sets.set_field
start = params.Cluster_Sets_BARCODE_cluster_sets.start
end = params.Cluster_Sets_BARCODE_cluster_sets.end
barcode_field = params.Cluster_Sets_BARCODE_cluster_sets.barcode_field
//* @style @condition:{method="set",set_field,start,end},{method="all",start,end},{method="barcode",barcode_field} @array:{method,failed,cluster_field,ident,length,prefix,cluster_tool,cluster_exec,set_field,start,end,barcode_field}  @multicolumn:{method,failed,nproc,cluster_field,ident,length,prefix,cluster_tool,cluster_exec},{set_field,start,end,barcode_field}

method = (method.size==2) ? method : [method[0],method[0]]
failed = (failed.size==2) ? failed : [failed[0],failed[0]]
cluster_field = (cluster_field.size==2) ? cluster_field : [cluster_field[0],cluster_field[0]]
ident = (ident.size==2) ? ident : [ident[0],ident[0]]
length = (length.size==2) ? length : [length[0],length[0]]
prefix = (prefix.size==2) ? prefix : [prefix[0],prefix[0]]
cluster_tool = (cluster_tool.size==2) ? cluster_tool : [cluster_tool[0],cluster_tool[0]]
cluster_exec = (cluster_exec.size==2) ? cluster_exec : [cluster_exec[0],cluster_exec[0]]
set_field = (set_field.size==2) ? set_field : [set_field[0],set_field[0]]
start = (start.size==2) ? start : [start[0],start[0]]
end = (end.size==2) ? end : [end[0],end[0]]
barcode_field = (barcode_field.size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]

def args_values = [];
[method, failed, cluster_field, ident, length, prefix, cluster_tool, cluster_exec, set_field, start, end, barcode_field].transpose().each { m, f, cf, i, l, p, ct, ce, sf, s, e, bf -> {
    f = (f=="true") ? "--failed" : ""
    p = (p=="") ? "" : "--prefix ${p}" 
    ce = (ce=="") ? "" : "--exec ${ce}" 
    sf = (m=="set") ? "-f ${sf}" : ""
    s = (m=="barcode") ? "" : "--start ${s}" 
    e = (m=="barcode") ? "" : (e=="") ? "" : "--end ${e}" 
    bf = (m=="barcode") ? "-f ${bf}" : ""
    args_values.add("${m} ${f} -k ${cf} --ident ${i} --length ${l} ${p} --cluster ${ct} ${ce} ${sf} ${s} ${e} ${bf}")
}}


if(mate=="pair"){
	// files
	readArray = reads.toString().split(' ')	
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	args_1 = args_values[0]
	args_2 = args_values[1]
	
	"""
	ClusterSets.py ${args_1} -s $R1  --nproc ${nproc}
	ClusterSets.py ${args_2} -s $R2  --nproc ${nproc}
	"""
}else{
	args_1 = args_values[0]
	"""
	ClusterSets.py ${args_1} -s $reads --nproc ${nproc}
	"""
}


}


process Cluster_Sets_VDJ_cluster_sets {

input:
 set val(name),file(reads) from g35_14_reads00_g32_14
 val mate from g_4_1_g32_14

output:
 set val(name),file("*_cluster-pass.fastq")  into g32_14_reads00_g33_15
 set val(name),file("*_cluster-fail.fastq") optional true  into g32_14_reads_failed11

script:
method = params.Cluster_Sets_VDJ_cluster_sets.method
failed = params.Cluster_Sets_VDJ_cluster_sets.failed
nproc = params.Cluster_Sets_VDJ_cluster_sets.nproc
cluster_field = params.Cluster_Sets_VDJ_cluster_sets.cluster_field
ident = params.Cluster_Sets_VDJ_cluster_sets.ident
length = params.Cluster_Sets_VDJ_cluster_sets.length
prefix = params.Cluster_Sets_VDJ_cluster_sets.prefix
cluster_tool = params.Cluster_Sets_VDJ_cluster_sets.cluster_tool
cluster_exec = params.Cluster_Sets_VDJ_cluster_sets.cluster_exec
set_field = params.Cluster_Sets_VDJ_cluster_sets.set_field
start = params.Cluster_Sets_VDJ_cluster_sets.start
end = params.Cluster_Sets_VDJ_cluster_sets.end
barcode_field = params.Cluster_Sets_VDJ_cluster_sets.barcode_field
//* @style @condition:{method="set",set_field,start,end},{method="all",start,end},{method="barcode",barcode_field} @array:{method,failed,cluster_field,ident,length,prefix,cluster_tool,cluster_exec,set_field,start,end,barcode_field}  @multicolumn:{method,failed,nproc,cluster_field,ident,length,prefix,cluster_tool,cluster_exec},{set_field,start,end,barcode_field}

method = (method.size==2) ? method : [method[0],method[0]]
failed = (failed.size==2) ? failed : [failed[0],failed[0]]
cluster_field = (cluster_field.size==2) ? cluster_field : [cluster_field[0],cluster_field[0]]
ident = (ident.size==2) ? ident : [ident[0],ident[0]]
length = (length.size==2) ? length : [length[0],length[0]]
prefix = (prefix.size==2) ? prefix : [prefix[0],prefix[0]]
cluster_tool = (cluster_tool.size==2) ? cluster_tool : [cluster_tool[0],cluster_tool[0]]
cluster_exec = (cluster_exec.size==2) ? cluster_exec : [cluster_exec[0],cluster_exec[0]]
set_field = (set_field.size==2) ? set_field : [set_field[0],set_field[0]]
start = (start.size==2) ? start : [start[0],start[0]]
end = (end.size==2) ? end : [end[0],end[0]]
barcode_field = (barcode_field.size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]

def args_values = [];
[method, failed, cluster_field, ident, length, prefix, cluster_tool, cluster_exec, set_field, start, end, barcode_field].transpose().each { m, f, cf, i, l, p, ct, ce, sf, s, e, bf -> {
    f = (f=="true") ? "--failed" : ""
    p = (p=="") ? "" : "--prefix ${p}" 
    ce = (ce=="") ? "" : "--exec ${ce}" 
    sf = (m=="set") ? "-f ${sf}" : ""
    s = (m=="barcode") ? "" : "--start ${s}" 
    e = (m=="barcode") ? "" : (e=="") ? "" : "--end ${e}" 
    bf = (m=="barcode") ? "-f ${bf}" : ""
    args_values.add("${m} ${f} -k ${cf} --ident ${i} --length ${l} ${p} --cluster ${ct} ${ce} ${sf} ${s} ${e} ${bf}")
}}


if(mate=="pair"){
	// files
	readArray = reads.toString().split(' ')	
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	args_1 = args_values[0]
	args_2 = args_values[1]
	
	"""
	ClusterSets.py ${args_1} -s $R1  --nproc ${nproc}
	ClusterSets.py ${args_2} -s $R2  --nproc ${nproc}
	"""
}else{
	args_1 = args_values[0]
	"""
	ClusterSets.py ${args_1} -s $reads --nproc ${nproc}
	"""
}


}


process Parse_header_parse_headers {

input:
 set val(name), file(reads) from g32_14_reads00_g33_15

output:
 set val(name),file("*${out}")  into g33_15_reads00_g15_10

script:
method = params.Parse_header_parse_headers.method
act = params.Parse_header_parse_headers.act
args = params.Parse_header_parse_headers.args


if(method=="collapse" || method=="copy" || method=="rename" || method=="merge"){
	out="_reheader.fastq"
	act = (act=="none") ? "" : "--act ${act}"
	"""
	ParseHeaders.py  ${method} -s ${reads} ${args} ${act}
	"""
}else{
	if(method=="table"){
			out=".tab"
			"""
			ParseHeaders.py ${method} -s ${reads} ${args}
			"""	
	}else{
		out="_reheader.fastq"
		"""
		ParseHeaders.py ${method} -s ${reads} ${args}
		"""		
	}
}


}

boolean isCollectionOrArray_bc(object) {    
    [Collection, Object[]].any { it.isAssignableFrom(object.getClass()) }
}

def args_creator_bc(barcode_field, primer_field, act, copy_field, mincount, minqual, minfreq, maxerror, prcons, maxgap, maxdiv, dep){
	def args_values;
    if(isCollectionOrArray_bc(barcode_field) || isCollectionOrArray_bc(primer_field) || isCollectionOrArray_bc(copy_field) || isCollectionOrArray_bc(mincount) || isCollectionOrArray_bc(minqual) || isCollectionOrArray_bc(minfreq) || isCollectionOrArray_bc(maxerror) || isCollectionOrArray_bc(prcons) || isCollectionOrArray_bc(maxgap) || isCollectionOrArray_bc(maxdiv) || isCollectionOrArray_bc(dep)){
    	primer_field = (isCollectionOrArray_bc(primer_field)) ? primer_field : [primer_field,primer_field]
    	act = (isCollectionOrArray_bc(act)) ? act : [act,act]
    	copy_field = (isCollectionOrArray_bc(copy_field)) ? copy_field : [copy_field,copy_field]
    	mincount = (isCollectionOrArray_bc(mincount)) ? mincount : [mincount,mincount]
    	minqual = (isCollectionOrArray_bc(minqual)) ? minqual : [minqual,minqual]
    	minfreq = (isCollectionOrArray_bc(minfreq)) ? minfreq : [minfreq,minfreq]
    	maxerror = (isCollectionOrArray_bc(maxerror)) ? maxerror : [maxerror,maxerror]
    	prcons = (isCollectionOrArray_bc(prcons)) ? prcons : [prcons,prcons]
    	maxgap = (isCollectionOrArray_bc(maxgap)) ? maxgap : [maxgap,maxgap]
    	maxdiv = (isCollectionOrArray_bc(maxdiv)) ? maxdiv : [maxdiv,maxdiv]
    	dep = (isCollectionOrArray_bc(dep)) ? dep : [dep,dep]
    	args_values = []
        [barcode_field,primer_field,act,copy_field,mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep].transpose().each { bf,pf,a,cf,mc,mq,mf,mr,pc,mg,md,d -> {
            bf = (bf=="") ? "" : "--bf ${bf}"
            pf = (pf=="") ? "" : "--pf ${pf}" 
            a = (a=="none") ? "" : "--act ${a}" 
            cf = (cf=="") ? "" : "--cf ${cf}" 
            mr = (mr=="none") ? "" : "--maxerror ${mr}" 
            pc = (pc=="none") ? "" : "--prcons ${pc}" 
            mg = (mg=="none") ? "" : "--maxgap ${mg}" 
            md = (md=="none") ? "" : "--maxdiv ${md}" 
            d = (d=="true") ? "--dep" : "" 
            args_values.add("${bf} ${pf} ${a} ${cf} -n ${mc} -q ${mq} --freq ${mf} ${mr} ${pc} ${mg} ${md} ${d}")
        }}
    }else{
        barcode_field = (barcode_field=="") ? "" : "--bf ${barcode_field}"
        primer_field = (primer_field=="") ? "" : "--pf ${primer_field}" 
        act = (act=="none") ? "" : "--act ${act}" 
        copy_field = (copy_field=="") ? "" : "--cf ${copy_field}" 
        maxerror = (maxerror=="none") ? "" : "--maxerror ${maxerror}" 
        prcons = (prcons=="none") ? "" : "--prcons ${prcons}" 
        maxgap = (maxgap=="none") ? "" : "--maxgap ${maxgap}" 
        maxdiv = (maxdiv=="none") ? "" : "--maxdiv ${maxdiv}" 
        dep = (dep=="true") ? "--dep" : "" 
        args_values = "${barcode_field} ${primer_field} ${act} ${copy_field} -n ${mincount} -q ${minqual} --freq ${minfreq} ${maxerror} ${prcons} ${maxgap} ${maxdiv} ${dep}"
    }
    return args_values
}


process Build_Consensus_build_consensus {

input:
 set val(name),file(reads) from g33_15_reads00_g15_10
 val mate from g_4_1_g15_10

output:
 set val(name),file("*_consensus-pass.fastq")  into g15_10_reads00_g21_16
 set val(name),file("BC*")  into g15_10_logFile10_g15_12
 set val(name),file("out*")  into g15_10_logFile22

script:
failed = params.Build_Consensus_build_consensus.failed
nproc = params.Build_Consensus_build_consensus.nproc
barcode_field = params.Build_Consensus_build_consensus.barcode_field
primer_field = params.Build_Consensus_build_consensus.primer_field
act = params.Build_Consensus_build_consensus.act
copy_field = params.Build_Consensus_build_consensus.copy_field
mincount = params.Build_Consensus_build_consensus.mincount
minqual = params.Build_Consensus_build_consensus.minqual
minfreq = params.Build_Consensus_build_consensus.minfreq
maxerror = params.Build_Consensus_build_consensus.maxerror
prcons = params.Build_Consensus_build_consensus.prcons
maxgap = params.Build_Consensus_build_consensus.maxgap
maxdiv = params.Build_Consensus_build_consensus.maxdiv
dep = params.Build_Consensus_build_consensus.dep
//* @style @condition:{act="none",},{act="min",copy_field},{act="max",copy_field},{act="sum",copy_field},{act="set",copy_field},{act="majority",copy_field} @array:{barcode_field,primer_field,act,copy_field,mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep} @multicolumn:{failed,nproc},{barcode_field,primer_field,act,copy_field}, {mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep}

args_values_bc = args_creator_bc(barcode_field, primer_field, act, copy_field, mincount, minqual, minfreq, maxerror, prcons, maxgap, maxdiv, dep)

// args 
if(isCollectionOrArray_bc(args_values_bc)){
	args_1 = args_values_bc[0]
	args_2 = args_values_bc[1]
}else{
	args_1 = args_values_bc
	args_2 = args_values_bc
}

failed = (failed=="true") ? "--failed" : "" 


if(mate=="pair"){
	// files
	readArray = reads.toString().split(' ')	
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	"""
	BuildConsensus.py -s $R1 ${args_1} --log BC_${name}_R1.log ${failed} --nproc ${nproc} >> out_${R1}_BC.log
	BuildConsensus.py -s $R2 ${args_2} --log BC_${name}_R2.log ${failed} --nproc ${nproc} >> out_${R1}_BC.log
	"""
}else{
	"""
	BuildConsensus.py -s $reads ${args_1} --outname ${name} --log BC_${name}.log ${failed} --nproc ${nproc} >> out_${R1}_BC.log
	"""
}


}


process Build_Consensus_parse_log_BC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "BC_log_table/$filename"}
input:
 set val(name),file(log_file) from g15_10_logFile10_g15_12
 val mate from g_4_1_g15_12

output:
 set val(name),file("*.tab")  into g15_12_logFile00_g15_14

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray} -f BARCODE SEQCOUNT CONSCOUNT PRCONS PRFREQ ERROR
"""

}


process Build_Consensus_report_Build_Consensus {

input:
 set val(name),file(log_files) from g15_12_logFile00_g15_14
 val matee from g_4_1_g15_14

output:
 file "*.rmd"  into g15_14_rMarkdown00_g15_16




shell:

if(matee=="pair"){
	readArray = log_files.toString().split(' ')	
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- c("Read 1", "Read 2")
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("cons_size", 
	        paste("Histogram of UMI read group sizes (reads per UMI) for",  
	              plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The x-axis indicates the number of reads in a UMI group and the y-axis is the 
	               number of UMI groups with that size. The Consensus and Total bars are overlayed 
	               (not stacked) histograms indicating whether the distribution has been calculated 
	               using the total number of reads (Total) or only those reads used for consensus 
	               generation (Consensus)."))
	figures("cons_prfreq", 
	        paste("Histograms showing the distribution of majority primer frequency for all UMI read groups for",
	              plot_titles[1], "(top) and", plot_titles[2], "(bottom)."))
	figures("cons_prsize", 
	        paste("Violin plots showing the distribution of UMI read group sizes by majority primer for",
	              plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "Only groups with majority primer frequency over the PRFREQ threshold set when running
	               BuildConsensus. Meaning, only retained UMI groups."))
	figures("cons_error", 
	        paste("Histogram showing the distribution of UMI read group error rates for",
	              plot_titles[1], "(top) and", plot_titles[2], "(bottom)."))
	figures("cons_prerror", 
	        paste("Violin plots showing the distribution of UMI read group error rates by majority primer for",
	              plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "Only groups with majority primer frequency over the PRFREQ threshold set when 
	               running BuildConsensus. Meaning, only retained UMI groups."))
	```
	
	```{r, echo=FALSE}
	consensus_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	consensus_log_2 <- loadLogTable(file.path(".", "!{R2}"))
	```
	
	# Generation of UMI Consensus Sequences
	
	Reads sharing the same UMI are collapsed into a single consensus sequence by
	the BuildConsensus tool. BuildConsensus considers several factors in determining
	the final consensus sequence, including the number of reads in a UMI group, 
	Phred quality scores (`Q`), primer annotations, and the number of mismatches 
	within a UMI group. Quality scores are used to resolve conflicting base calls in
	a UMI read group and the final consensus sequence is assigned consensus quality 
	scores derived from the individual base quality scores. The numbers of reads in a UMI
	group, number of matching primer annotations, and error rate (average base mismatches from 
	consensus) are used as strict cut-offs for exclusion of erroneous UMI read groups.
	Additionally, individual reads are excluded whose primer annotation differs from 
	the majority in cases where there are sufficient number of reads exceeding 
	the primer consensus cut-off.
	
	## Reads per UMI
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log_1, consensus_log_2, titles=plot_titles, 
	                   style="size", sizing="figure")
	```
	
	`r figures("cons_size")`
	
	## UMI read group primer frequencies
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log_1, consensus_log_2, titles=plot_titles, 
	                   style="prfreq", sizing="figure")
	```
	
	`r figures("cons_prfreq")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log_1, consensus_log_2, titles=plot_titles, 
	                   style="prsize", sizing="figure")
	```
	
	`r figures("cons_prsize")`
	
	## UMI read group error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log_1, consensus_log_2, titles=plot_titles, 
	                   style="error", sizing="figure")
	```
	
	`r figures("cons_error")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log_1, consensus_log_2, titles=plot_titles, 
	                   style="prerror", sizing="figure")
	```
	
	`r figures("cons_prerror")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{
	
	readArray = log_files.toString().split(' ')
	R1 = readArray[0]
	
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	
		
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("cons_size", "Histogram of UMI read group sizes (reads per UMI). 
	                      The x-axis indicates the number of reads 
	                      in a UMI group and the y-axis is the number of UMI groups 
	                      with that size. The Consensus and Total bars are overlayed
	                      (not stacked) histograms indicating whether the distribution
	                      has been calculated using the total number of reads (Total)
	                      or only those reads used for consensus generation (Consensus).")
	figures("cons_error", "Histogram showing the distribution of UMI read group error rates.")
	```
	
	```{r, echo=FALSE}
	consensus_log <- loadLogTable(file.path(".", "!{R1}"))
	```
	
	# Generation of UMI Consensus Sequences
	
	Reads sharing the same UMI are collapsed into a single consensus sequence by
	the BuildConsensus tool. BuildConsensus considers several factors in determining
	the final consensus sequence, including the number of reads in a UMI group, 
	Phred quality scores (`Q`), primer annotations, and the number of mismatches 
	within a UMI group. Quality scores are used to resolve conflicting base calls in
	a UMI read group and the final consensus sequence is assigned consensus quality 
	scores derived from the individual base quality scores. The numbers of reads in a UMI
	group, number of matching primer annotations, and error rate (average base mismatches from 
	consensus) are used as strict cut-offs for exclusion of erroneous UMI read groups.
	Additionally, individual reads are excluded whose primer annotation differs from 
	the majority in cases where there are sufficient number of reads exceeding 
	the primer consensus cut-off.
	
	## Reads per UMI
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log, style="size", sizing="figure")
	```
	
	`r figures("cons_size")`
	
	## UMI read group error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log, style="error", sizing="figure")
	```
	
	`r figures("cons_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}

}


process Build_Consensus_render_rmarkdown {

input:
 file rmk from g15_14_rMarkdown00_g15_16

output:
 file "*.html"  into g15_16_outputFileHTML00

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process collapse_sequences_collapse_seq {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_collapse-unique.fast.*$/) "reads_unique/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_collapse-duplicate.fast.*$/) "reads_duplicates/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_collapse-undetermined.fast.*$/) "reads_undetermined/$filename"}
input:
 set val(name), file(reads) from g15_10_reads00_g21_16
 val mate from g_4_1_g21_16

output:
 set val(name),  file("*_collapse-unique.fast*")  into g21_16_reads00
 set val(name),  file("*_collapse-duplicate.fast*") optional true  into g21_16_reads_duplicate11
 set val(name),  file("*_collapse-undetermined.fast*") optional true  into g21_16_reads_undetermined22
 file "CS_*"  into g21_16_logFile33
 set val(name),  file("out*")  into g21_16_logFile44

script:
max_missing = params.collapse_sequences_collapse_seq.max_missing
inner = params.collapse_sequences_collapse_seq.inner
fasta = params.collapse_sequences_collapse_seq.fasta
act = params.collapse_sequences_collapse_seq.act
uf = params.collapse_sequences_collapse_seq.uf
cf = params.collapse_sequences_collapse_seq.cf
nproc = params.collapse_sequences_collapse_seq.nproc
failed = params.collapse_sequences_collapse_seq.failed

inner = (inner=="true") ? "--inner" : ""
fasta = (fasta=="true") ? "--fasta" : ""
act = (act=="none") ? "" : "--act ${act}"
cf = (cf=="") ? "" : "--cf ${cf}"
uf = (uf=="") ? "" : "--uf ${uf}"
failed = (failed=="false") ? "" : "--failed"

readArray = reads.toString().split(' ')	
if(mate=="pair"){
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
}else{
	R1 = readArray[0]
}


"""
CollapseSeq.py -s ${reads} -n ${max_missing} ${fasta} ${inner} ${uf} ${cf} ${act} --log CS_${name}.log ${failed} >> out_${R1}_collapse.log
"""

}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
