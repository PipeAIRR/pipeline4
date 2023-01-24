$HOSTNAME = ""
params.outdir = 'results'  


if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 
if (!params.mate_single){params.mate_single = ""} 
if (!params.V_primers){params.V_primers = ""} 
if (!params.C_primers){params.C_primers = ""} 
if (!params.S_primer){params.S_primer = ""} 
// Stage empty file to be used as an optional input where required
ch_empty_file_1 = file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)
ch_empty_file_2 = file("$baseDir/.emptyfiles/NO_FILE_2", hidden:true)

if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_0_reads_g2_12}
 } else {  
	g_0_reads_g2_12 = Channel.empty()
 }

Channel.value(params.mate).into{g_1_mate_g2_15;g_1_mate_g2_12}
Channel.value(params.mate_single).into{g_4_mate_g3_0;g_4_mate_g3_5;g_4_mate_g5_9;g_4_mate_g5_11;g_4_mate_g6_9;g_4_mate_g6_11;g_4_mate_g9_12;g_4_mate_g9_14;g_4_mate_g15_10;g_4_mate_g15_12;g_4_mate_g26_9;g_4_mate_g26_11;g_4_mate_g32_12;g_4_mate_g32_14}
g_7_R1_primers_g5_11 = params.V_primers && file(params.V_primers, type: 'any').exists() ? file(params.V_primers, type: 'any') : ch_empty_file_1
g_8_R1_primers_g26_11 = params.C_primers && file(params.C_primers, type: 'any').exists() ? file(params.C_primers, type: 'any') : ch_empty_file_1
g_27_R1_primers_g6_11 = params.S_primer && file(params.S_primer, type: 'any').exists() ? file(params.S_primer, type: 'any') : ch_empty_file_1


process Assemble_pairs_assemble_pairs {

input:
 set val(name),file(reads) from g_0_reads_g2_12
 val mate from g_1_mate_g2_12

output:
 set val(name),file("*_assemble-pass.f*")  into g2_12_reads0_g3_0
 file "AP_*"  into g2_12_logFile1_g2_15
 set val(name),file("*_assemble-fail.f*") optional true  into g2_12_reads_failed22

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

def ref_file;
if(params.assemble_ref != null) ref_file = (assemble_reference.startsWith('NO_FILE')) ? "" : "-r ${assemble_reference}"

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
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
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
	
	echo \$align_exec
	echo \$dbexec
	
	
	
	AssemblePairs.py ${method} -1 ${R1} -2 ${R2} ${coord} ${rc} ${head_fields_R1} ${head_fields_R2} ${args} \$align_exec \$dbexec ${fasta} ${failed} --log AP_${name}.log ${nproc}
	"""

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}

}


process Assemble_pairs_parse_log_AP {

input:
 file log_file from g2_12_logFile1_g2_15
 val mate from g_1_mate_g2_15

output:
 file "*.tab"  into g2_15_logFile00

script:
field_to_parse = params.Assemble_pairs_parse_log_AP.field_to_parse
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ${field_to_parse}
"""


}


process Filter_Sequence_Quality_filter_seq_quality {

input:
 set val(name),file(reads) from g2_12_reads0_g3_0
 val mate from g_4_mate_g3_0

output:
 set val(name), file("*_${method}-pass.fastq")  into g3_0_reads0_g6_11
 file "FS_*"  into g3_0_logFile1_g3_5

script:
method = params.Filter_Sequence_Quality_filter_seq_quality.method
nproc = params.Filter_Sequence_Quality_filter_seq_quality.nproc
q = params.Filter_Sequence_Quality_filter_seq_quality.q
n_length = params.Filter_Sequence_Quality_filter_seq_quality.n_length
n_missing = params.Filter_Sequence_Quality_filter_seq_quality.n_missing
//* @style @condition:{method="quality",q}, {method="length",n_length}, {method="missing",n_missing} @multicolumn:{method,nproc}

if(method=="quality"){
	q = "-q ${q}"
	n_length = ""
	n_missing = ""
}else{
	if(method=="length"){
		q = ""
		n_length = "-n ${n_length}"
		n_missing = ""
	}else{
		q = ""
		n_length = ""
		n_missing = "-n ${n_missing}"
	}
}

readArray = reads.toString().split(' ')	

if(mate=="pair"){
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --outname ${name}_R1 --nproc ${nproc} --log FS_R1_${name}.log
	FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} --outname ${name}_R2 --nproc ${nproc} --log FS_R2_${name}.log
	"""
}else{
	R1 = readArray[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --outname ${name} --nproc ${nproc} --log FS_${name}.log
	"""
}


}


process Filter_Sequence_Quality_parse_log_FS {

input:
 file log_file from g3_0_logFile1_g3_5
 val mate from g_4_mate_g3_5

output:
 file "*.tab"  into g3_5_logFile00

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray}  -f ID QUALITY
"""

}



process Mask_Primer_CNU_Barcode_primer_S_MaskPrimers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_primers-fail.fastq$/) "MP_fail_reads_s_primer/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /MP_.*$/) "MP_log_C_table/$filename"}
input:
 val mate from g_4_mate_g6_11
 set val(name),file(reads) from g3_0_reads0_g6_11
 file R1_primers from g_27_R1_primers_g6_11

output:
 set val(name), file("*_primers-pass.fastq")  into g6_11_reads0_g26_11
 set val(name), file("*_primers-fail.fastq") optional true  into g6_11_reads_failed11
 file "MP_*"  into g6_11_logFile2_g6_9

script:
method = params.Mask_Primer_CNU_Barcode_primer_S_MaskPrimers.method
barcode_field = params.Mask_Primer_CNU_Barcode_primer_S_MaskPrimers.barcode_field
primer_field = params.Mask_Primer_CNU_Barcode_primer_S_MaskPrimers.primer_field
barcode = params.Mask_Primer_CNU_Barcode_primer_S_MaskPrimers.barcode
revpr = params.Mask_Primer_CNU_Barcode_primer_S_MaskPrimers.revpr
mode = params.Mask_Primer_CNU_Barcode_primer_S_MaskPrimers.mode
failed = params.Mask_Primer_CNU_Barcode_primer_S_MaskPrimers.failed
nproc = params.Mask_Primer_CNU_Barcode_primer_S_MaskPrimers.nproc
maxerror = params.Mask_Primer_CNU_Barcode_primer_S_MaskPrimers.maxerror
umi_length = params.Mask_Primer_CNU_Barcode_primer_S_MaskPrimers.umi_length
start = params.Mask_Primer_CNU_Barcode_primer_S_MaskPrimers.start
extract_length = params.Mask_Primer_CNU_Barcode_primer_S_MaskPrimers.extract_length
maxlen = params.Mask_Primer_CNU_Barcode_primer_S_MaskPrimers.maxlen
skiprc = params.Mask_Primer_CNU_Barcode_primer_S_MaskPrimers.skiprc
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
    
    // all
    bf = (bf=="") ? "" : "--bf ${bf}"
    pf = (pf=="") ? "" : "--pf ${pf}"
    bc = (bc=="false") ? "" : "--barcode"
    rp = (rp=="false") ? "" : "--revpr"
    args_values.add("${m} --mode ${md} ${bf} ${pf} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk}")
}}

if(mate=="pair"){
	args_1 = args_values[0]
	args_2 = args_values[1]
	
    readArray = reads.toString().split(' ')	
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
	
	"""
	MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed}
	MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed}
	"""
}else{
	args_1 = args_values[0]
	println args_1
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	
	"""
	echo -e "Assuming inputs for R1\n"
	
	MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed}
	"""
}

}


process Mask_Primer_CNU_Barcode_primer_S_parse_log_MP {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "MP_log_C/$filename"}
input:
 val mate from g_4_mate_g6_9
 file log_file from g6_11_logFile2_g6_9

output:
 file "*.tab"  into g6_9_logFile00

script:
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ID PRIMER BARCODE ERROR
"""

}



process Mask_Primer_C_primer_MaskPrimers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_primers-fail.fastq$/) "MP_fail_reads_c_primer/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /MP_.*$/) "MP_log_C_table/$filename"}
input:
 val mate from g_4_mate_g26_11
 set val(name),file(reads) from g6_11_reads0_g26_11
 file R1_primers from g_8_R1_primers_g26_11

output:
 set val(name), file("*_primers-pass.fastq")  into g26_11_reads0_g5_11
 set val(name), file("*_primers-fail.fastq") optional true  into g26_11_reads_failed11
 file "MP_*"  into g26_11_logFile2_g26_9

script:
method = params.Mask_Primer_C_primer_MaskPrimers.method
barcode_field = params.Mask_Primer_C_primer_MaskPrimers.barcode_field
primer_field = params.Mask_Primer_C_primer_MaskPrimers.primer_field
barcode = params.Mask_Primer_C_primer_MaskPrimers.barcode
revpr = params.Mask_Primer_C_primer_MaskPrimers.revpr
mode = params.Mask_Primer_C_primer_MaskPrimers.mode
failed = params.Mask_Primer_C_primer_MaskPrimers.failed
nproc = params.Mask_Primer_C_primer_MaskPrimers.nproc
maxerror = params.Mask_Primer_C_primer_MaskPrimers.maxerror
umi_length = params.Mask_Primer_C_primer_MaskPrimers.umi_length
start = params.Mask_Primer_C_primer_MaskPrimers.start
extract_length = params.Mask_Primer_C_primer_MaskPrimers.extract_length
maxlen = params.Mask_Primer_C_primer_MaskPrimers.maxlen
skiprc = params.Mask_Primer_C_primer_MaskPrimers.skiprc
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
    
    // all
    bf = (bf=="") ? "" : "--bf ${bf}"
    pf = (pf=="") ? "" : "--pf ${pf}"
    bc = (bc=="false") ? "" : "--barcode"
    rp = (rp=="false") ? "" : "--revpr"
    args_values.add("${m} --mode ${md} ${bf} ${pf} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk}")
}}

if(mate=="pair"){
	args_1 = args_values[0]
	args_2 = args_values[1]
	
    readArray = reads.toString().split(' ')	
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
	
	"""
	MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed}
	MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed}
	"""
}else{
	args_1 = args_values[0]
	println args_1
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	
	"""
	echo -e "Assuming inputs for R1\n"
	
	MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed}
	"""
}

}



process Mask_Primer_V_primer_MaskPrimers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /MP_.*$/) "MP_log_table/$filename"}
input:
 val mate from g_4_mate_g5_11
 set val(name),file(reads) from g26_11_reads0_g5_11
 file R1_primers from g_7_R1_primers_g5_11

output:
 set val(name), file("*_primers-pass.fastq")  into g5_11_reads0_g9_14
 set val(name), file("*_primers-fail.fastq") optional true  into g5_11_reads_failed11
 file "MP_*"  into g5_11_logFile2_g5_9

script:
method = params.Mask_Primer_V_primer_MaskPrimers.method
barcode_field = params.Mask_Primer_V_primer_MaskPrimers.barcode_field
primer_field = params.Mask_Primer_V_primer_MaskPrimers.primer_field
barcode = params.Mask_Primer_V_primer_MaskPrimers.barcode
revpr = params.Mask_Primer_V_primer_MaskPrimers.revpr
mode = params.Mask_Primer_V_primer_MaskPrimers.mode
failed = params.Mask_Primer_V_primer_MaskPrimers.failed
nproc = params.Mask_Primer_V_primer_MaskPrimers.nproc
maxerror = params.Mask_Primer_V_primer_MaskPrimers.maxerror
umi_length = params.Mask_Primer_V_primer_MaskPrimers.umi_length
start = params.Mask_Primer_V_primer_MaskPrimers.start
extract_length = params.Mask_Primer_V_primer_MaskPrimers.extract_length
maxlen = params.Mask_Primer_V_primer_MaskPrimers.maxlen
skiprc = params.Mask_Primer_V_primer_MaskPrimers.skiprc
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
    
    // all
    bf = (bf=="") ? "" : "--bf ${bf}"
    pf = (pf=="") ? "" : "--pf ${pf}"
    bc = (bc=="false") ? "" : "--barcode"
    rp = (rp=="false") ? "" : "--revpr"
    args_values.add("${m} --mode ${md} ${bf} ${pf} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk}")
}}

if(mate=="pair"){
	args_1 = args_values[0]
	args_2 = args_values[1]
	
    readArray = reads.toString().split(' ')	
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
	
	"""
	MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed}
	MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed}
	"""
}else{
	args_1 = args_values[0]
	println args_1
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	
	"""
	echo -e "Assuming inputs for R1\n"
	
	MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed}
	"""
}

}


process Mask_Primer_V_primer_parse_log_MP {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "MP_log/$filename"}
input:
 val mate from g_4_mate_g5_9
 file log_file from g5_11_logFile2_g5_9

output:
 file "*.tab"  into g5_9_logFile00

script:
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ID PRIMER BARCODE ERROR
"""

}


process Cluster_Sets_CPRIMER_cluster_sets {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /CS_.*$/) "CS_log_table/$filename"}
input:
 set val(name),file(reads) from g5_11_reads0_g9_14
 val mate from g_4_mate_g9_14

output:
 set val(name),file("*_cluster-pass.fastq")  into g9_14_reads0_g25_15
 file "CS_*"  into g9_14_logFile1_g9_12
 set val(name),file("*_cluster-fail.fastq") optional true  into g9_14_reads_failed22

script:
method = params.Cluster_Sets_CPRIMER_cluster_sets.method
failed = params.Cluster_Sets_CPRIMER_cluster_sets.failed
nproc = params.Cluster_Sets_CPRIMER_cluster_sets.nproc
cluster_field = params.Cluster_Sets_CPRIMER_cluster_sets.cluster_field
ident = params.Cluster_Sets_CPRIMER_cluster_sets.ident
length = params.Cluster_Sets_CPRIMER_cluster_sets.length
prefix = params.Cluster_Sets_CPRIMER_cluster_sets.prefix
cluster_tool = params.Cluster_Sets_CPRIMER_cluster_sets.cluster_tool
cluster_exec = params.Cluster_Sets_CPRIMER_cluster_sets.cluster_exec
set_field = params.Cluster_Sets_CPRIMER_cluster_sets.set_field
start = params.Cluster_Sets_CPRIMER_cluster_sets.start
end = params.Cluster_Sets_CPRIMER_cluster_sets.end
barcode_field = params.Cluster_Sets_CPRIMER_cluster_sets.barcode_field
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
	ClusterSets.py ${args_1} -s $R1  --log CS_${name}_R1.log --nproc ${nproc}
	ClusterSets.py ${args_2} -s $R2  --log CS_${name}_R2.log --nproc ${nproc}
	"""
}else{
	args_1 = args_values[0]
	"""
	ClusterSets.py ${args_1} -s $reads --log CS_${name}.log --nproc ${nproc}
	"""
}


}


process Parse_header_CPRIMER_parse_headers {

input:
 set val(name), file(reads) from g9_14_reads0_g25_15

output:
 set val(name),file("*${out}")  into g25_15_reads0_g32_14

script:
method = params.Parse_header_CPRIMER_parse_headers.method
act = params.Parse_header_CPRIMER_parse_headers.act
args = params.Parse_header_CPRIMER_parse_headers.args

out="_reheader.fastq"
if(method=="collapse" || method=="add" || method=="copy" || method=="rename" || method=="merge"){
	"""
	ParseHeaders.py  ${method} -s ${reads} ${args} --act ${act}
	"""
}else{
	if(method=="table"){
			out=".tab"
			"""
			ParseHeaders.py ${method} -s ${reads} -o ${name}.tab ${args}
			"""	
	}else{
		"""
		ParseHeaders.py ${method} -s ${reads} ${args}
		"""		
	}
}



}


process Cluster_Sets_CPRIMER_parse_log_BC_copy {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "CS_log/$filename"}
input:
 file log_file from g9_14_logFile1_g9_12
 val mate from g_4_mate_g9_12

output:
 file "*.tab"  into g9_12_logFile00

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray}  -f BARCODE SEQCOUNT CONSCOUNT PRCONS PRFREQ ERROR
"""

}


process Mask_Primer_C_primer_parse_log_MP {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "MP_log_C/$filename"}
input:
 val mate from g_4_mate_g26_9
 file log_file from g26_11_logFile2_g26_9

output:
 file "*.tab"  into g26_9_logFile00

script:
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ID PRIMER BARCODE ERROR
"""

}


process Cluster_Sets_BARCODE_cluster_sets {

input:
 set val(name),file(reads) from g25_15_reads0_g32_14
 val mate from g_4_mate_g32_14

output:
 set val(name),file("*_cluster-pass.fastq")  into g32_14_reads0_g33_15
 file "CS_*"  into g32_14_logFile1_g32_12
 set val(name),file("*_cluster-fail.fastq") optional true  into g32_14_reads_failed22

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
	ClusterSets.py ${args_1} -s $R1  --log CS_${name}_R1.log --nproc ${nproc}
	ClusterSets.py ${args_2} -s $R2  --log CS_${name}_R2.log --nproc ${nproc}
	"""
}else{
	args_1 = args_values[0]
	"""
	ClusterSets.py ${args_1} -s $reads --log CS_${name}.log --nproc ${nproc}
	"""
}


}


process Parse_header_BARCODE_parse_headers {

input:
 set val(name), file(reads) from g32_14_reads0_g33_15

output:
 set val(name),file("*${out}")  into g33_15_reads0_g15_10

script:
method = params.Parse_header_BARCODE_parse_headers.method
act = params.Parse_header_BARCODE_parse_headers.act
args = params.Parse_header_BARCODE_parse_headers.args

out="_reheader.fastq"
if(method=="collapse" || method=="add" || method=="copy" || method=="rename" || method=="merge"){
	"""
	ParseHeaders.py  ${method} -s ${reads} ${args} --act ${act}
	"""
}else{
	if(method=="table"){
			out=".tab"
			"""
			ParseHeaders.py ${method} -s ${reads} -o ${name}.tab ${args}
			"""	
	}else{
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

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /BC_.*$/) "BC_log/$filename"}
input:
 set val(name),file(reads) from g33_15_reads0_g15_10
 val mate from g_4_mate_g15_10

output:
 set val(name),file("*_consensus-pass.fastq")  into g15_10_reads0_g21_16
 file "BC_*"  into g15_10_logFile1_g15_12

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
	BuildConsensus.py -s $R1 ${args_1} --outname ${name}_R1 --log BC_${name}_R1.log ${failed} --nproc ${nproc}
	BuildConsensus.py -s $R2 ${args_2} --outname ${name}_R2 --log BC_${name}_R2.log ${failed} --nproc ${nproc}
	"""
}else{
	"""
	BuildConsensus.py -s $reads ${args_1} --outname ${name} --log BC_${name}.log ${failed} --nproc ${nproc}
	"""
}


}


process collapse_sequences_collapse_seq {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_collapse-unique.fast.*$/) "reads_unique/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_collapse-duplicate.fast.*$/) "reads_duplicates/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_collapse-undetermined.fast.*$/) "reads_undetermined/$filename"}
input:
 set val(name), file(reads) from g15_10_reads0_g21_16

output:
 set val(name),  file("*_collapse-unique.fast*")  into g21_16_reads00
 set val(name),  file("*_collapse-duplicate.fast*") optional true  into g21_16_reads_duplicate11
 set val(name),  file("*_collapse-undetermined.fast*") optional true  into g21_16_reads_undetermined22

script:
max_missing = params.collapse_sequences_collapse_seq.max_missing
inner = params.collapse_sequences_collapse_seq.inner
fasta = params.collapse_sequences_collapse_seq.fasta
act = params.collapse_sequences_collapse_seq.act
uf = params.collapse_sequences_collapse_seq.uf
cf = params.collapse_sequences_collapse_seq.cf
nproc = params.collapse_sequences_collapse_seq.nproc

inner = (inner=="true") ? "--inner" : ""
fasta = (fasta=="true") ? "--fasta" : ""
act = (act=="none") ? "" : "--act ${act}"
cf = (cf=="") ? "" : "--cf ${cf}"
uf = (uf=="") ? "" : "--uf ${uf}"

"""
CollapseSeq.py -s ${reads} -n ${max_missing} ${fasta} ${inner} ${uf} ${cf} ${act}
"""

}


process Build_Consensus_parse_log_BC_copy {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "BC_log_table/$filename"}
input:
 file log_file from g15_10_logFile1_g15_12
 val mate from g_4_mate_g15_12

output:
 file "*.tab"  into g15_12_logFile00

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray}  -f BARCODE SEQCOUNT CONSCOUNT PRCONS PRFREQ ERROR
"""

}


process Cluster_Sets_BARCODE_parse_log_BC_copy {

input:
 file log_file from g32_14_logFile1_g32_12
 val mate from g_4_mate_g32_12

output:
 file "*.tab"  into g32_12_logFile00

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray}  -f BARCODE SEQCOUNT CONSCOUNT PRCONS PRFREQ ERROR
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
