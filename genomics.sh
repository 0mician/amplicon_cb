#!/bin/bash

readonly FASTQ=$1   # path to raw fastq file
readonly OUTDIR=$2  # where the analysis will be saved
readonly BARCODE=$3 # barcode used (eg, "barcode01")
readonly FAST5=$4   # location of basecalled fast5 
readonly ID=$(basename "${FASTQ}" .fastq)
readonly FASTQDIR=$(dirname "${FASTQ}")

##########
# Assembly
##########

# clean reads of adapters
porechop -i "${FASTQ}" \
	 -o "${FASTQDIR}"/"${ID}".adaptrim.fastq

# removing reads that are too long (causes problem with assembly
# because of repeats). In this analysis, the amplicons are around 4
# kb, hence trimming everything above 5 kb. Note: very few long
# reads > 5kb. Probably arised during ligation with native barcodes
readonly awkcommand='
BEGIN {OFS = "\n"} 
{
    header = $0 ;
    getline seq ;
    getline qheader ;
    getline qseq ;
    
    if (length(seq) <= 5000)
    {
        print header, seq, qheader, qseq
    }
}'
awk "${awkcommand}" < "${FASTQDIR}"/"${ID}".adaptrim.fastq > "${FASTQDIR}"/"${ID}".adalentrim.fastq

# assembly with canu, the contigFilter was adjusted to allow
# single reads to span the entire assembly (short amplicon,
# so definitely a possibility). 
canu -p "${ID}" \
     -d "${HOME}"/"${ID}"_canu_asm_adalentrim \
     contigFilter="2 1000 1.0 1.0 2" \
     genomeSize=4.3k \
     -nanopore-raw "${FASTQDIR}"/"${ID}".adalentrim.fastq

################################################################
# verification of assembly by mapping reads back to draft genome
################################################################

readonly MAPDIR="${OUTDIR}"/"${ID}"_read_mapping
mkdir -p "${MAPDIR}"
mv "${FASTQDIR}"/"${ID}".adalentrim.fastq "${MAPDIR}"
mv "${HOME}"/"${ID}"_canu_asm_adalentrim "${OUTDIR}"
cp "${OUTDIR}"/"${ID}"_canu_asm_adalentrim/"${ID}".contigs.fasta "${MAPDIR}"

# indexing assembly
readonly REFERENCE="${MAPDIR}"/"${ID}".contigs.fasta
readonly PICARD="${HOME}"/Tools/picard.jar

bwa index -a bwtsw "${REFERENCE}"
samtools faidx "${REFERENCE}"
java -jar "${PICARD}" CreateSequenceDictionary \
     REFERENCE="${REFERENCE}" \
     OUTPUT="${REFERENCE}".dict

# mapping ont reads to assembly
bwa mem -x ont2d "${REFERENCE}" "${MAPDIR}"/"${ID}".adalentrim.fastq > "${MAPDIR}"/aln.sam
samtools sort "${MAPDIR}"/aln.sam -o "${MAPDIR}"/aln.sorted.bam
samtools index "${MAPDIR}"/aln.sorted.bam
samtools flagstat "${MAPDIR}"/aln.sorted.bam > "${MAPDIR}"/aln.sorted.stats

###################################################################
# polishing the assembly (obtention of a better consensus sequence)
###################################################################

mkdir "${OUTDIR}"/"${ID}"_polish
readonly POLISH_DIR="${OUTDIR}"/"${ID}"_polish

# extraction of fasta sequence from fast5
./nanopolish extract --type template "${FAST5}" > "${POLISH_DIR}"/pass.fasta

# mapping of reads against draft
bwa mem -x ont2d -t 8 "${REFERENCE}" \
    "${POLISH_DIR}"/pass.fasta > "${POLISH_DIR}"/aln.sam

samtools sort "${POLISH_DIR}"/aln.sam \
	 -o "${POLISH_DIR}"/aln.sort.bam

samtools index "${POLISH_DIR}"/aln.sort.bam

# nanopolish
python nanopolish_makerange.py --segment-length 1000 "${REFERENCE}" | \
    parallel --results "${POLISH_DIR}"/nanopolish.results -P 8 \
	     nanopolish variants --consensus "${POLISH_DIR}"/polished.{1}.fa -w {1} \
	     -r "${POLISH_DIR}"/"${ID}".adalentrim.fasta \
	     -b "${MAPDIR}"/aln.sorted.bam \
	     -g "${REFERENCE}" \
	     -t 4 --min-candidate-frequency 0.1

python nanopolish_merge.py "${POLISH_DIR}"/polished.*.fa > "${POLISH_DIR}"/am1_polished.fasta

######################
# MSA for verification
######################
cat "${REFERENCE}" "${POLISH_DIR}"/am1_polished.fasta > "${POLISH_DIR}"/polished_vs_non_polished.fasta
cat "${POLISH_DIR}"/am1_polished.fasta "${POLISH_DIR}"/sanger.fasta > "${POLISH_DIR}"/polished_vs_sanger.fasta

clustalo --infile= "${POLISH_DIR}"/polished_vs_non_polished.fasta \
	 --outfmt=clu --outfile="${POLISH_DIR}"/msa_polished_vs_non_polished.clustal --verbose

clustalo --infile= "${POLISH_DIR}"/polished_vs_sanger.fasta \
	 --outfmt=clu --outfile="${POLISH_DIR}"/msa_polished_vs_sanger.clustal --verbose
