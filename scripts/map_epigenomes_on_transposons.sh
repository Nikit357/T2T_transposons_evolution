#!/bin/bash
set -e

# Define some colors for logging
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${GREEN}Starting epigenome mapping process...${NC}\n"
# --- Processing file 1/39 ---
echo -e "${BLUE}Processing: 22Rv1.CTCF.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph 22Rv1.CTCF.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > 22Rv1.CTCF.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b 22Rv1.CTCF.chm13v2.0.bedGraph -c 4 -o mean > 22Rv1.CTCF.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm 22Rv1.CTCF.chm13v2.0.bedGraph
rm 22Rv1.CTCF.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped 22Rv1.CTCF.chm13v2.0.bw to 22Rv1.CTCF.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 2/39 ---
echo -e "${BLUE}Processing: 22Rv1.H3K27ac.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph 22Rv1.H3K27ac.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > 22Rv1.H3K27ac.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b 22Rv1.H3K27ac.chm13v2.0.bedGraph -c 4 -o mean > 22Rv1.H3K27ac.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm 22Rv1.H3K27ac.chm13v2.0.bedGraph
rm 22Rv1.H3K27ac.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped 22Rv1.H3K27ac.chm13v2.0.bw to 22Rv1.H3K27ac.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 3/39 ---
echo -e "${BLUE}Processing: BE2C.H3K4me1.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph BE2C.H3K4me1.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > BE2C.H3K4me1.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b BE2C.H3K4me1.chm13v2.0.bedGraph -c 4 -o mean > BE2C.H3K4me1.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm BE2C.H3K4me1.chm13v2.0.bedGraph
rm BE2C.H3K4me1.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped BE2C.H3K4me1.chm13v2.0.bw to BE2C.H3K4me1.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 4/39 ---
echo -e "${BLUE}Processing: BE2C.H3K9me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph BE2C.H3K9me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > BE2C.H3K9me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b BE2C.H3K9me3.chm13v2.0.bedGraph -c 4 -o mean > BE2C.H3K9me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm BE2C.H3K9me3.chm13v2.0.bedGraph
rm BE2C.H3K9me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped BE2C.H3K9me3.chm13v2.0.bw to BE2C.H3K9me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 5/39 ---
echo -e "${BLUE}Processing: BE2C.H3K27me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph BE2C.H3K27me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > BE2C.H3K27me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b BE2C.H3K27me3.chm13v2.0.bedGraph -c 4 -o mean > BE2C.H3K27me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm BE2C.H3K27me3.chm13v2.0.bedGraph
rm BE2C.H3K27me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped BE2C.H3K27me3.chm13v2.0.bw to BE2C.H3K27me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 6/39 ---
echo -e "${BLUE}Processing: BE2C.H3K36me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph BE2C.H3K36me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > BE2C.H3K36me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b BE2C.H3K36me3.chm13v2.0.bedGraph -c 4 -o mean > BE2C.H3K36me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm BE2C.H3K36me3.chm13v2.0.bedGraph
rm BE2C.H3K36me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped BE2C.H3K36me3.chm13v2.0.bw to BE2C.H3K36me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 7/39 ---
echo -e "${BLUE}Processing: C4-2B.CTCF.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph C4-2B.CTCF.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > C4-2B.CTCF.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b C4-2B.CTCF.chm13v2.0.bedGraph -c 4 -o mean > C4-2B.CTCF.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm C4-2B.CTCF.chm13v2.0.bedGraph
rm C4-2B.CTCF.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped C4-2B.CTCF.chm13v2.0.bw to C4-2B.CTCF.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 8/39 ---
echo -e "${BLUE}Processing: C4-2B.H3K27ac.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph C4-2B.H3K27ac.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > C4-2B.H3K27ac.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b C4-2B.H3K27ac.chm13v2.0.bedGraph -c 4 -o mean > C4-2B.H3K27ac.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm C4-2B.H3K27ac.chm13v2.0.bedGraph
rm C4-2B.H3K27ac.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped C4-2B.H3K27ac.chm13v2.0.bw to C4-2B.H3K27ac.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 9/39 ---
echo -e "${BLUE}Processing: Caco-2.H3K4me1.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph Caco-2.H3K4me1.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > Caco-2.H3K4me1.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b Caco-2.H3K4me1.chm13v2.0.bedGraph -c 4 -o mean > Caco-2.H3K4me1.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm Caco-2.H3K4me1.chm13v2.0.bedGraph
rm Caco-2.H3K4me1.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped Caco-2.H3K4me1.chm13v2.0.bw to Caco-2.H3K4me1.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 10/39 ---
echo -e "${BLUE}Processing: Caco-2.H3K9me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph Caco-2.H3K9me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > Caco-2.H3K9me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b Caco-2.H3K9me3.chm13v2.0.bedGraph -c 4 -o mean > Caco-2.H3K9me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm Caco-2.H3K9me3.chm13v2.0.bedGraph
rm Caco-2.H3K9me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped Caco-2.H3K9me3.chm13v2.0.bw to Caco-2.H3K9me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 11/39 ---
echo -e "${BLUE}Processing: HAP-1.H3K4me1.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph HAP-1.H3K4me1.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > HAP-1.H3K4me1.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b HAP-1.H3K4me1.chm13v2.0.bedGraph -c 4 -o mean > HAP-1.H3K4me1.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm HAP-1.H3K4me1.chm13v2.0.bedGraph
rm HAP-1.H3K4me1.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped HAP-1.H3K4me1.chm13v2.0.bw to HAP-1.H3K4me1.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 12/39 ---
echo -e "${BLUE}Processing: HAP-1.H3K4me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph HAP-1.H3K4me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > HAP-1.H3K4me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b HAP-1.H3K4me3.chm13v2.0.bedGraph -c 4 -o mean > HAP-1.H3K4me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm HAP-1.H3K4me3.chm13v2.0.bedGraph
rm HAP-1.H3K4me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped HAP-1.H3K4me3.chm13v2.0.bw to HAP-1.H3K4me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 13/39 ---
echo -e "${BLUE}Processing: HAP-1.H3K9me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph HAP-1.H3K9me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > HAP-1.H3K9me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b HAP-1.H3K9me3.chm13v2.0.bedGraph -c 4 -o mean > HAP-1.H3K9me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm HAP-1.H3K9me3.chm13v2.0.bedGraph
rm HAP-1.H3K9me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped HAP-1.H3K9me3.chm13v2.0.bw to HAP-1.H3K9me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 14/39 ---
echo -e "${BLUE}Processing: HAP-1.H3K27ac.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph HAP-1.H3K27ac.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > HAP-1.H3K27ac.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b HAP-1.H3K27ac.chm13v2.0.bedGraph -c 4 -o mean > HAP-1.H3K27ac.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm HAP-1.H3K27ac.chm13v2.0.bedGraph
rm HAP-1.H3K27ac.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped HAP-1.H3K27ac.chm13v2.0.bw to HAP-1.H3K27ac.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 15/39 ---
echo -e "${BLUE}Processing: HAP-1.H3K27me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph HAP-1.H3K27me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > HAP-1.H3K27me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b HAP-1.H3K27me3.chm13v2.0.bedGraph -c 4 -o mean > HAP-1.H3K27me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm HAP-1.H3K27me3.chm13v2.0.bedGraph
rm HAP-1.H3K27me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped HAP-1.H3K27me3.chm13v2.0.bw to HAP-1.H3K27me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 16/39 ---
echo -e "${BLUE}Processing: HAP-1.H3K36me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph HAP-1.H3K36me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > HAP-1.H3K36me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b HAP-1.H3K36me3.chm13v2.0.bedGraph -c 4 -o mean > HAP-1.H3K36me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm HAP-1.H3K36me3.chm13v2.0.bedGraph
rm HAP-1.H3K36me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped HAP-1.H3K36me3.chm13v2.0.bw to HAP-1.H3K36me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 17/39 ---
echo -e "${BLUE}Processing: HL-60.H3K27ac.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph HL-60.H3K27ac.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > HL-60.H3K27ac.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b HL-60.H3K27ac.chm13v2.0.bedGraph -c 4 -o mean > HL-60.H3K27ac.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm HL-60.H3K27ac.chm13v2.0.bedGraph
rm HL-60.H3K27ac.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped HL-60.H3K27ac.chm13v2.0.bw to HL-60.H3K27ac.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 18/39 ---
echo -e "${BLUE}Processing: MG63.H3K4me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph MG63.H3K4me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > MG63.H3K4me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b MG63.H3K4me3.chm13v2.0.bedGraph -c 4 -o mean > MG63.H3K4me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm MG63.H3K4me3.chm13v2.0.bedGraph
rm MG63.H3K4me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped MG63.H3K4me3.chm13v2.0.bw to MG63.H3K4me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 19/39 ---
echo -e "${BLUE}Processing: MG63.H3K9me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph MG63.H3K9me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > MG63.H3K9me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b MG63.H3K9me3.chm13v2.0.bedGraph -c 4 -o mean > MG63.H3K9me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm MG63.H3K9me3.chm13v2.0.bedGraph
rm MG63.H3K9me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped MG63.H3K9me3.chm13v2.0.bw to MG63.H3K9me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 20/39 ---
echo -e "${BLUE}Processing: MG63.H3K27me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph MG63.H3K27me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > MG63.H3K27me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b MG63.H3K27me3.chm13v2.0.bedGraph -c 4 -o mean > MG63.H3K27me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm MG63.H3K27me3.chm13v2.0.bedGraph
rm MG63.H3K27me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped MG63.H3K27me3.chm13v2.0.bw to MG63.H3K27me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 21/39 ---
echo -e "${BLUE}Processing: MG63.H3K36me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph MG63.H3K36me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > MG63.H3K36me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b MG63.H3K36me3.chm13v2.0.bedGraph -c 4 -o mean > MG63.H3K36me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm MG63.H3K36me3.chm13v2.0.bedGraph
rm MG63.H3K36me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped MG63.H3K36me3.chm13v2.0.bw to MG63.H3K36me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 22/39 ---
echo -e "${BLUE}Processing: RWPE1.CTCF.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph RWPE1.CTCF.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > RWPE1.CTCF.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b RWPE1.CTCF.chm13v2.0.bedGraph -c 4 -o mean > RWPE1.CTCF.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm RWPE1.CTCF.chm13v2.0.bedGraph
rm RWPE1.CTCF.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped RWPE1.CTCF.chm13v2.0.bw to RWPE1.CTCF.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 23/39 ---
echo -e "${BLUE}Processing: RWPE1.H3K27ac.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph RWPE1.H3K27ac.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > RWPE1.H3K27ac.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b RWPE1.H3K27ac.chm13v2.0.bedGraph -c 4 -o mean > RWPE1.H3K27ac.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm RWPE1.H3K27ac.chm13v2.0.bedGraph
rm RWPE1.H3K27ac.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped RWPE1.H3K27ac.chm13v2.0.bw to RWPE1.H3K27ac.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 24/39 ---
echo -e "${BLUE}Processing: RWPE2.CTCF.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph RWPE2.CTCF.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > RWPE2.CTCF.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b RWPE2.CTCF.chm13v2.0.bedGraph -c 4 -o mean > RWPE2.CTCF.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm RWPE2.CTCF.chm13v2.0.bedGraph
rm RWPE2.CTCF.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped RWPE2.CTCF.chm13v2.0.bw to RWPE2.CTCF.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 25/39 ---
echo -e "${BLUE}Processing: RWPE2.H3K27ac.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph RWPE2.H3K27ac.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > RWPE2.H3K27ac.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b RWPE2.H3K27ac.chm13v2.0.bedGraph -c 4 -o mean > RWPE2.H3K27ac.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm RWPE2.H3K27ac.chm13v2.0.bedGraph
rm RWPE2.H3K27ac.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped RWPE2.H3K27ac.chm13v2.0.bw to RWPE2.H3K27ac.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 26/39 ---
echo -e "${BLUE}Processing: RWPE2.H3K36me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph RWPE2.H3K36me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > RWPE2.H3K36me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b RWPE2.H3K36me3.chm13v2.0.bedGraph -c 4 -o mean > RWPE2.H3K36me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm RWPE2.H3K36me3.chm13v2.0.bedGraph
rm RWPE2.H3K36me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped RWPE2.H3K36me3.chm13v2.0.bw to RWPE2.H3K36me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 27/39 ---
echo -e "${BLUE}Processing: SJCRH30.H3K4me1.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph SJCRH30.H3K4me1.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > SJCRH30.H3K4me1.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b SJCRH30.H3K4me1.chm13v2.0.bedGraph -c 4 -o mean > SJCRH30.H3K4me1.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm SJCRH30.H3K4me1.chm13v2.0.bedGraph
rm SJCRH30.H3K4me1.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped SJCRH30.H3K4me1.chm13v2.0.bw to SJCRH30.H3K4me1.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 28/39 ---
echo -e "${BLUE}Processing: SJCRH30.H3K4me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph SJCRH30.H3K4me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > SJCRH30.H3K4me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b SJCRH30.H3K4me3.chm13v2.0.bedGraph -c 4 -o mean > SJCRH30.H3K4me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm SJCRH30.H3K4me3.chm13v2.0.bedGraph
rm SJCRH30.H3K4me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped SJCRH30.H3K4me3.chm13v2.0.bw to SJCRH30.H3K4me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 29/39 ---
echo -e "${BLUE}Processing: SJCRH30.H3K9me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph SJCRH30.H3K9me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > SJCRH30.H3K9me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b SJCRH30.H3K9me3.chm13v2.0.bedGraph -c 4 -o mean > SJCRH30.H3K9me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm SJCRH30.H3K9me3.chm13v2.0.bedGraph
rm SJCRH30.H3K9me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped SJCRH30.H3K9me3.chm13v2.0.bw to SJCRH30.H3K9me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 30/39 ---
echo -e "${BLUE}Processing: SJCRH30.H3K27ac.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph SJCRH30.H3K27ac.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > SJCRH30.H3K27ac.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b SJCRH30.H3K27ac.chm13v2.0.bedGraph -c 4 -o mean > SJCRH30.H3K27ac.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm SJCRH30.H3K27ac.chm13v2.0.bedGraph
rm SJCRH30.H3K27ac.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped SJCRH30.H3K27ac.chm13v2.0.bw to SJCRH30.H3K27ac.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 31/39 ---
echo -e "${BLUE}Processing: SJCRH30.H3K27me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph SJCRH30.H3K27me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > SJCRH30.H3K27me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b SJCRH30.H3K27me3.chm13v2.0.bedGraph -c 4 -o mean > SJCRH30.H3K27me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm SJCRH30.H3K27me3.chm13v2.0.bedGraph
rm SJCRH30.H3K27me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped SJCRH30.H3K27me3.chm13v2.0.bw to SJCRH30.H3K27me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 32/39 ---
echo -e "${BLUE}Processing: SJCRH30.H3K36me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph SJCRH30.H3K36me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > SJCRH30.H3K36me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b SJCRH30.H3K36me3.chm13v2.0.bedGraph -c 4 -o mean > SJCRH30.H3K36me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm SJCRH30.H3K36me3.chm13v2.0.bedGraph
rm SJCRH30.H3K36me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped SJCRH30.H3K36me3.chm13v2.0.bw to SJCRH30.H3K36me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 33/39 ---
echo -e "${BLUE}Processing: SJSA1.H3K4me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph SJSA1.H3K4me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > SJSA1.H3K4me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b SJSA1.H3K4me3.chm13v2.0.bedGraph -c 4 -o mean > SJSA1.H3K4me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm SJSA1.H3K4me3.chm13v2.0.bedGraph
rm SJSA1.H3K4me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped SJSA1.H3K4me3.chm13v2.0.bw to SJSA1.H3K4me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 34/39 ---
echo -e "${BLUE}Processing: SJSA1.H3K9me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph SJSA1.H3K9me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > SJSA1.H3K9me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b SJSA1.H3K9me3.chm13v2.0.bedGraph -c 4 -o mean > SJSA1.H3K9me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm SJSA1.H3K9me3.chm13v2.0.bedGraph
rm SJSA1.H3K9me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped SJSA1.H3K9me3.chm13v2.0.bw to SJSA1.H3K9me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 35/39 ---
echo -e "${BLUE}Processing: SJSA1.H3K27ac.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph SJSA1.H3K27ac.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > SJSA1.H3K27ac.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b SJSA1.H3K27ac.chm13v2.0.bedGraph -c 4 -o mean > SJSA1.H3K27ac.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm SJSA1.H3K27ac.chm13v2.0.bedGraph
rm SJSA1.H3K27ac.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped SJSA1.H3K27ac.chm13v2.0.bw to SJSA1.H3K27ac.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 36/39 ---
echo -e "${BLUE}Processing: SJSA1.H3K27me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph SJSA1.H3K27me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > SJSA1.H3K27me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b SJSA1.H3K27me3.chm13v2.0.bedGraph -c 4 -o mean > SJSA1.H3K27me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm SJSA1.H3K27me3.chm13v2.0.bedGraph
rm SJSA1.H3K27me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped SJSA1.H3K27me3.chm13v2.0.bw to SJSA1.H3K27me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 37/39 ---
echo -e "${BLUE}Processing: SJSA1.H3K36me3.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph SJSA1.H3K36me3.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > SJSA1.H3K36me3.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b SJSA1.H3K36me3.chm13v2.0.bedGraph -c 4 -o mean > SJSA1.H3K36me3.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm SJSA1.H3K36me3.chm13v2.0.bedGraph
rm SJSA1.H3K36me3.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped SJSA1.H3K36me3.chm13v2.0.bw to SJSA1.H3K36me3.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 38/39 ---
echo -e "${BLUE}Processing: VCaP.CTCF.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph VCaP.CTCF.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > VCaP.CTCF.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b VCaP.CTCF.chm13v2.0.bedGraph -c 4 -o mean > VCaP.CTCF.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm VCaP.CTCF.chm13v2.0.bedGraph
rm VCaP.CTCF.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped VCaP.CTCF.chm13v2.0.bw to VCaP.CTCF.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
# --- Processing file 39/39 ---
echo -e "${BLUE}Processing: VCaP.H3K27ac.chm13v2.0.bw${NC}"
echo '  Step 1/3: Converting bigWig to bedGraph and sorting...'
../bigWigToBedGraph VCaP.H3K27ac.chm13v2.0.bw stdout | sort -k1,1 -k2,2n > VCaP.H3K27ac.chm13v2.0.bedGraph
echo '  Step 2/3: Mapping features with bedtools...'
bedtools map -a ../T2T_repeat_masker_processed_sorted.bed -b VCaP.H3K27ac.chm13v2.0.bedGraph -c 4 -o mean > VCaP.H3K27ac.chm13v2.0.mapped_on_repeat_masker.bedGraph
echo '  Step 3/3: Cleaning up intermediate files...'
rm VCaP.H3K27ac.chm13v2.0.bedGraph
rm VCaP.H3K27ac.chm13v2.0.bw
echo -e "${GREEN}Successfully mapped VCaP.H3K27ac.chm13v2.0.bw to VCaP.H3K27ac.chm13v2.0.mapped_on_repeat_masker.bedGraph${NC}\n"
echo -e "${GREEN}All files processed successfully!${NC}"
