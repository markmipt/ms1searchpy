#!/bin/sh

TESTDIR="test_data"

if [ ! -d "$TESTDIR" ]; then
    echo "You must have the test dataset to run this test."
    exit 1
fi

cd "$TESTDIR"
rm -vf *.features* *.tsv *.txt

echo ""
echo "Starting ms1searchpy ..."
echo "------------------------"
time ms1searchpy -d sprot_ecoli_ups.fasta -ad 1 -nproc 6 -deeplc deeplc -deeplc_library deeplc.lib *.mzML \
    && echo "DirectMS1 run successful."

echo ""
echo "Starting ms1combine ..."
echo "-----------------------"
time ms1combine *_UPS_4_0?.features_PFMs_ML.tsv -out UPS_4 && echo "ms1combine run successful."
time ms1combine *_UPS_2_0?.features_PFMs_ML.tsv -out UPS_2 && echo "ms1combine run successful."

echo ""
echo "Starting ms1todiffacto ..."
echo "--------------------------"
time ms1todiffacto -dif diffacto -S1 *_UPS_4_0?.features_proteins.tsv -S2 *_UPS_2_0?.features_proteins.tsv \
    -norm median -out diffacto_output.tsv -min_samples 3 && echo "ms1todiffacto run successful."

echo ""
echo "Starting directms1quant ..."
echo "---------------------------"
time directms1quant -S1 *_UPS_4_0?.features_proteins_full.tsv -S2 *_UPS_2_0?.features_proteins_full.tsv -min_samples 3 \
    && echo "DirectMS1quant run successful."
