for i in {74..100}; do
    for c in {1..10}; do
        if ! test -f logs/dynamic-eqtl-permutecells/pb-cf-bin16-50k-0-5-F-$c-$i.out; then
            sbatch --mem=15G -J dp-$i-$c --error=logs/dynamic-eqtl-permutecells/pb-cf-bin16-50k-0-5-F-$c-$i.err --output=logs/dynamic-eqtl-permutecells/pb-cf-bin16-50k-0-5-F-$c-$i.out eqtl_dynamic_permutecells.sh "pseudobulk-cf" "50k" 0 5 "F" "bin16" $c 10 $i
        fi
        if ! test -f logs/dynamic-eqtl-permutecells/pb-cm-bin16-50k-0-5-F-$c-$i.out; then
            sbatch --mem=15G -J dp-$i-$c --error=logs/dynamic-eqtl-permutecells/pb-cm-bin16-50k-0-5-F-$c-$i.err --output=logs/dynamic-eqtl-permutecells/pb-cm-bin16-50k-0-5-F-$c-$i.out eqtl_dynamic_permutecells.sh "pseudobulk-cm" "50k" 0 5 "F" "bin16" $c 10 $i
        fi
    done
done
