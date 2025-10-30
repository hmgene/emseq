#!/usr/bin/env bash
diff <(metilene -X 8 -Y 8 -t 2 -m 8 -M 300 -a 2M -b 30M input.tsv 2> >(grep -v segmenting >&2) | sort -k1,1 -k2,2n -k3,3n) output.tsv
ex=$?
[[ $ex -ne 0 ]] && echo "Test failed" >&2 || echo "Test passed" >&2
exit $ex
