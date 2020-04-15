#!/bin/bash

set -euo pipefail

for test_type in extra gen_20 populations all_mask all_no_mask; do
    time nice ./acceptance_test.sh $test_type #&
done
wait
