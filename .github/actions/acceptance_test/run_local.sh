#!/bin/bash

set -euo pipefail

for test_type in extra gen_20 all_mask all_no_mask populations; do
    time nice ./acceptance_test.sh $test_type #&
done
wait
