#!/bin/sh
# properties = {"type": "single", "rule": "RepeatSeq_RepeatMasker_soft_01", "local": false, "input": ["/public/home/zpliu/Pan-genome/RepeatMasked/RepeatMasker/Jin668_allRepeats.lib", "/public/home/zpliu/Pan-genome/RepeatMasked/HC04-Repeat_TRF_softMasked.fa"], "output": ["/public/home/zpliu/Pan-genome/RepeatMasked/RepeatMasker/HC04-Repeat_TRF_softMasked.tbi"], "wildcards": {}, "params": {"soft_outdir": "/public/home/zpliu/Pan-genome/RepeatMasked/RepeatMasker"}, "log": ["logs/RepeatSeq_RepeatMasker_soft_01.log"], "threads": 20, "resources": {"mem_mb": 50000, "mem_mib": 47684, "disk_mb": 4565, "disk_mib": 4354, "tmpdir": "<TBD>"}, "jobid": 0, "cluster": {}}
cd '/public/home/zpliu/Pan-genome/RepeatMasked' && /public/home/zpliu/miniconda3/bin/python -m snakemake --snakefile '/public/home/zpliu/Pan-genome/RepeatMasked/mask.smk' --target-jobs 'RepeatSeq_RepeatMasker_soft_01:' --allowed-rules 'RepeatSeq_RepeatMasker_soft_01' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=50000' 'mem_mib=47684' 'disk_mb=4565' 'disk_mib=4354' --wait-for-files '/public/home/zpliu/Pan-genome/RepeatMasked/.snakemake/tmp.l12fcbnb' '/public/home/zpliu/Pan-genome/RepeatMasked/RepeatMasker/Jin668_allRepeats.lib' '/public/home/zpliu/Pan-genome/RepeatMasked/HC04-Repeat_TRF_softMasked.fa' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'input' 'code' 'software-env' 'mtime' 'params' --skip-script-cleanup  --conda-frontend 'mamba' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --configfiles '/public/home/zpliu/Pan-genome/RepeatMasked/inputConfig.json' --latency-wait 5 --scheduler 'greedy' --scheduler-solver-path '/public/home/zpliu/miniconda3/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && exit 0 || exit 1
