import gzip

cpg_species = ['panTro2', 'ponAbe2', 'rheMac2']

rule build_seqbility:
    output: directory(paths['seqbility_root'])
    shell:
        'wget --quiet -O {output}.tar.bz2 {config[urls][seqbility]} \n'
        'bzip2 -d {output}.tar.bz2 \n'
        'tar -xf {output}.tar \n'
        'rm {output}.tar \n'
        'mv {output}-20091110 {output} \n'
        'cd {output} \n'
        'make \n'

rule samtools_index:
    input: ancient(paths['reference_genome'])
    output: paths['reference_genome'] + '.fai'
    singularity:
        config['containers']['samtools']
    shell:
        'samtools faidx {input}'

rule chrom_sizes:
    # keep only digits, make bed file
    input: paths['reference_genome'] + '.fai'
    output: paths['chrom_sizes']
    shell:
        'sed -n -e '
             "'s/^\([1-9][0-9]\?\\t\)\([0-9]\+\).*/\\10\\t\\2/p' "
             '{input} '
        '| sort -n -k1 '
        '> {output}'

checkpoint split_fa:
    input:
        seqbility=rules.build_seqbility.output,
        reference=paths['reference_genome']
    output:
        directory(Path(paths['split_fa']).parents[0])
    shell:
        'mkdir -p {output} \n'
        '{input.seqbility}/splitfa '
            '{input.reference} 35 '
            '| split '
                '--lines 20000000 '
                '--filter=\'gzip > $FILE\' '
                '- '  # read from stdin
                '{output}/split_ '  # suffix

bwa_outputs = multiext(paths['reference_genome'],
                       '.bwt', '.pac', '.ann', '.amb', '.sa')

rule bwa_index:
    input: ancient(paths['reference_genome'])
    output: bwa_outputs
    singularity:
        config['containers']['bwa']
    shell:
        'bwa index '
            '-a bwtsw '
            '{input} '
            '2> /dev/null '

rule bwa_aln:
    input:
        fa=paths['split_fa'],
        reference=paths['reference_genome'],
        ind=rules.bwa_index.output
    output:
        temp(paths['mapped_sai'])
    singularity:
        config['containers']['bwa']
    shell:
        'bwa aln '
            '-R 1000000 '
            '-O 3 '
            '-E 3 '
            '{input.reference} '
            '<(zcat {input.fa}) '
        '> {output} '
        '2> /dev/null '

rule bwa_samse:
    input:
        sai=paths['mapped_sai'],
        fa=paths['split_fa'],
        reference=paths['reference_genome']
    output:
        temp(paths['mapped_sam'])
    singularity:
        config['containers']['bwa']
    priority: 1
    shell:
        'bwa samse '
            '{input.reference} '
            '{input.sai} '
            '<(zcat {input.fa}) '
        '2> /dev/null '
        '| gzip > {output} '

def gen_raw_mask_input(wildcards):
    checkpoints.split_fa.get(**wildcards)
    parts = glob_wildcards(paths['split_fa']).part
    return {
            'sams': expand(paths['mapped_sam'], part=sorted(parts)),
            'seqbility': rules.build_seqbility.output,
            }

rule gen_raw_mask:
    input:
        unpack(gen_raw_mask_input)
    output:
        paths['raw_map_mask']
    shell:
        'gzip -dc {input.sams} '
            '| {input.seqbility}/gen_raw_mask.pl '
            '2> /dev/null '
            '> {output}'

rule gen_mapability_mask:
    input:
        seqbility=rules.build_seqbility.output,
        raw=paths['raw_map_mask']
    output:
        paths['map_mask']
    shell:
        '{input.seqbility}/gen_mask '
            '-l 35 '
            '-r 0.5 '
            '{input.raw} '
        '> {output}'

rule map_fa_to_bed:
    input:
        paths['map_mask']
    output:
        paths['map_bed']
    run:
        threshold = '3'  # only retain unique
        chrom = ''
        chrom_pos = 0
        region_start = 0
        in_region = False

        with open(output[0], 'w') as outfile, \
                open(input[0], 'r') as infile:
            for line in infile:
                line = line.strip()
                if line[0] == '>':
                    if in_region:  # emit region from last chrom
                        outfile.write(f'{chrom}\t{region_start}\t{chrom_pos}\n')
                    chrom = line.split()[0][1:]
                    chrom_pos = 0
                    region_start = 0
                    in_region = False
                    continue
                
                for i,b in enumerate(line):
                    # start new region
                    if not in_region and b == threshold:
                        region_start = chrom_pos + i
                        in_region = True
                    elif in_region and b != threshold:
                        outfile.write(f'{chrom}\t{region_start}\t{chrom_pos+i}\n')
                        in_region = False

                chrom_pos += len(line)
            if in_region:  # emit region from last chrom
                outfile.write(f'{chrom}\t{region_start}\t{chrom_pos}\n')

def get_autosome_input(wildcards):
    if Path(paths['mapability_mask']).exists():
        return {}
    return {
            'bed': paths['map_bed'],
            'chrom_sizes': paths['chrom_sizes']
            }

rule get_autosomes:  # invert genome, and clean up seq temp
    input:
        unpack(get_autosome_input)
    output:
        paths['mapability_mask']
    singularity:
        config['containers']['bedtools']
    shell:
        "sed -n '/^[0-9]\{{1,2\}}[^_0-9]/p' {input.bed} "
            '| sort -k1,1n -k2,2n -k3,3n '
        '| bedtools subtract '
            '-a {input.chrom_sizes} '
            '-b - '
        '> {output} \n'
        'rm -rf {paths[seqbility_temp]}'

rule generate_indel_mask:
    input:
        paths['modern_vcf']
    output:
        temp(paths['indel_mask'])
    run:
        WINDOW = 5
        with open(output[0], 'w') as outfile:
            for line in gzip.open(input[0], 'rt'):
                if line.startswith('#'):
                    continue
                chrom, pos, _, ref, alt, *rest = line.split()
                if ',' in ref:
                    print('Unexpected input for indel mask, contains comma')
                    print(chrom, pos, ref, alt)
                    continue
                pos = int(pos)
                if len(ref) > 1:
                    outfile.write(f'{chrom}\t{pos-WINDOW-1}\t{pos+WINDOW+len(ref)}\n')
                if len(alt) > 1:
                    outfile.write(f'{chrom}\t{pos-WINDOW-1}\t{pos+WINDOW}\n')

def download_axt_params(wildcards):
    return config['urls']['syn_net_url'].format(
            species_upper=wildcards.species[0].upper() + wildcards.species[1:],
            species=wildcards.species,
            chrom=wildcards.chrom)

rule download_axt:
    output:
        paths['syn_net']
    params:
        url=download_axt_params
    resources:
        wget_instances=1
    shell:
        'wget --quiet '
            '-O {output} '
            '{params.url} '

rule generate_simple_snps:
    output:
        paths['simple_snps']
    input:
        paths['syn_net']
    run:
        with gzip.open(output[0], 'wt') as outfile, \
                gzip.open(input[0], 'rt') as infile:
            while True:
                header = infile.readline()
                if header == '':
                    return
                human = infile.readline().upper()
                other = infile.readline().upper()
                infile.readline()  # blank line
                _, chrom, start, *rest = header.split()

                chrom = chrom[3:]  # split 'chr'
                start = int(start)

                ind = 0
                for hum, oth in zip(human, other):
                    if hum == '-':  # skip del
                        continue
                    if hum != oth and oth != '-':
                        outfile.write(f'{chrom}\t{start+ind}\t{hum}\t{oth}\n')
                    ind += 1

def generate_cpg_mask_input(wildcards):
    inputs = {species: expand(paths['simple_snps'],
        species=species,
        chrom=chromosomes) for species in cpg_species}
    inputs['modern'] = str(Path(workflow.snakefile).parents[0] / 'abridged_variants.gz')
    inputs['reference'] = paths['reference_genome']
    return inputs

rule generate_cpg_mask:
    output:
        paths['cpg_temp_mask']
    input:
        unpack(generate_cpg_mask_input)
    script:
        'generate_cpg_mask.py'

rule move_cpg_mask:
    output:
        temp(paths['cpg_mask'])
    input:
        paths['cpg_temp_mask']
    shell:
        'sort -k1,1n -k2,2n -k3,3n '
            '{input} > {output}\n'
        'rm -rf {paths[cpg_temp]}'

rule seg_dups_to_bed:
    output:
        paths['seg_dups_bed']
    input:
        paths['seg_dups']
    singularity:
        config['containers']['bedtools']
    shell:
        'zcat {input} '
        '| cut -f 2-4 '
        "| sed -nr 's/^chr([1-9][0-9]?\\t)/\\1/p'"
        '| bedtools merge -i - '
        '| sort -n -k1 -k2'
        '> {output} '

rule process_1kg_bed:
    output:
        paths['bed_1kg_processed']
    input:
        paths['bed_1kg']
    shell:
        'cut -f1-3 {input} '
        "| sed 's/^chr//' "
        '> {output}'

def intersect_included_beds_input(wildcards):
    if wildcards.sample_name == 'altai':
        return [paths['altai_filter'],
                paths['bed_1kg_processed']]
    elif wildcards.sample_name == 'denisovan':
        return [paths['denisovan_filter'],
                paths['bed_1kg_processed']]
    raise ValueError(f'Unknown sample type {wildcards.sample_name}')

rule intersect_included_beds:
    output:
        temp(paths['included_regions'])
    input:
        intersect_included_beds_input
    singularity:
        config['containers']['bedtools']
    shell:
        'bedtools intersect '
            '-a {input[0]} '
            '-b {input[1]} '
        '| sort -n -k1 -k2 -k3 '
            '> {output}'

def union_excluded_beds_input(wildcards):
    if os.path.exists(paths['excluded_regions']):
        return ""
    else:
        return [
                expand(paths['indel_mask'], chrom=chromosomes),
                paths['mapability_mask'],
                paths['cpg_mask'],
                paths['seg_dups_bed']
                ]
rule union_excluded_beds:
    output:
        paths['excluded_regions']
    input: union_excluded_beds_input
    singularity:
        config['containers']['bedtools']
    shell:
        'sort '
            '-n -k1 -k2 -k3 '
            '--merge {input} '
            '| bedtools merge -i - '
            '> {output} '

rule subtract_excluded:
    output:
        temp(paths['included_mask'])
    input:
        inc=paths['included_regions'],
        exc=paths['excluded_regions'],
    singularity:
        config['containers']['bedtools']
    shell:
        'bedtools subtract '
            '-a {input.inc} '
            '-b {input.exc} '
        '| sort -n -k1 -k2 -k3 '
            '> {output} '

rule generate_mask:
    output:
        paths['mask_file']
    input:
        mask=paths['included_mask'],
        sizes=paths['chrom_sizes']
    params:
        chrom=lambda wildcards: f'^{wildcards.chrom}\\t'
    singularity:
        config['containers']['bedtools']
    shell:
        'bedtools subtract '
            "-a <(sed -n '/{params.chrom}/p' {input.sizes}) "
            "-b <(sed -n '/{params.chrom}/p' {input.mask}) "
        '| sort -n -k1 -k2 -k3 '
            '> {output} '
