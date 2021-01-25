exes = ['generate_gt', 'ibdmix']

rule make_executables:
    input:
        [os.path.join(paths['source_root'], 'src', f)
                for f in os.listdir(os.path.join(paths['source_root'], 'src'))
                if (os.path.join(paths['source_root'], 'src', f)
                    and not f.endswith('.sh'))
                ],
        [os.path.join(paths['source_root'], 'include/IBDmix', f)
                for f in os.listdir(os.path.join(paths['source_root'],
                    'include/IBDmix'))
                if os.path.join(paths['source_root'], 'include/IBDmix', f)]

    output:
        expand(paths['exe'], exe=exes)

    params:
        targets=expand(
                str(Path(paths['source_root']) / 'snake_build' /
                    'src' / '{exe}'),
                exe=exes),
        root=str(Path(paths['exe']).parents[0]),
        build_dir=str(Path(paths['source_root']) / 'snake_build')

    envmodules:
        config['cmake_module']

    threads: 4

    shell:
        '''
        mkdir -p {params.build_dir}
        cd {params.build_dir}
        cmake .. -DCMAKE_BUILD_TYPE=Release
        cmake --build . -j {threads}
        cd -
        for f in {params.targets}; do mv $f {params.root} ; done
        rm -rf {params.build_dir}
        '''

rule generate_gt:
    input:
        modern=paths['modern_vcf'],
        archaic=paths['archaic_vcf'],
        exe=paths['exe'].format(exe='generate_gt')

    output:
        paths['genotype_file']

    threads: 2

    shell:
        '{input.exe} '
            '--archaic <(zcat -f {input.archaic}) '
            '--modern <(zcat {input.modern}) '
            '--output >(gzip > {output}) '

def get_ibd_input(wildcards):
    result = {'genotype': paths['genotype_file'].format(**wildcards),
              'exe': paths['exe'].format(exe='ibdmix')}

    if 'population' in wildcards.keys() and 'sample_file' in paths:
        result['samples'] = paths['sample_file'].format(**wildcards)

    if 'mask_file' in paths:
        result['mask'] = paths['mask_file'].format(**wildcards)

    return result

def get_ibd_options(wildcards, input):
    result = ''
    if 'samples' in input.keys():
        result += f'--sample {input.samples} '
    if 'mask' in input.keys():
        result += f'--mask {input.mask} '

    return result

rule ibdmix:
    input:
        unpack(get_ibd_input)

    output:
        temp(paths['ibd_output'])

    params: options=get_ibd_options

    shell:
        '{input.exe} '
            '--genotype <(zcat {input.genotype}) '
            '--output >(gzip > {output}) '
            '{params.options} '
            '{config[IBDmix][options]} '

