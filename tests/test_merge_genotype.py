from IBDmix.merge_genotype import VCF_File, File_Merger
from io import StringIO
import pytest
import numpy as np
from numpy.testing import assert_array_equal as aae


def test_file_merger_run(mocker):
    # mock most of the function away, just test the right functions will be
    # called
    mock_open = mocker.patch('IBDmix.merge_genotype.open')
    mock_gzip = mocker.patch('IBDmix.merge_genotype.gzip.open')
    mock_vcf = mocker.patch.object(VCF_File, '__init__', return_value=None)
    mock_merge = mocker.patch.object(File_Merger, 'merge')

    merger = File_Merger('dir/test.vcf', 'test.vcf.gz')
    merger.run('test_*.txt')
    mock_open.assert_called_once_with('dir/test.vcf', 'r')
    mock_gzip.assert_called_once_with('test.vcf.gz', 'rt')
    assert mock_vcf.call_args_list == [
        mocker.call(mocker.ANY, 'dir/test.vcf'),
        mocker.call(mocker.ANY, 'test.vcf.gz'),
    ]
    mock_merge.return_value.to_csv.assert_called_once_with(
        'test_*.txt', sep='\t')

    mocker.resetall()
    merger = File_Merger('dir/test.vcf.gz', 'test.vcf')
    merger.run('test_*.txt')
    mock_open.assert_called_once_with('test.vcf', 'r')
    mock_gzip.assert_called_once_with('dir/test.vcf.gz', 'rt')
    assert mock_vcf.call_args_list == [
        mocker.call(mocker.ANY, 'dir/test.vcf.gz'),
        mocker.call(mocker.ANY, 'test.vcf'),
    ]
    mock_merge.return_value.to_csv.assert_called_once_with(
        'test_*.txt', sep='\t')


@pytest.fixture
def sample_vcf(tmp_path):
    vcf = (
        '##fileformat=VCFv4.1\n'
        '##fileDate=10122015_22h01m13s\n'
        '##source=SHAPEIT2.v837\n'
        '##FORMAT=<ID=GT,Type=String,Description="Phased\tGenotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t'
        'INFO\tFORMAT\tI1\tI2\tI3\tI4\tI5\n'
        '1\t846687\t1:846687\tC\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0\n'
        '1\t846688\t1:846688\tG\tA\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0\n'
        '1\t846742\t1:846742\tC\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0\n'
        '1\t846758\t1:846758\tG\tA\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0\n'
        '1\t846808\t1:846808\tC\tT\t.\tPASS\t.\tGT'
        '\t0|0\t0|0\t0|1\t0|0\t0/0:249:21.05:0,21,2\n'
        '1\t846824\t1:846824\tC\tT\t.\tPASS\t.\tGT\t0a0\t1|0\t1|1\t.|.\t./.\n'
        '1\t846824\t1:846\tTC\tT\t.\tPASS\t.\tGT\t0a0\t1|0\t1|1\t.|.\t./.\n'
        '1\t846824\t1:846\tT\tAT\t.\tPASS\t.\tGT\t0a0\t1|0\t1|1\t.|.\t./.\n'
    )
    path = tmp_path / 'data.vcf'
    path.write_text(vcf)
    return VCF_File(StringIO(vcf), path)


@pytest.fixture
def sample_arch_vcf(tmp_path):
    vcf = (
        '##reference=file:///mnt/solexa/Genomes/hg19_1000g/whole_genome.fa\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tAltaiNea\n'
        '1\t10001\t.\tT\tC\t51.05\t.\t.\t.\t0/1:249:21.05:0,21,265:1,2:0,1:\n'
        '1\t10002\t.\tA\tT\t36.01\t.\t.\t.\t1/0:250:6.01:0,6,69:171,184:0,9\n'
        '1\t10003\t.\tA\tT\t36.01\t.\t.\t.\t1/1:250:6.02:0,6,72:222,201:1,0\n'
        '1\t846808\t.\tC\t.\t39.01\t.\t.\t.\t0/0:250:9.02:0,9,97:0,1:268,211\n'
        '1\t846824\t.\tC\tT\t39.01\t.\t.\t.\t0/0:250:9.02:0,9,97:0,1:268,211\n'
        '1\t10005\t.\tC\t.\t41.99\t.\t.\t.\t0/0:250:12:0,12,119:0,0:341,287\n'
        '1\t10006\t.\tC\t.\t36.01\t.\t.\t.\t0/0:249:6.02:0,6,72:0,1:408,347\n'
        '1\t10622\t.\tT\t.\t.\t.\t.\tGT:A:C:G:T:IR\t./.:0,0:0,0:0,1:1,0:0\n'
    )
    path = tmp_path / 'data2.vcf'
    path.write_text(vcf)
    return VCF_File(StringIO(vcf), path)


@pytest.fixture
def sample_arch_vcf_2(tmp_path):
    vcf = (
        '##reference=file:///mnt/solexa/Genomes/hg19_1000g/whole_genome.fa\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tAltaiNea\tT\n'
        '1\t10001\t.\tT\tC\t51.05\t.\t.\t.\t0/1\t0/0\n'
        '1\t10622\t.\tT\t.\t.\t.\t.\tGT:A:C:G:T:IR\t1/1\t0/1\n'
        '1\t10623\t.\tT\t.\t.\t.\t.\tGT:A:C:G:T:IR\t./.\t0/0\n'
        '1\t10624\t.\tT\t.\t.\t.\t.\tGT:A:C:G:T:IR\t0/0\t./.\n'
        '1\t10625\t.\tT\t.\t.\t.\t.\tGT:A:C:G:T:IR\t./.\t./.\n'
    )
    path = tmp_path / 'data2.vcf'
    path.write_text(vcf)
    return VCF_File(StringIO(vcf), path)


def test_file_merger_merge(sample_vcf, sample_arch_vcf):
    merger = File_Merger('arch', 'mod')
    result = merger.merge(sample_arch_vcf, sample_vcf).compute()

    aae(result.values,
        np.array([
            [1, 'T', 'C', 1, 0, 0, 0, 0, 0],
            [1, 'A', 'T', 1, 0, 0, 0, 0, 0],
            [1, 'A', 'T', 2, 0, 0, 0, 0, 0],
            [1, 'C', 'T', 0, 0, 0, 1, 0, 0],
            [1, 'C', 'T', 0, 0, 1, 2, 9, 9],
        ], dtype=object))

    aae(result.index.values,
        np.array([10001, 10002, 10003, 846808, 846824]))


def test_vcf_file_init(sample_vcf):
    assert sample_vcf.header == (
        'CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t'
        'INFO\tFORMAT\tI1\tI2\tI3\tI4\tI5').split()


def test_vcf_get_dataframe(sample_vcf):
    result = sample_vcf.get_dataframe()

    # chunk of 4 lines
    aae(result.values,
        np.array([
            [1, 'C', 'T', 0, 0, 0, 0, 0],
            [1, 'G', 'A', 0, 0, 0, 0, 0],
            [1, 'C', 'T', 0, 0, 0, 0, 0],
            [1, 'G', 'A', 0, 0, 0, 0, 0],
            [1, 'C', 'T', 0, 0, 1, 0, 0],
            [1, 'C', 'T', 0, 1, 2, 9, 9],
        ], dtype=object))


def test_vcf_get_dataframe_arch(sample_arch_vcf):
    result = sample_arch_vcf.get_dataframe(skip_non_informative=True)

    # chunk of 4 lines
    aae(result.values,
        np.array([
            [1, 'T', 'C', 1],
            [1, 'A', 'T', 1],
            [1, 'A', 'T', 2],
            [1, 'C', '.', 0],
            [1, 'C', 'T', 0],
            [1, 'C', '.', 0],
            [1, 'C', '.', 0],
        ], dtype=object))


def test_vcf_get_dataframe_mult_arch(sample_arch_vcf_2):
    result = sample_arch_vcf_2.get_dataframe(skip_non_informative=True)
    print(result.compute().values)

    # chunk of 4 lines
    aae(result.values,
        np.array([
            [1, 'T', 'C', 1, 0],
            [1, 'T', '.', 2, 1],
            [1, 'T', '.', 9, 0],
            [1, 'T', '.', 0, 9],
        ], dtype=object))


def test_vcf_get_sample_header(sample_vcf):
    assert sample_vcf.get_sample_header() == 'I1 I2 I3 I4 I5'.split()


def test_vcf_get_initial_header(sample_vcf):
    assert sample_vcf.get_initial_header() == 'CHROM POS REF ALT'.split()
