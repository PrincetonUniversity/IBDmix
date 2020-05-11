from pytest import approx
import translate_distance
from io import StringIO
import random


'''
The following tests set up a linear range of some sort over input (bp) 0-100
then test the outputs of all values in the range.
Later tests randomize the slopes, start-end points and order of testing
'''


def test_singleton():
    infile = StringIO(
        "chrom\tpos\trate\tmap\n"
        ".\t0\t.\t3\n"
    )
    mapper = translate_distance.gen_mapper(infile)
    for i in range(100):
        assert mapper.predict(i) == 3

    assert mapper.predict(-1) == 3
    assert mapper.predict(101) == 3


def test_single_slope():
    infile = StringIO(
        "chrom\tpos\trate\tmap\n"
        ".\t0\t.\t0\n"
        ".\t100\t.\t10\n"
    )
    mapper = translate_distance.gen_mapper(infile)
    for i in range(100):
        assert mapper.predict(i) == approx(i/10)

    assert mapper.predict(-1) == approx(0)
    assert mapper.predict(101) == approx(10)


def test_piecewise():
    infile = StringIO(
        "chrom\tpos\trate\tmap\n"
        ".\t0\t.\t0\n"
        ".\t50\t.\t5\n"
        ".\t100\t.\t15\n"
    )
    mapper = translate_distance.gen_mapper(infile)
    for i in range(50):
        assert mapper.predict(i) == approx(i/10)

    for i in range(50, 100):
        assert mapper.predict(i) == approx(5 + (i-50)/5)

    assert mapper.predict(-1) == approx(0)
    assert mapper.predict(101) == approx(15)


def test_piecewise_random():
    infile = StringIO(
        "chrom\tpos\trate\tmap\n"
        ".\t0\t.\t0\n"
        ".\t50\t.\t5\n"
        ".\t100\t.\t15\n"
    )
    mapper = translate_distance.gen_mapper(infile)
    inputs = list(range(50))
    random.shuffle(inputs)
    for i in inputs:
        assert mapper.predict(i) == approx(i/10)

    inputs = list(range(50, 100))
    random.shuffle(inputs)
    for i in inputs:
        assert mapper.predict(i) == approx(5 + (i-50)/5)

    assert mapper.predict(-1) == approx(0)
    assert mapper.predict(101) == approx(15)


def test_random():
    random.seed(71873)
    rand_sample()
    seeds = random.sample(range(0, 100_000), 100)
    for seed in seeds:
        print(seed)
        random.seed(seed)
        rand_sample()


def rand_sample():
    # build up StringIO object
    num_samples = random.randrange(1, 25)
    slopes = [random.random() for _ in range(num_samples)]
    starts = sorted(random.sample(range(100), num_samples+1))

    infile = "header\n"
    values = [0]
    for i in range(num_samples):
        infile += f'.\t{starts[i]}\t.\t{values[-1]:0.5f}\n'
        values.append(values[-1] + slopes[i] * (starts[i+1] - starts[i]))
    infile += f'.\t{starts[-1]}\t.\t{values[-1]:0.5f}\n'

    mapper = translate_distance.gen_mapper(StringIO(infile))
    inputs = list(range(50))
    random.shuffle(inputs)
    for i in inputs:
        if i < starts[0]:
            assert mapper.predict(i) == 0, mapper.show()
        elif i >= starts[-1]:
            assert mapper.predict(i) == approx(values[-1], rel=1e-3)
        else:
            found = False
            for ind in range(len(starts)):
                if starts[ind] <= i < starts[ind+1]:
                    found = True
                    value = values[ind] + slopes[ind] * (i - starts[ind])
                    assert mapper.predict(i) == approx(value, rel=1e-3)
            assert found
