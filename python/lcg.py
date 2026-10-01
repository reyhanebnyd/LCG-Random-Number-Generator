"""Linear congruential generator (LCG): x_{n+1} = (a * x_n + c) mod m."""


def lcg_sequence(seed, a, c, m):
    """Generate the LCG cycle starting at seed, stopping when the seed recurs.

    Returns the list [seed, x_1, x_2, ..., x_period] where x_period == seed.
    len(result) - 1 is the period of the generator for these parameters.
    """
    sequence = [seed]
    x = seed
    while True:
        x = (a * x + c) % m
        sequence.append(x)
        if x == seed:
            break
    return sequence


if __name__ == "__main__":
    m, a, c, seed = 17, 7, 3, 6
    sequence = lcg_sequence(seed, a, c, m)
    print(sequence)
    print(f"period: {len(sequence) - 1}")
