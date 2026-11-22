#!/usr/bin/env python3
"""Check the depth-six obstruction to triangular Vardy layer stacks."""

from collections import Counter, deque


WIDTH = 14
MASK = (1 << WIDTH) - 1
ROWS = (471, 1787, 2543, 2687, 2715, 2727)


def image_table(rows):
    image = [0] * (1 << WIDTH)
    for word in range(1, len(image)):
        low = (word & -word).bit_length() - 1
        image[word] = image[word & (word - 1)] ^ rows[low]
    return image


def coset_depths(rows):
    columns = (*rows, *(1 << i for i in range(WIDTH)))
    distance = [-1] * (1 << WIDTH)
    distance[0] = 0
    queue = deque([0])
    while queue:
        syndrome = queue.popleft()
        for column in columns:
            reached = syndrome ^ column
            if distance[reached] == -1:
                distance[reached] = distance[syndrome] + 1
                queue.append(reached)
    return distance


def main():
    expected = {0: 1, 1: 28, 2: 378, 3: 3276, 4: 7735, 5: 4888, 6: 78}
    for initial in ROWS:
        a = [((initial << i) | (initial >> (WIDTH - i))) & MASK for i in range(WIDTH)]
        transpose = [sum(((row >> i) & 1) << j for j, row in enumerate(a)) for i in range(WIDTH)]
        for rows in (a, transpose):
            image = image_table(rows)
            assert all((x.bit_count() + y.bit_count()) % 2 == 0 for x, y in enumerate(image))
            assert min(x.bit_count() + image[x].bit_count() for x in range(1, len(image))) == 8
            depths = coset_depths(rows)
            assert Counter(depths) == expected
            deep = [syndrome for syndrome, depth in enumerate(depths) if depth == 6]
            assert all(depths[x ^ y] < 6 for x in deep for y in deep)
    print("Verified the depth-six obstruction for six Vardy circulants and their duals.")


if __name__ == "__main__":
    main()
