import math


def radians_to_degrees(radians):
    """Convert radians to degrees."""
    return radians * (180.0 / math.pi)


class IterateTractAndRing:
    """An iterator to traverse through tract ids and ring ids in a skymap.

    This iterator yields tuples of (tract_index, ring_index) for each tract in the skymap.

    Note that the first tract is the south pole, and is considered to be in ring -1, but is
    assumed to be excluded from ``ring_nums``. If ``add_poles`` is True, the iterator will
    include the south pole (tract 0) and the north pole (last tract) in the iteration.

    Parameters
    ----------
    ring_nums : list of int
        A list where each element represents the number of tracts in each ring (eg, [5, 10, 15] will
        be interpreted as 5 tracts in ring 0, 10 in ring 1, and 15 in ring 2).
    add_poles : bool, optional
        If True, the iterator will include the south pole (tract 0) and the north pole (last tract).
        If False, it will only iterate through the rings. Default is True.
    """

    def __init__(self, ring_nums, add_poles=True):
        self.ring_nums = ring_nums
        if add_poles:
            self.total_tracts = sum(ring_nums) + 2
            self.current_tract = 0
            self.current_ring = -1
        else:
            self.total_tracts = sum(ring_nums)
            self.current_tract = 1
            self.current_ring = 0

    def __iter__(self):
        return self

    def __next__(self):
        # End iteration if we have processed all tracts.
        if self.current_tract >= self.total_tracts:
            raise StopIteration
        tract_and_ring = (self.current_tract, self.current_ring)

        # Increase tract.
        self.current_tract += 1

        # Check if we need to move to the next ring.clear
        if self.current_ring == -1:
            self.current_ring += 1
        elif self.current_tract > sum(self.ring_nums[: self.current_ring + 1]):
            self.current_ring += 1

        return tract_and_ring
