import itertools, numbers, operator, math

def partition(it, n):
    return itertools.takewhile(bool, (list(itertools.islice(it, n)) for _ in itertools.count(0)))

def pairwise_2(iterable):
    previt = next(iterable)
    while True:
        nextit = next(iterable)
        yield (previt, nextit)
        previt = nextit

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)

def norm(v):
    x, y = v
    return (float(x)**2+float(y)**2)**0.5

def IsEqual(a, b, delta=0.0000001):
    return abs(b[0] - a[0]) < delta and abs(b[1] - a[1]) < delta

def readvicsekdata(fp, start=0, stop=-1, step=1, smallsteps=(0,)):
    it = itertools.islice(fp, 1, None)

    for i in itertools.count(0):
        
        try:
            positions_line = next(it)
            angles_line = next(it)
        except StopIteration:
            break

        if i < start:
            continue
        elif i >= stop and stop >= 0:
            break

        if i % step not in smallsteps:
            continue

        positions = list(partition(itertools.imap(float, positions_line.strip().split()), 2))
        angles = map(float, angles_line.split())

        yield (i, positions, angles)

def readlendata(fp, start=0, stop=-1, step=1, smallsteps=(0,)):
    it = itertools.islice(fp, 1, None)

    for i in itertools.count(0):
        
        try:
            lengths_line = next(it)
        except StopIteration:
            break

        if i < start:
            continue
        elif i >= stop and stop >= 0:
            break

        if i % step not in smallsteps:
            continue

        lengths = list(partition(itertools.imap(float, lengths_line.strip().split()), 2))

        yield (i, lengths)


def writevicsekdata(fw, positions, angles):
    fw.write(" ".join("{} {}".format(p[0], p[1]) for p in positions))
    fw.write("\n")
    fw.write(" ".join(itertools.imap(str, angles)))
    fw.write("\n")

def writelendata(fw, lengths):
    fw.write(" ".join("{} {}".format(l[0], l[1]) for l in lengths))
    fw.write("\n")

def readiddata(fr, start=0, stop=-1):
    next(fr)
    for line in fr:
        yield map(int, line.split(" "))

def readtrackingdata(fr, start=0, stop=-1):
    it = readvicsekdata(fr, start=start, stop=stop, step=1)
    yield next(it)
    try:
        while True:
            yield (next(it), next(it))
    except StopIteration:
        pass

def tracking_ids(fr, start=0, stop=-1):
    it = readtrackingdata(fr, start, stop)
    _, positions,_ = next(it)
    line = 4
    n = len(positions)
    current_id = n
    prev_ids = range(n)
    yield (prev_ids, positions)
    for (_, positions1, _), (_, positions2, _) in it:
        new_ids = []
        assert(len(prev_ids) == len(positions1))
        for p2 in positions2:
            found_particle = False
            for i, p1 in enumerate(positions1):
                if IsEqual(p1, p2, 0.0001):
                    assert prev_ids[i] not in new_ids, "Duplicate position {} found in line {}".format(p1, line+2)
                    new_ids.append(prev_ids[i])
                    found_particle = True
                    break
            if not found_particle:
                new_ids.append(current_id)
                current_id += 1
        yield (new_ids, positions2)
        prev_ids = new_ids
        line += 4

class Vector2D:

    def __init__(self, x, y, L = None):
        self.x = x
        self.y = y
        self.L = L

    def __add__(self, rhs):
        return Vector2D(self.x+rhs.x, self.y+rhs.y, self.L)

    def __sub__(self, rhs):
        return Vector2D(self.x-rhs.x, self.y-rhs.y, self.L)

    def __mul__(self, rhs):
        if isinstance(rhs, Vector2D):
            return self.x*rhs.x + self.y*rhs.y
        elif isinstance(rhs, numbers.Number):
            return Vector2D(self.x * rhs, self.y * rhs, self.L)
        else:
            raise TypeError()

    def __div__(self, rhs):
        if isinstance(rhs, numbers.Number):
            return Vector2D(self.x / rhs, self.y / rhs, self.L)
        else:
            raise TypeError()

    def Norm(self):
        return (self.x**2+self.y**2)**0.5

    def Normalized(self):
        n = self.Norm()
        if n == 0:
            raise RuntimeError("Vector has norm 0")
        return Vector2D(self.x / n, self.y / n, self.L)

    def __str__(self):
        return "{} {}".format(self.x, self.y)

    def __repr__(self):
        return "Vector2D({},{})".format(self.x, self.y)

    def centered(self, L=None):
        L = L or self.L
        return Vector2D(math.fmod(self.x, L), math.fmod(self.y, L), L)

    def distance_periodic(self, rhs, L=None):
        L = L or self.L
        absx = abs(rhs.x-self.x)
        absy = abs(rhs.y-self.y)
        minx = min(absx, L-absx)
        miny = min(absy, L-absy)
        return math.sqrt(minx*minx + miny*miny)

def combine_dicts(a, b, op=operator.add):
    return dict((k, op(a[k], b[k])) for k in (set(b) & set(a)))
